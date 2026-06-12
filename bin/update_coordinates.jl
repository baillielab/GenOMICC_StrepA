using CSV
using DataFrames

function _temp_analysis()
    combine(groupby(pvar, [:KGP_QC_RESULT]), nrow)
    nsamples = 31
    transform!(pvar,
        :ALT_FREQS => (x -> sqrt.(x .* (1 .- x) / nsamples)) => "ALT_FREQS_SE"
    )
    transform!(pvar,
        [:ALT_FREQS, :ALT_FREQS_SE, :KGP_ALT_FREQ, :KGP_QC_RESULT] =>
            ByRow((f, se, kf, qc) -> qc ∈ ("DROP", "FLIP", "SWITCH_FLIPPED")  ? missing : qc == "MATCH" ? f - 1.96se < kf < f + 1.96se : f - 1.96se < 1 - kf < f + 1.96se) 
                => "KGP_IN_BOUND"
    )
    combine(groupby(pvar, [:KGP_QC_RESULT, :KGP_IN_BOUND]), nrow)
    
    switches = subset(pvar, :KGP_QC_RESULT => x -> x .== "SWITCH", skipmissing=true)
    matches = subset(pvar, :KGP_QC_RESULT => x -> x .== "MATCH", skipmissing=true)
    freqs_0 = collect(0:0.001:1)
    ci_up = [x + 1.96sqrt(x*(1-x)/31) for x in freqs_0]
    ci_low = [x - 1.96sqrt(x*(1-x)/31) for x in freqs_0]
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="KGP", ylabel="HC24 Dataset", title="Matched Alleles")
    scatter!(ax, 1 .- switches.KGP_ALT_FREQ, switches.ALT_FREQS, color=(:blue, 0.5), label="Switches")
    scatter!(ax, matches.KGP_ALT_FREQ, matches.ALT_FREQS, color=(:orange, 0.5), label="Matches")
    lines!(ax, freqs_0, ci_up, color=:red)
    lines!(ax, freqs_0, ci_low, color=:red, label="Confidence Band")
    fig[1,2] = axislegend(ax)
    fig
    save("frequency_comparison_kgp_hc24.png", fig)

end

function load_kgp_freqs(kgp_freq_file)
    kgp_freqs = CSV.read(kgp_freq_file, DataFrame)
    select!(kgp_freqs, 
        "#CHROM" => (x -> string.("chr", x)) => "CHROM",
        "ID" => ByRow(x -> parse(Int, split(x, ":")[2])) => "POS",
        "ID" => "KGP_ID",
        ["REF", "ALT", "ALT_FREQS"] => ByRow((ref, alt, alt_freq) -> alt_freq <= 0.5 ? (ref, alt, alt_freq) : (alt, ref, 1 - alt_freq)) => ["KGP_REF", "KGP_ALT", "KGP_ALT_FREQ"]
    )
end

function kgp_qc(row; threshold=0.2)
    if ismissing(row.KGP_ID)
        return "DROP"

    elseif length(row.ALT) > 1 || length(row.REF) > 1 || length(row.KGP_ALT) > 1 || length(row.KGP_REF) > 1
        return "DROP"

    elseif row.REF == row.KGP_REF && row.ALT == row.KGP_ALT
        return "MATCH"

    elseif row.REF == row.KGP_REF && row.ALT == "."
        return "MATCH"

    elseif row.REF == row.KGP_ALT && (row.ALT == "." || row.ALT == row.KGP_REF)
        kgp_alleles = Set([row.KGP_ALT, row.KGP_REF])
        if kgp_alleles == Set(["A", "T"]) || kgp_alleles == Set(["C", "G"])
            if abs(row.ALT_FREQS - 0.5) < threshold
                return "DROP"
            else
                return "FLIP"
            end
        else
            return "SWITCH"
        end
    else
        complement = Dict("A" => "T", "T" => "A", "C" => "G", "G" => "C")
        if row.REF == complement[row.KGP_REF] && (row.ALT == "." || row.ALT == complement[row.KGP_ALT]) 
            return "FLIP"
        elseif row.REF == complement[row.KGP_ALT] && (row.ALT == "." || row.ALT == complement[row.KGP_REF])
            return "SWITCH_FLIPPED"
        else
            return "DROP"
        end
    end
end

function load_pvar_with_freqs(pvar_file, freqs_file)
    pvar = CSV.read(pvar_file, DataFrame)
    freqs = CSV.read(freqs_file, DataFrame, select=["ID", "ALT_FREQS", "OBS_CT"])
    return leftjoin!(pvar, freqs, on=["ID"])
end

function write_list_to_file(file, list)
    open(file, "w") do io
        for variant_id in list
            println(io, variant_id)
        end
    end
end

function main(bed_prefix, grch38_coords_file, kgp_freq_file)
    bim_file = "$bed_prefix.bim"
    # Load GRCh37 bim and GRCh38 mapping
    bim = CSV.read(bim_file, DataFrame; header=[:CHROM, :ID, :COORD, :POS, :A1, :A2])
    new_coords = CSV.read(grch38_coords_file, DataFrame; header=[:NEW_CHROM, :START, :NEW_POS, :ID, :SUCCESS])
    # Update coordinates
    leftjoin!(bim, select(new_coords, :ID, :NEW_CHROM, :NEW_POS), on=:ID)
    bim.CHROM = map(zip(bim.CHROM, bim.NEW_CHROM)) do (old_chr, new_chr)
        ismissing(new_chr) ? old_chr : new_chr
    end
    bim.POS = map(zip(bim.POS, bim.NEW_POS)) do (old_pos, new_pos)
        ismissing(new_pos) ? old_pos : new_pos
    end
    # Make tempdir
    tmp_dir = mktempdir()
    tmp_prefix = joinpath(tmp_dir, basename(bed_prefix))
    tmp_bim_file = joinpath(tmp_dir, "$tmp_prefix.bim")
    # Write new bim file
    CSV.write(tmp_bim_file, select(bim, :CHROM, :ID, :COORD, :POS, :A1, :A2), header=false, delim="\t")
    # Filter unmapped variants from bim
    remove_list_file = joinpath(tmp_dir, "remove.txt")
    write_list_to_file(remove_list_file, bim[ismissing.(bim.NEW_POS), :ID])
    run(`plink2 \
        --bed $bed_prefix.bed \
        --fam $bed_prefix.fam \
        --bim $tmp_bim_file \
        --output-chr chrMT \
        --sort-vars \
        --make-pgen \
        --out $tmp_prefix \
        --exclude $remove_list_file`
    )
    # KGP-based QC
    run(`plink2 \
        --pfile $tmp_prefix \
        --output-chr chrMT \
        --freq \
        --out $tmp_prefix`
    )
    pvar = load_pvar_with_freqs(string(tmp_prefix, ".pvar"), string(tmp_prefix, ".afreq"))
    kgp_freqs = load_kgp_freqs(kgp_freq_file)
    leftjoin!(pvar, kgp_freqs, on=["#CHROM" => "CHROM", "POS"])
    pvar.KGP_QC_RESULT = map(eachrow(pvar)) do row
        kgp_qc(row; threshold=0.2)
    end
    # Update unknown ALT in pvar when possible
    allele_map = [ismissing(ref) || ismissing(alt) ? missing : Dict(ref => alt, alt => ref) for (ref, alt) in zip(pvar.KGP_REF, pvar.KGP_ALT)]
    pvar.ALT = map(enumerate(zip(pvar.REF, pvar.ALT))) do (index, (ref, alt))
        if alt == "."
            index_allele_map = allele_map[index]
            if ismissing(index_allele_map)
                alt
            elseif haskey(index_allele_map, ref)
                index_allele_map[ref]
            else
                alt
            end
        else
            alt
        end
    end
    CSV.write(string(tmp_prefix, ".pvar"), pvar[!, ["#CHROM", "POS", "ID", "REF", "ALT"]], delim="\t")
    # Make drop and flip lists
    snps_to_drop_file = joinpath(tmp_dir, "variants.qc_drop.txt")
    write_list_to_file(snps_to_drop_file, pvar[pvar.KGP_QC_RESULT .== "DROP", :ID])
    snps_to_flip_file = joinpath(tmp_dir, "variants.qc_flip.txt")
    write_list_to_file(snps_to_flip_file, pvar[pvar.KGP_QC_RESULT .== "FLIP", :ID])
    # Make bed because plink2 does not support --flip
    run(`plink2 \
        --pfile $tmp_prefix \
        --make-bed \
        --out $tmp_prefix
    `)

    run(`plink \
        --bfile $tmp_prefix \
        --flip $snps_to_flip_file \
        --output-chr chrMT \
        --exclude $snps_to_drop_file \
        --make-bed \
        --allow-extra-chr \
        --freq \
        --out $bed_prefix.qced
    `)
    return 0
end

main(ARGS[1], ARGS[2], ARGS[3])