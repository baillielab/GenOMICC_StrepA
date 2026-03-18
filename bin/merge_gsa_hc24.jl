using DataFrames
using CSV

function make_filtered_bcf(input_prefix, output_prefix, shared_variants_file; keep_sample_file=nothing)
    cmd_vec = ["plink2", 
        "--pfile", input_prefix, 
        "--export", "bcf", "vcf-dosage=DS", "id-delim=_",
        "--out", output_prefix,
        "--extract", shared_variants_file
    ]
    if keep_sample_file !== nothing
        append!(cmd_vec, ["--keep", keep_sample_file])
    end
    run(Cmd(cmd_vec))
    run(Cmd(["bcftools", "index",  string(output_prefix, ".bcf")]))
end

function main(gsa_dir, hc24_dir; output_dir="merged_HC24_GSA")
    isdir(output_dir) || mkdir(output_dir)
    for chr in 1:22
        @info "Merging GSA/HC24 chromosome: $chr"
        # Filter files to shared variants
        ## Identify shared variants
        gsa_input_prefix = joinpath(gsa_dir, "chr$chr")
        gsa_variants = CSV.read(string(gsa_input_prefix, ".pvar"), DataFrame, comment="##")
        hc24_input_prefix = joinpath(hc24_dir, "chr$chr")
        hc24_variants = CSV.read(string(hc24_input_prefix, ".pvar"), DataFrame, comment="##")
        shared_variants = intersect(gsa_variants.ID, hc24_variants.ID)
        workdir = mktempdir()
        shared_variants_file = joinpath(workdir, "shared_variants.txt")
        open(shared_variants_file, "w") do io
            for variant in shared_variants
                println(io, variant)
            end
        end
        ## GSA has duplicated samples (to trick TOPMed), we filter them.
        gsa_samples = CSV.read(string(gsa_input_prefix, ".psam"), DataFrame)
        gsa_samples_to_keep = gsa_samples[.!endswith.(gsa_samples[!, "#IID"], "copy"), "#IID"]
        gsa_samples_to_keep_file = joinpath(workdir, "gsa_keep_samples.txt")
        open(gsa_samples_to_keep_file, "w") do io
            for sample in gsa_samples_to_keep
                println(io, sample)
            end
        end
        gsa_shared_prefix = joinpath(workdir, "chr$chr.gsa")
        make_filtered_bcf(gsa_input_prefix, gsa_shared_prefix, shared_variants_file; keep_sample_file=gsa_samples_to_keep_file)
        ## HC24 samples are all retained
        hc24_shared_prefix = joinpath(workdir, "chr$chr.hc24")
        make_filtered_bcf(hc24_input_prefix, hc24_shared_prefix, shared_variants_file)
        # Merge the two files with bcftools
        merged_bcf_file = joinpath(output_dir, "chr$chr.bcf")
        run(Cmd([
            "bcftools", "merge", 
            string(hc24_shared_prefix, ".bcf"), 
            string(gsa_shared_prefix, ".bcf"),
            "-Ob", "-o", merged_bcf_file,
            "--write-index"
        ]))
    end
end

main(ARGS[1], ARGS[2])