# GenOMICC_StrepA

## Install the project dependencies

- Install [Julia](https://docs.julialang.org/en/v1/manual/installation/)

- Install Julia dependencies:

```bash
julia --project=. --startup-file=no -e'using Pkg; Pkg.instantiate()'
```

## Imputing EGAD00001004558

We integrate data from [this study](https://ega-archive.org/studies/EGAS00001003421). In principle, this has already been done and the imputed genotypes can be downloaded from ODAP (ask admin member) or accessed from datastore at `/exports/cmvm/datastore/eb/groups/baillielab2/strepA_EGAD00001004558` (). 

Here are the steps that were performed to obtain these imputed genotypes. We assume the source data is located in the following directory:

```bash
EGAD_DATA_DIR=/gpfs/igmmfs01/eddie/ISARIC4C/olivier/data/EGAD00001004558
KGP_DATA_DIR=/home/olabayle/isaric/olivier/data/kgp-merged-unrelated-or3
mkdir imputation_work
mkdir imputation_work/mappings
# Get sexes
awk -F',' 'NR>1 {print $6 "\t" $6 "\t" $5}' $EGAD_DATA_DIR/GSA2024_1140_DATA/samples.csv > imputation_work/sexes.tsv
# Make BED file for: STREPGENE_HG19_GSA
plink2 \
--vcf $EGAD_DATA_DIR/EGAF00002318977/STREPGENE_HG19_GSA.vcf \
--id-delim _ \
--update-sex imputation_work/sexes.tsv \
--make-bed \
--out imputation_work/STREPGENE_HG19_GSA
plink2 \
--bfile imputation_work/STREPGENE_HG19_GSA \
--freq \
--out imputation_work/STREPGENE_HG19_GSA
# Make file to map to GRCh38
awk '{
  chr=$1
  pos=$4
  start=pos-1
  end=pos
  print "chr"chr"\t"start"\t"end"\t"$2
}' imputation_work/STREPGENE_HG19_GSA.bim > imputation_work/STREPGENE_HG19_GSA.variants.GRCh37.bed

# Make BED file for: STREPGENE_HG19_GSA
plink2 \
--vcf $EGAD_DATA_DIR/EGAF00002318978/STREPGENE_HG19_HC24.vcf \
--id-delim _ \
--update-sex imputation_work/sexes.tsv \
--make-bed \
--out imputation_work/STREPGENE_HG19_HC24
plink2 \
--bfile imputation_work/STREPGENE_HG19_HC24 \
--freq \
--out imputation_work/STREPGENE_HG19_HC24
# Make file to map to GRCh38
awk '{
  chr=$1
  pos=$4
  start=pos-1
  end=pos
  print "chr"chr"\t"start"\t"end"\t"$2
}' imputation_work/STREPGENE_HG19_HC24.bim > imputation_work/STREPGENE_HG19_HC24.variants.GRCh37.bed
```

Now go on [ucsc](https://genome.ucsc.edu/cgi-bin/hgLiftOver) and lift over the `variants.GRCh37.bed` files to GRCh38. We can update the coordinates in the bim file using the `bin/update_coordinates.jl`:

```bash
julia --project=. --startup-file=no bin/update_coordinates.jl imputation_work/STREPGENE_HG19_GSA imputation_work/mappings/STREPGENE_HG19_GSA.variants.GRCh38.bed $KGP_DATA_DIR/kgp.merged.unrelated.afreq
plink2 --bfile imputation_work/STREPGENE_HG19_GSA.qced --make-bed --set-all-var-ids @:#:\$r:\$a --chr 1-22,X --freq --out imputation_work/STREPGENE_HG19_GSA.qced.ready
# This is because TOPMed requires at least 20 samples to proceed
cp imputation_work/STREPGENE_HG19_GSA.qced.ready.bed imputation_work/STREPGENE_HG19_GSA.qced.ready.bis.bed
cp imputation_work/STREPGENE_HG19_GSA.qced.ready.bim imputation_work/STREPGENE_HG19_GSA.qced.ready.bis.bim
cp imputation_work/STREPGENE_HG19_GSA.qced.ready.fam imputation_work/STREPGENE_HG19_GSA.qced.ready.bis.fam
awk 'BEGIN{OFS="\t"} {$1 = $1 "_copy"; $2 = $2 "_copy"; print}' imputation_work/STREPGENE_HG19_GSA.qced.ready.bis.fam > imputation_work/STREPGENE_HG19_GSA.qced.ready.bis.temp.fam
mv imputation_work/STREPGENE_HG19_GSA.qced.ready.bis.temp.fam imputation_work/STREPGENE_HG19_GSA.qced.ready.bis.fam
plink -bfile imputation_work/STREPGENE_HG19_GSA.qced.ready -bmerge imputation_work/STREPGENE_HG19_GSA.qced.ready.bis --make-bed --out imputation_work/STREPGENE_HG19_GSA.qced.ready.duplicated
```

```bash
julia --project=. --startup-file=no bin/update_coordinates.jl imputation_work/STREPGENE_HG19_HC24 imputation_work/mappings/STREPGENE_HG19_HC24.variants.GRCh38.bed $KGP_DATA_DIR/kgp.merged.unrelated.afreq
plink2 --bfile imputation_work/STREPGENE_HG19_HC24.qced --make-bed --set-all-var-ids @:#:\$r:\$a --chr 1-22,X --freq --out imputation_work/STREPGENE_HG19_HC24.qced.ready
```

Then open a persistent session and run:

```
nextflow run nf-topmed-imputation/main.nf -profile eddie -resume -with-report -with-trace -c config/imputation_STREPGENE_HG19_GSA.config
```

and

```
nextflow run nf-topmed-imputation/main.nf -profile eddie -resume -with-report -with-trace -c config/imputation_STREPGENE_HG19_HC24.config
```

to impute the genotypes.

Finally, the resulting genotypes can be merged via:

```bash
julia --project=. --startup-file=no bin/merge_gsa_hc24.jl results_GSA/ results_HC24/
```
## Merging EGAD00001004558 with GenOMICC

## Reproducing The GWAS

To run the GWAS, we only need the [WDL-GWAS](https://github.com/olivierlabayle/WDL-GWAS) pipeline, make sure you have followed the installation isntructions and run the following:

```bash
input_config=config/genomicc-inputs.strepA.meta.plink2.json
compiled_input_config=config/genomicc-inputs.strepA.meta.plink2.dx.json
output_dir=/or3_strepa_plink2
commit=98c732469053100dff96cf3b2a4ebb9c70fafa43

# Clone WDL-GWAS pipeline
git clone https://github.com/olivierlabayle/WDL-GWAS
cd WDL-GWAS && git checkout $commit && cd ..

# Compile the workflow
java -jar $DX_COMPILER_PATH compile WDL-GWAS/workflows/gwas.wdl \
    -f -project $RAP_PROJECT_ID \
    -extras WDL-GWAS/config/extras.json \
    -reorg \
    -folder /workflows/gwas \
    -inputs ${input_config}

# Run the workflow
dx run -y \
    -f ${compiled_input_config} \
    --priority high \
    --preserve-job-outputs \
    --destination ${output_dir} \
    /workflows/gwas/gwas
```

## Debugging a Run With A DNA Nexus Instance

To connect to a DNA Nexus instance, run:

```bash
instance_type=mem2_ssd1_v2_x8
dx run \
    --instance-type $instance_type \
    -imax_session_length="48h" \
    -y \
    --ssh app-cloud_workstation
```

You can download files via:

```bash
dx download file_id_1 file_id_2 ... 
```

And run the container with:

```bash
docker_tag=optional_fp
docker run -it --rm -v $PWD:/mnt/data olivierlabayle/wdl-gwas:$docker_tag /bin/bash
```

Then you can run code as normal, for instance to enter the Julia REPL:

```bash
julia --project=/opt/PopGen --startup-file=no --sysimage=/opt/PopGen/sysimage.so --threads=auto
```