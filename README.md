# GenOMICC_StrepA

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