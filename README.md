
# TFXcan
This pipeline tests TF binding-GWAS trait associations using SNP-based predictors of TF binding.

## Version: 
TFXcan v3.0

## Usage/Command:

1. conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
2. snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ --resources load=45

 The `--resources load=45` flag makes sure that the PredictDB part of the pipeline does not run more than 9 jobs at a time on midway3 i.e 9*5. Any number could have been used but I chose multiples of 5. If your cluster allows you to run more than 100 jobs at a time, you can up this number.

### To use screen [preferred]:

1. screen
2. conda activate << conda environment >>  (see software section)
3. export PATH=$PATH:/project2/haky/temi/software/homer/bin
4. snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ --resources load=45

## Software: 

This pipeline depends on a number of software to do the following:

1. Finemap GWAS SNPs (SuSie); optional because you can decide not to finemap and just use the top SNPs per locus [default]
2. Predict with Enformer (this dependency is optional) (Enformer, GPUs, pytorch)
3. Train models of TF binding that is linear on SNPs (Nextflow, predictDB)
4. Test TF binding-GWAS trait association (PrediXcan, Summary-PrediXcan, MetaXcan)

We suggest the following to have a hitch-free environment:

1. Use conda to create an environment and install the software with the [environment file](/beagle3/haky/users/shared_software/TFXcan-pipeline-tools)

All of these software are self-contained in this repository. You only need to install the conda environment. 

## Input:

In general, the pipeline expects:

1. A yaml config or parameters file. Details are [here](./minimal/pipeline_minimal.yaml)
2. A metadata sheet of the GWAS summary statistics. Details are [here](./minimal/minimal_gwas.txt)

The GWAS summary statistics file should have the following columns
(others headers are allowed but will be ignored): 

    |chrom|pos|variant_id|ref|alt|pval|zscore|beta|se|
    |---|---|---|---|---|---|---|---|---|
    |1|134|1_134_A_G|A|G|0.0001|0.1|0.7|0.1|

    - chrom: (character or string) 1,2,3, e.t.c (No chromosomes X, Y, or M e.t.c)

    - pos: (numeric) 134 (bp coordinates)

    - variant_id: chrom_pos_ref_alt

    - pval: (numeric) GWAS pvalues

    - zscore: GWAS zscores; you can pre-calculate this from the beta and standard errors (beta/se)

* The framework assumes that genomic coordinates are in hg38 coordinates

## Output:
The output of the pipeline is the association results of the GWAS trait with the TF binding, and it can be found in the `data/.../output` folder. The output is a summary ***.TFXcan.csv file of the association results file with the following columns:


#### Notes: 


## Updates:

[X] To predict TF/tissue binding, the pipeline takes in a dataframe of weights. 

|feature|TF/tissue1|TF/tissue2|...|TF/tissueN|
|---|---|---|---|---|
|f1|0.1|0.2|...|0.1|
|f2|0.1|0.2|...|0.1|
|...|...|...|...|...|
|f5313|0.1|0.2|...|0.1|

[X] The pipeline now matches SNPs with the reference panel and uses the matched SNPs for the PredictDB training. This is to ensure that the SNPs used for the PredictDB training are the same as the SNPs used for the GWAS.

[X] All software necessary for TFXcan are shipped with the pipeline. You only need to install the conda environment.
