
# TFXcan
TFXcan was developed to test transcription factor (TF) binding-GWAS trait associations using SNP-based predictors of TF binding.

These SNP-based predictors are developed using Enformer, a sequence-to-function deep learning model. These predictors are called Enpact predictors. Enpact weights are linear combinations of epigenomic features and are stored [here](./weights). 

## Version: 
TFXcan v4.0

## Usage/Command:

### Using Colab
This is a version that runs a minimal TFXcan. You can edit and input your own GWAS. 
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/hakyimlab/TFXcan-snakemake/blob/main/colab/TFXcan_colab.ipynb)


### Access to a computing cluster

If you have access to a computing cluster, you could download the colab notebook, connect to your GPUs (if you have that resource), and run. 

Otherwise, you can submit a snakemake job as below:
1. conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
2. snakemake -s snakefile.smk --configfile config/pipeline.asthma.children.yaml --profile profiles/simple/ --resources load=45

 The `--resources load=45` flag makes sure that the PredictDB part of the pipeline does not run more than 9 jobs at a time on midway3 i.e 9*5. Any number could have been used but I chose multiples of 5. If your cluster allows you to run more than 100 jobs at a time, you can up this number.

### To use screen [preferred]:

1. screen
2. conda activate << conda environment >>  (see software section)
3. export PATH=$PATH:/project2/haky/temi/software/homer/bin
4. snakemake -s snakefile.smk --configfile config/pipeline.asthma.children.yaml --profile profiles/simple/ --resources load=45

## Software: 

This pipeline depends on a number of software to do the following:

1. Predict with Enformer (this dependency is optional) (Enformer, GPUs, pytorch)
2. Train models of TF binding that is linear on SNPs (Nextflow, predictDB)
3. Test TF binding-GWAS trait association (PrediXcan, Summary-PrediXcan, MetaXcan)

We suggest the following to have a hitch-free environment:

1. Use conda to create an environment and install the software with the [environment file](/beagle3/haky/users/shared_software/TFXcan-pipeline-tools)

All of these software are self-contained in this repository. You only need to install the conda environment. 

## Input:

In general, the pipeline expects:

1. A yaml config or parameters file. Details are [here](./minimal/pipeline_minimal.yaml)
2. A metadata sheet of the GWAS summary statistics. Details are [here](./minimal/minimal_gwas.txt)
3. A number of files needed in the yaml config file. These can be downloaded from [here](https://uchicago.box.com/shared/static/kffo3k9zl16irrveysbnr05ww14qq1vd.gz). This is a direct download link. You will need to decompress this archive. 

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

- The framework assumes that genomic coordinates are in hg38 coordinates

- Weights file: You will need this in the config yaml file (see `enpact_weights`). The weights file is a dataframe of TF binding predictors. You can find and use examples from [here](./weights/).

## Output:
The output of the pipeline is the association results of the GWAS trait with the TF binding, and it can be found in the `data/.../output` folder. The output is a summary ***.TFXcan.csv file of the association results file with the following columns:

    |tfbs|zscore|effect_size|pvalue|var_g|pred_perf_r2|pred_perf_pval|pred_perf_qval|n_snps_used|n_snps_in_cov|n_snps_in_model|


#### Notes:


## Updates:

[X] To predict TF/tissue binding, the pipeline takes in a dataframe of weights. 

|feature|TF/tissue1|TF/tissue2|...|TF/tissueN|
|---|---|---|---|---|
|f1|0.1|0.2|...|0.1|
|f2|0.1|0.2|...|0.1|
|...|...|...|...|...|
|f5313|0.1|0.2|...|0.1|

[X] SNPs are matched with the reference panel and uses the matched SNPs for the PredictDB training. This is to ensure that the SNPs used for the PredictDB training are the same as the SNPs used for the GWAS.

[X] All software necessary for TFXcan are shipped with the pipeline. You only need to install the conda environment.
