
## TFXcan: 

This pipeline tests TF binding-GWAS trait associations.

## Usage:

1. conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
2. snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/

## Software: 

This pipeline depends on a number of software to do the following:

1. Finemap GWAS SNPs (SuSie)

2. Predict with Enformer (this dependency is optional) (Enformer, GPUs, pytorch)

3. Train models of TF binding that is linear on SNPs (Nextflow, predictDB)

4. Test TF binding-GWAS trait association (PrediXcan, Summary-PrediXcan, MetaXcan)

We suggest the following to have a hitch-free environment:

1. Use conda to create an environment and install the software with the [environment file](/beagle3/haky/users/shared_software/TFXcan-pipeline-tools)

## Main inputs:

In general, the pipeline expects:

1. A yaml config or parameters file. Details are [here]()

2. A metadata sheet of the GWAS summary statistics. Details are [here]()

#### Notes: 

* GWAS summary stats with the following headers: 

    - chrom: 1,2,3, e.t.c (No chromosomes X, Y, or M e.t.c)

    - pos: 134 (bp coordinates)

    - variant_id: 1_345_A_G

    - pval: GWAS pvalues

    - zscore: GWAS zscores; you can pre-calculate this from the beta and standard errors (beta/se)

* The framework assumes that genomic coordinates are in hg38 coordinates



## Main output:
The final output is a summary ***.TFXcan.csv file of the association results.

