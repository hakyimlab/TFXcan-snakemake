
# TFXcan

This pipeline tests TF binding-GWAS trait associations using SNP-based predictors of TF binding.

## Usage:

1. conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
2. snakemake -s snakefile.smk --configfile config/pipeline.yaml --profile profiles/simple/ --resources load=45

## Software: 

This pipeline depends on a number of software to do the following:

1. Finemap GWAS SNPs (SuSie); optional because you can decide not to finemap and just use the top SNPs per locus

2. Predict with Enformer (this dependency is optional) (Enformer, GPUs, pytorch)

3. Train models of TF binding that is linear on SNPs (Nextflow, predictDB)

4. Test TF binding-GWAS trait association (PrediXcan, Summary-PrediXcan, MetaXcan)

We suggest the following to have a hitch-free environment:

1. Use conda to create an environment and install the software with the [environment file](/beagle3/haky/users/shared_software/TFXcan-pipeline-tools)

## Input:

In general, the pipeline expects:

1. A yaml config or parameters file. Details are [here]()

2. A metadata sheet of the GWAS summary statistics. Details are [here]()

## Command:


## Output:
The output of the pipeline is the association results of the GWAS trait with the TF binding, and it can be found in the `data/.../output` folder. The output is a summary ***.TFXcan.csv file of the association results file with the following columns:


#### Notes: 

* GWAS summary stats with the following headers (others headers are allowed but will be ignored): 

    - chrom: 1,2,3, e.t.c (No chromosomes X, Y, or M e.t.c)

    - pos: 134 (bp coordinates)

    - variant_id: 1_345_A_G i.e. chrom_pos_ref_alt

    - pval: GWAS pvalues

    - zscore: GWAS zscores; you can pre-calculate this from the beta and standard errors (beta/se)

* The framework assumes that genomic coordinates are in hg38 coordinates
