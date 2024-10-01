#!/bin/bash

# Author: Temi
# Description: S-PrediXcan on TFPred Scores for any phenotype
# Usage: sbatch --export=gwas_file=phenotype_gwas_sumstats.txt.gz,phenotype=phenotype spredixcan_ADT_tfscores.sbatch
# Date: Mon Aug 7 2023
# Variables: gwas_file, phenotype

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=summary_TFXcan
#SBATCH --account=pi-haky
#SBATCH --output=/project/haky/users/temi/projects/TFXcan-snakemake/output/logs/summary_TFXcan.out
#SBATCH --error=/project/haky/users/temi/projects/TFXcan-snakemake/output/logs/summary_TFXcan.err
#SBATCH --time=06:00:00
#SBATCH --partition=bigmem

date

slurm_workdir=${SLURM_SUBMIT_DIR}
# SLURM_O_WORKDIR=${slurm_workdir}/run
# mkdir -p ${SLURM_O_WORKDIR}
# cd $SLURM_O_WORKDIR
echo "INFO - Working directory is $SLURM_SUBMIT_DIR"
echo "INFO - Jobid is $SLURM_JOBID"
echo "INFO - Running on host `hostname`"

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/imlabtools 

# variables
phenotype=${1}
model_db=${2}
outputfile=${3}
gwas_folder=${4}
gwas_pattern=${5}
covariances=${6}

printf "INFO - Running S-TFXcan for %s\n" ${phenotype}

exec_file=/beagle3/haky/users/temi/software/MetaXcan/software/SPrediXcan.py

# note that the ...
python3 ${exec_file} \
    --model_db_path ${model_db} \
    --model_db_snp_key varID \
    --covariance ${covariances} \
    --gwas_folder ${gwas_folder} \
    --gwas_file_pattern ${gwas_pattern} \
    --chromosome_column chrom \
    --snp_column SNP \
    --effect_allele_column alt \
    --non_effect_allele_column ref \
    --zscore_column zscore \
    --keep_non_rsid \
    --output_file ${outputfile} \
    --verbosity 9 \
    --throw