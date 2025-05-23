#!/bin/bash

# Author: Temi
# Description: S-PrediXcan on TFPred Scores for any phenotype
# Usage: sbatch --export=gwas_file=phenotype_gwas_sumstats.txt.gz,phenotype=phenotype spredixcan_ADT_tfscores.sbatch
# Date: Mon Aug 7 2023
# Variables: gwas_file, phenotype


date

slurm_workdir=${SLURM_SUBMIT_DIR}
# SLURM_O_WORKDIR=${slurm_workdir}/run
# mkdir -p ${SLURM_O_WORKDIR}
# cd $SLURM_O_WORKDIR
echo "INFO - Working directory is $SLURM_SUBMIT_DIR"
echo "INFO - Jobid is $SLURM_JOBID"
echo "INFO - Running on host `hostname`"

source ~/.bashrc

# variables
phenotype=${1}
model_db=${2}
outputfile=${3}
gwas_folder=${4}
gwas_pattern=${5}
covariances=${6}
exec_file=${7}
conda_environment=${8}

conda activate ${conda_environment} 

printf "INFO - Running S-TFXcan for %s\n" ${phenotype}

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



    #     --gwas_folder ${base_dir}/data/GWAS \
    # --gwas_file_pattern ".*_Schumacher.gwas_ss.txt.gz" \


    # awk '{gsub(/:/, "_"); print}' Covariances.varID.txt > Covariances.varID.txt

    # sed -i -e 's/:/_/g' Covariances.varID.txt

    #sed -e 's/:/_/g' Covariances.txt > Covariances.varID.txt

    # gzip -d predict_db_EUR_AR_Prostate_logistic_filtered.txt.gz
    # sed -e 's/:/_/g' predict_db_EUR_AR_Prostate_logistic_filtered.txt > Covariances.varID.txt