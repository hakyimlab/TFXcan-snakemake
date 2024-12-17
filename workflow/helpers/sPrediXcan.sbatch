#!/bin/bash

# Author: Temi
# Description: S-PrediXcan on TFPred Scores
# Usage: sbatch s-predixcan.sbatch
# Date: Mon Aug 7 2023
# Dependencies: 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=spredixcan_tfscores
#SBATCH --account=pi-haky
#SBATCH --output=/project2/haky/temi/projects/l-TFPred/logs/spredixcan_tfscores.out
#SBATCH --error=/project2/haky/temi/projects/l-TFPred/logs/spredixcan_tfscores.err
#SBATCH --time=24:00:00	
#SBATCH --partition=bigmem

# module load openmpi
# module load parallel

slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}/run
mkdir -p ${SLURM_O_WORKDIR}
echo Working directory is $SLURM_O_WORKDIR
cd $SLURM_O_WORKDIR

echo Jobid: $SLURM_JOBID
echo Running on host `hostname`

printf "Starting to run\n"

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/imlabtools 

# variables
model_name="EUR_AR_Prostate_linear"
base_dir="/project2/haky/temi/projects/l-TFPred/experiments/prostate_cancer"
model_dir="/project2/haky/temi/projects/l-TFPred/models"
model_db=${model_dir}/${model_name}/filtered_db/gtex_v7_${model_name}_filtered_signif.db
exec_file=/beagle3/haky/users/temi/software/MetaXcan/software/SPrediXcan.py

# note that the ...
${exec_file} \
    --model_db_path ${model_db} \
    --model_db_snp_key varID \
    --covariance ${model_dir}/${model_name}/database/Covariances.txt \
    --gwas_folder ${base_dir}/data/GWAS \
    --gwas_file_pattern ".*_Schumacher.gwas_ss.txt.gz" \
    --chromosome_column chrom \
    --snp_column variant_id \
    --effect_allele_column A1 \
    --non_effect_allele_column A2 \
    --position_column pos \
    --freq_column maf \
    --beta_column beta \
    --se_column se \
    --pvalue_column p_value \
    --zscore_column zscore \
    --keep_non_rsid \
    --output_file ${base_dir}/output/${model_name}.TFScores.Prostate.spredixcan.csv \
    --verbosity 9 \
    --throw

    #     --gwas_folder ${base_dir}/data/GWAS \
    # --gwas_file_pattern ".*_Schumacher.gwas_ss.txt.gz" \


    # awk '{gsub(/:/, "_"); print}' Covariances.varID.txt > Covariances.varID.txt

    # sed -i -e 's/:/_/g' Covariances.varID.txt

    #sed -i -e 's/:/_/g' Covariances.txt > Covariances.varID.txt