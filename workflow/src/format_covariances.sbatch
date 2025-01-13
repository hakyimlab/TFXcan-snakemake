#!/bin/bash

# Author: Temi
# Description: Format covariances to be predixcan compatible
# Usage: sbatch --export=gwas_file=phenotype_gwas_sumstats.txt.gz,phenotype=phenotype spredixcan_ADT_tfscores.sbatch
# Date: Mon Sep 30 2024
# Variables: gwas_file, phenotype

date

slurm_workdir=${SLURM_SUBMIT_DIR}
# SLURM_O_WORKDIR=${slurm_workdir}/run
# mkdir -p ${SLURM_O_WORKDIR}
# cd $SLURM_O_WORKDIR
echo "INFO - Working directory is $SLURM_SUBMIT_DIR"
echo "INFO - Jobid is $SLURM_JOBID"
echo "INFO - Running on host `hostname`"

#source ~/.bashrc
#conda activate /beagle3/haky/users/temi/software/conda_envs/imlabtools 

# variables
compressed_covariance_file=${1}
sed_covariance_file=${2} # Covariances.varID.txt.gz

# unpigz -c input.vcf.gz | sed 's/^chr//' | pigz > output.vcf.gz
unpigz -ck ${compressed_covariance_file} | sed -e 's/:/_/g' | pigz > ${sed_covariance_file}

# decompress old file 
#  ${decompressed_covariance_file} > ${sed_covariance_file} # replace : with _ in the file
# tar -cz ${sed_covariance_file}.gz ${sed_covariance_file} # compress the file
# rm ${sed_covariance_file} ${decompressed_covariance_file} # remove the old files