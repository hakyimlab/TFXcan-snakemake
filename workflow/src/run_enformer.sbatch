#!/bin/bash

#SBATCH 
#SBATCH --partition=beagle3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=enformer_predict
#SBATCH --account=pi-haky
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:2
#SBATCH --output=/project/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/slurm-%j.out
#SBATCH --error=/project/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/slurm-%j.err

date
source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib
echo Found `nvidia-smi --list-gpus | wc -l` gpus

# make sure that in your config, provider is set to "local" before running this script
python3 ${1} --parameters ${2}

#  sbatch workflow/src/run_enformer.sbatch /beagle3/haky/users/shared_pipelines/enformer_pipeline_aggregate/scripts/enformer_predict.py data/prostate_cancer_risk_2024-09-23/enformer_parameters/enformer_parameters_prostate_cancer_risk_pc_risk.json