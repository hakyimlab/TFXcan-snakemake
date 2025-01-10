#!/bin/bash

#SBATCH --job-name=enformer_predict
#SBATCH --partition=beagle3
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/enformer_predict.err
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/enformer_predict.out
#SBATCH --account=pi-haky

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib

nvidia-smi

# make sure that in your config, provider is set to "local" before running this script
python3 ${1} --parameters ${2}