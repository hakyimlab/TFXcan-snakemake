#!/bin/bash

#SBATCH --job-name=enformer_predict
#SBATCH --partition=beagle3
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=8
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --error=slurm_logs/enformer_predict.err
#SBATCH --output=slurm_logs/enformer_predict.out
#SBATCH --account=beagle3-exusers

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib

nvidia-smi

# make sure that in your config, provider is set to "local" before running this script
python3 ${1} --parameters ${2}