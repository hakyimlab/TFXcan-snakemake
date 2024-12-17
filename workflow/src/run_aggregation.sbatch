#!/bin/bash

#SBATCH --partition=beagle3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=enformer_predict
#SBATCH --account=pi-haky
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:2
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/slurm-%j.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/slurm-%j.err

date
source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib
echo Found `nvidia-smi --list-gpus | wc -l` gpus

# make sure that in your config, provider is set to "local" before running this script
python3 ${1} --parameters ${2}

#  python3 workflow/src/aggregate.py --metadata_file data/breast_cancer_2024-12-11/enformer_parameters/aggregation_config_breast_cancer_breast_cancer.json --agg_types aggByCollect --output_directory data/breast_cancer_2024-12-11/aggregated_predictions/breast_cancer --hpc caslake --parsl_executor local