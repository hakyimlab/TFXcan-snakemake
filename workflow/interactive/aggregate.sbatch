#!/bin/bash

#SBATCH --job-name=aggregate_predictions
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --error=/project/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/aggregate_predictions.err
#SBATCH --output=/project/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/aggregate_predictions.log
#SBATCH --account=pi-haky

date
source ~/.bashrc

conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools

output_directory=${1}
mkdir -p ${output_directory}
metadata_file=${2}
agg_types="aggByCollect"
hpc="caslake"
parsl_executor="local"

python3 workflow/src/aggregate.py \
    --metadata_file ${metadata_file} \
    --agg_types ${agg_types} \
    --output_directory ${output_directory} \
    --hpc ${hpc} \
    --parsl_executor ${parsl_executor}


# sbatch workflow/interactive/aggregate.sbatch /scratch/midway3/temi/predictions_folder/predictions_folder/prostate_cancer_risk_pc_risk/predictions_2024-09-23/aggregation_folder data/prostate_cancer_risk_2024-09-23/enformer_parameters/aggregation_config_prostate_cancer_risk_pc_risk.json