#!/bin/bash

#SBATCH --job-name=enformer_merge
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/enformer_merge.err
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/enformer_merge.out
#SBATCH --account=pi-haky

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib

python3 ${1} \
    --config ${2} \
    --output_directory ${3} \
    --output_filename ${4}

# python3 src/enformer_merge.py \
#     --config /beagle3/haky/users/temi/projects/enformer-predict/enformer_test/enformer_minimal/predictions_2024-10-11/metadata/enformer_minimal.aggregation_config.yaml \
#     --output_directory /beagle3/haky/users/temi/projects/enformer-predict/enformer_test/enformer_minimal/predictions_2024-10-11/ \
#     --output_filename enformer_minimal.predictions.h5

# python3 src/enformer_merge.py \
#     --config /beagle3/haky/users/temi/projects/chrom-enpact/output/aracena200/predictions_2024-10-15/metadata/aracena200.aggregation_config.yaml \
#     --output_directory /beagle3/haky/users/temi/projects/chrom-enpact/output/aracena200/predictions_2024-10-15/ \
#     --output_filename aracena200.predictions.h5

# make sure that in your config, provider is set to "local" before running this script
# python3 src/enformer_merge.py \
#     --predictions_logfile /beagle3/haky/users/temi/projects/chrom-enpact/output/aracena2000/predictions_2024-10-11/predictions_log/aracena2000.prediction_log.txt \
#     --predictions_directory /beagle3/haky/users/temi/projects/chrom-enpact/output/aracena2000/predictions_2024-10-11/predictions \
#     --expected_shape 1,5313 \
#     --output_directory /beagle3/haky/users/temi/projects/chrom-enpact/output/aracena2000/predictions_2024-10-11 \
#     --output_filename aracena2000.predictions.h5

echo "INFO - Done!"