#!/bin/bash

#SBATCH --job-name=enformer_process
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=8G
#SBATCH --error=/beagle3/haky/users/temi/projects/enformer-predict/logs/enformer_process.err
#SBATCH --output=/beagle3/haky/users/temi/projects/enformer-predict/logs/enformer_process.out
#SBATCH --account=pi-haky

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib

python3 src/enformer_process.py \
    --merged_h5_file /beagle3/haky/users/temi/projects/enformer-predict/enformer_test/enformer_minimal/predictions_2024-10-11/enformer_minimal.predictions.h5 \
    --process_by_haplotype \
    --process_function 'sum' \
    --output_basename /beagle3/haky/users/temi/projects/enformer-predict/enformer_test/enformer_minimal/enformer_minimal.2024-10-11.processed

echo "INFO - Done!"