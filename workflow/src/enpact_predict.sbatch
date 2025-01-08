#!/bin/bash

#SBATCH --job-name=enpact_predict
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=8G
#SBATCH --error=/beagle3/haky/users/temi/projects/enformer-predict/logs/enpact_predict.err
#SBATCH --output=/beagle3/haky/users/temi/projects/enformer-predict/logs/enpact_predict.out
#SBATCH --account=pi-haky

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/beagle3/haky/users/shared_software/TFXcan-pipeline-tools/lib

python3 src/enpact_predict.py \
    --matrix /beagle3/haky/users/temi/projects/enformer-predict/enformer_test/enformer_minimal/enformer_minimal.2024-10-11.processed.matrix.tsv.gz \
    --weights /beagle3/haky/users/temi/projects/Enpact/files/ENPACT_734_2024-07-26.compiled_weights.lambda.1se.txt.gz \
    --metadata /beagle3/haky/users/temi/projects/enformer-predict/enformer_test/enformer_minimal/enformer_minimal.2024-10-11.processed.metadata.tsv \
    --split \
    --output_basename /beagle3/haky/users/temi/projects/enformer-predict/enformer_test/enpact_predictions/enformer_minimal.2024-10-11 \

echo "INFO - Done!"