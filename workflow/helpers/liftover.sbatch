#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=liftover
#SBATCH --account=pi-haky
#SBATCH --output=/project/haky/users/temi/projects/TFXcan-snakemake/output/logs/liftover.out
#SBATCH --error=/project/haky/users/temi/projects/TFXcan-snakemake/output/logs/liftover.err
#SBATCH --time=06:00:00 
#SBATCH --partition=beagle3

# module load openmpi
# module load parallel

date

slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}/run
mkdir -p ${SLURM_O_WORKDIR}
cd $SLURM_O_WORKDIR
echo "INFO - Working directory is $SLURM_O_WORKDIR"
echo "INFO - Jobid is $SLURM_JOBID"
echo "INFO - Running on host `hostname`"

source ~/.bashrc
conda activate /project/haky/users/temi/software/conda_envs/r-env 


Rscript /project/haky/users/temi/projects/TFXcan-snakemake/workflow/helpers/liftover.R --input_file ${1} --chain_file ${2} --output_file ${3}


printf "INFO - Liftover completed\n"