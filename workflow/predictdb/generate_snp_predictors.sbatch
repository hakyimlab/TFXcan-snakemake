#!/bin/bash
#SBATCH --job-name=lEnpact_training
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --error=slurm_logs/lEnpact_training.err
#SBATCH --output=slurm_logs/lEnpact_training.log
#SBATCH --account=pi-haky

date
source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/predictdb-env

export JAVA_HOME=/beagle3/haky/users/temi/software/conda_envs/predictdb-env/

phenotype=${1}
s_path=${2}
loci_annotation=${3}
enpact_scores=${4}
ref_genotypes=${5}
ref_annotations=${6}
nextflow_executable=${7}

echo "Running ${phenotype}"
# outdir
mkdir -p ${s_path} && cd ${s_path}

echo $(which nextflow)

export NXF_DISABLE_CHECK_LATEST=true

# run nextflow
nextflow run ${nextflow_executable} \
    --gene_annotation "${loci_annotation}" \
    --snp_annotation "${ref_annotations}" \
    --genotype "${ref_genotypes}" \
    --gene_exp "${enpact_scores}" \
    --prefix ${phenotype} \
    --outdir ${s_path}/models/ \
    -resume \
    -profile midway3 \
    --nested_cv
    #--keepIntermediate 