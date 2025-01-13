#!/bin/bash

#SBATCH --job-name=lEnpact_training
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/lEnpact_training.err
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan-snakemake/slurm_logs/lEnpact_training.log
#SBATCH --account=pi-haky

# source ~/.bashrc
# conda activate /project/haky/users/temi/projects/TFXcan-snakemake/software/predictdb_nf
# module load python
# module load R
# # Rscript -e 'install.packages("/project2/haky/temi/projects/TFXcan-snakemake/workflow/software/optparse_1.7.4.tar.gz", lib="/home/temi/R/x86_64-pc-linux-gnu-library/4.2", repos = NULL, type = "source")'

# #install.packages(pkgname.tar.gz, repos = NULL, type =‘‘source’’)
# #R CMD INSTALL /project2/haky/temi/projects/TFXcan-snakemake/workflow/helpers/optparse_1.7.4.tar.gz
# #export JAVA_HOME=/usr/lib/jvm/jre-1.8.0-openjdk-1.8.0.352.b08-2.el8_6.x86_64

# export JAVA_HOME=/project/haky/users/temi/projects/TFXcan-snakemake/software/predictdb_nf/

date
source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/predictdb-env
# module load python
# module load R
#Rscript -e 'install.packages("/project2/haky/temi/projects/TFXcan-snakemake/workflow/software/optparse_1.7.4.tar.gz", lib="/home/temi/R/x86_64-pc-linux-gnu-library/4.2", repos = NULL, type = "source")'
#export JAVA_HOME=/usr/lib/jvm/jre-1.8.0-openjdk-1.8.0.352.b08-2.el8_6.x86_64

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