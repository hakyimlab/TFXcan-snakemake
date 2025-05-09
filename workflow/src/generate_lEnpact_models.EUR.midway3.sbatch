#!/bin/bash

#SBATCH --job-name=lEnpact_training
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --error=/project/haky/users/temi/projects/TFXcan/lEnpact_training.err
#SBATCH --output=/project/haky/users/temi/projects/TFXcan/lEnpact_training.log
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
conda activate /project/haky/users/temi/software/conda_envs/predictdb-env
module load python
module load R
#Rscript -e 'install.packages("/project2/haky/temi/projects/TFXcan-snakemake/workflow/software/optparse_1.7.4.tar.gz", lib="/home/temi/R/x86_64-pc-linux-gnu-library/4.2", repos = NULL, type = "source")'
#export JAVA_HOME=/usr/lib/jvm/jre-1.8.0-openjdk-1.8.0.352.b08-2.el8_6.x86_64

export JAVA_HOME=/project/haky/users/temi/software/conda_envs/predictdb-env/

phenotype=${1}
s_path=${2}
loci_annotation=${3}
enpact_scores=${4}
geno_dir=/project/haky/data/1000G/population_data/EUR/annot_files

echo "Running ${phenotype}"
# outdir
mkdir -p ${s_path}/${phenotype} && cd ${s_path}/${phenotype}

echo $(which nextflow)

export NXF_DISABLE_CHECK_LATEST=true

# run nextflow
nextflow run /project/haky/users/temi/projects/TFXcan-snakemake/workflow/PredictDb-nextflow/main.nf \
    --gene_annotation "${loci_annotation}" \
    --snp_annotation "${geno_dir}/EUR.snp_annot.txt" \
    --genotype "${geno_dir}/EUR.geno.txt" \
    --gene_exp "${enpact_scores}" \
    --prefix ${phenotype} \
    --outdir ${s_path}/${phenotype}/models/ \
    --keepIntermediate \
    -resume \
    -profile midway3 \
    --nested_cv

#for TISSUE in adipose_0yr adipose_2yr adipose_avg liver_0yr liver_2yr liver_avg muscle_0yr muscle_2yr muscle_avg; do sbatch --job-name=${TISSUE} --export=TISSUE=${TISSUE} gen_models.sbatch; done

# SBATCH --error=/project/haky/users/temi/projects/TFXcan-snakemake/output/logs/lEnpact_training.err
# SBATCH --output=/project/haky/users/temi/projects/TFXcan-snakemake/output/logs/lEnpact_training.log


# sbatch workflow/helpers/generate_lEnpact_models.EUR.midway3.sbatch prostate_cancer_risk /project/haky/users/temi/projects/TFXcan-snakemake/output/lEnpact /project/haky/users/temi/projects/TFXcan-snakemake/data/PCR_GWAS_2024-06-04/predictdb/prostate_cancer_risk.tf_tissue_annot.txt /project/haky/users/temi/projects/TFXcan-snakemake/data/PCR_GWAS_2024-06-04/predictdb/prostate_cancer_risk.enpact_scores.txt

# sbatch workflow/helpers/generate_lEnpact_models.EUR.midway3.sbatch prostate_cancer_risk /project/haky/users/temi/projects/TFXcan-snakemake/output/lEnpact /project/haky/users/temi/projects/TFXcan-snakemake/data/PCR_GWAS_2024-06-04/predictdb/prostate_cancer_risk.tf_tissue_annot.txt /project/haky/users/temi/projects/TFXcan-snakemake/data/PCR_GWAS_2024-06-04/predictdb/prostate_cancer_risk.noChr6.enpact_scores.txt


# sbatch /project/haky/users/temi/projects/TFXcan-snakemake/workflow/helpers/generate_lEnpact_models.EUR.midway3.sbatch prostate_cancer_risk /project/haky/users/temi/projects/TFXcan/lEnpact /project/haky/users/temi/projects/TFXcan-snakemake/data/PCR_GWAS_2024-06-04/predictdb/prostate_cancer_risk.FT.tf_tissue_annot.txt /project/haky/users/temi/projects/TFXcan-snakemake/data/PCR_GWAS_2024-06-04/predictdb/prostate_cancer_risk.FT.enpact_scores.txt

# /project/haky/users/temi/projects/TFXcan-snakemake/data/archives/predictdb/prostate_cancer_risk

# sbatch /project/haky/users/temi/projects/TFXcan-snakemake/workflow/helpers/generate_lEnpact_models.EUR.midway3.sbatch prostate_cancer_risk /project/haky/users/temi/projects/TFXcan/lEnpact /project/haky/users/temi/projects/TFXcan-snakemake/data/archives/predictdb/prostate_cancer_risk/prostate_cancer_risk.tf_tissue_annot.txt /project/haky/users/temi/projects/TFXcan-snakemake/data/archives/predictdb/prostate_cancer_risk/prostate_cancer_risk.enpact_scores.txt

# sbatch /project/haky/users/temi/projects/TFXcan-snakemake/workflow/helpers/generate_lEnpact_models.EUR.midway3.sbatch PcR /project/haky/users/temi/projects/TFXcan/lEnpact /project/haky/users/temi/projects/TFXcan-snakemake/data/PCR_GWAS_2024-06-04/predictdb/prostate_cancer_risk.FT.tf_tissue_annot.txt /project/haky/users/temi/projects/TFXcan-snakemake/data/PCR_GWAS_2024-06-04/predictdb/prostate_cancer_risk.FT.enpact_scores.txt



# sbatch /project/haky/users/temi/projects/TFXcan-snakemake/workflow/helpers/generate_lEnpact_models.EUR.midway3.sbatch \
#     PcR \
#     /project/haky/users/temi/projects/TFXcan/lEnpact \
#     /project/haky/users/temi/projects/TFXcan-snakemake/data/PCR_GWAS_2024-06-04/predictdb/prostate_cancer_risk.FT.tf_tissue_annot.txt \
#     /project/haky/users/temi/projects/TFXcan-snakemake/data/PCR_GWAS_2024-06-04/predictdb/prostate_cancer_risk.FT.enpact_scores.txt