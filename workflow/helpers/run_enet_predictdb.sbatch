#!/bin/bash

#SBATCH --job-name=lEnpact_training_unnested
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --error=/project2/haky/temi/projects/l-TFPred/logs/lEnpact_training_unnested.err
#SBATCH --output=/project2/haky/temi/projects/l-TFPred/logs/lEnpact_training_unnested.out
#SBATCH --account=pi-haky

source ~/.bashrc
conda activate /project2/haky/festus/miniconda3/envs/predictdb_nf
module load python
module load R
Rscript -e 'install.packages("/project2/haky/temi/projects/TFXcan-snakemake/workflow/software/optparse_1.7.4.tar.gz", lib="/home/temi/R/x86_64-pc-linux-gnu-library/4.2", repos = NULL, type = "source")'

#install.packages(pkgname.tar.gz, repos = NULL, type =‘‘source’’)
#R CMD INSTALL /project2/haky/temi/projects/TFXcan-snakemake/workflow/helpers/optparse_1.7.4.tar.gz
#export JAVA_HOME=/usr/lib/jvm/jre-1.8.0-openjdk-1.8.0.352.b08-2.el8_6.x86_64

export JAVA_HOME=/project2/haky/festus/miniconda3/envs/predictdb_nf/

phenotype="asthma_children" 
s_path="/project2/haky/temi/projects/l-TFPred/models" 
loci_annotation="/project2/haky/temi/projects/TFXcan-snakemake/data/predictdb/asthma_children/asthma_children.tf_tissue_annot.txt"
enpact_scores="/project2/haky/temi/projects/TFXcan-snakemake/data/predictdb/asthma_children/asthma_children.enpact_scores.txt"
geno_dir=/project2/haky/Data/1000G/population_data/EUR/annot_files
export PATH=$PATH:/project2/haky/festus/github/peer/build/peertool/bin

echo "Running ${phenotype}"

# outdir
mkdir -p ${s_path}/${phenotype} && cd ${s_path}/${phenotype}

# run nextflow
nextflow run /project2/haky/festus/github/PredictDb-nextflow/main.nf \
    --gene_annotation ${loci_annotation} \
    --snp_annotation ${geno_dir}/EUR.snp_annot.txt \
    --genotype ${geno_dir}/EUR.geno.txt \
    --gene_exp ${enpact_scores} \
    --prefix ${phenotype} \
    --outdir ${s_path}/${phenotype}/models/ \
    --keepIntermediate \
    -resume \
    -profile midway3

    


#for TISSUE in adipose_0yr adipose_2yr adipose_avg liver_0yr liver_2yr liver_avg muscle_0yr muscle_2yr muscle_avg; do sbatch --job-name=${TISSUE} --export=TISSUE=${TISSUE} gen_models.sbatch; done

