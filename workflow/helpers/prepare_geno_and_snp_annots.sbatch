#!/bin/bash
# Author: Temi
# Description: Convert vcfs to dosages.txt file
# Usage: sbatch calculate_variant_correlations.sbatch ...
# Date: Fri Dec 15 2023
# Dependencies: 
# parameters: phenotype, 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=12G
#SBATCH --job-name=QC_filter_vcfs
#SBATCH --account=pi-haky
#SBATCH --output=/project2/haky/temi/projects/TFPred/logs/process_vcf.out
#SBATCH --error=/project2/haky/temi/projects/TFPred/logs/process_vcf.err
#SBATCH --time=02:00:00	
#SBATCH --partition=bigmem

slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}/run
mkdir -p ${SLURM_O_WORKDIR}
echo Working directory is $SLURM_O_WORKDIR
cd $SLURM_O_WORKDIR

echo Jobid: $SLURM_JOBID
echo Running on host `hostname`

source ~/.bashrc
conda activate compbio-tools
module load parallel

#ref_panel_directory="/project2/haky/Data/1000G/plink_ref_panel"
vcfs_dir='/project2/haky/Data/1000G/vcf_snps_only'
metadata_dir='/project2/haky/Data/1000G/metadata'
population_data_dir=$"/project2/haky/Data/1000G/population_data"
rscript=/project2/haky/temi/software/miniconda3/envs/r-env/bin/Rscript

## declare an array variable
declare -a popcode=("AFR" "AMR" "EAS" "EUR" "SAS")
declare -a chromosomes=( {1..22} X )
# use plink to filter and make bed files 
for pop in "${popcode[@]}"; do
    echo "INFO - processing ${pop} samples"
    pop_dir=${population_data_dir}/${pop}/bfiles
    out_dir=${population_data_dir}/${pop}/annot_files
    mkdir -p ${pop_dir} ${out_dir}

    if [[ ! -f "${out_dir}/${pop}.geno.txt.gz" && ! -f "${out_dir}${pop}.snp_annot.txt.gz" ]]; then
        echo ${rscript} /project2/haky/temi/projects/TFXcan-snakemake/workflow/helpers/prepare.R --pop ${pop} --snp_annot_pattern ${vcfs_dir}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz --geno_pattern ${vcfs_dir}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz --output_folder ${out_dir}
    fi
done