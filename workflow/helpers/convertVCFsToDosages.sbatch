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


## CHANGE THESE ======
vcfs_pattern="chr{}_CWAS_SNPs_2023-05-15_phased.vcf.gz" # "merged_cwas_genotypes_chrRenamed.vcf.gz" # must end with vcf.gz

vcfs_pattern_no_basename=`basename ${vcfs_pattern} .vcf.gz` # ../vcfs/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased
pop="BACA"

#ref_panel_directory="/project2/haky/Data/1000G/plink_ref_panel"
base_dir='/project2/haky/Data/baca_cwas/vcfs/hg19'


vcfs_dir='/project2/haky/Data/baca_cwas/vcfs/hg19/split_vcf'

#metadata_dir='/project2/haky/Data/1000G/metadata'
rscript=/project2/haky/temi/software/miniconda3/envs/r-env/bin/Rscript
plink_dir="${base_dir}/plink"
dosages_dir="${base_dir}/dosages"
mkdir -p ${plink_dir}
mkdir -p ${dosages_dir}


# first check that the plink files exist or not
famExistingFiles=$(printf "%s\n" {1..22} | xargs -I {} bash -c 'if [ ! -f '${plink_dir}/${vcfs_pattern_no_basename}.fam' ]; then echo {} ; fi')
bedExistingFiles=$(printf "%s\n" {1..22} | xargs -I {} bash -c 'if [ ! -f '${plink_dir}/${vcfs_pattern_no_basename}.bed' ]; then echo {} ; fi')
bimExistingFiles=$(printf "%s\n" {1..22} | xargs -I {} bash -c 'if [ ! -f '${plink_dir}/${vcfs_pattern_no_basename}.bim' ]; then echo {} ; fi')
trawExistingFiles=$(printf "%s\n" {1..22} | xargs -I {} bash -c 'if [ ! -f '${plink_dir}/${vcfs_pattern_no_basename}.traw' ]; then echo {} ; fi')

if [[ -n ${famExistingFiles} && -n ${bedExistingFiles} && -n ${bimExistingFiles} && -n ${trawExistingFiles} ]]; then 
    echo "INFO - Creating plink files"
    printf "%s\n" {1..22} X | parallel -j 23 "plink --make-bed --snps-only --vcf ${vcfs_dir}/${vcfs_pattern} --maf 0.01 --geno 0.01 --hwe 1e-5 --set-missing-var-ids '@_#_\$1_\$2' --out ${plink_dir}/${vcfs_pattern_no_basename} --recode A-transpose --keep-allele-order --double-id"
else 
    echo "INFO - All necessary bfiles exist" #&& exit 
fi

# --keep ${metadata_dir}/${pop}.samples.txt

echo "INFO - Creating dosage files"
printf "%s\n" {1..22} X | parallel -j 10 "${rscript} /project2/haky/temi/projects/TFXcan-snakemake/workflow/helpers/encode_alleles.R --bfile ${plink_dir}/${vcfs_pattern_no_basename} --output_prefix ${dosages_dir}/${vcfs_pattern_no_basename}"

## Merge into one giant dosage file
if [[ ! -f "${dosages_dir}/${pop}.geno.txt.gz" && ! -f "${dosages_dir}/${pop}.snp_annot.txt.gz" ]]; then
    ${rscript} /project2/haky/temi/projects/TFXcan-snakemake/workflow/helpers/prepare.R --pop ${pop} --snp_annot_pattern ${dosages_dir}/${vcfs_pattern_no_basename}.snp_annot.txt.gz --geno_pattern ${dosages_dir}/${vcfs_pattern_no_basename}.geno.txt.gz --output_folder ${dosages_dir}
fi



# ## declare an array variable
# declare -a popcode=("AFR" "AMR" "EAS" "EUR" "SAS")
# declare -a chromosomes=( {1..22} X )






# # use plink to filter and make bed files 
# for pop in "${popcode[@]}"; do
#     echo "INFO - processing ${pop} samples"
#     pop_dir=${population_data_dir}/${pop}/bfiles
#     mkdir -p ${pop_dir}

#     #${pop_dir}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased
#     # echo $(printf "%s\n" {1..22} X | xargs -I {} bash -c 'echo {}')
#     # test if bfiles files exists
#     famExistingFiles=$(printf "%s\n" {1..22} X | xargs -I {} bash -c 'if [ ! -f '${pop_dir}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.fam' ]; then echo {} ; fi')
#     bedExistingFiles=$(printf "%s\n" {1..22} X | xargs -I {} bash -c 'if [ ! -f '${pop_dir}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.bed' ]; then echo {} ; fi')
#     bimExistingFiles=$(printf "%s\n" {1..22} X | xargs -I {} bash -c 'if [ ! -f '${pop_dir}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.bim' ]; then echo {} ; fi')
#     trawExistingFiles=$(printf "%s\n" {1..22} X | xargs -I {} bash -c 'if [ ! -f '${pop_dir}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.traw' ]; then echo {} ; fi')

    
# done

# echo "INFO - Creating dosage files for all population-chromosome combinations"
# parallel -j 10 "${rscript} /project2/haky/temi/projects/TFXcan-snakemake/workflow/helpers/encode_alleles.R --bfile ${population_data_dir}/{1}/bfiles/ALL.chr{2}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased --output_prefix ${population_data_dir}/{1}/bfiles/ALL.chr{2}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased" ::: "${popcode[@]}" ::: "${chromosomes[@]}"


# for pop in "${popcode[@]}"; do
#     echo "INFO - creating megafile for ${pop} samples"
#     pop_dir=${population_data_dir}/${pop}/bfiles
#     out_dir=${population_data_dir}/${pop}/annot_files
#     mkdir -p ${pop_dir} ${out_dir}

#     if [[ ! -f "${out_dir}/${pop}.geno.txt.gz" && ! -f "${out_dir}${pop}.snp_annot.txt.gz" ]]; then
#         ${rscript} /project2/haky/temi/projects/TFXcan-snakemake/workflow/helpers/prepare.R --pop ${pop} --snp_annot_pattern ${pop_dir}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.snp_annot.txt.gz --geno_pattern ${pop_dir}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.geno.txt.gz --output_folder ${out_dir}
#     fi
# done

# printf "%s\n" {1..22} X Y | parallel -j 23 "ls -1 ./ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.fam" 

# tt=$(printf "%s\n" {1..22} X Y | xargs -I {} bash -c 'if [ ! -f ./ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.fam ]; then echo {} ; fi')

# parallel -j 23 "ls -1 ./ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.bim"

# ls ALL.chrY.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.fam && echo $?
# if [ test -f "./ALL.chrX.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.fam" ]; then echo "File exists" else echo "File does not exist" fi

# printf "%s\n" {1..22} X Y | parallel -j 23 "find . ! -name . -prune -name 'ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.bim' -type f | read" # | read



# find . ! -name . -prune -name '*.fam' -type f | read

# for pop in "${popcode[@]}"; do
#     for chromosome in {1..22} X; do

#     done
#     echo "INFO - processing ${pop} samples"
#     pop_dir=${population_data_dir}/${pop}/bfiles
#     path_pattern="${pop_dir}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased"

#     if [ ! -f ${path_pattern}.geno.txt.gz ] && [ ! -f ${path_pattern}.snp_annot.txt.gz ]; then
#         printf '\n%s\n' "INFO - preparing ${pattern}: R portion"
#         ${rscript} ${utility_functions} --file_prefix ${plink_results_pattern} --output_prefix ${formatted_files_folder}/${pattern} --command "create_genotype_dosage_and_snp_annot_files"
#     fi
# done