#!/bin/bash
# Author: Temi
# Description: Convert vcfs to dosages.txt file
# Usage: sbatch calculate_variant_correlations.sbatch ...
# Date: Fri Dec 15 2023
# Dependencies: 
# parameters: phenotype, 

#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=calculate_variant_correlations
#SBATCH --account=pi-haky
#SBATCH --output=/project2/haky/temi/projects/TFPred/logs/calculate_variant_correlations.out
#SBATCH --error=/project2/haky/temi/projects/TFPred/logs/calculate_variant_correlations.err
#SBATCH --time=04:00:00	
#SBATCH --partition=caslake

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
module load mpich/

# parameter definitions
population_data_dir="/project2/haky/Data/1000G/population_data"
LD_variants_dir="/project2/haky/Data/1000G/LD/LD_variants"
LD_matrices_dir="/project2/haky/Data/1000G/LD/LD_matrices"
LD_dir='/project2/haky/Data/1000G/LD'
#output_directory=${ref_panel_directory}/LD
mkdir -p ${LD_variants_dir} ${LD_matrices_dir} ${LD_dir}

# first calculate the number of nodes available 
# obtain the list of SNPs in that LD window
# --> I can do this with a R script
# printf "%s\n" {1..22} X | parallel -j 23 "/project2/haky/temi/software/miniconda3/envs/r-env/bin/Rscript /project2/haky/temi/projects/TFXcan-snakemake/workflow/helpers/obtain_snps_in_block.R --chromosome {} --partition_file /project2/haky/Data/LD_blocks/hg38/EUR/hg38_fourier_ls-all.bed --plink_bim_file /project2/haky/Data/1000G/plink_ref_panel/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.bim --LD_folder ${LD_DIR}"

ls ${LD_dir}/chr*.SNPsLDBlock.txt | xargs -n 1 basename > blocks.txt
blocks_list=blocks.txt
x=`wc -l ${blocks_list} | awk '{print $1}'` # count the number of lines (samples or individuals)
y=${SLURM_NNODES} # >> should be 4 # count the number of nodes allocated
ll=$(( ($x + $y - 1) / $y )) # do the math such that you have maximum ll equally split

split --lines="${ll}" --numeric-suffixes=1 --suffix-length=1 "${blocks_list}" "block."

# prepare nodes
nodelist=$(scontrol show hostname $SLURM_NODELIST)
printf "%s\n" "${nodelist[@]}" > local_hostfiles.txt
split --lines=1 --numeric-suffixes=1 --suffix-length=1 local_hostfiles.txt "local_hostfile."

function gnuParallelComputeLDByBlock () {
  blockNode=${1} # chr12.SNPsLDBlock.txt, chr18.SNPsLDBlock.txt, chr1.SNPsLDBlock.txt
  blockFolder=${2} # 
  outputFolder=${3}
  populationDataFolder=${4}

  declare -a popcode=("AFR" "AMR" "EAS" "EUR" "SAS")

  while read blockFile; do # each p is chr12.SNPsLDBlock.txt
    chromosome=$(echo ${blockFile} | cut -d "." -f1)
      
    for population in "${popcode[@]}"; do

      echo "INFO - processing ${chromosome} and ${population}\n"

      LDMatrixFolder=${outputFolder}/${population}/${chromosome}
      mkdir -p ${LDMatrixFolder}

      parallel -j 10 "plink --r --ld-window-r2 0 --bfile ${populationDataFolder}/${population}/bfiles/ALL.${chromosome}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased --out ${LDMatrixFolder}/{}.LDMatrix --ld-snp-list ${blockFolder}/LD_variants/${chromosome}/{}.SNPsLDBlock.txt" < ${blockFolder}/${blockFile}
    done

  done < ${blockNode}
}

export -f gnuParallelComputeLDByBlock

for suf in `seq 1 ${SLURM_NNODES}`; do
  (
    mpirun -np 1 --hostfile "local_hostfile.${suf}" bash -c "gnuParallelComputeLDByBlock block.${suf} ${LD_dir} ${LD_matrices_dir}/ ${population_data_dir}"
  ) & sleep 1
  done 
wait