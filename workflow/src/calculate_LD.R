# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--output_pattern", help='A list of files to combine'),
    make_option("--bfile_pattern", help='reference panel correlation matrix'),
    make_option("--snplist_pattern", help='A list of files to combine')
)

opt <- parse_args(OptionParser(option_list=option_list))  

library(data.table)
library(tidyverse)
library(susieR)
library(glue)

# opt <- list()
# opt$output_pattern <- '/project2/haky/temi/projects/TFXcan-snakemake/data/LD/asthma_adults/asthma_adults.chr{}'
# opt$snplist_pattern <- '/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/asthma_adults/asthma_adults.chr{}.SNPsForLD'
# opt$bfile_pattern <- '/project2/haky/Data/1000G/plink_ref_panel/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased'

chrom_filter <- paste(c(1:22, 'X', 'Y'), sep='')

ld_files <- sapply(chrom_filter, function(f){
    ff <- gsub('\\{\\}', f, opt$snplist_pattern)
    if(file.exists(ff)){
        return(ff)
    }
}) 

ld_files <- Filter(Negate(is.null), ld_files)
valid_chrom <- names(ld_files)


cmd <- print(
    glue(
        "printf '%s\\n' {paste(valid_chrom, collapse=' ')} | parallel \"plink --r --ld-window-r2 0 --bfile {opt$bfile_pattern} --out {opt$output_pattern} --ld-snp-list {opt$snplist_pattern}\""
    )
)

system(cmd)



# sapply(valid_chrom, function(vc){
#     nn <- names(ld_file)
#     print(nn)
#     # if(file.exists(ld_file)){
#     #     glue('module load parallel; \n
#     #         plink printf "%s\n" {1..22} X | parallel -j 23 "plink --r --ld-window-r2 0 --bfile ${ref_panel_directory}/plink_ref_panel/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased --out ${output_directory}/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased"')
#     # }
# })




# snplist <- data.table::fread(opt$snp_list, header=F) %>%
#     dplyr::mutate(chrom=as.character(gsub('chr', '', chrom))) %>%
#     dplyr::filter(chrom == opt$chromosome)

# # QC filter
# data.table::fwrite(as.data.frame(sumstats$SNP), file=glue('{opt$snp_list}'), quote=F, row.names=F, col.names=F)

#  {params.rscript} workflow/src/create_correlation_matrix_from_reference_genome.R --output_file {output} --reference_genome_txt_file {params.reference_genome_txt_file}