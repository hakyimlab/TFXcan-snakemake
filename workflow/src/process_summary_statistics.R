# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--summary_stats_file", help='A transcription factor e.g. AR'),
    make_option("--output_folder", help='the output folder'),
    make_option("--phenotype", help = 'a GWAS phenotype')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glue)

print(opt)

# opt <- list()
# opt$summary_stats_file <- '/project2/haky/temi/projects/TFXcan-snakemake/data/sumstats/testosterone.gwas_sumstats.ALL.filtered.txt.gz'
# opt$output_folder <- '/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/testosterone'
# opt$phenotype <- 'testosterone'

chrom_filter <- c(1:22, 'X', 'Y')
if(!dir.exists(opt$output_folder)){dir.create(opt$output_folder)}

dt <- data.table::fread(opt$summary_stats_file) %>%
    dplyr::mutate(chrom = dplyr::case_when(
        grepl('chr', chrom, fixed=TRUE) ~ gsub('chr', '', chrom),
        .default = as.character(chrom)
    )) %>%
    dplyr::filter(chrom %in% chrom_filter) %>%
    dplyr::mutate(SNP = paste(chrom, pos, alt, ref, sep='_')) %>%
    dplyr::group_by(chrom)
print(head(dt))

split_dt <- dt %>%
    base::split(f=.$chrom)

nn <- names(split_dt)
lapply(seq_along(split_dt), function(i){
    data.table::fwrite(split_dt[[i]], file=glue('{opt$output_folder}/chr{nn[i]}.sumstats.txt.gz'), 
        quote=F, col.names=T, row.names=F, sep='\t', compress='gzip')
})
    
# dt %>%
#     dplyr::pull(SNP) %>% 
#     as.data.frame() %>%
#     data.table::fwrite(., file=glue('{opt$output_file_basename}.SNPsForLD'), quote=F, col.names=F, row.names=F)

# sumstats <- data.table::fread(dt) %>%
#     dplyr::mutate(chrom=as.character(gsub('chr', '', chrom))) %>%
#     dplyr::filter(chrom == opt$chromosome)

# # QC filter
# data.table::fwrite(as.data.frame(sumstats$SNP), file=glue('{opt$snp_list}'), quote=F, row.names=F, col.names=F)
