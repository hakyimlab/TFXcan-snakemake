# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--summary_stats_file", help='A GWAS summary statistics file; should be a tsv file with columns: chrom, pos, ref, alt, pval, beta, se, zscore'),
    make_option("--output_folder", help='the output folder'),
    make_option("--phenotype", help = 'a GWAS phenotype'),
    make_option("--pvalue_threshold", default=5e-8, type='numeric', help = 'the pvalue threshold for significance; default is 5e-8'),
    make_option('--diagnostics_file', type='character', default=NULL, help='A file to write diagnostics to; default is NULL i.e no diagnostics file will be written')
)

#source('/project/haky/users/temi/projects/TFXcan-snakemake/workflow/src/modules.R')

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glue)

print(opt)

# opt <- list()
# opt$summary_stats_file <- '/project2/haky/temi/projects/TFXcan-snakemake/data/sumstats/asthma_children.liftover.logistic.assoc.tsv.gz'
# opt$pvalue_threshold <- 5e-8
# opt$output_folder <- '/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/testosterone'
# opt$phenotype <- 'asthma_children'
# opt$diagnostics_folder <- '/project2/haky/temi/projects/TFXcan-snakemake/data/diagnostics'

chrom_filter <- c(1:22)
if(!dir.exists(opt$output_folder)){dir.create(opt$output_folder)}

dt <- data.table::fread(opt$summary_stats_file) 

# important checks
if(!all(c('chrom', 'pos', 'ref', 'alt', 'pval', 'beta', 'se') %in% colnames(dt))){
    stop('ERROR - The summary statistics file must have columns: chrom, pos, ref, alt, pval, beta, se')
}

dt <- dt %>%
    dplyr::mutate(chrom = dplyr::case_when(
        grepl('chr', chrom, fixed=TRUE) ~ gsub('chr', '', chrom),
        .default = as.character(chrom)
    )) %>%
    dplyr::filter(chrom %in% chrom_filter) %>%
    dplyr::mutate(SNP = paste(chrom, pos, ref, alt, sep='_'), adjP = pval) %>%
    dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1) %>% # still need to properly resolve this
    dplyr::mutate(gwas_significant=ifelse(pval <= opt$pvalue_threshold, 'YES', 'NO')) 
    
if(!'rsid' %in% colnames(dt)){
    dt <- dt %>%
        dplyr::mutate(rsid = paste(chrom, pos, ref, alt, sep=':'))
}

if(!'zscore' %in% colnames(dt)){
    dt <- dt %>%
        dplyr::mutate(zscore = beta/se)
}

# split and save 

split_dt <- dt %>%
    base::split(f=.$chrom)

nn <- names(split_dt)
lapply(seq_along(split_dt), function(i){
    data.table::fwrite(split_dt[[i]], file=glue('{opt$output_folder}/chr{nn[i]}.sumstats.txt.gz'), 
        quote=F, col.names=T, row.names=F, sep='\t', compress='gzip')
})