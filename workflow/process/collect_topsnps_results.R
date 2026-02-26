# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript collect_topsnps_results.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--selection_dir", help='[Input] A list of files to combine'),
    make_option("--filtered_sumstats", help='[Output] Processed summary statistics file'),
    make_option("--enformer_loci", help='[Input] Processed list of loci for Enformer inference of epigenomic features'),
    make_option("--phenotype", help='[Input] Name of the phenotype'),
    make_option("--limit_number_of_loci", default = NULL, help='[Input] How many top loci should TFXcan be run at? Default is everything, but > 200 may be computationally expensive.', type='integer')
)

opt <- parse_args(OptionParser(option_list=option_list))  
print(opt)

library(data.table)
library(tidyverse)
library(glue)

# opt <- list()
# opt$selection_dir <- '/project/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2024-09-23/filtering/prostate_cancer_risk'
# opt$phenotype <- 'prostate_cancer_risk'
# opt$filtered_sumstats <- '/project/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2024-09-23/collection/prostate_cancer_risk.filteredGWAS.topSNPs.txt.gz'
# opt$enformer_loci <- '/project/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2024-09-23/collection/prostate_cancer_risk.EnformerLoci.topSNPs.txt'

if(!dir.exists(dirname(opt$filtered_sumstats))){
    dir.create(dirname(opt$filtered_sumstats), recursive = T)
}

chrom_filter <- c(1:22)

dir.create(dirname(opt$filtered_sumstats), recursive = T, showWarnings = F)
dir.create(dirname(opt$enformer_loci), recursive = T, showWarnings = F)

finemapped_pattern <- file.path(opt$selection_dir, glue::glue('{opt$phenotype}.chr{chrom_filter}.filteredGWAS.topSNPs.txt'))

fg <- lapply(finemapped_pattern, function(each_file){
    if(file.exists(each_file)){
        dt <- data.table::fread(each_file)
        return(dt)
    }
}) %>%
    Filter(Negate(is.null), .) %>%
    do.call('rbind', .) %>%
    as.data.frame()

data.table::fwrite(fg, file=opt$filtered_sumstats, compress='gzip', quote=F, row.names=F, sep = '\t')

if(!is.null(opt$limit_number_of_loci)){
    subfg <- fg %>%
        dplyr::arrange(pval) %>%
        dplyr::slice(1:opt$limit_number_of_loci)
} else {
    subfg <- as.data.frame(fg)
}


finemapped_pattern <- file.path(opt$selection_dir, glue::glue('{opt$phenotype}.chr{chrom_filter}.EnformerLoci.topSNPs.txt'))

lg <- lapply(finemapped_pattern, function(each_file){
    if(file.exists(each_file)){
        dt <- data.table::fread(each_file, header=F)
        return(dt)
    }
}) %>%
    Filter(Negate(is.null), .) %>%
    do.call('rbind', .) %>%
    as.data.frame() 

loci_subfg <- paste0('chr', subfg$chr, '_', subfg$pos, '_', subfg$pos + 1)

lg %>% dplyr::filter(V1 %in% loci_subfg) %>%
    data.table::fwrite(., file=opt$enformer_loci, quote=F, row.names=F, sep = '\t', col.names=F)