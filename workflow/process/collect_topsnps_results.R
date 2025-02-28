# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--selection_dir", help='A list of files to combine'),
    make_option("--filtered_sumstats", help='reference panel correlation matrix'),
    make_option("--enformer_loci", help='reference panel correlation matrix'),
    make_option("--phenotype", help='summary statistics sample size')
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

if(nrow(fg) > 100){
    subfg <- fg %>%
        dplyr::arrange(pval) %>%
        dplyr::slice(1:100)
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