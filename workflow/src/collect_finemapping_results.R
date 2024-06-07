# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--finemapping_dir", help='A list of files to combine'),
    make_option("--filtered_sumstats", help='reference panel correlation matrix'),
    make_option("--enformer_loci", help='reference panel correlation matrix'),
    make_option("--phenotype", help='summary statistics sample size')
)

opt <- parse_args(OptionParser(option_list=option_list))  
print(opt)

library(data.table)
library(tidyverse)
library(susieR)
library(glue)

# opt <- list()
# opt$finemapping_dir <- '/project2/haky/temi/projects/TFXcan-snakemake/data/finemapping/asthma_children'
# opt$phenotype <- 'asthma_children'
# opt$filtered_sumstats <- '/project2/haky/temi/projects/TFXcan-snakemake/data/collection/asthma_children/asthma_children.filteredGWAS.txt.gz'
# opt$enformer_loci <- '/project2/haky/temi/projects/TFXcan-snakemake/data/collection/asthma_children/asthma_children.EnformerLoci.txt'

if(!dir.exists(dirname(opt$filtered_sumstats))){
    dir.create(dirname(opt$filtered_sumstats), recursive = T)
}

chrom_filter <- c(1:22)

finemapped_pattern <- file.path(opt$finemapping_dir, glue::glue('{opt$phenotype}.chr{chrom_filter}.filteredGWAS.txt.gz'))

lapply(finemapped_pattern, function(each_file){
    if(file.exists(each_file)){
        dt <- data.table::fread(each_file)
        return(dt)
    }
}) %>%
    Filter(Negate(is.null), .) %>%
    do.call('rbind', .) %>%
    as.data.frame() %>%
    data.table::fwrite(., file=opt$filtered_sumstats, compress='gzip', quote=F, row.names=F, sep = '\t')


finemapped_pattern <- file.path(opt$finemapping_dir, glue::glue('{opt$phenotype}.chr{chrom_filter}.EnformerLoci.txt'))

lapply(finemapped_pattern, function(each_file){
    if(file.exists(each_file)){
        dt <- data.table::fread(each_file, header=F)
        return(dt)
    }
}) %>%
    Filter(Negate(is.null), .) %>%
    do.call('rbind', .) %>%
    as.data.frame() %>%
    data.table::fwrite(., file=opt$enformer_loci, quote=F, row.names=F, sep = '\t', col.names=F)