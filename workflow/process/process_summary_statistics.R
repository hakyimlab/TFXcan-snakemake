# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--summary_stats_file", help='A GWAS summary statistics file; should be a tsv file with columns: chrom, pos, ref, alt, pval, beta, se, zscore'),
    make_option("--output_folder", help='the output folder'),
    make_option("--annotation_file", help = 'a GWAS phenotype'),
    make_option("--pvalue_threshold", default=5e-8, type='numeric', help = 'the pvalue threshold for significance; default is 5e-8'),
    make_option('--diagnostics_file', type='character', default=NULL, help='A file to write diagnostics to; default is NULL i.e no diagnostics file will be written')
)

#source('/project/haky/users/temi/projects/TFXcan-snakemake/workflow/src/modules.R')

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glue)
library(bigsnpr)

print(opt)


# setwd('/beagle3/haky/users/temi/projects/TFXcan-snakemake')
# opt <- list()
# opt$summary_stats_file <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/sumstats/standing_height.gwas_sumstats.processed.txt.gz'
# opt$pvalue_threshold <- 5e-8
# opt$output_folder <- '/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/testosterone'
# opt$annotation_file <- '/project2/haky/Data/1000G/population_data/EUR/annot_files/EUR.snp_annot.txt'
# opt$diagnostics_folder <- '/project2/haky/temi/projects/TFXcan-snakemake/data/diagnostics'

chrom_filter <- c(1:22)
if(!dir.exists(opt$output_folder)){dir.create(opt$output_folder)}

if(!file.exists(opt$annotation_file)){
    stop('ERROR - The reference annotation file does not exist.')
}

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

dt <- dt %>% dplyr::mutate(chrom = as.numeric(chrom))

# split and match the snps 
split_dt <- dt %>%
    base::split(f=.$chrom)

## use purrr to match snps
annotation <- data.table::fread(opt$annotation_file) %>% 
    data.table::setDT() %>%
    dplyr::rename(a0=ref_vcf, a1=alt_vcf)

nn <- names(split_dt)

matched_stats <- purrr::map(nn, function(ch){
    x <- split_dt[[ch]]
    x <- x %>% dplyr::rename(chr=chrom, a0=ref, a1=alt)
    annot <- annotation %>% dplyr::filter(chr == ch) #%>% dplyr::mutate(chr = as.numeric(chr))
    # match snps
    matched_stats <- tryCatch({
        mm <- bigsnpr::snp_match(x, annot, return_flip_and_rev = TRUE) |> data.table::setDT() %>% 
            dplyr::select(chrom=chr, pos, ref=a0, alt=a1, rsid, varID, beta, se, zscore, pval, adjP, gwas_significant)
        return(mm)
    }, error = function(e){
        print(glue('ERROR -  Not enough variants have been matched for chromosome {ch}.'))
        return(NULL)
    })
    return(matched_stats)
}, .progress = TRUE)

names(matched_stats) <- nn
matched_stats <- Filter(Negate(is.null), matched_stats) 

lapply(seq_along(matched_stats), function(i){
    data.table::fwrite(split_dt[[i]], file=glue('{opt$output_folder}/chr{nn[i]}.sumstats.txt.gz'), 
        quote=F, col.names=T, row.names=F, sep='\t', compress='gzip')
})