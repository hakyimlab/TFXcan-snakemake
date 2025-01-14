# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--summary_stats_pattern", help='A GWAS summary statistics file; should be a tsv file with columns: chrom, pos, ref, alt, pval, beta, se, zscore'),
    make_option("--annotation_file", help='the output folder'),
    make_option("--chromosome", help = 'a GWAS phenotype'),
    make_option("--output_file", default=5e-8, type='numeric', help = 'the pvalue threshold for significance; default is 5e-8')
)

#source('/project/haky/users/temi/projects/TFXcan-snakemake/workflow/src/modules.R')

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glue)
library(bigsnpr)

print(opt)

opt <- list()
opt$summary_stats_pattern <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/standing_height_2025-01-10/processed_sumstats/standing_height/chr{}.sumstats.txt.gz'
opt$annotation_file <- '/project2/haky/Data/1000G/population_data/EUR/annot_files/EUR.snp_annot.txt'
opt$chromosome <- '1'
opt$output_file <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/standing_height_2025-01-10/processed_sumstats/standing_height/chr{}.matched.sumstats.txt.gz'

# read in the summary statistics file
ss_file <- gsub('\\{\\}', opt$chromosome, opt$summary_stats_pattern)
aa_file <- opt$annotation_file
if(!(file.exists(ss_file) & file.exists(aa_file))){
    stop('ERROR - Either or both of the summary statistics and/or the reference annotation files does not exist')
} else {
    sumstats <- data.table::fread(ss_file) %>% 
        data.table::setDT() %>% 
        dplyr::rename(chr=chrom, a0=ref, a1=alt)
    annotation <- data.table::fread(aa_file) %>% 
        data.table::setDT() %>%
        dplyr::filter(chr == opt$chromosome) %>%
        dplyr::rename(a0=ref_vcf, a1=alt_vcf)
}

# match snps
matched_stats <- bigsnpr::snp_match(sumstats, annotation, return_flip_and_rev = TRUE) |> data.table::setDT() %>% 
    dplyr::select(chrom=chr, pos, ref=a0, alt=a1, rsid, varID, beta, se, zscore, pval, adjP, gwas_significant)

# write out
output_file <- gsub('\\{\\}', opt$chromosome, opt$output_file)
data.table::fwrite(matched_stats, output_file, sep='\t', quote=F, row.names=F, col.names=T, compress = 'gzip')

