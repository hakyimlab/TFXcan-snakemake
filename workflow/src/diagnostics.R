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

opt <- list()

opt$sumstats_file_pattern <- '/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/asthma_children/chr{}.sumstats.txt.gz'