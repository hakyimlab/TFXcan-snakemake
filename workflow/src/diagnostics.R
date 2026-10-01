# Author: Temi
# Date: Thursday August 10 2023
# Description: INCOMPLETE/SCRATCH -- not called by any rule. As checked in, it parses CLI
#   options below but then immediately discards them (opt gets overwritten by a hardcoded
#   list), sets one hardcoded sumstats file-pattern variable, and does nothing further.
#   Looks like the start of a diagnostics script (e.g. for inspecting processed GWAS
#   summary statistics per chromosome) that was never finished.
# Usage: Rscript diagnostics.R [options] (currently a no-op beyond the two lines below)

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--summary_stats_file", help='A GWAS summary statistics file'),
    make_option("--output_folder", help='the output folder'),
    make_option("--phenotype", help = 'a GWAS phenotype')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glue)

print(opt)

# hardcoded scratch values below override the CLI opt above -- see header note
opt <- list()

opt$sumstats_file_pattern <- '/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/asthma_children/chr{}.sumstats.txt.gz'