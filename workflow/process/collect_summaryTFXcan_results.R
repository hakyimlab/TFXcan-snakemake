# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library(optparse))
#library(optparse)
option_list <- list(
    make_option("--input_files_pattern", help='a pattern for the predicted enpact scores'), 
    make_option("--phenotype", default = NULL, help='the phenotype name'),
    make_option("--output_file")
)

opt <- optparse::parse_args(OptionParser(option_list=option_list))
print(opt)

library(data.table)
library(tidyverse)
library(glue)

# opt <- list()
# opt$input_files_pattern <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/T2D_2025-01-24/summaryTFXcan/t2d_suzuki/.*-t2d_suzuki.enpactScores.spredixcan.csv'
# opt$phenotype <- 't2d_suzuki'

# grep
sFiles <- list.files(path = dirname(opt$input_files_pattern), pattern = basename(opt$input_files_pattern), full.names = T)

if(length(sFiles) == 0){
    message('ERROR - No summary TFXcan files found.')
    quit(-1)
}

if(is.null(opt$phenotype)){
    message(glue("INFO - Found {length(sFiles)} summary TFXcan files to collect."))
} else {
    message(glue("INFO - Found {length(sFiles)} summary TFXcan files for {opt$phenotype} to collect."))
}

# split the input files
#sFiles <- base::strsplit(opt$input_files, split = ',')[[1]]
dFiles <- purrr::map(.x=sFiles, .f=function(each_file){
    if(file.exists(each_file)){
        dt <- data.table::fread(each_file)
        return(dt)
    } else {
        return(NULL)
    }
})
# remove the NULLs
dFiles <- Filter(Negate(is.null), dFiles)
# bind the rows
dt <- dplyr::bind_rows(dFiles)
# rename the columns
dt <- dt %>% dplyr::select(-c(gene_name)) %>% dplyr::rename(tfbs = gene)
# write the file
data.table::fwrite(dt, opt$output_file, sep = "\t", quote = F, row.names = F)