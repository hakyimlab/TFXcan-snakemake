# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library(optparse))
#library(optparse)
option_list <- list(
    make_option("--input_files_pattern", help='[Input] A headerless file listing all the summary TFXcan results to be collected'), 
    make_option("--phenotype", default = NULL, help='[Input] the phenotype name'),
    make_option("--output_file", help='[Output] The output file')
)

opt <- optparse::parse_args(OptionParser(option_list=option_list))
print(opt)

library(data.table)
library(tidyverse)
library(glue)

# opt <- list()
# opt$input_files_pattern <- "/scratch/beagle3/temi/summary_tfxcan/prostate_cancer_risk/*.prostate_cancer_risk.summaryTFXcan.csv"
# opt$phenotype <- "prostate_cancer_risk"

# grep
#sFiles <- list.files(path = dirname(opt$input_files_pattern), pattern = basename(opt$input_files_pattern), full.names = T)
sFiles <- data.table::fread(opt$input_files_pattern, header = FALSE)$V1
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
        if(nrow(dt) == 0){
            return(NULL)
        }
        return(dt)
    } else {
        return(NULL)
    }
}, .progress = TRUE)

# bn <- lapply(dFiles, ncol) |> unlist(); which(bn != 12)
# remove the NULLs
dFiles <- Filter(Negate(is.null), dFiles)
# print(head(dFiles[[3]]))
# bind the rows
dt_summary <- dplyr::bind_rows(dFiles)
# rename the columns
# print(head(dFiles[[1]]))
# print(head(dt))
dt_summary <- dt_summary %>% dplyr::select(-gene_name) %>% dplyr::rename(tfbs = gene)
#colnames(dt)[colnames(dt) == 'gene'] <- 'tfbs'
print(head(dt_summary))
# write the file
data.table::fwrite(dt_summary, opt$output_file, sep = "\t", quote = F, row.names = F)
