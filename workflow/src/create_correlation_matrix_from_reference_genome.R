# this runs one locus at a time

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--reference_genome_txt_file", help='A list of files to combine'),
    make_option("--output_file", help='A list of files to combine')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(susieR)
library(glue)
library(Matrix)

### collect and prepare covariances
print(glue('INFO - collecting and preparing covariance .rds file'))
data.table::fread(opt$reference_genome_txt_file) %>%
    tibble::column_to_rownames('varID') %>%
    as.matrix() %>%
    t() %>%
    scale() %>%
    cor() %>%
    as(., 'sparseMatrix') %>%
    saveRDS(file = glue('{opt$output_file}'))

# saveRDS(R, file = glue('{opt$output_file}'))