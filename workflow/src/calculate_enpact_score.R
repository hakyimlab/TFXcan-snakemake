# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--input_file", help='A transcription factor e.g. AR'),
    make_option("--output_file", help='An output file'),
    make_option("--enpact_models_directory"),
    make_option("--filters_date"),
    make_option("--filters_type")
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glmnet)
library(glue)

# setwd("/project2/haky/temi/projects/TFXcan-snakemake/")
# opt <- list()
# opt$input_file <- "data/aggregated_predictions/test_aggByCollect_testosterone.csv.gz"
# opt$output_file <- "data/enpact_predictions/test_testosterone_2023-12-07.csv.gz"
# opt$enpact_models_directory <-  "/project2/haky/temi/projects/TFPred-snakemake/output/models"
# opt$filters_date <- "2023-12-01"
# opt$filters_type <- "logistic"

# filter and get the models
models_list <- list.files(opt$enpact_models_directory, pattern = ifelse(
        typeof(opt$filters_date) == 'character', opt$filters_date, NULL
        ), 
    full.names = TRUE) %>%
        list.files(., pattern = ifelse(
        typeof(opt$filters_type) == 'character', opt$filters_type, NULL
        ),
    full.names = TRUE)

models_names <- base::strsplit(models_list, split = '/') %>%
    sapply(., function(each_name){
        each_name[length(each_name)]
    }) %>%
    base::strsplit(., split = '_|\\.') %>%
    sapply(., function(each_name){
        paste(each_name[2], each_name[3], sep='_')
    })

if(length(models_list) != length(models_names)){
    stop("Something is wrong")
} else {
    names(models_list) <- models_names
}

# predict
dt <- data.table::fread(opt$input_file)
X <- as.matrix(dt[, -c(1)])

models_list['AR_Prostate'] <- "/project2/haky/temi/projects/TFPred-snakemake/output/models/cistrome_AR_Prostate_2023-11-14/aggByCollect_AR_Prostate.logistic.rds"

# models_list <- c(models_list[1:4], models_list['AR_Prostate'])

predictions <- purrr::map(.x=models_list, .f = function(each_model){
    #print(each_model)
    if(file.exists(each_model)){
        model <- readRDS(each_model)
        tryCatch({
            link_pred <- predict(model, X, s = "lambda.1se", type = 'link') |> as.vector()
            return(link_pred)
        }, error = function(e){
            return(NULL)
        })
    } else {
        return(NULL)
    }
}) 

predictions <- Filter(Negate(is.null), predictions) %>%
    do.call('cbind', .)
rownames(predictions) <- dt$id
predictions %>%
    as.data.frame() %>%
    tibble::rownames_to_column('locus') %>%
    data.table::fwrite(., file=opt$output_file, compress='gzip', quote=F, sep=',', col.names = TRUE, row.names = FALSE)
print(glue('INFO - {opt$output_file} has been saved.'))
