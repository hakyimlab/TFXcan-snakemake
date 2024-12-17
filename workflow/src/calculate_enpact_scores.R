# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--input_file", help='A transcription factor e.g. AR'),
    make_option("--output_file", help='An output file'),
    make_option("--enpact_models_directory"),
    make_option("--enpact_models_metadata"),
    make_option("--personalized", default=FALSE, action="store_true"),
    make_option("--verbose", default=FALSE, action="store_true")
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glmnet)
library(glue)

# setwd("/project/haky/users/temi/projects/TFXcan-snakemake/")
# opt <- list()
# opt$input_file <- "data/prostate_cancer_risk_2024-09-30/aggregated_predictions/prostate_cancer_risk/NA20827_aggByCollect_prostate_cancer_risk.csv.gz"
# opt$output_file <- "/scratch/midway3/temi/enpact_predictions/prostate_cancer_risk/NA20827.prostate_cancer_risk.aggByCollect.2024-09-30.csv.gz"
# opt$enpact_models_directory <-  "/beagle3/haky/users/temi/projects/TFPred-snakemake"
# opt$enpact_models_metadata <-  "metadata/models734.prostate.txt"


mt <- data.table::fread(opt$enpact_models_metadata, na.strings = '') %>%
    dplyr::select(model, path) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(path))

# mt <- data.table::fread('/project/haky/users/temi/projects/TFXcan-snakemake/metadata/models734.prostate.txt')

# filter and get the models

# nn <- c()
# models_list <- apply(mt, 1, function(each_row){
#     name <- each_row[1]
#     model <- each_row[2]
#     mpath <- file.path(opt$enpact_models_directory, model)
#     if(file.exists(mpath)){
#         nn <<- append(nn, name)
#         return(mpath)
#     }
# })


mpaths <- file.path(opt$enpact_models_directory, mt$path)
names(mpaths) <- mt$model
mpaths <- mpaths[file.exists(mpaths)]
models_list <- as.list(mpaths)

# names(models_list) <- nn
# models_list <- as.list(models_list)

# print(models_list[1:3])

# predict
dt <- data.table::fread(opt$input_file)
X <- as.matrix(dt[, -c(1)])

#models_list['AR_Prostate'] <- "/project2/haky/temi/projects/TFPred-snakemake/output/models/cistrome_AR_Prostate_2023-11-14/aggByCollect_AR_Prostate.logistic.rds"

# models_list <- c(models_list[1:4], models_list['AR_Prostate'])

predictions <- purrr::map(.x=models_list, .f = function(each_model){
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

# mm <- readRDS("/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2024-07-26/models/ARNT_BoneMarrow/ARNT_BoneMarrow_2024-07-26.logistic.rds")
# predict(mm, X, s = "lambda.1se", type = 'link') |> as.vector()

predictions <- Filter(Negate(is.null), predictions) %>%
    do.call('cbind', .)
rownames(predictions) <- dt$id
predictions %>%
    as.data.frame() %>%
    tibble::rownames_to_column('locus') %>%
    data.table::fwrite(., file=opt$output_file, compress='gzip', quote=F, sep=',', col.names = TRUE, row.names = FALSE)
print(glue('INFO - {opt$output_file} has been saved.'))