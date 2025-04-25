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

setwd("/beagle3/haky/users/temi/projects/TFXcan-snakemake/")
opt <- list()
opt$input_file <- "data/prostate_cancer_risk_2024-09-30/aggregated_predictions/prostate_cancer_risk/NA20827_aggByCollect_prostate_cancer_risk.csv.gz"
opt$output_file <- "/scratch/midway3/temi/enpact_predictions/prostate_cancer_risk/NA20827.prostate_cancer_risk.aggByCollect.2024-09-30.csv.gz"
opt$enpact_weights <- "/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/statistics/ENPACT_734_2025-04-24.compiled_weights.lambda.1se.txt.gz"
# opt$enpact_models_directory <-  "/beagle3/haky/users/temi/projects/TFPred-snakemake"
# opt$enpact_models_metadata <-  "metadata/models734.prostate.txt"


#
opt$data_directory <- "/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2024-09-30/aggregated_predictions/prostate_cancer_risk"
opt$individuals <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/metadata/EUR_individuals.1KG.txt'
opt$output_file <- '/beagle3/haky/users/temi/projects/Enpact/misc/reruns/enpact_predictions/ENPACT_48.predictions.2025-04-24.rds.gz'


eur_individuals <- data.table::fread(opt$individuals, header = FALSE)
#enpact_features_directory <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2024-09-30/aggregated_predictions/prostate_cancer_risk'
individual_enpact_features <- glue::glue("{opt$data_directory}/{eur_individuals$V1}_aggByCollect_prostate_cancer_risk.csv.gz") |> as.vector()
names(individual_enpact_features) <- eur_individuals$V1
individual_enpact_features <- individual_enpact_features[file.exists(individual_enpact_features)]

 # weights 
weights <- data.table::fread(opt$enpact_weights) %>% 
    dplyr::select(-feature) %>% as.matrix()

Y_hats <- purrr:::map(.x=individual_enpact_features, .f = function(each_file){
    dt <- data.table::fread(each_file)
    X <- as.matrix(dt[, -c(1)])

    stopifnot(dim(X)[2] == dim(weights)[1])

    # prediction
    y_hat <- X %*% weights
    rownames(y_hat) <- dt$id

    return(y_hat)

}, .progress = TRUE)

library(abind)

# collect the names of the locus and individuals
loci <- purrr::map(.x=Y_hats, .f=rownames) %>%
    base::Reduce(intersect, .) |> unique()
models <- purrr::map(.x=Y_hats, .f=colnames) %>%
    base::Reduce(intersect, .) |> unique()
individuals <- names(individual_enpact_features)

print(models)
print(length(individuals))

# ensure the dimensions are matched
out <- purrr::map(.x=Y_hats, function(each_dt){
    X <- as.matrix(each_dt)
    return(X[loci, models])
})

names(out) <- individuals

# combine into an array
myarray <- tryCatch({
    abind::abind(out, along=3)
}, error = function(e){
    abind::abind(out, along=2)
})

dimarray <- dim(myarray)
if(length(dimarray) == 2){
    #myarray <- array(myarray, dim=c(dimarray[1], 1, dimarray[2]))
    saveRDS(myarray, file = opt$output_file, compress = "gzip")
} else if(length(dimarray) == 3){
    reshapedarray <- aperm(myarray, c(1, 3, 2), resize=TRUE)
    # print(dimarray)
    # print(dimnames(reshapedarray))
    saveRDS(reshapedarray, file = opt$output_file, compress = "gzip")
}



reshapedarray <- readRDS(opt$output_file)

models <- dimnames(reshapedarray)[[3]]

# create the predictdb input

X_pdb <- purrr::map(.x=models, .f = function(each_file){
    X <- reshapedarray[, , each_file]
    rn <- rownames(X)
    # append the model to the rownames
    rownames(X) <- paste0(each_file, '_', rn)
    X <- as.data.table(X, keep.rownames = 'NAME')
    return(X)
}, .progress = TRUE)

names(X_pdb) <- models

# save the data
outdir <- '/beagle3/haky/users/temi/projects/Enpact/misc/reruns/predictdb_data'
purrr::map(names(X_pdb), .f = function(each_name){
    data.table::fwrite(X_pdb[[each_name]], file=glue::glue("{outdir}/prostate_cancer_risk.{each_name}.enpact_scores.tsv"), quote=F, sep='\t', col.names = TRUE, row.names = FALSE)
})


# copy the metadata

from_dir <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2024-09-30/predictdb'
to_dir <- '/beagle3/haky/users/temi/projects/Enpact/misc/reruns/predictdb_data'

files_to_copy <- glue::glue("{from_dir}/{models}/prostate_cancer_risk.tf_tissue_annot.txt") |> as.vector()
files_rename <- glue::glue("{to_dir}/prostate_cancer_risk.{models}.annotations.txt") |> as.vector()
sapply(seq_along(files_to_copy), function(i){
    file.copy(from = files_to_copy[i], to = files_rename[i], overwrite = FALSE)
})

# generate lenpact models

scratch_dir <- '/scratch/midway3/temi/lEnpact/pcr_retrain'

generate_sbatch <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/workflow/predictdb/generate_snp_predictors.sbatch'
cmds <- purrr::map(.x=models[1:3], .f = function(each_model){
    phenotype <- 'prostate_cancer_risk'
    s_path <- file.path(scratch_dir, each_model)
    if(!dir.exists(s_path)){
        dir.create(s_path, recursive = TRUE)
    }
    loci_annotation <- glue::glue("{to_dir}/{each_model}.prostate_cancer_risk.annotations.txt") |> as.vector()
    enpact_scores <- glue::glue("{to_dir}/{each_model}.prostate_cancer_risk.enpact_scores.tsv") |> as.vector()
    ref_genotypes <- '/project2/haky/Data/1000G/population_data/EUR/annot_files/EUR.geno.txt'
    ref_annotations <- '/project2/haky/Data/1000G/population_data/EUR/annot_files/EUR.snp_annot.txt'
    nextflow_executable <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/workflow/PredictDb-nextflow/main.nf'

    return(glue::glue("cd {s_path} && sbatch {generate_sbatch} {phenotype} {s_path} {loci_annotation} {enpact_scores} {ref_genotypes} {ref_annotations} {nextflow_executable}"))
    
}, .progress = TRUE)


sapply(cmds, function(each_cmd){
    system(each_cmd)
})












# data
dt <- data.table::fread(opt$input_file)
X <- as.matrix(dt[, -c(1)])

# weights 
weights <- data.table::fread(opt$enpact_weights) %>% 
    dplyr::select(-feature) %>% as.matrix()

stopifnot(dim(X)[2] == dim(weights)[1])

# prediction
y_hat <- X %*% weights
rownames(y_hat) <- dt$id

y_hat %>%
    as.data.frame() %>%
    tibble::rownames_to_column('locus') %>%
    data.table::fwrite(., file=opt$output_file, compress='gzip', quote=F, sep=',', col.names = TRUE, row.names = FALSE)
print(glue('INFO - {opt$output_file} has been saved.'))

#models_list['AR_Prostate'] <- "/project2/haky/temi/projects/TFPred-snakemake/output/models/cistrome_AR_Prostate_2023-11-14/aggByCollect_AR_Prostate.logistic.rds"

# models_list <- c(models_list[1:4], models_list['AR_Prostate'])

# predictions <- purrr::map(.x=models_list, .f = function(each_model){
#     if(file.exists(each_model)){
#         model <- readRDS(each_model)
#         tryCatch({
#             link_pred <- predict(model, X, s = "lambda.1se", type = 'link') |> as.vector()
#             return(link_pred)
#         }, error = function(e){
#             return(NULL)
#         })
#     } else {
#         return(NULL)
#     }
# }) 

# # mm <- readRDS("/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2024-07-26/models/ARNT_BoneMarrow/ARNT_BoneMarrow_2024-07-26.logistic.rds")
# # predict(mm, X, s = "lambda.1se", type = 'link') |> as.vector()

# predictions <- Filter(Negate(is.null), predictions) %>%
#     do.call('cbind', .)
# rownames(predictions) <- dt$id
# predictions %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column('locus') %>%
#     data.table::fwrite(., file=opt$output_file, compress='gzip', quote=F, sep=',', col.names = TRUE, row.names = FALSE)
# print(glue('INFO - {opt$output_file} has been saved.'))

















# mt <- data.table::fread(opt$enpact_models_metadata, na.strings = '') %>%
#     dplyr::select(model, path) %>%
#     dplyr::distinct() %>%
#     dplyr::filter(!is.na(path))

# # mt <- data.table::fread('/project/haky/users/temi/projects/TFXcan-snakemake/metadata/models734.prostate.txt')

# # filter and get the models

# # nn <- c()
# # models_list <- apply(mt, 1, function(each_row){
# #     name <- each_row[1]
# #     model <- each_row[2]
# #     mpath <- file.path(opt$enpact_models_directory, model)
# #     if(file.exists(mpath)){
# #         nn <<- append(nn, name)
# #         return(mpath)
# #     }
# # })


# mpaths <- file.path(opt$enpact_models_directory, mt$path)
# names(mpaths) <- mt$model
# mpaths <- mpaths[file.exists(mpaths)]
# models_list <- as.list(mpaths)

# names(models_list) <- nn
# models_list <- as.list(models_list)

# print(models_list[1:3])