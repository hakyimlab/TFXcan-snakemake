# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--input_files", help='a pattern for the predicted enpact scores'),
    make_option("--output_db"),
    make_option("--output_file"),
    make_option("--individuals")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

library(data.table)
library(tidyverse)
library(glue)
library(abind)

# setwd("/project2/haky/temi/projects/enpact-predict-snakemake")
# opt <- list()
# opt$input_files <- "data/BC_GWAS/enpact_predictions/breast_cancer/{}.breast_cancer.aggByCollect.2024-06-06.csv.gz"
# opt$individuals <- "HG00353,HG01613,NA12750,HG00180,NA20513"

# split the individuals by the separator
individuals <- base::strsplit(opt$individuals, split = ',')[[1]]

# create a pattern given the input pattern
scores_files <- sapply(individuals, function(idi){
    gsub(pattern = "\\{\\}", replacement = idi, x = opt$input_files)
})

# 
retained_names <- c()

dt_files <- purrr::map(.x=seq_along(scores_files), function(i){
    each_file <- scores_files[i]
    if(file.exists(each_file)){
        dt <- data.table::fread(each_file) %>%
            tibble::column_to_rownames('locus')
        retained_names <<- append(retained_names, names(scores_files)[i])
        return(dt)
    }   
})

# collect the names of the locus and individuals
loci <- purrr::map(.x=dt_files, .f=rownames) %>%
    base::Reduce(intersect, .) |> unique()
models <- purrr::map(.x=dt_files, .f=colnames) %>%
    base::Reduce(intersect, .) |> unique()

print(models)

# ensure the dimensions are matched
dt_files <- purrr::map(.x=dt_files, function(each_dt){
    X <- as.matrix(each_dt)
    return(X[loci, models])
})

names(dt_files) <- retained_names

# combine into an array
myarray <- abind::abind(dt_files, along=3)
dimarray <- dim(myarray)

reshapedarray <- aperm(myarray, c(1, 3, 2), resize=TRUE)

print(dimarray)
print(dimnames(reshapedarray))
saveRDS(reshapedarray, file = opt$output_db, compress = "gzip")

# process and prepare for predictDB; create a predictDB input EnpactScores file
tarray <- as.array(reshapedarray)
tarray <- aperm(tarray, c(3, 1, 2), resize=TRUE)
dta <- dim(tarray)
farray <- as.array(tarray)
dim(farray) <- c(dta[1] * dta[2], dta[3])

dnames <- dimnames(reshapedarray)

mnames <- dnames[[3]]
lnames <- dnames[[1]]

rnames <- expand.grid(mnames, lnames, KEEP.OUT.ATTRS = TRUE, stringsAsFactors = FALSE) %>%
    tidyr::unite(mn, Var1:Var2, sep = ':') %>% dplyr::pull(mn)

rownames(farray) <- rnames
colnames(farray) <- dnames[[2]]
fdt <- as.data.frame(farray) %>% tibble::rownames_to_column('NAME')
data.table::fwrite(fdt, file = opt$output_file, compress = "gzip", sep = '\t', quote = F, col.names = T, row.names = F)