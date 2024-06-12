# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--enpact_scores_file", help='input: A file containing the Enpact scores for the models'),
    make_option("--formatted_escores_file", help='output: The formatted Enpact scores file'),
    make_option("--formatted_annot_file", help='output: The formatted annotation file'),
    make_option("--include_models", default = NULL),
    make_option("--filtered_GWAS_file", help='input: The filtered GWAS file'),
    make_option("--blacklist", default = NULL, help='input: Remove loci within these blacklist regions')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glue)

# opt <- list()
# opt$include_models <- NULL
# opt$blacklist <- '/project/haky/users/temi/projects/TFXcan-snakemake/metadata/HLA_blocks.txt'
# opt$filtered_GWAS_file <- '/project/haky/users/temi/projects/TFXcan-snakemake/data/PCRISK_GWAS_2024-06-07/collection/pc_risk.filteredGWAS.txt.gz'
# opt$enpact_scores_file <- '/project/haky/users/temi/projects/TFXcan-snakemake/data/PCRISK_GWAS_2024-06-07/enpactdb/pc_risk.enpact_scores.txt.gz'
# dt$NAME[startsWith(dt$NAME, 'HSF1_Mammary-Gland')] |> print()


print(opt)
if(!is.null(opt$include_models)){
    include_models <- data.table::fread(include_models, header = F) %>%
        pull(V1)
} else {
    include_models <- NULL
}


enpact_scores_dt <- data.table::fread(opt$enpact_scores_file)

if(!is.null(opt$blacklist)){
    blacklist <- data.table::fread(opt$blacklist, header = T) %>%
        dplyr::select(chrom, start, end)

    # read in the annotation file too and filter out the blacklist regions
    annot_dt_1 <- data.table::fread(opt$filtered_GWAS_file) %>%
        dplyr::mutate(chr = gsub('chr', '', chrom), start = pos, end = pos + 1, gene_type='protein_coding') %>%
        dplyr::filter(!dplyr::between(pos, blacklist$start, blacklist$end))
} else {
    blacklist <- NULL
    # read in the annotation file too
    annot_dt_1 <- data.table::fread(opt$filtered_GWAS_file) %>%
        dplyr::mutate(chr = gsub('chr', '', chrom), start = pos, end = pos + 1, gene_type='protein_coding')
}

annot_dt_2 <- enpact_scores_dt %>%
    dplyr::select(NAME) %>%
    tidyr::separate(col = NAME, into=c('model', 'locus'), sep=':', remove=F) %>%
    tidyr::separate(col = locus, into = c('chrom', 'start', 'end'), sep= '_', remove=F) %>%
    dplyr::mutate(chr = gsub('chr', '', chrom), start = as.integer(start), end= as.integer(end)) %>% 
    dplyr::select(chr, start, end, NAME)

annot_dt <- dplyr::inner_join(annot_dt_1, annot_dt_2, by = c('chr' = 'chr', 'start' = 'start', 'end' = 'end')) %>%
    dplyr::select(chr, start, end, gene_id=NAME, gene_name=NAME, gene_type)

annot_dt <- annot_dt %>%
    dplyr::mutate(gene_name = gsub(':', '_',gene_name), gene_id = gene_name) 

enpact_scores_dt <- enpact_scores_dt %>%
    dplyr::mutate(NAME = gsub(':', '_', NAME))

annot_dt <- annot_dt %>% dplyr::mutate(gene_id = gsub('-', '', gene_id), gene_name = gsub('-', '', gene_name))
enpact_scores_dt <- enpact_scores_dt %>% dplyr::mutate(NAME = gsub('-', '', NAME))

enpact_scores_dt <- enpact_scores_dt %>%
    dplyr::filter(NAME %in% annot_dt$gene_id)

data.table::fwrite(annot_dt, file = opt$formatted_annot_file, quote=F, row.names = F, col.names = T, sep = '\t')
data.table::fwrite(enpact_scores_dt, file = opt$formatted_escores_file, quote=F, row.names = F, col.names = T, sep = '\t')




# dt <- data.table::fread('/project2/haky/Data/1000G/metadata/EUR.samples.txt', header=F) %>%
#     dplyr::select(V1) %>%
#     data.table::fwrite('/project2/haky/temi/projects/TFXcan-snakemake/metadata/EUR_individuals.1KG.txt', col.names=F, row.names=F, quote = F)


#names(enpact_scores_dt) <- enpact_models

# calculate the variance and all that ========
# enpact_scores_var <- enpact_scores_dt %>%
#     tibble::column_to_rownames('NAME') %>%
#     as.matrix() %>%
#     base::apply(., 1, var) %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column('model') %>%
#     tidyr::separate('model', into=c('model', 'locus'), sep=':', remove=T) %>%
#     dplyr::rename(varEnpactScore=3)

# enpact_scores_var %>%
#     dplyr::mutate(model = factor(model)) %>%
#     ggplot(aes(model, log(varEnpactScore))) +
#     geom_point(shape='.', size=200) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# enpact_scores_var %>%
#     dplyr::mutate(locus = factor(locus)) %>%
#     ggplot(aes(locus, log(varEnpactScore))) +
#     geom_point(shape='.') +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# individuals <- base::strsplit(opt$individuals, split=',')[[1]]
# enpact_scores <- lapply(individuals, function(each_individual){
#     ff <- gsub('\\{\\}', each_individual, opt$input_file_pattern)
#     print(ff)
#     if(!file.exists(ff)){
#         return(NULL)
#     } else {
#             dt <- data.table::fread(ff)
            
#             if(!is.null(include_models)){
#                 return(dt %>%dplyr::select(all_of('locus', include_models))
#                 )
#             } else {
#                 return(dt)
                
#             }
#     }
# })

# names(enpact_scores) <- individuals

# # get the TF tissue models:
# enpact_models <- colnames(enpact_scores[[1]])[-1]
# #print(enpact_models)

# enpact_scores_dt <- lapply(seq_along(enpact_models), function(i){
#     each_model <- enpact_models[[i]]
#     ot <- lapply(seq_along(enpact_scores), function(i){
#         each_pred <- enpact_scores[[i]]
#         dt <- each_pred %>% dplyr::select(all_of(c('locus', each_model)))
#         return(dt)
#     }) %>%
#         Filter(Negate(is.null), .) %>%
#         purrr::reduce(., dplyr::inner_join, by=c('locus' = 'locus')) %>%
#         data.table::setnames(., colnames(.), new=c('locus', individuals)) %>%
#         dplyr::mutate(model = enpact_models[i]) %>%
#         tidyr::unite(col = 'NAME',  c('model', 'locus'), sep=':')

#     return(ot)
# }) %>%
#     do.call('rbind', .)