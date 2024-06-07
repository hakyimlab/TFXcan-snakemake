# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--summary_stats_file", help='A transcription factor e.g. AR'),
    make_option("--output_folder", help='the output folder'),
    make_option("--phenotype", help = 'a GWAS phenotype'),
    make_option("--pvalue_threshold", default=5e-8, type='numeric', help = 'a GWAS phenotype'),
    make_option('--diagnostics_file', type='character', default=NULL, help='')
)

source('/project/haky/users/temi/projects/TFXcan-snakemake/workflow/src/modules.R')

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glue)

print(opt)

# opt <- list()
# opt$summary_stats_file <- '/project2/haky/temi/projects/TFXcan-snakemake/data/sumstats/asthma_children.liftover.logistic.assoc.tsv.gz'
# opt$pvalue_threshold <- 5e-8
# opt$output_folder <- '/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/testosterone'
# opt$phenotype <- 'asthma_children'
# opt$diagnostics_folder <- '/project2/haky/temi/projects/TFXcan-snakemake/data/diagnostics'

chrom_filter <- c(1:22)
if(!dir.exists(opt$output_folder)){dir.create(opt$output_folder)}

dt <- data.table::fread(opt$summary_stats_file) %>%
    dplyr::mutate(chrom = dplyr::case_when(
        grepl('chr', chrom, fixed=TRUE) ~ gsub('chr', '', chrom),
        .default = as.character(chrom)
    )) %>%
    dplyr::filter(chrom %in% chrom_filter) %>%
    dplyr::mutate(SNP = paste(chrom, pos, ref, alt, sep='_'), adjP = pval) %>%
    dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1) %>% # still need to properly resolve this
    dplyr::mutate(gwas_significant=ifelse(pval <= opt$pvalue_threshold, 'YES', 'NO')) 
    
if(!'rsid' %in% colnames(dt)){
    dt <- dt %>%
        dplyr::mutate(rsid = paste(chrom, pos, ref, alt, sep=':'))
}

dt <- dt %>%
    dplyr::group_by(chrom)

# dt <- data.table::fread(opt$summary_stats_file) %>%
#     dplyr::mutate(chrom = dplyr::case_when(
#         grepl('chr', chrom, fixed=TRUE) ~ gsub('chr', '', chrom),
#         .default = as.character(chrom)
#     )) %>%
#     dplyr::filter(chrom %in% chrom_filter) %>%
#     dplyr::mutate(SNP = paste(chrom, pos, ref, alt, sep='_')) %>%
#     dplyr::mutate(gwas_significant="YES") %>%
#     dplyr::group_by(chrom) 

split_dt <- dt %>%
    base::split(f=.$chrom)

nn <- names(split_dt)
lapply(seq_along(split_dt), function(i){
    data.table::fwrite(split_dt[[i]], file=glue('{opt$output_folder}/chr{nn[i]}.sumstats.txt.gz'), 
        quote=F, col.names=T, row.names=F, sep='\t', compress='gzip')
})


# if(!is.null(opt$diagnostics_file)){
#     if(!dir.exists(dirname(opt$diagnostics_file))){
#         dir.create(dirname(opt$diagnostics_file), recursive = TRUE)
#     }

#     print(glue('INFO - writing diagnostics for {opt$phenotype}'))

#     ah <- dt %>%
#         dplyr::group_by(chrom, gwas_significant) %>%
#         dplyr::summarise(n=n(), .groups = "drop") %>%
#         tidyr::complete(chrom, gwas_significant, fill = list(n=0)) %>%
#         tidyr::pivot_wider(id_cols = chrom, names_from=gwas_significant, values_from = n)

#     fil <- file(opt$diagnostics_file, open = "a")
#     cat("#### GWAS diagnostics", file = fil, sep = '\n')
#     cat(glue("## Number of GWAS significant loci at {opt$pvalue_threshold}"), file = fil, sep = '\n\n', append=T)
#     cat(paste0(colnames(ah), collapse = '\t'), file = fil, append = T, sep = '\n')
#     cat(apply(ah, 1, paste0, collapse='\t'), file = fil, append = T, sep = '\n')
#     close(fil)

#     out <- dt %>%
#         dplyr::select(chrom, pos, rsid, pval) %>%
#         setNames(c("chrom", "locus", 'id', 'pvalue')) %>%
#         prepare_manhattan_dt()

#     png(glue("{dirname(opt$diagnostics_file)}/{opt$phenotype}.gwas_manhattan.png"))
#     plot_manhattan_dt(out, signif= opt$pvalue_threshold, plot_id = FALSE)
#     dev.off()

#     png(glue("{dirname(opt$diagnostics_file)}/{opt$phenotype}.gwas_qqunif.png"))
#     qqunif(dt$pval)
#     dev.off()

# }
    
# dt %>%
#     dplyr::pull(SNP) %>% 
#     as.data.frame() %>%
#     data.table::fwrite(., file=glue('{opt$output_file_basename}.SNPsForLD'), quote=F, col.names=F, row.names=F)

# sumstats <- data.table::fread(dt) %>%
#     dplyr::mutate(chrom=as.character(gsub('chr', '', chrom))) %>%
#     dplyr::filter(chrom == opt$chromosome)

# # QC filter
# data.table::fwrite(as.data.frame(sumstats$SNP), file=glue('{opt$snp_list}'), quote=F, row.names=F, col.names=F)
