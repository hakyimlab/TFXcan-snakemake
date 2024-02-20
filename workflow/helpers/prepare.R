


suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--pop", help='A list of files to combine'),
    make_option("--snp_annot_pattern", help='reference panel correlation matrix'),
    make_option("--geno_pattern", help='reference panel correlation matrix'),
    make_option("--output_folder", help='summary statistics sample size')
)

opt <- parse_args(OptionParser(option_list=option_list))  

library(data.table)
library(tidyverse)
library(susieR)
library(glue)

valid_chr <- c(1:22)
# geno_files <- file.path(files_folder, glue('*{valid_chr}*.geno.txt.gz')) |> Sys.glob() |> unique()
# snpannot_files <- file.path(files_folder, glue('*{valid_chr}*.snp_annot.txt.gz')) |> Sys.glob() |> unique()
# print(glue('INFO - Found {length(geno_files)} genotype .txt.gz files to merge.')) 
# print(glue('INFO - Found {length(snpannot_files)} snp annotation .txt.gz files to merge.'))

#print(geno_files) ; print(snpannot_files)

# prepare and save geno patterns
geno_file <- purrr::map(valid_chr, function(each_chr){
    ff <- gsub('\\{\\}', each_chr, opt$geno_pattern)
    if(file.exists(ff)){
        geno_dt <- data.table::fread(ff)
    }
}, .progress = T) %>% base::Filter(Negate(is.null), .) %>% 
    dplyr::bind_rows(.)

data.table::fwrite(geno_file, file=glue("{opt$output_folder}/{opt$pop}.geno.txt"), sep='\t', quote=F, row.names=F)
data.table::fwrite(geno_file, file=glue("{opt$output_folder}/{opt$pop}.geno.txt.gz"), sep='\t', quote=F, row.names=F, compress='gzip')

rm('geno_file')

# prepare and save snp files
snpannot_file <- purrr::map(valid_chr, function(each_chr){
    ff <- gsub('\\{\\}', each_chr, opt$snp_annot_pattern)
    if(file.exists(ff)){
        dt <- data.table::fread(ff) %>% dplyr::mutate(varID=as.character(varID))
        return(dt)
    }
}) %>% base::Filter(Negate(is.null), .) %>%
    dplyr::bind_rows(.)


data.table::fwrite(snpannot_file, file=glue("{opt$output_folder}/{opt$pop}.snp_annot.txt"), sep='\t', quote=F, row.names=F)
data.table::fwrite(snpannot_file, file=glue("{opt$output_folder}/{opt$pop}.snp_annot.txt.gz"), sep='\t', quote=F, row.names=F, compress='gzip')

# if(length(geno_files) == 0 | length(snpannot_files) == 0){
#     stop(glue('ERROR - Either there are no geno .txt.gz files snp_annot .txt.gz files or both in {files_folder}.'))
# } else {

#     all_geno_files <- purrr::map(geno_files, function(gfile){
#         #ff <- glue('{files_folder}/ALL.{chrom}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.geno.txt')
#         if(file.exists(gfile)){
#             dt <- data.table::fread(gfile) %>% dplyr::mutate(varID=as.character(varID))
#             return(dt)
#         }
#     }, .progress=TRUE)
#     geno_file <- dplyr::bind_rows(all_geno_files)
#     data.table::fwrite(geno_file, file=glue("{output_folder}/all_chrs.geno.txt"), sep='\t', quote=F, row.names=F)
#     data.table::fwrite(geno_file, file=glue("{output_folder}/all_chrs.geno.txt.gz"), sep='\t', quote=F, row.names=F, compress='gzip')

#     rm('all_geno_files')

#     all_snpannot_files <- purrr::map(snpannot_files, function(sfile){
#         #ff <- glue('{files_folder}/ALL.{chrom}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.snp_annot.txt')
#         if(file.exists(sfile)){
#             dt <- data.table::fread(sfile) %>% dplyr::mutate(chr=as.character(chr))
#             return(dt)
#         }
#     }, .progress=TRUE)
#     snpannot_file <- dplyr::bind_rows(all_snpannot_files)
#     data.table::fwrite(snpannot_file, file=glue("{output_folder}/all_chrs.snp_annot.txt"), sep='\t', quote=F, row.names=F)
#     data.table::fwrite(snpannot_file, file=glue("{output_folder}/all_chrs.snp_annot.txt.gz"), sep='\t', quote=F, row.names=F, compress='gzip')
    
#     rm('all_snpannot_files')

#     if(combine == 'yes'){
#         metadata='all_chrs'
#         print("INFO - combining snp and geno files into one.")
#         merged_dt <- snpannot_file %>% 
#             dplyr::select(chr, varID, pos, ref_vcf, alt_vcf, maf) %>% 
#             dplyr::inner_join(geno_file, by=c('varID' = 'varID'))

#         ff1 <- glue("{output_folder}/{metadata}.text_dosages.txt.gz")
#         data.table::fwrite(merged_dt, file=ff1, sep='\t', quote=F, row.names=F, col.names=F, compress = 'gzip')

#         print("INFO - writing out samples file.")
#         data.frame(colnames(geno_file)[-1], colnames(geno_file)[-1]) %>%
#             data.table::fwrite(file=glue("{output_folder}/samples.text_dosages.txt"), sep='\t', quote=F, row.names=F, col.names=F)

#     }
# }