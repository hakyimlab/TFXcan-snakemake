# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript process_summary_statistics.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--summary_stats_file", help='[Input] A GWAS summary statistics file; should be a tsv file with columns: chrom, pos, ref, alt, pval, beta, se, zscore'),
    make_option("--output_folder", help='[Output] The output folder'),
    make_option("--annotation_file", help = '[Input] A file that contains allele information on the reference population data'),
    make_option("--[Input] pvalue_threshold", default=5e-8, type='numeric', help = 'the pvalue threshold for significance; default is 5e-8'),
    make_option('--[Output] diagnostics_file', type='character', default=NULL, help='A file to write diagnostics to; default is NULL i.e no diagnostics file will be written')
)


opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glue)
library(bigsnpr)

print(opt)

chrom_filter <- c(1:22)
if(!dir.exists(opt$output_folder)){dir.create(opt$output_folder)}

if(!file.exists(opt$annotation_file)){
    stop('ERROR - The reference annotation file does not exist.')
}

dt <- data.table::fread(opt$summary_stats_file) 

# important checks
if(!all(c('chrom', 'pos', 'ref', 'alt', 'pval', 'beta', 'se') %in% colnames(dt))){
    stop('ERROR - The summary statistics file must have columns: chrom, pos, ref, alt, pval, beta, se')
}

# set columns as types
dt <- dt %>%
    dplyr::mutate(chrom = as.numeric(chrom), pos = as.numeric(pos), ref = as.character(ref), alt = as.character(alt), pval = as.numeric(pval), beta = as.numeric(beta), se = as.numeric(se))

dt <- dt %>%
    dplyr::mutate(chrom = dplyr::case_when(
        grepl('chr', chrom, fixed=TRUE) ~ gsub('chr', '', chrom),
        .default = as.character(chrom)
    )) %>%
    dplyr::filter(chrom %in% chrom_filter) %>%
    dplyr::mutate(SNP = paste(chrom, pos, ref, alt, sep='_')) %>% 
    dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1) %>% # still need to properly resolve this
    dplyr::mutate(gwas_significant=ifelse(pval <= opt$pvalue_threshold, 'YES', 'NO'), chrom = as.numeric(chrom))
    
if(!'rsid' %in% colnames(dt)){
    dt <- dt %>%
        dplyr::mutate(rsid = paste(chrom, pos, ref, alt, sep=':'))
}

if(!'zscore' %in% colnames(dt)){
    dt <- dt %>%
        dplyr::mutate(zscore = beta/se)
}

# dt <- dt %>% dplyr::mutate(chrom = as.numeric(chrom))
## use purrr to match snps
annotation <- data.table::fread(opt$annotation_file) %>% 
    data.table::setDT() %>%
    dplyr::rename(a0=ref_vcf, a1=alt_vcf)

# split and match the snps 
split_dt <- dt %>%
    base::split(f=.$chrom)
sumstat_chroms <- names(split_dt)

# split and match the snps 
split_annot <- annotation %>%
    base::split(f=.$chr)
annot_chroms <- names(split_annot)

matched_stats <- purrr::map(sumstat_chroms, function(ch){
    x <- split_dt[[ch]] %>% dplyr::rename(chr=chrom, a0=ref, a1=alt)
    annot <- split_annot[[ch]]
    # match snps
    dt.matched_stats <- tryCatch({
        bigsnpr::snp_match(x, annot, return_flip_and_rev = TRUE, match.min.prop = 0) |> data.table::setDT() %>% 
            dplyr::select(chrom=chr, pos, ref=a0, alt=a1, rsid, varID, beta, se, zscore, pval, gwas_significant)
    }, error = function(e){
        message(glue('WARNING - Not enough variants have been matched for chromosome {ch}.'))
        print(e)
        return(NULL)
    })
    return(dt.matched_stats)
}, .progress = TRUE)

names(matched_stats) <- sumstat_chroms
matched_stats <- Filter(Negate(is.null), matched_stats) 

# print(names(matched_stats))

purrr::map(names(matched_stats), function(ns){
    data.table::fwrite(matched_stats[[ns]], file=glue('{opt$output_folder}/chr{ns}.sumstats.txt.gz'), 
        quote=F, col.names=T, row.names=F, sep='\t', compress='gzip')
}, .progress = TRUE)















# opt <- list()
# opt$summary_stats_file <- "/beagle3/haky/users/shared_data/from_david/asthma_children_hg38_chr17_35000000_45000000_Format.tsv.gz"
# opt$annotation_file <- '/project2/haky/Data/1000G/population_data/EUR/annot_files/EUR.snp_annot.txt'
# opt$pvalue_threshold <- 5e-18
# opt$diagnostics_file <- NULL


# setwd('/beagle3/haky/users/temi/projects/TFXcan-snakemake')
# opt <- list()
# opt$summary_stats_file <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/sumstats/type_2_diabetes.suzuki.gwas_sumstats.hg38.processed.txt.gz'
# opt$pvalue_threshold <- 5e-8
# opt$output_folder <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/T2D_2025-01-24/processed_sumstats/t2d_suzuki'
# opt$annotation_file <- '/project2/haky/Data/1000G/population_data/EUR/annot_files/EUR.snp_annot.txt'
# opt$diagnostics_folder <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/T2D_2025-01-24/diagnostics'


# x <- split_dt[['1']]
# x <- x %>% dplyr::rename(chr=chrom, a0=ref, a1=alt)
# annot <- annotation %>% dplyr::filter(chr == '1') #%>% dplyr::mutate(chr = as.numeric(chr))
# # match snps
# matched_stats <- tryCatch({
#     mm <- bigsnpr::snp_match(x, annot, return_flip_and_rev = TRUE, strand_flip = TRUE) |> data.table::setDT() %>% 
#         dplyr::select(chrom=chr, pos, ref=a1, alt=a0, rsid, varID, beta, se, zscore, pval, qval, gwas_significant)
#     return(mm)
# }, error = function(e){
#     print(glue('ERROR -  Not enough variants have been matched for chromosome {ch}.'))
#     return(NULL)
# })

# bx <- x[1:100000, ] %>% dplyr::mutate(chr = paste('chr', chr, sep=''))
# xx <- snp_modifyBuild(
#   bx,
#   '/project2/haky/Data/liftover/liftOver2',
#   from = "hg19",
#   to = "hg38",
#   check_reverse = TRUE,
#   #local_chain = '/project2/haky/Data/liftover/ucsc_chainfiles/hg19ToHg38.over.chain.gz',
#   base_url = "https://hgdownload.soe.ucsc.edu/goldenPath/"
# )

# nn <- names(split_dt)[1:3





# x <- split_dt[['17']]
# x <- x %>% dplyr::rename(chr=chrom, a0=ref, a1=alt)
# annot <- annotation %>% dplyr::filter(chr == 17) #%>% dplyr::mutate(chr = as.numeric(chr))
# # match snps
# mm <- bigsnpr::snp_match(x, annot, return_flip_and_rev = TRUE, match.min.prop = 0) %>% 
#     data.table::setDT() %>% 
#     dplyr::select(chrom=chr, pos, ref=a0, alt=a1, rsid, varID, beta, se, zscore, pval, gwas_significant)

# head(mm)