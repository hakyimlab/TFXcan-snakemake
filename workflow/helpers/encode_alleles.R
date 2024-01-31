# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--bfile", help='A list of files to combine'),
    make_option("--output_prefix", help='reference panel correlation matrix')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(glue)
library(tidyverse)
library(data.table)

glue('INFO - processing {opt$bfile}')

snp_annot <- fread(glue::glue("{opt$bfile}.bim")) %>% 
    setnames(.,names(.), c("chr", "snp", "CM", "pos", "alt_vcf", "ref_vcf")) %>%
    dplyr::mutate(chr = gsub('chr', '', chr)) %>%
    dplyr::mutate(rsid = paste(chr, pos, ref_vcf, alt_vcf, sep=':')) %>%
    dplyr::mutate(maf = 0.01) %>%  
    dplyr::mutate(varID = str_replace_all(rsid,":","_")) %>%
    dplyr::select(chr, pos, varID, ref_vcf, alt_vcf, maf, rsid)

data.table::fwrite(snp_annot, file=glue("{opt$output_prefix}.snp_annot.txt.gz"), sep='\t', quote=F, row.names=F, compress='gzip')


genotype <- readr::read_table(glue::glue("{opt$bfile}.traw")) %>% 
    dplyr::mutate(CHR = gsub('chr', '', CHR)) %>%
    tidyr::unite('varID', CHR, POS, COUNTED, ALT, sep = '_', remove=FALSE) %>%
    dplyr::mutate(varID = paste(varID, sep='_')) %>%
    dplyr::select(-c(CHR,`(C)M`,POS, COUNTED, ALT, SNP)) %>% 
    setnames(.,names(.),gsub("0_", "", colnames(.)))

data.table::fwrite(genotype, file=glue("{opt$output_prefix}.geno.txt.gz"), sep='\t', quote=F, row.names=F, compress='gzip')





# 1_51479_T_A_b37 1       2       2       2       2       2       1       1       1       2       1       1       2       1       2       1       1       2       1       2    >
# 1_54490_G_A_b37 1       2       2       2       2       2       1       1       1       2       1       1       2       1       2       1       1       2       1       2    >
# 1_54708_G_C_b37 1       0       1       1       1       0       2       1       1       1       1       2       2       2       0       1       1       1       1       2    >
# 1_54716_C_T_b37 1       0       1       1       1       1       2       1       1       1       1       2       2       2       1       1       1       1       1       2    >
# 1_54753_T_G_b37 2