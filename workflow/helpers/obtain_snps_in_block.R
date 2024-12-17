# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--chromosome", help='A list of files to combine'),
    make_option("--partition_file", help='reference panel correlation matrix'),
    make_option("--plink_bim_file", help='A list of files to combine'),
    make_option("--LD_folder", help='A list of files to combine')
)

opt <- parse_args(OptionParser(option_list=option_list))  

library(data.table)
library(tidyverse)
library(susieR)
library(glue)

print(opt)

# opt <- list()
# opt$chromosome <- '1'
# opt$partition_file <- '/project2/haky/Data/LD_blocks/hg38/EUR/hg38_fourier_ls-all.bed'
# opt$plink_bim_file <- '/project2/haky/Data/1000G/plink_ref_panel/ALL.chr1.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.bim'
# opt$LD_folder <- '/project2/haky/Data/1000G/LD'

bfile <- data.table::fread(opt$plink_bim_file) 
colnames(bfile) <- c('chrom', 'SNP', 'pos_cm', 'bp_coord', 'A1', 'A2')

if(!opt$chromosome %in% bfile$chrom){
    stop(glue('ERROR - {opt$chromosome} does not exist'))
}

chr_folder <- file.path(opt$LD_folder, paste0('chr', opt$chromosome))
if(!dir.exists(chr_folder)){
    dir.create(chr_folder)
}
    
pfile <- data.table::fread(opt$partition_file, header=F) %>%
    dplyr::filter(V1 == paste0('chr', opt$chromosome)) %>%
    dplyr::mutate(block=1:nrow(.))
colnames(pfile) <- c('chrom', 'start', 'stop', 'block')


pfile %>%
    base::split(., f=.$block) %>%
    purrr::map(.x = ., .f=function(block_i){
        fname <- paste(block_i$chrom, block_i$start, block_i$stop, sep='_')
        out_null <- bfile %>%
            dplyr::filter(dplyr::between(bp_coord, block_i$start, block_i$stop)) %>%
            dplyr::mutate(SNP=paste0('chr', SNP)) %>%
            dplyr::pull(SNP) %>%
            as.data.frame() %>%
            data.table::fwrite(., file=glue('{chr_folder}/{fname}.SNPsLDBlock.txt'), col.names=F, row.names=F, quote=F)
        
        return(fname)
    }) %>%
    unlist() %>%
    as.data.frame() %>%
    data.table::fwrite(., file=glue('{opt$LD_folder}/chr{opt$chromosome}.SNPsLDBlock.txt'), col.names=F, row.names=F, quote=F)

# if(!is.null(opt$node)){
#     node <- as.numeric(opt$node)
#     data.table::fwrite(., file=glue('{opt$LD_folder}/chr{opt$chromosome}.SNPsLDBlock.txt.{node}'), col.names=F, row.names=F, quote=F)
# } else {
    
# }
    

