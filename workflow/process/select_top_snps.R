# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))


option_list <- list(
    make_option("--chromosome", help='Chromsome number e.g. 1, 2, 3, ..., 22'),
    make_option("--sumstats", help='A GWAS summary statistics file for the chromsome; should be a tsv file with columns: chrom, pos, ref, alt, pval, beta, se, zscore'),
    make_option("--LDBlocks_info", help='A file for the LD blocks of where to run SuSie'),
    make_option("--output_folder", help='The folder to put SuSie results in'),
    make_option("--phenotype", help='A GWAS phenotype'),
    make_option('--diagnostics_file', type='character', default=NULL, help='A file to write diagnostics to; default is NULL i.e no diagnostics file will be written')
)

opt <- parse_args(OptionParser(option_list=option_list))  
print(opt)

library(data.table) |> suppressPackageStartupMessages()
library(tidyverse) |> suppressPackageStartupMessages()
library(glue) |> suppressPackageStartupMessages()

# opt <- list()
# opt$chromosome <- '3'
# opt$sumstats <- "/project/haky/users/temi/projects/TFXcan-snakemake/data/PCRISK_GWAS_2024-06-07/processed_sumstats/pc_risk/chr3.sumstats.txt.gz"
# #"/project2/haky/temi/projects/TFPred/data/asthma/asthma_children.logistic.assoc.tsv.gz" #"/project2/haky/temi/projects/TFXcan-snakemake/data/sumstats/prostate_cancer_risk.gwas_sumstats.ALL.filtered.txt.gz"
# #opt$LD_window <- 200000
# opt$LDBlocks_info <- '/project/haky/users/temi/projects/TFXcan-snakemake/metadata/EUR_LDBlocks.hg38.updated.bed'
# opt$diagnostics_file <- NULL
# opt$phenotype <- 'pc_risk'

if(!dir.exists(opt$output_folder)){
    dir.create(opt$output_folder, recursive = T)
}

sumstats <- data.table::fread(opt$sumstats) %>%
    dplyr::filter(chrom == opt$chromosome) %>%
    dplyr::rename(a0 = ref, a1=alt, chr=chrom) %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    dplyr::filter(gwas_significant == 'YES')

if(nrow(sumstats) == 0){
    message(glue("INFO - No significant GWAS hits for chromosome {opt$chromosome}"))
    # data.table::fwrite(glue("{opt$output_folder}/{opt$phenotype}.chr{opt$chromosome}.filteredGWAS.topSNPs.txt"), sep='\t', quote=FALSE, row.names=FALSE)
    # data.table::fwrite(file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.EnformerLoci.topSNPs.txt')), quote=F, row.names=F, sep = '\t', col.names=F)

    opt <- options(show.error.messages=FALSE)
    on.exit(options(opt))
    quit()
}

LD_block <- data.table::fread(opt$LDBlocks_info) %>%
    dplyr::rename(chrom=V1, start=V2, stop=V3) %>%
    dplyr::filter(chrom == paste0('chr', opt$chromosome)) %>%
    dplyr::mutate(split=seq_len(nrow(.)))

if(!is.null(opt$diagnostics_file)){
    if(!file.exists(dirname(opt$diagnostics_file))){
        if(!dir.exists(dirname(opt$diagnostics_file))){
            dir.create(dirname(opt$diagnostics_file), recursive = TRUE)
        }
    }

    diagfile <- file(opt$diagnostics_file, open = "a")
    cat("#### Susie diagnostics", file = diagfile, sep = '\n')
} else {
    diagfile <- NULL
}

findTopSNPPerLDBlock <- function(ld_window, summary_stats, diagnostics_file=NULL){
    # get the LD block
    ft <- sumstats %>%
        dplyr::filter(dplyr::between(pos, ld_window$start, ld_window$stop)) 
    if(nrow(ft) == 0){
        return(NULL)
    } else {
        return(ft[which.min(ft$pval),])
    }
}

# for each LD block, run find top SNP

LD_block_split <- base::split(LD_block, LD_block$split)

topsnps <- base::lapply(LD_block_split, findTopSNPPerLDBlock, sumstats, opt$diagnostics_file)

# close the connection if necessary
if(!is.null(opt$diagnostics_file)){
    close(diagfile)
}

topsnps <- dplyr::bind_rows(topsnps) %>%
    dplyr::mutate(phenotype = opt$phenotype) %>%
    dplyr::select(chr, pos, a0, a1, pval, beta, se, zscore, phenotype)
topsnps %>%
    data.table::fwrite(glue("{opt$output_folder}/{opt$phenotype}.chr{opt$chromosome}.filteredGWAS.topSNPs.txt"), sep='\t', quote=FALSE, row.names=FALSE)

topsnps %>%
    dplyr::mutate(locus = paste(chr, pos, pos + 1, sep='_')) %>%
    dplyr::mutate(locus = dplyr::case_when(
        !grepl('chr', locus, fixed=TRUE) ~ paste('chr', locus, sep=''),
        .default = as.character(locus)
    )) %>%
    dplyr::pull(locus) %>%
    base::data.frame(.) %>%
    data.table::fwrite(file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.EnformerLoci.topSNPs.txt')), quote=F, row.names=F, sep = '\t', col.names=F)

