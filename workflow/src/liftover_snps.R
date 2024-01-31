# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--sumstats", help='A list of files to combine'),
    make_option("--chainfile", help='reference panel correlation matrix'),
    make_option("--zscores_column", help='A list of files to combine'),
    make_option("--n", help='summary statistics sample size', default = 200000L),
    make_option("--L", help='summary statistics sample size', default = 10L),
    make_option("--chromosome", help='A list of files to combine'),
    make_option("--LD_window", default=200000, help='A list of files to combine'),
    make_option("--output_file", help='A list of files to combine')
)

opt <- parse_args(OptionParser(option_list=option_list)) 

# install and load rtracklayer
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

if (!requireNamespace("rtracklayer", quietly = TRUE)){
    BiocManager::install("rtracklayer", version = "3.8")
}

library(rtracklayer)

opt <- list()
opt$chromosome <- '1'
opt$chainfile <- '/project2/haky/Data/liftover/ensembl_chainfiles/GRCh37_to_GRCh38_tabs.chain'
opt$sumstats <- "/project2/haky/temi/projects/TFPred/data/asthma/asthma_children.logistic.assoc.tsv.gz" #"/project2/haky/temi/projects/TFXcan-snakemake/data/sumstats/prostate_cancer_risk.gwas_sumstats.ALL.filtered.txt.gz"

sumstats <- data.table::fread(opt$sumstats) %>%
    #dplyr::filter(chr == opt$chromosome) %>%
    dplyr::mutate(chr=as.character(gsub('chr', '', chr)), SNP=paste(chr, pos, alt, ref, sep='_')) %>%
    dplyr::rename(zscore=zstat, chrom=chr)


# specify coordinates to liftover
grObject <- GenomicRanges::GRanges(seqnames=sumstats$chr, ranges=IRanges(start=sumstats$pos, end=sumstats$pos), variant_id=sumstats$SNP)

# import the chain file
chainObject <- rtracklayer::import.chain(opt$chainfile)
# run liftOver
results <- rtracklayer::liftOver(grObject, chainObject) %>% 
    as.data.frame() %>%
    dplyr::select(seqnames, hg38_pos=start, variant_id) 

head(results)

# join and save
hg38_sumstats <- sumstats %>%
    dplyr::inner_join(results, by=c('SNP' = 'variant_id')) %>%
    dplyr::mutate(SNP=paste(chrom, hg38_pos, alt, ref, sep='_')) %>%
    dplyr::select(chrom, pos=hg38_pos, rsid, ref, alt, beta, se, zscore, pval, SNP)


data.table::fwrite(hg38_sumstats, file='/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/asthma_children.liftover.logistic.assoc.tsv.gz', compress='gzip', row.names=F, quote=F)