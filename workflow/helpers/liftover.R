
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--input_file", help='A list of files to combine'),
    make_option("--chain_file", help='reference panel correlation matrix'),
    make_option("--output_file", help='reference panel correlation matrix')
)

opt <- parse_args(OptionParser(option_list=option_list)) 

library(data.table)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

if(!dir.exists(dirname(opt$output_file))){
    dir.create(dirname(opt$output_file), recursive = T)
}

# load chain file
chain <- rtracklayer::import.chain(opt$chain_file)

# load input file
input <- data.table::fread(opt$input_file)

# convert to GRanges
gr <- with(input, GenomicRanges::GRanges(seqnames = chr, IRanges(pos, pos + 1), strand = '*', chrom, pos, alt, ref, variant_id, beta, se, pval, zscore))

# liftover
lo <- rtracklayer::liftOver(gr, chain)

# convert to dataframe and format
lf <- as.data.frame(lo) %>%
    dplyr::select(-c(group, group_name, width, strand, pos, seqnames)) %>%
    dplyr::select(chrom, pos = start, alt, ref, variant_id, beta, se, pval, zscore)

# save output
data.table::fwrite(lf, file=opt$output_file, sep='\t', quote=F, row.names=F, compress = 'gzip')