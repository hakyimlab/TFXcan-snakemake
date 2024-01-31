# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--sumstats_file", help='A list of files to combine'),
    make_option("--correlation_matrix", help='reference panel correlation matrix'),
    make_option("--zscores_column", help='A list of files to combine'),
    make_option("--n", help='summary statistics sample size', default = 200000L),
    make_option("--L", help='summary statistics sample size', default = 10L),
    make_option("--chromosome", help='A list of files to combine'),
    make_option("--LD_window", default=200000, help='A list of files to combine'),
    make_option("--output_file", help='A list of files to combine')
)

opt <- parse_args(OptionParser(option_list=option_list))  

library(data.table)
library(tidyverse)
library(susieR)
library(glue)

opt <- list()
opt$LDfile_pattern <- '/project2/haky/temi/projects/TFXcan-snakemake/data/LD/asthma_adults/asthma_adults.chr{}.ld'
opt$snplist_pattern <- '/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/asthma_adults/asthma_adults.chr{}.SNPsForLD'
opt$output_matrix <- '/project2/haky/Data/1000G/plink_ref_panel/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased'
opt$sumstats <- "/project2/haky/temi/projects/TFXcan-snakemake/data/sumstats/asthma_adults.gwas_sumstats.ALL.filtered.txt.gz"
opt$n <- 200000
opt$L <- 10

chrom_filter <- paste(c(1:22, 'X', 'Y'), sep='')

ld_files <- sapply(chrom_filter, function(f){
    ff <- gsub('\\{\\}', f, opt$snplist_pattern)
    if(file.exists(ff)){
        return(ff)
    }
}) 

ld_files <- Filter(Negate(is.null), ld_files)
valid_chrom <- names(ld_files)

sumstats <- data.table::fread(opt$sumstats) %>%
    dplyr::filter(chrom=='6')

for(f in chrom_filter[6]){
    ss <- gsub('\\{\\}', f, opt$snplist_pattern)
    ll <- gsub('\\{\\}', f, opt$LDfile_pattern)

    if(file.exists(ss) && file.exists(ll)){
        ss <- data.table::fread(ss, header=F) 
        ll <- data.table::fread(ll)

    }
}

# filter where the zscores are avaialable
sumstats_ld <- ll %>%
    dplyr::filter(SNP_A %in% sumstats$SNP & SNP_B %in% sumstats$SNP) %>%
    dplyr::select(SNP_A, SNP_B, R)

sumstats$SNP %in% ll$SNP_A

# https://stackoverflow.com/questions/57903948/creating-a-correlation-matrix-from-a-data-frame-in-r
sumstats_ld1 <- data.frame(SNP_A = sumstats_ld$SNP_B, SNP_B = sumstats_ld$SNP_A, R = sumstats_ld$R)
df <- rbind(sumstats_ld, sumstats_ld1)
df1 <- as.data.frame.matrix(xtabs(R ~ ., data=df)) |> as.matrix()
diag(df1) <- 1

# ensure rows and columns are appropriately ordered
R <- df1[sumstats$SNP, sumstats$SNP]

# ensure summary statistics rows match the correlation matrix
fitted_rss <- susieR::susie_rss(z=sumstats$zstat, n = opt$n, R = R, L = opt$L)
susie_plot(fitted_rss, y="PIP")


# data_file <- tempfile(fileext = ".RData")
# data_url <- paste0("https://raw.githubusercontent.com/stephenslab/susieR/",
#                    "master/inst/datafiles/SummaryConsistency1k.RData")

# curl_download(data_url,data_file)
# load(data_file)
# zflip = SummaryConsistency$z
# ld = SummaryConsistency$ldref