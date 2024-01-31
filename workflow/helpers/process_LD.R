# this runs one locus at a time

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--LD_folder", help='A list of files to combine'),
    make_option("--pattern", help='A list of files to combine'),
    make_option("--output_file", help='A list of files to combine')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(susieR)
library(glue)
library(Matrix)
library(bigmemory)
print(opt)

# opt <- list()
# opt$LD_folder <- '/project2/haky/Data/1000G/LD'
# opt$pattern <- 'ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased'

chromosomes <- c(1:22, 'X')

out <- purrr::map(.x=chromosomes, .f=function(x){
    ff <- gsub('\\{\\}', x, opt$pattern)
    if(file.exists(glue::glue('{opt$LD_folder}/{ff}.ld'))){
        print(glue::glue('{opt$LD_folder}/{ff}.ld'))
        dt <- data.table::fread(glue::glue('{opt$LD_folder}/{ff}.ld')) %>%
            dplyr::select(SNP_A, SNP_B, R) %>%
            dplyr::mutate(SNP_A = as.character(SNP_A),
                          SNP_B = as.character(SNP_B)) %>% 
            tidyr::spread(SNP_B, R) %>% 
            tibble::column_to_rownames("SNP_A") %>% 
            as.big.matrix(x=., type='double', backingfile=glue::glue('{opt$LD_folder}/{ff}.bigmatrix.bin'), descriptorfile=glue::glue('{opt$LD_folder}/{ff}.bigmatrix.desc'))
        return(dt)
        # dt1 <- data.frame(SNP_A = dt$SNP_B, SNP_B = dt$SNP_A, R = dt$R)
        # df <- rbind(dt, dt)
        # df1 <- as.data.frame.matrix(xtabs(R ~ ., data=df))
        # diag(df1) <- 1
    }
})

# data <- do.call('rbind', out)

# print(str(data))

# # https://stackoverflow.com/questions/57903948/creating-a-correlation-matrix-from-a-data-frame-in-r
# data1 <- data.frame(SNP_A = data$SNP_B, SNP_B = data$SNP_A, R = data$R)
# df <- rbind(data, data1)
# df1 <- as.data.frame.matrix(xtabs(R ~ ., data=df))
# diag(df1) <- 1

saveRDS(out, file = glue('{opt$output_file}'))