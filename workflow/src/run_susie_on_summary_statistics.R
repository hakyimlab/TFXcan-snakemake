# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--chromosome", help='A list of files to combine'),
    make_option("--sumstats", help='reference panel correlation matrix'),
    make_option("--n", help='summary statistics sample size', default = 1000000L),
    make_option("--L", help='summary statistics sample size', default = 10L),
    make_option("--LDBlocks_info", help='A list of files to combine'),
    make_option("--pip_threshold", default=0.5, help='threshold to filter SNPs by PIP', type='numeric'),
    make_option("--genotypes_dosages_pattern", default='/project2/haky/Data/1000G/population_data/EUR/bfiles/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.geno.txt.gz', help=''),
    make_option("--output_folder", help='threshold to filter SNPs by PIP'),
    make_option("--phenotype", help='threshold to filter SNPs by PIP'),
    make_option('--diagnostics_file', type='character', default=NULL, help='')
)

opt <- parse_args(OptionParser(option_list=option_list))  

print(opt)

library(data.table) |> suppressPackageStartupMessages()
library(tidyverse) |> suppressPackageStartupMessages()
library(susieR) |> suppressPackageStartupMessages()
library(glue) |> suppressPackageStartupMessages()
library(bigsnpr) |> suppressPackageStartupMessages()

if(!dir.exists(opt$output_folder)){
    dir.create(opt$output_folder)
}

# opt <- list()
# opt$chromosome <- '3'
# opt$sumstats <- "/project/haky/users/temi/projects/TFXcan-snakemake/data/processed_sumstats/prostate_cancer_risk/chr3.sumstats.txt.gz"
# #"/project2/haky/temi/projects/TFPred/data/asthma/asthma_children.logistic.assoc.tsv.gz" #"/project2/haky/temi/projects/TFXcan-snakemake/data/sumstats/prostate_cancer_risk.gwas_sumstats.ALL.filtered.txt.gz"
# opt$LDMatrix_folder <- "/project/haky/data/1000G/LD/LD_matrices/EUR"
# #opt$LD_window <- 200000
# opt$LDBlocks_info <- '/project/haky/users/temi/projects/TFXcan-snakemake/metadata/hg38_fourier_ls-all.bed'
# opt$n <- 1000000L
# opt$L <- 10L
# opt$pip_threshold <- 0.5
# opt$genotypes_dosages_pattern <- '/project/haky/data/1000G/population_data/EUR/bfiles/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.geno.txt.gz'
# opt$diagnostics_file <- NULL

# LDMat_chr <- file.path(opt$LDMatrix_folder, paste0('chr', opt$chromosome))
sumstats <- data.table::fread(opt$sumstats) %>%
    dplyr::filter(chrom == opt$chromosome) %>%
    dplyr::rename(a0 = ref, a1=alt, chr=chrom) %>%
    dplyr::mutate(chr = as.character(chr))
    #dplyr::mutate(chrom=as.character(paste('chr', chrom, sep='')))
    #dplyr::filter(chr == opt$chromosome) %>%
    # dplyr::mutate(chr=as.character(gsub('chr', '', chr))) %>% #, SNP=paste(chr, pos, alt, ref, sep='_')) %>%
    # dplyr::rename(zscore=zstat, chrom=chr)

LD_block <- data.table::fread(opt$LDBlocks_info) %>%
    dplyr::rename(chrom=V1, start=V2, stop=V3) %>%
    dplyr::filter(chrom == paste0('chr', opt$chromosome)) %>%
    dplyr::mutate(split=1:nrow(.))

genotypes_file <- gsub('\\{\\}', opt$chromosome, opt$genotypes_dosages_pattern)
if(!file.exists(genotypes_file)){
    stop('ERROR - File does not exist')
} else {
    genotypes_chr <- data.table::fread(genotypes_file) %>% 
        data.table::setDT()
}



# use bigsnpr to correct the alleles
gdt <- genotypes_chr %>%
    dplyr::select(varID) %>%
    tidyr::separate(varID, into=c('chr', 'pos', 'a0', 'a1'), remove=F) %>%
    dplyr::mutate(pos = as.integer(pos))

matched_stats <- bigsnpr::snp_match(sumstats, gdt) |> data.table::setDT()
matched_stats <- matched_stats %>%
    dplyr::select(chrom=chr, pos, a0, a1, rsid, beta, se, zscore, pval, adjP, gwas_significant, varID)

# matched_stats <- matched_stats %>%
#     dplyr::select(chrom=chr, pos, a0, a1, beta, rsid, zscore, varID)

data.table::fwrite(matched_stats, file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.bigsnpr.txt.gz')), compress='gzip', quote=F, row.names=F, sep = '\t')


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

#sum(sumstats$SNP %in% genotypes_chr$varID)

# ld_window <- split_ld_blocks[[64]]
# summary_stat <- matched_stats
# susie_fit <- fitted_rss
# threshold <- opt$pip_threshold

getChosenLoci <- function(susie_fit, threshold, ldblock=NULL){
    # take the summary of fit and summarize

    summary_dt <- base::summary(susie_fit)[['vars']] %>% 
        dplyr::group_by(cs) %>%
        dplyr::summarise(sum_pip = sum(variable_prob, na.rm=TRUE)) %>%
        dplyr::filter(cs >= 0.95)

    which_cs <- summary_dt$cs

    which_vars <- base::summary(susie_fit)[['vars']] %>%
        dplyr::group_by(cs) %>%
        dplyr::filter(variable_prob > threshold) %>%
        dplyr::pull(variable) %>%
        as.numeric()

    if(length(which_vars > 0)){
        # ensure the credible sets
        cs_info <- base::summary(susie_fit)[['vars']] %>%
            as.data.frame() %>%
            dplyr::filter(variable %in% which_vars)

        vars_info <- susie_fit$pip[which_vars] %>%
            as.data.frame() %>%
            tibble::rownames_to_column('SNP') %>%
            setNames(., c('SNP', 'PIP'))

        dt <- cbind(cs_info, vars_info)
        dt$ldBlock <- ldblock
        return(dt)

    } else {
        return(NULL)
    }
    
}

runSusiePerLDBlock <- function(ld_window, summary_stat, conn=NULL){
    ldfile_basename <- paste(ld_window$chrom, ld_window$start, ld_window$stop, sep='_')
    vqueries <- summary_stat %>%
        dplyr::filter(dplyr::between(pos, ld_window$start, ld_window$stop))

    #print(vqueries[1:5, ])

    if(!"YES" %in% vqueries$gwas_significant){
        #msg <- glue('INFO - No GWAS significant variants in the LD block {ldfile_basename}')
        # print(msg)
        # if(!is.null(conn)){
        #     cat(msg, file = conn, sep = '\n')
        # }
        return(list(loci=NULL, fit=NULL))
    } else {
        msg <- glue('INFO - GWAS significant variant(s) found in the LD block {ldfile_basename}')
        print(msg)
        if(!is.null(conn)){
            cat(msg, file = conn, sep = '\n')
        }

        X <- genotypes_chr %>%
            dplyr::filter(varID %in% vqueries$varID) %>%
            tibble::column_to_rownames('varID') %>%
            as.matrix() %>% 
            t() %>% 
            scale()
        Rmat <- cor(X)

        xt <- vqueries %>%
            dplyr::filter(varID %in% row.names(Rmat)) %>%
            dplyr::arrange(match(varID, row.names(Rmat)))
        zscores <- xt %>% dplyr::pull(zscore)

        result <- tryCatch(
            withCallingHandlers({   
                    fitted_rss <- susieR::susie_rss(z=zscores, n = opt$n, R = Rmat, L=opt$L)
                    #chosen_loci <- fitted_rss$pip[fitted_rss$pip >= opt$pip_threshold]

                    chosen_loci <- getChosenLoci(fitted_rss, threshold = opt$pip_threshold, ldblock=ldfile_basename)
                    if(is.null(chosen_loci)){
                        msg <- glue('INFO - No Susie significant variant(s) found in the LD block {ldfile_basename}')
                        print(msg)
                        if(!is.null(conn)){
                            cat(msg, file = conn, sep = '\n')
                        }
                    }
                    
                    #print(head(chosen_loci))
                    return(list(loci=chosen_loci, fit=fitted_rss))
                },
                warning = function(w){
                    message_txt <<- trimws(paste0("WARNING: ", w))
                    invokeRestart("muffleWarning")
                }),
            error=function(e){
                print(glue('ERROR - Susie found errors'))
                return(list(loci=NULL, fit=NULL))
            }
        )
        return(result)
    }
}


# ?bigsnpr::snp_match  # strand information is crucial in bioinformatic analysis
# matched.gwas <- as_tibble(bigsnpr::snp_match(sumstats = gwas[["1452"]], 
#                           info_snp = bigSNP[["1452"]]$map) %>% 
#                           dplyr::rename(og_index = `_NUM_ID_.ss`) %>% 
#                           dplyr::rename(bigSNP_index = `_NUM_ID_`))

# chr8_125398675_127646866

split_ld_blocks <- LD_block %>%
        base::split(., f=.$split)

# mstats <- sumstats %>% dplyr::rename(varID = SNP)

susieRun <- list()
for(i in 1:length(split_ld_blocks)){
    blockName <- paste(split_ld_blocks[[i]]$chrom, split_ld_blocks[[i]]$start, split_ld_blocks[[i]]$stop, sep='_')
    susieRun[[blockName]] <- runSusiePerLDBlock(split_ld_blocks[[i]], matched_stats, conn=diagfile)
}

# saveRDS(susieRun, file = file.path(glue('/project/haky/users/temi/projects/TFXcan-snakemake/misc/chr{opt$chromosome}.susie.RDS')))

# # chr8_23039544_24817205
# # chr8_1213245_2095102
# # chr3_85533081_87360582

# ld_window <- split_ld_blocks[[56]]
# summary_stat <- matched_stats
# rt <- runSusiePerLDBlock(ld_window, mstats, conn=diagfile)

# susie_plot(susieRun$chr8_125398675_127646866$fit, 'PIP')

# png('/project/haky/users/temi/projects/TFXcan-snakemake/plt.png')
# susie_plot(rt$fit, 'PIP')
# plot(1:50)
# dev.off()

# close the connection if necessary
if(!is.null(opt$diagnostics_file)){
    close(diagfile)
}


# t1 <- list(a=list(a=1, b=2), c=list(a=NULL, b=2), d=list(a=NULL, b=2))
# t2 <- list(a=list(a=1, b=NULL), c=list(a=NULL, b=NULL), d=list(a=NULL, b=NULL))
# all(rapply(t2, length) == 0)

# check that not everything is NULL
isNullList <- all(rapply(susieRun, length) == 0)

if(isNullList == FALSE){
    # filter lCausalSnps 
    filteredSNPs <- lapply(susieRun, function(x){
        if(is.null(x$loci)){
            return(NULL)
        }
        dt <- as.data.frame(x$loci)
        if(nrow(dt) == 0){
            return(NULL)
        }
        return(dt)
    }) 
    filteredSNPs <- Filter(Negate(is.null), filteredSNPs) %>%
        do.call(rbind, .)
    row.names(filteredSNPs) <- NULL

    tryCatch({
        if(!is.null(filteredSNPs)){
            data.table::fwrite(filteredSNPs, file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.finemapped.txt.gz')), compress='gzip', quote=F, row.names=F, sep = '\t')
        }
    }, error=function(e){
        print(glue('ERROR - Cannot save the results from chr{opt$chromosome}'))
    })
    
    

    #isNullList <- function(x){all(!lengths(x))}
    # filter the convergence
    convergence <- sapply(susieRun, function(x){
        if(is.null(x$fit)){
            return(NULL)
        }
        return(x$fit$converged)
    }) %>% Filter(Negate(is.null), .) %>%
        as.data.frame() %>% t() %>% as.data.frame() %>%
        tibble::rownames_to_column('ldBlock') %>%
        dplyr::rename(convergence=2)
    data.table::fwrite(convergence, file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.convergence.txt.gz')), compress='gzip', quote=F, row.names=F, sep = '\t')

    if(!is.null(filteredSNPs)){
        filteredGWAS <- matched_stats %>%
            dplyr::filter(varID %in% filteredSNPs$SNP)
        
        if(!nrow(filteredGWAS) == 0){
            data.table::fwrite(filteredGWAS, file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.filteredGWAS.txt.gz')), compress='gzip', quote=F, row.names=F, sep = '\t')

            ft <- filteredGWAS %>%
                dplyr::mutate(locus = paste(chrom, pos, pos + 1, sep='_')) %>%
                dplyr::mutate(locus = dplyr::case_when(
                    !grepl('chr', locus, fixed=TRUE) ~ paste('chr', locus, sep=''),
                    .default = as.character(locus)
                )) %>%
                dplyr::pull(locus) %>%
                base::data.frame(.)
            data.table::fwrite(as.data.frame(ft), file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.EnformerLoci.txt')), quote=F, row.names=F, sep = '\t', col.names=F)
        }
    }
    # write out everything
    saveRDS(susieRun, file = file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.susie.RDS')))
}








# getChosenLoci <- function(susie_fit, threshold){
#     # take the summary of fit and summarize

#     summary_dt <- base::summary(susie_fit)[['vars']] %>% 
#         dplyr::group_by(cs) %>%
#         dplyr::summarise(sum_pip = sum(variable_prob, na.rm=TRUE)) %>%
#         dplyr::filter(cs >= 1)

#     which_cs <- summary_dt %>%
#         dplyr::filter(sum_pip >= threshold) %>%
#         dplyr::pull(cs)

#     which_vars <- summary(susie_fit)[["cs"]] %>%
#         as.data.frame() %>%
#         dplyr::filter(cs %in% which_cs) %>%
#         dplyr::pull(variable) %>%
#         base::strsplit(., split=',') %>%
#         unlist() %>% as.numeric()

#     vars_dt <- susie_fit$pip[which_vars] %>%
#         as.data.frame() %>%
#         tibble::rownames_to_column('SNP') %>%
#         setNames(., c('SNP', 'PIP'))

#     vars_dt <- cbind(which_vars, vars_dt)

#     return(vars_dt)
# }


# runSusiePerLDBlock(split_ld_blocks[[33]], sumstats, conn=diagfile)
# susieRun[["chr12_55272053_57155077"]]$loci

# gg <- susieRun[["chr12_55272053_57155077"]]$fit

# getChosenLoci(fitted_rss, 0.5)



# # take the summary of fit and summarize
# summary_dt <- base::summary(gg)[['vars']] %>% 
#     dplyr::group_by(cs) %>%
#     dplyr::summarise(sum_pip = sum(variable_prob)) %>%
#     dplyr::filter(cs >= 1)

# which_cs <- summary_dt %>%
#     dplyr::filter(sum_pip >= opt$pip_threshold) %>%
#     dplyr::pull(cs)

# which_vars <- summary(gg)[["cs"]] %>%
#     as.data.frame() %>%
#     dplyr::filter(cs %in% which_cs) %>%
#     dplyr::pull(variable) %>%
#     base::strsplit(., split=',') %>%
#     setNames(., which_cs) %>%
#     purrr::map_df(., ~ as.data.frame(.x), .id="id") %>%
#     setNames(., c('cs', 'var_id')) %>%
#     dplyr::mutate(var_id = as.numeric(var_id))

# vars_dt <- gg$pip[which_vars$var_id] %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column('SNP') %>%
#     setNames(., c('SNP', 'PIP'))

# vars_dt <- cbind(which_vars, vars_dt)

# susie_plot(gg, y="PIP")

# tryCatch(
#       withCallingHandlers(
#         {
#           error_text <- "No error."
#           list(value = hurz(x), error_text = error_text)
#         }, 
#         warning = function(e) {
#           error_text <<- trimws(paste0("WARNING: ", e))
#           invokeRestart("muffleWarning")
#         }
#       )

# nums <- c(1, 0, -1, 2, 3)

# for(i in nums){
#     result <- tryCatch(
#         withCallingHandlers({
#                 lognums <- log(i)
#                 print(lognums)
#                 message_txt <- 'success'
#                 return(list(result=lognums, message=message_txt))
#             },
#             warning = function(w){
#                 message_txt <<- trimws(paste0("WARNING: ", w))
#                 #invokeRestart("muffleWarning")
#             }),
#         error=function(e){
#             print(glue('ERROR - Susie found errors'))
#             return(list(result=NULL, message='Errors'))
#         }
#     )
#     print(result)
# }



#     # sumstats %>% dplyr::filter(rsid == 'rs28407950')
# ld_window <- LD_block %>% dplyr::filter(start <= 32658571 & stop >= 32658571)
# vqueries <- sumstats %>%
#     dplyr::filter(dplyr::between(pos, ld_window$start, ld_window$stop)) %>%
#     dplyr::pull(SNP)

# variants_query <- paste0(vqueries, collapse = '\\|')
# variants_query <- paste0("'", variants_query, "'")
# variants_cmd <- paste0("zgrep ", variants_query, glue(" {genotypes_file}"))

#     dplyr::filter(varID %in% vqueries) %>%
#     tibble::column_to_rownames('varID') %>%
#     as.matrix() %>% 
#     t() %>% 
#     scale()
# Rmat <- cor(X)

# zscores <- sumstats %>%
#     dplyr::filter(SNP %in% row.names(Rmat)) %>%
#     dplyr::pull(zscore)

# # jj <- susieR::kriging_rss(z=zscores, R = Rmat, n=opt$n)
# # jj$plot

# fitted_rss <- susieR::susie_rss(z=zscores, n = opt$n, R = Rmat, L=opt$L)
# susie_plot(fitted_rss, y="PIP")
# susie_summary <- summary(fitted_rss)
# run_susie(split_ld_blocks[[2]], sumstats)



    # setNames(c('id', individuals)) %>% 
    #  
    


# run_susie <- function(ld_window, summary_stat){
#     ldfile_basename <- paste(ld_window$chrom, ld_window$start, ld_window$stop, sep='_')
#     print(ldfile_basename)

#     # filter for SNPs within the LD block
#     sumstats_dt <- summary_stat %>%
#         dplyr::filter(chrom == ld_window$chrom) %>%
#         dplyr::filter(dplyr::between(pos, ld_window$start, ld_window$stop))
    
#     # if there are no SNPs, return NULL else, continue
#     if(nrow(sumstats_dt) == 0){
#         print(glue('INFO - No variants in the GWAS is in the LD block'))
#         return(NULL) # nolint
#     }

#     # read and process the LD matrix corresponding to this region
#     ldmat <- file.path(opt$LDMatrix_folder, ld_window$chrom, paste0(ldfile_basename, '.LDMatrix.ld')) |> 
#         data.table::fread() %>% 
#         dplyr::filter(SNP_A %in% sumstats_dt$SNP & SNP_B %in% sumstats_dt$SNP) %>%
#         dplyr::select(SNP_A, SNP_B, R) # select only the snps and corr. coefficient

#     df <- rbind(ldmat, data.frame(SNP_A = ldmat$SNP_B, SNP_B = ldmat$SNP_A, R = ldmat$R))
#     Rmat <- as.data.frame.matrix(xtabs(R ~ ., data=df)) |> as.matrix() # this creates a matrix of corr. coefficients
#     diag(Rmat) <- 1

#     if(!identical(rownames(Rmat), colnames(Rmat))){
#         print(glue('ERROR - Dimensions in R matrix do not match'))
#         return(NULL)
#         if(!identical(sumstats_dt$SNP, colnames(Rmat))){
#             print(glue('ERROR - Dimensions in summary statistics and R matrix do not match'))
#             return(NULL)
#         }
#     }
#     #Rmat <- df1[valid_snps$SNP, valid_snps$SNP]
#     dim_cond <- length(unique(c(nrow(sumstats_dt), nrow(Rmat), ncol(Rmat)))) == 1
#     if( dim_cond ){
#         run <- tryCatch({
#             fitted_rss <- susieR::susie_rss(z=sumstats_dt$zscore, n = opt$n, R = Rmat, L = opt$L)
#             if(is.null(dim(summary(fitted_rss)))){
#                 print(glue('INFO - Susie found nothing'))
#                 return(NULL)
#             } else {
#                 return(summary(fitted_rss))
#             }
#             }, 
#             error=function(e){
#                 print(glue('ERROR - Susie found errors'))
#                 return(NULL)
#             }
#         )
#         return(run)
#     }
# }

# split_ld_blocks <- LD_block %>%
#         base::split(., f=.$split)

# susieRun <- list()
# for(i in 1:length(split_ld_blocks)){
#     susieRun[[i]] <- run_susie(split_ld_blocks[[i]], sumstats)
# }

# print(glue('INFO - saving susieR results'))
# saveRDS(susieRun, file=glue("/project2/haky/temi/projects/TFXcan-snakemake/data/susie_result/chr{opt$chromosome}_asthma_children.susie.RData"))


# sumstats %>% dplyr::filter(rsid == 'rs28407950')
# ld_window <- LD_block %>% dplyr::filter(start <= 32658571 & stop >= 32658571)
# ldfile_basename <- paste(ld_window$chrom, ld_window$start, ld_window$stop, sep='_')

# ldmat_dt <- file.path(opt$LDMatrix_folder, ld_window$chrom, paste0(ldfile_basename, '.LDMatrix.ld')) |> 
#     data.table::fread() %>% dplyr::select(SNP_A, SNP_B, R)

# ss_dt <- sumstats %>% dplyr::filter(dplyr::between(pos, ld_window$start, ld_window$stop)) %>%
#     dplyr::mutate(SNP=paste0('chr', SNP))

# ldmat_filt <- ldmat_dt %>%
#     dplyr::filter(SNP_A %in% ss_dt$SNP & SNP_B %in% ss_dt$SNP)
# ldmat_dt1 <- data.frame(SNP_A = ldmat_filt$SNP_B, SNP_B = ldmat_filt$SNP_A, R = ldmat_filt$R)
# df <- rbind(ldmat_filt, ldmat_dt1)
# df1 <- as.data.frame.matrix(xtabs(R ~ ., data=df)) |> as.matrix() # this creates a matrix of corr. coefficients
# diag(df1) <- 1

# ss_dt_use <- ss_dt %>%
#     dplyr::filter(SNP %in% rownames(df1))

# identical(rownames(df1), colnames(df1))
# identical(ss_dt_use$SNP, colnames(df1))

# rmat <- df1[,]

# jj <- susieR::kriging_rss(z=ss_dt_use$zscore, R = rmat, n=opt$n)

# fitted_rss <- susieR::susie_rss(z=ss_dt_use$zscore, n = opt$n, R = rmat)
# susie_plot(fitted_rss, y="PIP")

# summary(fitted_rss)
# run_susie(split_ld_blocks[[2]], sumstats)

# #  Q <- cov(mData)
# Q <- cov(rmat)
# eigen(Q)$values |> sort() |> plot()

# Q <- cov(rmat) + diag(1e-6, nrow(rmat))
# eigen(Q)$values |> sort() |> plot()

# Qadj <- Matrix::nearPD(cov(rmat))$mat

# ff <- susieR::kriging_rss(z=ss_dt_use$zscore, R = Qadj, n=opt$n)

# fitted_rss <- susieR::susie_rss(z=ss_dt_use$zscore, n = opt$n, R = as.matrix(Qadj))
# susie_plot(fitted_rss, y="PIP")

# all(rmat[lower.tri(rmat)] == 1, rmat[upper.tri(rmat)] == 1)

# diag(rmat)[diag(rmat) != 1]


# out <- LD_block[1:2] %>%
#         base::split(., f=.$split) %>%
#         purrr::map(.f=function(each_row){
#             ldfile_basename <- paste(each_row$chrom, each_row$start, each_row$stop, sep='_')
#             print(ldfile_basename)
#             sumstats_dt <- sumstats %>%
#                 dplyr::filter(dplyr::between(pos, each_row$start, each_row$stop))

#             if(nrow(sumstats_dt) == 0){
#                 return(NULL)
#             } else {
#                 # read and process the LD matrix
#                 ldmat_dt <- file.path(opt$LDMatrix_folder, each_row$chrom, paste0(ldfile_basename, '.LDMatrix.ld')) |> 
#                     data.table::fread() %>%
#                     dplyr::select(SNP_A, SNP_B, R)

#                 ldmat_dt1 <- data.frame(SNP_A = ldmat_dt$SNP_B, SNP_B = ldmat_dt$SNP_A, R = ldmat_dt$R)
#                 df <- rbind(ldmat_dt, ldmat_dt1)
#                 df1 <- as.data.frame.matrix(xtabs(R ~ ., data=df)) |> as.matrix()
#                 diag(df1) <- 1

#                 # select and processthe SNPs in the LD matrix 
#                 zscores_dt <- data.frame(SNP=row.names(df1))
#                 valid_snps <- sumstats_dt %>%
#                     dplyr::select(SNP, zscore)
#                 invalid_snps <- zscores_dt %>%
#                     dplyr::anti_join(valid_snps, by='SNP') %>%
#                     dplyr::mutate(zscore=0)

#                 # print(head(valid_snps))
#                 # print(head(invalid_snps))
            
#                 dt <- rbind(valid_snps, invalid_snps)
#                 df1 <- df1[dt$SNP, dt$SNP]

#                 run <- tryCatch({
#                     fitted_rss <- susieR::susie_rss(z=dt$zscore, n = opt$n, R = df1, L = opt$L)
#                     if(is.null(dim(summary(fitted_rss)))){
#                         return(NULL)
#                     } else {
#                         return(summary(fitted_rss))
#                     }
#                 }, 
#                 error=function(e){
#                     # print(length(dt$SNP))
#                     # print(dim(df1))
#                     # print(glue('ERROR - {ldfile_basename}'))
#                     # print(dt)
#                     # print(df1[1:5, 1:5])
#                     return(NULL)
                    
#                 })
                
#                 return(run)
#             }
#         })

# out <- out[!purrr::map_lgl(out, is.null)]
# summ_out <- lapply(out, summary)

# susie_plot(out[[4]], 'PIP')
# summary(out[[1]])
# coef(out[[1]])

# purrr::map()

# block1 <- LD_block[3, ]

# sumstats_window <- sumstats %>%
#     dplyr::filter(
#             dplyr::between(pos, left=block1$start, right=block1$stop)
#     )

# ff <- gsub('\\{\\}', opt$chromosome, opt$correlation_matrix)
# if(file.exists(ff)){
#     LD <- data.table::fread(ff) 
# } else {
#     stop(glue::glue('ERROR - {ff} does not exist'))
# }

# LD_window <- LD %>%
#     dplyr::filter(SNP_A %in% sumstats_window$SNP |
#         SNP_B %in% sumstats_window$SNP) %>%
#     dplyr::select(SNP_A, SNP_B, R)


# LD %>%
#     dplyr::filter(BP_A == 2291684)

# ll <- sumstats[1, ]

# pad_window <- ceiling(opt$LD_window/2)


# sumstats_window <- sumstats %>%
#     dplyr::filter(
#             dplyr::between(pos, left=ll$pos - pad_window, right=ll$pos + pad_window)
#     )

# LD %>%
#     dplyr::filter(SNP_A %in% sumstats_window$SNP)


# LD_window <- LD %>%
#     dplyr::filter(SNP_A %in% sumstats_window$SNP &
#         SNP_B %in% sumstats_window$SNP) %>%
#     dplyr::select(SNP_A, SNP_B, R)

# # https://stackoverflow.com/questions/57903948/creating-a-correlation-matrix-from-a-data-frame-in-r
# LD_window1 <- data.frame(SNP_A = LD_window$SNP_B, SNP_B = LD_window$SNP_A, R = LD_window$R)
# df <- rbind(LD_window, LD_window1)
# df1 <- as.data.frame.matrix(xtabs(R ~ ., data=df))
# diag(df1) <- 1











# dt <- data.table::fread(opt$processed_sumstats_file)
# R_mat <- readRDS(opt$correlation_matrix)

# # ensure summary statistics rows match the correlation matrix
# fitted_rss <- susieR::susie_rss(z=dt[, opt$zscores_column] n = opt$n, R = R_mat, L = opt$L)

# print(glue('INFO - saving susieR results'))
# saveRDS(fitted_rss, file = glue('{opt$output_file}'))