
# Description:
# Author: 
# Date: 











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


#sum(sumstats$SNP %in% genotypes_chr$varID)

# ld_window <- split_ld_blocks[[64]]
# summary_stat <- matched_stats
# susie_fit <- fitted_rss
# threshold <- opt$pip_threshold



# ?bigsnpr::snp_match  # strand information is crucial in bioinformatic analysis
# matched.gwas <- as_tibble(bigsnpr::snp_match(sumstats = gwas[["1452"]], 
#                           info_snp = bigSNP[["1452"]]$map) %>% 
#                           dplyr::rename(og_index = `_NUM_ID_.ss`) %>% 
#                           dplyr::rename(bigSNP_index = `_NUM_ID_`))

# chr8_125398675_127646866