# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))


option_list <- list(
    make_option("--chromosome", help='Chromsome number e.g. 1, 2, 3, ..., 22'),
    make_option("--sumstats", help='A GWAS summary statistics file for the chromsome; should be a tsv file with columns: chrom, pos, ref, alt, pval, beta, se, zscore'),
    make_option("--n", help='summary statistics sample size', default = 1000000L),
    make_option("--L", help='summary statistics L parameter used by SuSie', default = 10L),
    make_option("--LDBlocks_info", help='A file for the LD blocks of where to run SuSie'),
    make_option("--pip_threshold", default=0.5, help='threshold to filter SNPs by PIP; default is 0.5', type='numeric'),
    make_option("--genotypes_dosages_pattern", default='/project2/haky/Data/1000G/population_data/EUR/bfiles/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.geno.txt.gz', help='a pattern to find genotype files for the chromosome'),
    make_option("--output_folder", help='The folder to put SuSie results in'),
    make_option("--phenotype", help='A GWAS phenotype'),
    make_option('--diagnostics_file', type='character', default=NULL, help='A file to write diagnostics to; default is NULL i.e no diagnostics file will be written')
)

opt <- parse_args(OptionParser(option_list=option_list))  

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
# opt$sumstats <- "/project/haky/users/temi/projects/TFXcan-snakemake/data/PCRISK_GWAS_2024-06-07/processed_sumstats/pc_risk/chr3.sumstats.txt.gz"
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

LD_block <- data.table::fread(opt$LDBlocks_info) %>%
    dplyr::rename(chrom=V1, start=V2, stop=V3) %>%
    dplyr::filter(chrom == paste0('chr', opt$chromosome)) %>%
    dplyr::mutate(split=seq_len(nrow(.)))

genotypes_file <- gsub('\\{\\}', opt$chromosome, opt$genotypes_dosages_pattern)
if(!file.exists(genotypes_file)){
    stop('ERROR - File does not exist')
} else {
    genotypes_chr <- data.table::fread(genotypes_file) %>% 
        data.table::setDT()
}

runSusiePerLDBlock <- function(ld_window, summary_stat, genotypes_chr, susie_options = list(n=1e6, L=10, pip_threshold = 0.5), conn=NULL){

    stopifnot(is.list(susie_options), is.data.frame(ld_window) | is.data.table(ld_window))

    ldfile_basename <- paste(ld_window$chrom, ld_window$start, ld_window$stop, sep='_')
    vqueries <- summary_stat %>%
        dplyr::filter(dplyr::between(pos, ld_window$start, ld_window$stop))

    if(!"YES" %in% vqueries$gwas_significant){
        return(list(loci=NULL, fit=NULL))
    } else {
        msg <- glue('INFO - GWAS significant variant(s) found in the LD block {ldfile_basename}')
        print(msg)
        if(!is.null(conn)){
            cat(msg, file = conn, sep = '\n')
        }

        # prepare the LD matrix from the genotypes
        X <- genotypes_chr %>%
            dplyr::filter(varID %in% vqueries$varID) %>%
            tibble::column_to_rownames('varID') %>%
            as.matrix() %>% 
            t() %>% 
            scale()
        Rmat <- cor(X)

        # match and arrange and prepare zscores
        xt <- vqueries %>%
            dplyr::filter(varID %in% row.names(Rmat)) %>%
            dplyr::arrange(match(varID, row.names(Rmat)))
        zscores <- xt %>% dplyr::pull(zscore)

        # set susie options
        


        # run susieR
        result <- tryCatch(
            withCallingHandlers({   
                    fitted_rss <- susieR::susie_rss(z=zscores, n = susie_options$n, R = Rmat, L=susie_options$L)
                    #chosen_loci <- fitted_rss$pip[fitted_rss$pip >= opt$pip_threshold]

                    chosen_loci <- getChosenLoci(fitted_rss, threshold = susie_options$pip_threshold, ldblock=ldfile_basename)
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

getChosenLoci <- function(susie_fit, threshold, ldblock=NULL){
    # take the summary of fit and summarize

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

# use bigsnpr to correct the alleles
print(glue("INFO - Using bigsnpr to correct and match summary statistics alleles with reference genotype alleles"))
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



# split into approx. independent loci
split_ld_blocks <- LD_block %>%
        base::split(., f=.$split)

# mstats <- sumstats %>% dplyr::rename(varID = SNP)

susieRun <- list()
for(i in seq_along(split_ld_blocks)){
    blockName <- paste(split_ld_blocks[[i]]$chrom, split_ld_blocks[[i]]$start, split_ld_blocks[[i]]$stop, sep='_')
    susieRun[[blockName]] <- runSusiePerLDBlock(split_ld_blocks[[i]], matched_stats, genotypes_chr = genotypes_chr, conn=diagfile)
}

# saveRDS(susieRun, file = file.path(glue('/project/haky/users/temi/projects/TFXcan-snakemake/misc/chr{opt$chromosome}.susie.RDS')))

# # chr8_23039544_24817205
# # chr8_1213245_2095102
# # chr3_85533081_87360582
# susieRun <- list()
# ld_window <- split_ld_blocks[[56]]
# mstats <- matched_stats
# blockName <- paste(ld_window$chrom, ld_window$start, ld_window$stop, sep='_')
# susieRun[[blockName]] <- runSusiePerLDBlock(ld_window, mstats, genotypes_chr = genotypes_chr,conn=diagfile)
# all(rapply(susieRun, length) == 0)

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

# check that not everything is NULL and write out files appropriately
isNullList <- all(rapply(susieRun, length) == 0)

if(isNullList == TRUE){
    stop(glue("INFO - All SuSie runs on chromosome {opt$chromosome} are empty; this usually means that SuSie could not finemap any variant on this chromosome."))
} else if(isNullList == FALSE) {

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

    # filter for top GWAS SNPs only if present
    fSNP <- as.data.frame(filteredSNPs)
    ySNP <- dplyr::filter(summary_stat, varID %in% fSNP$SNP)
    if("YES" %in% ySNP$gwas_significant){
        ySNP <- ySNP[ySNP$gwas_significant == 'YES', ]
    }

    ySNP %>%
        dplyr::mutate(locus = paste(chrom, pos, pos + 1, sep='_')) %>%
        dplyr::mutate(locus = dplyr::case_when(
            !grepl('chr', locus, fixed=TRUE) ~ paste('chr', locus, sep=''),
            .default = as.character(locus)
        )) %>%
        dplyr::pull(locus) %>%
        base::data.frame(.) %>%
        data.table::fwrite(as.data.frame(ft), file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.EnformerLoci.topSNPs.txt')), quote=F, row.names=F, sep = '\t', col.names=F)


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
            data.table::fwrite(as.data.frame(ft), file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.EnformerLoci.allSNPs.txt')), quote=F, row.names=F, sep = '\t', col.names=F)
        }
    }
    # write out everything
    saveRDS(susieRun, file = file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.susie.RDS')))
}