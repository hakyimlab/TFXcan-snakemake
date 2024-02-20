# this runs one locus at a time
# change to calculate LD
# https://merrimanlab.github.io/docs/locuszooms/multiple_lz_plots/

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--chromosome", help='chromosome e.g. 1,2,3,4,5,21 e.t.c.'),
    make_option("--sumstats", help='GWAS summary statistics with columns: chrom, pos, rsid, ref, alt, beta, se, zscore, pval, adjP, gwas_significant column{NO, YES}'),
    make_option("--n", help='GWAS sample size if available. Else, 1 million is assumed.', default = 1000000L),
    make_option("--L", help='L parameter in susie', default = 10L),
    make_option("--LDBlocks_info", default='/project2/haky/Data/LD_blocks/hg38/EUR/hg38_fourier_ls-all.bed', help='LD block information'),
    make_option("--pip_threshold", default=0.5, help='threshold to filter SNPs by Credible sets and PIPs', type='numeric'),
    make_option("--genotypes_dosages_pattern", default='/project2/haky/Data/1000G/population_data/EUR/bfiles/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.geno.txt.gz', help='Genotype dosages'),
    make_option("--output_folder", help='An output folder for Susie runs'),
    make_option("--phenotype", help='GWAS phenotype'),
    make_option('--diagnostics_file', type='character', default=NULL, help='Name for a file that will be created and that will hold information about Susie diagnostics')
)

opt <- parse_args(OptionParser(option_list=option_list))  

library(data.table) |> suppressPackageStartupMessages()
library(tidyverse) |> suppressPackageStartupMessages()
library(susieR) |> suppressPackageStartupMessages()
library(glue) |> suppressPackageStartupMessages()
library(bigsnpr) |> suppressPackageStartupMessages()

# opt <- list()
# opt$chromosome <- '12'
# opt$sumstats <- "/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/asthma_children/chr12.sumstats.txt.gz"
# #"/project2/haky/temi/projects/TFPred/data/asthma/asthma_children.logistic.assoc.tsv.gz" #"/project2/haky/temi/projects/TFXcan-snakemake/data/sumstats/prostate_cancer_risk.gwas_sumstats.ALL.filtered.txt.gz"
# opt$LDBlocks_info <- '/project2/haky/Data/LD_blocks/hg38/EUR/hg38_fourier_ls-all.bed'
# opt$n <- 1000000L
# opt$L <- 10L
# opt$pip_threshold <- 0.5
# opt$genotypes_dosages_pattern <- '/project2/haky/Data/1000G/population_data/EUR/bfiles/ALL.chr{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.geno.txt.gz'
# opt$diagnostics_file <- NULL

# chrom   pos     rsid    ref     alt     beta    se      zscore  pval    SNP     adjP    gwas_significant
# 10      15394   rs183305313     G       A       0.165417        0.110305        1.49964 0.133708        10_15394_G_A    1       NO
# 10      23147   rs35418599      C       T       -0.0180976      0.0187856       -0.963376       0.335359        10_23147_C_T    1       NO
# 10      44187   rs185642176     C       T       0.0273441       0.0264028       1.03565 0.300365        10_44187_C_T    1       NO
# 10      48323   rs184120752     C       A       -0.0333634      0.0484602       -0.68847        0.491157        10_48323_C_A    1       NO
# 10      48486   rs10904045      C       T       -0.0252516      0.0152247       -1.6586 0.0971968       10_48486_C_T    1       NO
# 10      48598   rs189409193     C       T       -0.0173221      0.111979        -0.154691       0.877065        10_48598_C_T    1       NO
# 10      48601   rs11251906      C       A       -0.0986587      0.045683        -2.15964        0.0308007       10_48601_C_A    1       NO
# 10      48971   rs565286582     G       A       -0.253552       0.204812        -1.23798        0.215725        10_48971_G_A    1       NO
# 10      49134   rs6560828       G       A       -0.0206479      0.0149594       -1.38026        0.167505        10_49134_G_A    1       NO
# 10      49156   rs11251929      A       G       0.00446223      0.0155045       0.287802        0.773498        10_49156_A_G    1       NO
# 10      49413   rs4320883       G       A       -0.0134982      0.0170535       -0.791523       0.428639        10_49413_G_A    1       NO


sumstats <- data.table::fread(opt$sumstats) %>%
    dplyr::filter(chrom == opt$chromosome) %>%
    dplyr::rename(a0 = ref, a1=alt, chr=chrom) %>%
    dplyr::mutate(chr = as.character(chr))

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

if(!dir.exists(opt$output_folder)){
    dir.create(opt$output_folder)
}

# use bigsnpr to correct the alleles
gdt <- genotypes_chr %>%
    dplyr::select(varID) %>%
    tidyr::separate(varID, into=c('chr', 'pos', 'a0', 'a1'), remove=F) %>%
    dplyr::mutate(pos = as.integer(pos))

matched_stats <- bigsnpr::snp_match(sumstats, gdt) |> data.table::setDT()
matched_stats <- matched_stats %>%
    dplyr::select(chrom=chr, pos, a0, a1, rsid, beta, se, zscore, pval, adjP, gwas_significant, varID)

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

getChosenLoci <- function(susie_fit, threshold){
    # take the summary of fit and summarize

    summary_dt <- base::summary(susie_fit)[['vars']] %>% 
        dplyr::group_by(cs) %>%
        dplyr::summarise(sum_pip = sum(variable_prob, na.rm=TRUE)) %>%
        dplyr::filter(cs >= 1)

    which_cs <- summary_dt %>%
        dplyr::filter(sum_pip >= threshold) %>%
        dplyr::pull(cs)

    which_vars <- summary(susie_fit)[["cs"]] %>%
        as.data.frame() %>%
        dplyr::filter(cs %in% which_cs) %>%
        dplyr::pull(variable) %>%
        base::strsplit(., split=',') %>%
        unlist() %>% as.numeric()

    vars_dt <- susie_fit$pip[which_vars] %>%
        as.data.frame() %>%
        tibble::rownames_to_column('SNP') %>%
        setNames(., c('SNP', 'PIP'))

    vars_dt <- cbind(which_vars, vars_dt)

    return(vars_dt)
}

runSusiePerLDBlock <- function(ld_window, summary_stat, conn=NULL){
    ldfile_basename <- paste(ld_window$chrom, ld_window$start, ld_window$stop, sep='_')
    vqueries <- summary_stat %>%
        dplyr::filter(dplyr::between(pos, ld_window$start, ld_window$stop))

    #print(vqueries[1:5, ])

    if(!"YES" %in% vqueries$gwas_significant){
        msg <- glue('INFO - No GWAS significant variants in the LD block {ldfile_basename}')
        print(msg)
        if(!is.null(conn)){
            cat(msg, file = conn, sep = '\n')
        }
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

        zscores <- vqueries %>%
            dplyr::filter(varID %in% row.names(Rmat)) %>%
            dplyr::pull(zscore)

        result <- tryCatch(
            withCallingHandlers({   
                    fitted_rss <- susieR::susie_rss(z=zscores, n = opt$n, R = Rmat, L=opt$L)
                    #chosen_loci <- fitted_rss$pip[fitted_rss$pip >= opt$pip_threshold]

                    chosen_loci <- getChosenLoci(fitted_rss, threshold = opt$pip_threshold)
                    if(!is.null(chosen_loci)){
                        chosen_loci$ldBlock <- ldfile_basename
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

split_ld_blocks <- LD_block %>%
        base::split(., f=.$split)

susieRun <- list()
for(i in 1:length(split_ld_blocks)){
    blockName <- paste(split_ld_blocks[[i]]$chrom, split_ld_blocks[[i]]$start, split_ld_blocks[[i]]$stop, sep='_')
    susieRun[[blockName]] <- runSusiePerLDBlock(split_ld_blocks[[i]], matched_stats, conn=diagfile)
}

# close the connection if necessary
if(!is.null(opt$diagnostics_file)){
    close(diagfile)
}


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

    # write out everything
    saveRDS(susieRun, file = file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.susie.RDS')))
    data.table::fwrite(filteredSNPs, file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.finemapped.txt.gz')), compress='gzip', quote=F, row.names=F, sep = '\t')
    data.table::fwrite(convergence, file=file.path(opt$output_folder, glue('{opt$phenotype}.chr{opt$chromosome}.convergence.txt.gz')), compress='gzip', quote=F, row.names=F, sep = '\t')

}