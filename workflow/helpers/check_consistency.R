

snp_annots <- data.table::fread('ALL.chr6.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.snp_annot.txt.gz')
sumstats <- data.table::fread("/project2/haky/temi/projects/TFXcan-snakemake/data/processed_sumstats/asthma_children.liftover.logistic.assoc.tsv.gz") %>%
    dplyr::filter(chrom == '6')




