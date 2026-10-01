#!/bin/bash
# Converts the GWAS Catalog harmonised file (GCST90428427, hg38) in minimal/data/ into the columns
# workflow/process/process_summary_statistics.R expects: chrom pos variant_id ref alt pval zscore beta se rsid
# effect_allele -> alt, other_allele -> ref. Drops rows with missing p-value, beta or se.
# Usage (from the repo root): bash minimal/prepare_minimal_gwas.sh

set -euo pipefail
in=minimal/data/GCST90428427.h.tsv.gz
out=minimal/data/minimal.gwas_sumstats.hg38.processed.txt.gz

zcat "$in" | awk -F'\t' -v OFS='\t' '
    NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i
              print "chrom", "pos", "variant_id", "ref", "alt", "pval", "zscore", "beta", "se", "rsid"; next }
    {
        chrom = $col["chromosome"]; pos = $col["base_pair_location"]
        ref = $col["other_allele"]; alt = $col["effect_allele"]
        beta = $col["beta"]; se = $col["standard_error"]; p = $col["p_value"]; rsid = $col["rsid"]
        if (p == "NA" || beta == "NA" || se == "NA" || se + 0 == 0) next
        print chrom, pos, chrom "_" pos "_" ref "_" alt, ref, alt, p, beta / se, beta, se, rsid
    }' | gzip > "$out"

echo "wrote $out: $(zcat "$out" | tail -n +2 | wc -l) SNPs"
