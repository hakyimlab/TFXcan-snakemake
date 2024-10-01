


rule all:
    input:
        expand(os.path.join(INPUT_SUMSTATS, '{phenotype}.liftover.logistic.assoc.tsv.gz'), phenotype = run_list.keys()),
        expand(os.path.join(PROCESSED_SUMSTATS, '{phenotype}'), phenotype = run_list.keys()),
        expand(os.path.join(COLLECTION_DIR, '{phenotype}.EnformerLoci.topSNPs.3.txt'), phenotype = run_list.keys()),
        expand(os.path.join(ENFORMER_PARAMETERS, f'enformer_parameters_{runname}_{{phenotype}}.json'), phenotype = run_list.keys()),
        expand(os.path.join(ENFORMER_PARAMETERS, f'aggregation_config_{runname}_{{phenotype}}.json'), phenotype = run_list.keys()),
        # expand(os.path.join(CHECKPOINTS_DIR, '{phenotype}.checkpoint'), phenotype = run_list.keys()),
        expand(os.path.join(AGGREGATED_PREDICTIONS, '{phenotype}'), phenotype = run_list.keys()),
        expand(os.path.join(ENPACT_PREDICTIONS, '{phenotype}'), phenotype = run_list.keys()),
        expand(os.path.join(ENPACT_DB, "{phenotype}.enpact_scores.array.rds.gz"), phenotype = run_list.keys()),
        expand(os.path.join(ENPACT_DB, "{phenotype}.enpact_scores.txt.gz"), phenotype = run_list.keys()),
        expand(os.path.join(PREDICTDB_DATA, '{phenotype}.enpact_scores.txt'), phenotype = run_list.keys()),
        expand(os.path.join(PREDICTDB_DATA, '{phenotype}.tf_tissue_annot.txt'), phenotype = run_list.keys()),
        expand(os.path.join('output', 'lEnpact', '{phenotype}', 'models/filtered_db/predict_db_{phenotype}_filtered.db'), phenotype = run_list.keys()),
        expand(os.path.join('output', 'lEnpact', '{phenotype}', 'models/database/predict_db_{phenotype}.txt.gz'), phenotype = run_list.keys()),
        expand(os.path.join('output', 'lEnpact', '{phenotype}', 'models/database/Covariances.varID.txt'), phenotype = run_list.keys()),
        expand(os.path.join('output', 'summary', '{phenotype}', '{phenotype}.enpactScores.spredixcan.csv'), phenotype = run_list.keys())