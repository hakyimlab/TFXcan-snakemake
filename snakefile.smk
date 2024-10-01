# Description: This snakemake file runs the full TFXcan pipeline
# Author: Temi
# Date: Wed Mar 29 2023
# Usage: --> see README

import pandas as pd
import os, glob, sys, re
from snakemake.io import glob_wildcards, Wildcards, expand
import numpy as np
import collections
from datetime import datetime

sys.path.append('workflow/src')

#import module
print_progress = False
# configfile: '/project2/haky/temi/projects/TFXcan-snakemake/config/parameters.yaml'
# from snakemake.utils import validate
# validate(config, "config.schema.yaml")

rundate = config['date'] if isinstance(config['date'], str) else datetime.today().strftime('%Y-%m-%d')
runname = config['runname']
runmeta = f"{runname}_{rundate}"

# directories
INPUT_SUMSTATS = os.path.join('data', 'sumstats')
DATA_DIR = os.path.join('data', f"{runname}_{rundate}")
PROCESSED_SUMSTATS = os.path.join(DATA_DIR, 'processed_sumstats')
FILTERING_DIR = os.path.join(DATA_DIR, 'filtering')
COLLECTION_DIR = os.path.join(DATA_DIR, 'collection')
ENFORMER_PARAMETERS = os.path.join(DATA_DIR, 'enformer_parameters')
ENFORMER_PREDICTIONS = os.path.join(config['scratch_dir'], 'predictions_folder') if os.path.exists(config['scratch_dir']) else os.path.join(DATA_DIR, 'predictions_folder')
AGGREGATED_PREDICTIONS = os.path.join(DATA_DIR, 'aggregated_predictions')
ENPACT_PREDICTIONS = os.path.join(config['scratch_dir'], 'enpact_predictions') if os.path.exists(config['scratch_dir']) else os.path.join(DATA_DIR, 'enpact_predictions')
CHECKPOINTS_DIR = os.path.join(DATA_DIR, 'checkpoints')
PREDICTDB_DATA = os.path.join(DATA_DIR, 'predictdb')
ENPACT_DB = os.path.join(DATA_DIR, 'enpactdb')
LENPACT_DIR = os.path.join(DATA_DIR, 'lEnpact')
SUMMARYTFXCAN_DIR = os.path.join(DATA_DIR, 'summaryTFXcan')
BENCHMARK_DIR = os.path.join(DATA_DIR, 'benchmark')

def read_metadata(mtdt_file):
    dd = pd.read_csv(mtdt_file)
    return(dict(zip(dd.phenotype.tolist(), dd.sumstat.tolist())))

run_list = read_metadata(config["metadata"])

# checkpoint functions
def collect_processed_summary_statistics(wildcards):
    chromosomes = glob_wildcards(os.path.join(PROCESSED_SUMSTATS, '{phenotype}', 'chr{chrom}.sumstats.txt.gz')).chrom
    output = expand(os.path.join(PROCESSED_SUMSTATS, "{phenotype}", 'chr{chrom}.sumstats.txt.gz'), chrom=chromosomes, phenotype = run_list.keys())
    return output

def collect_filtered_summary_statistics(wildcards):
    #checkpoint_output = checkpoints.process_summary_statistics.get(**wildcards).output[0]
    chromosomes = glob_wildcards(os.path.join(FILTERING_DIR, '{phenotype}', 'chr{chrom}.sumstats.txt.gz')).chrom #glob_wildcards(os.path.join(checkpoint_output, 'chr{chrom}.sumstats.txt.gz')).chrom
    output = expand(os.path.join(PROCESSED_SUMSTATS, "{phenotype}", 'chr{chrom}.sumstats.txt.gz'), chrom=chromosomes, phenotype = run_list.keys())
    # print(output)
    return output

def collect_chromosomes(wildcards):
    #checkpoint_output = checkpoints.process_summary_statistics.get(**wildcards).output[0]
    chromosomes = glob_wildcards(os.path.join(PROCESSED_SUMSTATS, "{phenotype}", 'chr{chrom}.sumstats.txt.gz')).chrom
    return(chromosomes)

def collect_aggregated_individuals():
    phenotype, individuals, _ = glob_wildcards(os.path.join(AGGREGATED_PREDICTIONS, '{phenotype}', f'{{individual}}_{config["enformer"]["aggtype"]}_{{pp}}.csv.gz'))
    return(individuals)

def collect_enpact_individuals():
    #ids = glob_wildcards(os.path.join(ENPACT_PREDICTIONS, '{phenotype}', '{idi}.*.csv.gz')).idi
    #print(ENPACT_PREDICTIONS)
    _, ids,_ = glob_wildcards(os.path.join(ENPACT_PREDICTIONS, '{phenotype}', f'{{idi}}.{{phenop}}.{config["enformer"]["aggtype"]}.{rundate}.csv.gz'))

    #print(ids)
    return(ids)
# rule all:
#     input:
#         expand(os.path.join(INPUT_SUMSTATS, '{phenotype}.liftover.logistic.assoc.tsv.gz'), phenotype = run_list.keys()),
#         expand(os.path.join(PROCESSED_SUMSTATS, '{phenotype}'), phenotype = run_list.keys()),
#         expand(os.path.join(COLLECTION_DIR, '{phenotype}.EnformerLoci.topSNPs.txt'), phenotype = run_list.keys()),
#         expand(os.path.join(ENFORMER_PARAMETERS, f'enformer_parameters_{runname}_{{phenotype}}.json'), phenotype = run_list.keys()),
#         expand(os.path.join(ENFORMER_PARAMETERS, f'aggregation_config_{runname}_{{phenotype}}.json'), phenotype = run_list.keys()),
#         # expand(os.path.join(CHECKPOINTS_DIR, '{phenotype}.checkpoint'), phenotype = run_list.keys()),
#         expand(os.path.join(AGGREGATED_PREDICTIONS, '{phenotype}'), phenotype = run_list.keys()),


rule all:
    input:
        # expand(os.path.join(INPUT_SUMSTATS, '{phenotype}.liftover.logistic.assoc.tsv.gz'), phenotype = run_list.keys()),
        # expand(os.path.join(PROCESSED_SUMSTATS, '{phenotype}'), phenotype = run_list.keys()),
        # #expand(os.path.join(COLLECTION_DIR, '{phenotype}.EnformerLoci.topSNPs.3.txt'), phenotype = run_list.keys()),
        # expand(os.path.join(ENFORMER_PARAMETERS, f'enformer_parameters_{runname}_{{phenotype}}.json'), phenotype = run_list.keys()),
        # expand(os.path.join(ENFORMER_PARAMETERS, f'aggregation_config_{runname}_{{phenotype}}.json'), phenotype = run_list.keys()),
        # expand(os.path.join(CHECKPOINTS_DIR, '{phenotype}.checkpoint'), phenotype = run_list.keys()),
        expand(os.path.join(AGGREGATED_PREDICTIONS, '{phenotype}'), phenotype = run_list.keys()),
        expand(os.path.join(ENPACT_PREDICTIONS, '{phenotype}'), phenotype = run_list.keys()),
        expand(os.path.join(ENPACT_DB, "{phenotype}.enpact_scores.array.rds.gz"), phenotype = run_list.keys()),
        expand(os.path.join(ENPACT_DB, "{phenotype}.enpact_scores.txt.gz"), phenotype = run_list.keys()),
        expand(os.path.join(PREDICTDB_DATA, '{phenotype}.enpact_scores.txt'), phenotype = run_list.keys()),
        expand(os.path.join(PREDICTDB_DATA, '{phenotype}.tf_tissue_annot.txt'), phenotype = run_list.keys()),
        expand(os.path.join(LENPACT_DIR, '{phenotype}', 'models/filtered_db/predict_db_{phenotype}_filtered.db'), phenotype = run_list.keys()),
        expand(os.path.join(LENPACT_DIR, '{phenotype}', 'models/database/predict_db_{phenotype}.txt.gz'), phenotype = run_list.keys()),
        expand(os.path.join(LENPACT_DIR, '{phenotype}', 'models/database/predict_db_{phenotype}.txt'), phenotype = run_list.keys()),
        expand(os.path.join(LENPACT_DIR, '{phenotype}', 'models/database/Covariances.varID.txt'), phenotype = run_list.keys()),
        expand(os.path.join(SUMMARYTFXCAN_DIR, '{phenotype}', '{phenotype}.enpactScores.spredixcan.csv'), phenotype = run_list.keys())
        # expand(os.path.join('output', 'lEnpact', '{phenotype}', 'models/filtered_db/predict_db_{phenotype}_filtered.db'), phenotype = run_list.keys()),
        # expand(os.path.join('output', 'lEnpact', '{phenotype}', 'models/database/predict_db_{phenotype}.txt.gz'), phenotype = run_list.keys()),
        # expand(os.path.join('output', 'lEnpact', '{phenotype}', 'models/database/Covariances.varID.txt'), phenotype = run_list.keys()),
        # expand(os.path.join('output', 'summary', '{phenotype}', '{phenotype}.enpactScores.spredixcan.csv'), phenotype = run_list.keys())

if config['runSusie'] == True:
    print(f'INFO - Running ENFORMER is set to True. Running the pipeline with ENFORMER...')
    include: 'workflow/rules/ruleAll.Susie.smk'
    include: 'workflow/rules/TFXcan.Susie.smk'
elif config['runSusie'] == False:
    print(f'INFO - Running Susie is set to False. Choosing top SNPs per significant GWAS loci...')
    #include: 'workflow/rules/ruleAll.noSusie.smk'
    include: 'workflow/rules/TFXcan.noSusie.smk'
else:
    print(f"ERROR - [FATAL] Please specify whether to run ENFORMER or not. Exiting...")
    sys.exit(1)