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
PREDICTDB_DATA = os.path.join(DATA_DIR, 'predictdb.temp')
ENPACT_DB = os.path.join(DATA_DIR, 'enpactdb')
LENPACT_DIR = os.path.join(config['scratch_dir'], 'lEnpact', runmeta) if os.path.exists(config['scratch_dir']) else os.path.join(DATA_DIR, 'lEnpact')
SUMMARYTFXCAN_DIR = os.path.join(DATA_DIR, 'summaryTFXcan')
SUMMARY_OUTPUT = os.path.join(DATA_DIR, 'output')
BENCHMARK_DIR = os.path.join(DATA_DIR, 'benchmark')

def read_metadata(mtdt_file):
    dd = pd.read_csv(mtdt_file)
    return(dict(zip(dd.phenotype.tolist(), dd.sumstat.tolist())))

run_list = read_metadata(config["metadata"])
# read in the list of models
enpact_models_list = pd.read_table(config["enpact"]["weights"]).columns[1:].tolist() #.model.tolist()#[0:5]

# ensure that there is the reference panel annotation file
if 'processing' in config.keys() and 'reference_annotations' in config['processing'].keys(): 
    REFERENCE_ANNOTATIONS = config['processing']['reference_annotations']
elif 'predictdb' in config.keys() and 'reference_annotations' in config['predictdb'].keys():
    REFERENCE_ANNOTATIONS = config['predictdb']['reference_annotations']

# checkpoint functions
def collect_processed_summary_statistics(wildcards):
    checkpoint_output = checkpoints.process_summary_statistics.get(phenotype = wildcards.phenotype).output[0]
    chromosomes = glob_wildcards(os.path.join(checkpoint_output, 'chr{chrom}.sumstats.txt.gz')).chrom
    # output = expand(os.path.join(PROCESSED_SUMSTATS, f'{wildcards.phenotype}', 'chr{chrom}.sumstats.txt.gz'), chrom=chromosomes)
    output = expand(os.path.join(checkpoint_output, 'chr{chrom}.sumstats.txt.gz'), chrom=chromosomes)
    return output

def collect_chromosomes(wildcards):
    #checkpoint_output = checkpoints.process_summary_statistics.get(**wildcards).output[0]
    #checkpoint_output = #checkpoints.process_summary_statistics.get(**wildcards).output[0]
    checkpoint_output = checkpoints.process_summary_statistics.get(phenotype = wildcards.phenotype).output[0]
    #chromosomes = glob_wildcards(os.path.join(checkpoint_output, 'chr{chrom}.sumstats.txt.gz')).chrom #glob_wildcards(os.path.join(checkpoint_output, 'chr{chrom}.sumstats.txt.gz')).chrom
    #chromosomes = glob_wildcards(os.path.join(PROCESSED_SUMSTATS, f'{wildcards.phenotype}', 'chr{chrom}.sumstats.txt.gz')).chrom
    chromosomes = glob_wildcards(os.path.join(checkpoint_output, 'chr{chrom}.sumstats.txt.gz')).chrom
    return(chromosomes)


# def collect_completed_models(model, phenotype):
#     checkpoint_output = checkpoints.summary_TFXcan.get(model = model, phenotype = phenotype).output[0]
#     mlist = list()
#     for m in enpact_models_list:
#         checkpoint_output = checkpoints.summary_TFXcan.get(model = m, **wc).output[0]
#         m = glob_wildcards(checkpoint_output).model
#     return(m)

def collect_completed_summary_tfxcan(wc):
    # for m in enpact_models_list:
    #     checkpoint_output = checkpoints.summary_TFXcan.get(model = m, **wc).output[0]
    checkpoint_output = checkpoints.summary_TFXcan.get(model = m, **wc).output[0]
    #checkpoint_output = os.path.dirname(checkpoints.summary_TFXcan.get(phenotype = wc.phenotype).output[0])
    #exp = expand(os.path.join(SUMMARYTFXCAN_DIR, f'{wc.phenotype}', f'{{model}}-{wc.phenotype}.enpactScores.spredixcan.csv'), model = enpact_models_list)

    sxcan_files = glob_wildcards(os.path.join(SUMMARYTFXCAN_DIR, f"{wc.phenotype}", '{sxcan}')).sxcan
    exp = [os.path.join(SUMMARYTFXCAN_DIR, f"{wc.phenotype}", f"{sxcan}") for sxcan in sxcan_files]

    print(exp)
    return(exp)

    #lambda wildcards: checkpoints.summary_TFXcan.get(model = enpact_models_list[0], **wildcards).output[0]

    # checkpoint_output = checkpoints.summary_TFXcan.get(**wildcards).output[0]
    # print(checkpoint_output)
    # #sxcan_files = glob_wildcards(os.path.join(SUMMARYTFXCAN_DIR, f"{wc.phenotype}", '{sxcan}')).sxcan

    # expand(os.path.join(SUMMARYTFXCAN_DIR, '{phenotype}', '{model}-{phenotype}.enpactScores.spredixcan.csv'), phenotype = run_list.keys(), model = enpact_models_list)

    # print(sxcan_files)

    # data/T2D_2025-01-24/summaryTFXcan/t2d_suzuki/AR_Breast-t2d_suzuki.enpactScores.spredixcan.csv
    # checkpoint_output = checkpoints.summary_TFXcan.get(**wildcards).output[0]
    # mm, ph, _ = glob_wildcards(checkpoint_output)
    # return(os.path.join(SUMMARYTF`CAN_DIR, "{mm}", "{ph}", '{pp}.enpactScores.spredixcan.csv'))
    #print(f"Wildcards: {wildcards}")
    # checkpoint_output = checkpoints.summary_TFXcan.get(**wildcards).output[0]
    # print(**wildcards)
    #exp = expand(os.path.join(SUMMARYTFXCAN_DIR, '{phenotype}', '{model}-{phenotype}.enpactScores.spredixcan.csv'), zip, model=mm, phenotype=ph)
    #exp = [os.path.join(SUMMARYTFXCAN_DIR, f"{wc.phenotype}", f"{sxcan}") for sxcan in sxcan_files]
    # return(exp)

print(f"INFO - Found {len(run_list)} phenotypes for TFXcan analysis.")
print(run_list)

onstart:
    print(run_list)
    print(f"INFO - Found {len(enpact_models_list)} models for TF-tissue:GWAS association tests.")
    if config['runSusie'] == True:
        print(f'INFO - Running ENFORMER is set to True. Running the pipeline with ENFORMER...')
    elif config['runSusie'] == False:
        print(f'INFO - Running Susie is set to False. Choosing top SNPs per significant GWAS loci...')


rule all:
    input:
        #[os.path.join(INPUT_SUMSTATS, f'{sumstat}') for sumstat in run_list.values()],
        #expand(os.path.join(PROCESSED_SUMSTATS, '{phenotype}'), phenotype = run_list.keys()),
        # expand(os.path.join(ENFORMER_PARAMETERS, f'enformer_parameters_{runname}_{{phenotype}}.yaml'), phenotype = run_list.keys()),
        # expand(os.path.join(ENFORMER_PARAMETERS, f'aggregation_config_{runname}_{{phenotype}}.yaml'), phenotype = run_list.keys()),
        # expand(os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.h5'), phenotype = run_list.keys()),
        # expand(os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.processed.metadata.tsv'), phenotype = run_list.keys()),
        # expand(os.path.join(AGGREGATED_PREDICTIONS, f'{{phenotype}}.{runmeta}.processed.matrix.h5.gz'), phenotype = run_list.keys()),
        expand(os.path.join(PREDICTDB_DATA, "{phenotype}", '{phenotype}.{model}.enpact_scores.txt'), phenotype = run_list.keys(), model = enpact_models_list),
        expand(os.path.join(PREDICTDB_DATA, "{phenotype}", '{phenotype}.{model}.annotation.txt'), phenotype = run_list.keys(), model = enpact_models_list),
        expand(os.path.join(LENPACT_DIR, '{phenotype}', "{model}", 'models/filtered_db/predict_db_{phenotype}_filtered.db'), phenotype = run_list.keys(), model = enpact_models_list),
        expand(os.path.join(LENPACT_DIR, '{phenotype}', "{model}", 'models/filtered_db/predict_db_{phenotype}_filtered.txt.gz'), phenotype = run_list.keys(), model = enpact_models_list),
        expand(os.path.join(LENPACT_DIR, '{phenotype}', "{model}", 'models/filtered_db/Covariances.varID.txt.gz'), phenotype = run_list.keys(), model = enpact_models_list),
        expand(os.path.join(SUMMARYTFXCAN_DIR, '{phenotype}', "{model}-{phenotype}.enpactScores.spredixcan.csv"), phenotype = run_list.keys(), model = enpact_models_list),
        expand(os.path.join(SUMMARY_OUTPUT, '{phenotype}.enpactScores.{rundate}.spredixcan.txt'), phenotype = run_list.keys(), rundate = [rundate])

if config['runSusie'] == True:
    include: 'workflow/rules/ruleAll.Susie.smk'
    include: 'workflow/rules/TFXcan.Susie.smk'
elif config['runSusie'] == False:
    #include: 'workflow/rules/ruleAll.noSusie.smk'
    include: 'workflow/rules/TFXcan.noSusie.smk'
else:
    print(f"ERROR - [FATAL] Please specify whether to run ENFORMER or not. Exiting...")
    sys.exit(1)