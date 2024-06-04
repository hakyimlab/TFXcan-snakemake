

import pandas as pd
import os, glob, sys, re
from snakemake.io import glob_wildcards, Wildcards, expand
import numpy as np
import collections

sys.path.append('workflow/src')

#import module
print_progress = False
# configfile: '/project2/haky/temi/projects/TFXcan-snakemake/config/parameters.yaml'
# from snakemake.utils import validate
# validate(config, "config.schema.yaml")

# directories
DATA_DIR = 'data' 
INPUT_SUMSTATS = os.path.join(DATA_DIR, 'sumstats')
PROCESSED_SUMSTATS = os.path.join(DATA_DIR, 'processed_sumstats')
FINEMAPPING_DIR = os.path.join(DATA_DIR, 'finemapping')
COLLECTION_DIR = os.path.join(DATA_DIR, 'collection')
ENFORMER_PARAMETERS = os.path.join(DATA_DIR, 'enformer_parameters')
ENFORMER_PREDICTIONS = os.path.join(DATA_DIR, 'enformer_predictions') #os.path.join(config['scratch_dir'], f"{config['runname']}_{config['date']}", 'predictions_folder') if os.path.exists(config['scratch_dir']) else os.path.join(DATA_DIR, 'enformer_predictions')
AGGREGATED_PREDICTIONS = os.path.join(DATA_DIR, 'aggregated_predictions')
ENPACT_PREDICTIONS = os.path.join(DATA_DIR, 'enpact_predictions')
CHECKPOINTS_DIR = os.path.join(DATA_DIR, 'checkpoints')
PREDICTDB_DATA = os.path.join(DATA_DIR, 'predictdb')

def read_metadata(mtdt_file):
    dd = pd.read_csv(mtdt_file)
    #return({'phenotypes': dd.phenotype.tolist(), 'sumstats': dd.sumstat.tolist()})
    return(dict(zip(dd.phenotype.tolist(), dd.sumstat.tolist())))
    #return(dict(zip(dd.phenotype.tolist(), dd.sumstat.tolist())))

run_list = read_metadata(config["metadata"])





rule prepare_files_for_predictDB:
    input: os.path.join(ENPACT_PREDICTIONS, '{phenotype}')
    output: 
        enpact_scores = os.path.join(PREDICTDB_DATA, '{phenotype}', '{phenotype}.enpact_scores.txt'),
        annot_file = os.path.join(PREDICTDB_DATA, '{phenotype}', '{phenotype}.tf_tissue_annot.txt')
    params:
        rscript = config['rscript'],
        input_files = lambda wildcards, output: os.path.join(ENPACT_PREDICTIONS, wildcards.phenotype, f'{{}}.{wildcards.phenotype}.{config["enformer"]["aggtype"]}.{config["date"]}.csv.gz'),
        individuals = lambda wildcards: ','.join(collect_aggregated_individuals()),
        output_dir = os.path.join(PREDICTDB_DATA, '{phenotype}'),
        jobname = '{phenotype}',
        filtered_gwas = rules.collect_finemapping_results.output.finemapped_sumstats
    resources:
        partition="caslake"
    shell: "{params.rscript} workflow/src/prepare_files_for_predictDB.R --individuals {params.individuals} --input_file_pattern {params.input_files} --output_directory {params.output_dir} --filtered_GWAS_file {params.filtered_gwas} --phenotype {params.jobname}"


rule generate_lEnpact_models:
    input:
        enpact_scores = rules.prepare_files_for_predictDB.output.enpact_scores,
        annot_file = rules.prepare_files_for_predictDB.output.annot_file
    output:
        lEnpact_model = os.path.join('output', 'lEnpact', '{phenotype}/models/filtered_db/predict_db_{phenotype}_filtered.db')
    params:
        jobname = '{phenotype}',
        output_dir = os.path.abspath(os.path.join('output', 'lEnpact')),
        annot_file = os.path.abspath(rules.prepare_files_for_predictDB.output.annot_file),
        enpact_scores = os.path.abspath(rules.prepare_files_for_predictDB.output.enpact_scores)
    resources:
        partition="caslake",
        time="36:00:00"
    shell: "workflow/helpers/generate_lEnpact_models.EUR.midway3.sbatch {wildcards.phenotype} {params.output_dir} {params.annot_file} {params.enpact_scores}"
    