




import os, sys, re
import pandas as pd, numpy as np
import argparse
import multiprocessing
import itertools

global split_and_save

# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--matrix", help="[Input] The matrix of epigenomic features", type=str)
parser.add_argument("--weights", help="[Input] The weights", type=str)
parser.add_argument("--metadata", help="[Input] Metadata corresponding to the input matrix", type=str)
parser.add_argument("--split", help="[Input] Split output by TF-tissue", action=argparse.BooleanOptionalAction)
parser.add_argument("--output_basename", help="[Output] Basename for output files", type=str)
parser.add_argument("--subset_of_loci", help="[Input] Subset of loci to process", type=str, default=None)
args = parser.parse_args()

import os, re

# xt = pd.read_hdf('/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/childhood_asthma_2025-12-12/aggregated_predictions/childhood_asthma.childhood_asthma_2025-12-12.processed.matrix.h5.gz').to_numpy()

# df_metadata = pd.read_table("/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/childhood_asthma_2025-12-12/aggregated_predictions/childhood_asthma.childhood_asthma_2025-12-12.processed.metadata.tsv")

# xt.shape; xt[0:5, 0:5]

aggr_directory = "/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2024-09-30/aggregated_predictions/prostate_cancer_risk"

contents = os.listdir(aggr_directory)
# new_text = re.sub(r'\d+', 'NUM', text)
inds = [re.sub("_aggByCollect_prostate_cancer_risk.csv.gz", "", f) for f in contents]
# preds = {ind: pd.read_csv(os.path.join(aggr_directory, f"{ind}_aggByCollect_prostate_cancer_risk.csv.gz")) for ind in inds}

preds = list()
for ind in inds:
    pred = pd.read_csv(os.path.join(aggr_directory, f"{ind}_aggByCollect_prostate_cancer_risk.csv.gz"))
    pred['individual'] = ind
    preds.append(pred[['id', 'individual'] + [c for c in pred if c not in ['id', 'individual']]])

dt_preds = pd.concat(preds)
dt_preds.rename(columns={'id': 'locus'}, inplace=True)
dt_preds = dt_preds.sort_values(by=['locus', 'individual'])

dt_preds.iloc[0:5, 0:5]

# write out the metadata
metadata_file = f'/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2025-03-13/aggregated_predictions/prostate_cancer_risk.prostate_cancer_risk_2025-03-13.processed.metadata.tsv'
dt_preds[['locus', 'individual']].to_csv(metadata_file, sep = '\t', index = False)

xt = dt_preds.drop(['locus', 'individual'], axis = 1)
xt.shape
xt.iloc[0:5, 0:5]
# data_file = f'{args.output_basename}.matrix.tsv.gz'
data_file = f'/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2025-03-13/aggregated_predictions/prostate_cancer_risk.prostate_cancer_risk_2025-03-13.processed.matrix.h5.gz'
xt.to_hdf(data_file, key = 'matrix', mode='w', complevel=9)
#xt.to_csv(data_file, sep = '\t', index = False, compression =