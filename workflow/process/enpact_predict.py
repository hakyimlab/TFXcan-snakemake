

import os, sys, re
import pandas as pd, numpy as np
import argparse
import multiprocessing
import itertools

global split_and_save

# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--matrix", help="Path to file", type=str)
parser.add_argument("--weights", help="Path to file", type=str)
parser.add_argument("--metadata", help="Path to file", type=str)
parser.add_argument("--split", help="", action=argparse.BooleanOptionalAction)
parser.add_argument("--output_basename", help="", type=str)
args = parser.parse_args()

# '/beagle3/haky/users/temi/projects/Enpact/files/ENPACT_734_2024-07-26.compiled_weights.lambda.1se.txt.gz'
# enpact weights
dt_weights = pd.read_table(args.weights).drop(columns = ['feature'])
mnames = dt_weights.columns.tolist()
weights = dt_weights.to_numpy()

xt = pd.read_table(args.matrix)
enpact_predictions = np.matmul(xt, weights)
enpact_predictions = pd.DataFrame(enpact_predictions)
enpact_predictions.columns = mnames

# metadata
df_metadata = pd.read_table(args.metadata)

# merge
df_predictions = pd.concat([df_metadata, enpact_predictions], axis = 1)

# create the directory
if not os.path.exists(os.path.dirname(args.output_basename)):
    os.makedirs(os.path.dirname(args.output_basename), exist_ok=True)


def split_and_save(tf_tissue, df, output_basename):
    df = df[['locus', 'individual', tf_tissue]]
    df = df[['locus', 'individual', tf_tissue]].pivot(index='locus', columns='individual', values=tf_tissue)
    new_indices = [f'{tf_tissue}_{dd}' for dd in df.index]
    df.insert(0, 'NAME', new_indices)

    bdt = [dd.split('_') for dd in new_indices]
    cdt = [(re.sub("chr", "", gg[2]), gg[3], gg[4], '_'.join(gg), '_'.join(gg), 'protein_coding') for gg in bdt]
    mdt = pd.DataFrame(cdt, columns=['chr', 'start', 'end', 'gene_id', 'gene_name', 'gene_type'])

    df.to_csv(f'{output_basename}.{tf_tissue}.enpact_scores.txt', sep = '\t', index = False)
    mdt.to_csv(f'{output_basename}.{tf_tissue}.annotation.txt', sep = '\t', index = False)

    return(0)

if args.split:
    pool = multiprocessing.Pool(16)
    outputs_list = pool.starmap(split_and_save, itertools.product(mnames, [df_predictions], [args.output_basename]))

    # for tf_tissue in mnames:
    #     df_tf_tissue = df[['locus', 'individual', tf_tissue]]
    #     df_tf_tissue = df_tf_tissue[['locus', 'individual', tf_tissue]].pivot(index='locus', columns='individual', values=tf_tissue)
    #     new_indices = [f'{tf_tissue}_{dd}' for dd in df_tf_tissue.index]
    #     #df_tf_tissue['NAME'] = new_indices
    #     df_tf_tissue.insert(0, 'NAME', new_indices)

    #     bdt = [dd.split('_') for dd in new_indices]
    #     cdt = [(re.sub("chr", "", gg[2]), gg[3], gg[4], '_'.join(gg), '_'.join(gg), 'protein_coding') for gg in bdt]
    #     mdt = pd.DataFrame(cdt, columns=['chr', 'start', 'end', 'gene_id', 'gene_name', 'gene_type'])

    #     df_tf_tissue.to_csv(f'{args.output_basename}.{tf_tissue}.enpact_scores.txt', sep = '\t', index = False)
    #     mdt.to_csv(f'{args.output_basename}.{tf_tissue}.annotation.txt', sep = '\t', index = False)
