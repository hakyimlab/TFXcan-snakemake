

import os, sys, re, glob
import dataclasses
from dataclasses import dataclass, fields
from typing import Union, Any
import pandas as pd, numpy as np
import argparse

def list_of_ints(arg):
    return list(map(int, arg.split(',')))

# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--config", help="Path to file", type=str, default=None)
parser.add_argument("--predictions_logfile", help="Path to file", type=str)
parser.add_argument("--predictions_directory", help="Path", type=str)
parser.add_argument("--expected_shape", help="", type=list_of_ints)
parser.add_argument("--output_directory", help="", type=str, required=True)
parser.add_argument("--output_filename", help="", type=str, required=True)
args = parser.parse_args()

# some locations and folders
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
batch_utils_path = os.path.join(script_path, 'modules')
sys.path.append(batch_utils_path)

import mergeUtils

if args.config is not None and os.path.exists(args.config):
    import yaml
    print(f'INFO - Merging with config file from {args.config}')
    with open(f'{args.config}') as f:
        parameters = yaml.load(f, Loader=yaml.FullLoader)
    args = parser.set_defaults(**parameters)
    args, unknown = parser.parse_known_args()

print(args)

# merge the hdf5 files
nps = pd.read_csv(args.predictions_logfile, sep='\t').shape[0]
xshape = (nps, args.expected_shape[0], args.expected_shape[1])

print(f'INFO - Expected shape of combined hdf5 file is: {xshape}')

# list the files in the directory
output_files_list = glob.glob(os.path.join(args.predictions_directory, '*.h5'))
print(f'INFO - Found {len(output_files_list)} files in the directory')

# create the output directory
if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

ofile = os.path.join(args.output_directory, args.output_filename)

# print(ofile)

out = mergeUtils.merge_hdf5_files(output_files_list, xshape = xshape, output_file = ofile)
print(f"INFO - Combined HDF5 is at {out}")

# print(f'INFO - Combined HDF5 is at {ofile}')
# python3 src/enformer_merge.py --config /beagle3/haky/users/temi/projects/chrom-enpact/output/aracena200/predictions_2024-10-15/metadata/aracena200.aggregation_config.yaml
# python3 src/enformer_merge.py --output_directory /beagle3/haky/users/temi/projects/chrom-enpact/output/aracena200/predictions_2024-10-15/ --output_filename pp.h5 --predictions_logfile '/beagle3/haky/users/temi/projects/chrom-enpact/output/aracena200/predictions_2024-10-15/predictions_log/aracena200.prediction_log.txt' --predictions_directory='/beagle3/haky/users/temi/projects/chrom-enpact/output/aracena200/predictions_2024-10-15/predictions' --expected_shape="1,5313"




# # I want to modify arguments as appropriate so I will be using dataclasses
# @dataclass
# class Arguments:
#     # must supply
#     expected_shape: Union[tuple, list]
#     output_directory: str
#     output_filename: str
#     predictions_logfile: str
#     predictions_directory: str

# arguments = Arguments(expected_shape = tuple(args.expected_shape), 
#                     output_directory = args.output_directory, 
#                     output_filename = args.output_filename, 
#                     predictions_logfile = args.predictions_logfile, 
#                     predictions_directory = args.predictions_directory)

# print(arguments)