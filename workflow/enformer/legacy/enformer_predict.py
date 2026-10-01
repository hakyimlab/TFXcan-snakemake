# Description: This script is used for ENFORMER inference
# Author: Temi
# Date: Wed 25 Jan 2023
# Usage: 

import os, sys, time, argparse

# needed arguments 
parser = argparse.ArgumentParser()
parser.add_argument("--parameters", help="Path to YAML file of parameters and directives to be used by ENFORMER", type=str, required=True)
args = parser.parse_args()

# some locations and folders
whereis_script = os.path.dirname(__file__) #os.path.dirname(sys.argv[0]) # or os.path.dirname(__file__)
script_path = os.path.abspath(whereis_script)
print(script_path)
batch_utils_path = os.path.join(script_path, 'modules')
sys.path.append(batch_utils_path)

exec(open(os.path.join(batch_utils_path, 'main.py')).read(), globals())

# (__name__ == '__main__') or (
if (__name__ == '__main__') or (__name__ == 'enformer_predict'):
    #check_input_parameters.check_inputs(args.param_config)

    # job stats will always be written
    import time
    job_start = time.perf_counter()
    main(args, script_path)
    job_end = time.perf_counter()

    job_runtime = job_end - job_start
    print(f'INFO - Completed job in {job_runtime} seconds.')



# parsl==2023.4.17