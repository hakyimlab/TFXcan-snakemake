# Utilities to predict using enformer
# Author: Temi
# Date: Thurs Feb 2 2023
# This has the following functions:
# - get_model
# - enformer_predict_on_sequence
# - batch
# - run_batch_predictions
# - return_prediction_function
# - generate_n_batches
# - generate_batch_n_elems
# - make_h5_db
# - make_h5_db_parsl
# - check_predictions_and_logs
# - return_check_function
# - enformer_predict_on_batch


import functools
import numpy as np, time, h5py, os, sys
import tensorflow as tf
import sequencesUtils, checksUtils, collectUtils, loggerUtils, saveUtils, batchUtils

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

@functools.lru_cache(5)
def get_model(model_path):
    """
    Return a tensorflow model

    Parameters:
        model_path: str
            A path to where the tensorflow model exists
    Returns: 
        a tensorflow model
    """
    import tensorflow as tf
    return tf.saved_model.load(model_path).model

def enformer_predict_on_sequence(model, sample_input):
    """
    given a compatible sequence that has been one-hot encoded, predict on ENFORMER

    Parameters:
        model: a tensorflow model
        sample_input: a (1, 393216, 4) np.array that is a one-hot encoding of a sequence

    Returns: A dictionary
        of the form {'haplotype': _predictions_}
        _predictions_ is a numpy array of shape (17, 5313) numpy array of predictions
    """
    
    prediction_output = {}
    for haplotype, sequence_encoding in sample_input.items():
        if not sequence_encoding.shape == (1, 393216, 4):
            raise Exception(f'[ERROR] Fatal. Input sequence shape is not appropriate')
        # prediction = model.predict_on_batch(sequence_encoding)['human'].numpy()[: , range(448 - 8, (448 + 8 + 1)), : ]
        prediction = model.predict_on_batch(sequence_encoding)['human'].numpy()

        prediction_output[haplotype] = prediction
        
    return(prediction_output)

#@python_app
def check_for_batches_to_run(batch_dictionary, module_directives): #
    """
    Predict and save on a given batch of regions in the genome

    This function also filters the regions in the batch for those that have been predicted and logged

    Parameters:
        batch_regions: list
            A list of regions with each element (region) in the form `chr_start_end`.
        batch_num: num
            The number of the batch e.g. batch 1, batch 2, e.t.c.
        samples: list
            a list of samples: [a, b, c]
        output_dir: str (path)
            Where should the predictions be saved? Predictions are saved as `{sample}/{haplotype0, haplotype1, haplotype2}/{region}_predictions.h5` 
        path_to_vcf: str (path)
            The path to the vcf file
        prediction_logfiles_folder: str (path)
            When predictions are made, they are logged in this folder and saved as {sample}_log.csv. Also useful for checking the log to prevent re-predicting an existing prediction. 
        script_path: str (path), default is the path to where this file is.
            The path to this module.
        sequence_source: one of 'personalized' or 'reference'

    Returns: num
        A single value of either 0 (if predictions were successful) or 1 (if there was a error).
        Check the call logs or stacks for the source of the error. 
    """

    import sys, os, faulthandler, time, importlib, yaml

    #print(batch_dictionary)
    sys.path.append(module_directives['path_to_modules'])
    exec(open(os.path.join(module_directives['path_to_modules'], 'datatypeUtils.py')).read(), globals())
    

    # spec = importlib.util.spec_from_file_location("predictionUtils", module_directives.predictionUtils)
    # predictionUtils = importlib.util.module_from_spec(spec)
    # #predictUtils_two.tmp_config_path = tmp_config_path
    # spec.loader.exec_module(predictionUtils)

    spec = importlib.util.spec_from_file_location("checksUtils", module_directives['checksUtils'])
    checksUtils = importlib.util.module_from_spec(spec)
    #predictUtils_two.tmp_config_path = tmp_config_path
    spec.loader.exec_module(checksUtils)

    # spec = importlib.util.spec_from_file_location("datatypeUtils", module_directives.checksUtils)
    # datatypeUtils = importlib.util.module_from_spec(spec)
    # #predictUtils_two.tmp_config_path = tmp_config_path
    # spec.loader.exec_module(datatypeUtils)
    # exec()


    batch_directives = batch_dictionary #AttrDict(batch_dictionary)

    # mpath = os.path.join(batch_directives.script_path, 'modules') #os.path.dirname(__file__) #
    # sys.path.append(mpath)
    
    faulthandler.enable() # to figure out where segmentation error is coming from
    

    # try:
    #     import predictionRunUtils
    # except ModuleNotFoundError as merr:
    #     raise Exception(f'ERROR - {type(merr).__name__} at run_batch_predictions. Cannot locate either of predictionRunUtils.')

    # I want to import predictUtils_two but I need it to be dynamic and imported with the path to the config file
    # try:
    #     import checksUtils, predictionUtils
    #     #import predictUtils_two
    # except ModuleNotFoundError as merr:
    #     raise Exception(f'ERROR - {type(merr).__name__} at run_batch_predictions. Cannot locate either of `checkUtils` or `predictionUtils` modules.')

    # check_results = {sample: checksUtils.check_queries(sample=sample, queries=batch_directives.batch_regions, output_dir=batch_directives.output_directory, prediction_logfiles_folder=batch_directives.prediction_logfiles_folder, sequence_source=batch_directives.sequence_source) for sample in batch_directives.samples}

    check_results = {sample: checksUtils.check_queries(sample=sample, queries=batch_directives["batch_regions"], output_dir=batch_directives['output_directory'], prediction_logfiles_folder=batch_directives['prediction_logfiles_folder'], sequence_source=batch_directives['sequence_source']) for sample in batch_directives['samples']}

    filtered_check_result = {k: v for k, v in check_results.items() if k is not None}

    # remove empty dict values
    for sample in list(filtered_check_result.keys()):
        if not filtered_check_result[sample]:
            del filtered_check_result[sample]

    return((filtered_check_result, batch_directives))

# def make_h5_db(h5_file, csv_file, files_list, files_path, dataset):
#     import h5py
#     import pandas as pd

#     with h5py.File(f"{h5_file}", "w") as f_dst:
#         #h5files = [f for f in os.listdir(f'{project_dir}/') if f.endswith(".h5")]

#         dset = f_dst.create_dataset(f'{dataset}_dataset', shape=(len(files_list), 17, 5313), dtype='f4')
#         for i, filename in enumerate(files_list):
#             with h5py.File(files_path[i]) as f_src:
#                 dset[i] = f_src[filename]

#     pd.DataFrame(files_list, columns=['region']).to_csv(f"{csv_file}")

#     return(0)

# def make_h5_db_parsl(use_parsl, fxn=make_h5_db):
#     '''
#     Decorate or not the `make_h5_db` function based on whether `use_parsl` is true or false

#     Returns: 
#         function object
#         The function if parsl is not used
#         The parsl decorated function if parsl is meant to be used
    
#     '''
#     from parsl.app.app import python_app
#     if use_parsl == True:
#         return python_app(fxn)
#     elif use_parsl == False:
#         return fxn

#         # if batch_directives.mainRunScript is not None:
#         #     spec = importlib.util.spec_from_file_location("mainRun", batch_directives.mainRunScript)
#         # else:
#         #     raise Exception('ERROR - The mainRun.py script is missing')

#         # mainRun = importlib.util.module_from_spec(spec)
#         # mainRun.tmp_config_path = batch_directives.tmp_config_path
#         # spec.loader.exec_module(mainRun)