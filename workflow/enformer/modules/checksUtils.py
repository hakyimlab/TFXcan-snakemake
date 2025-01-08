# Author: Temi
# DAte: Wed Sept 25 2024
# Description: This script is used for ENFORMER inference
# The following functions are available:
# - return_samples_to_predict_on
# - return_sample_logging_type
# - check_queries
# - check_predictions_and_logs
# - return_check_function

import os, sys, parsl

def return_samples_to_predict_on(query, logging_list_per_sample):
    ss = [sample for sample in logging_list_per_sample.keys() if query in logging_list_per_sample[sample]]
    return(ss)

def return_sample_logging_type(sample, query_region, logging_dictonary):
    result = [elem['logtype'] for elem in logging_dictonary[sample] if elem['query'] == query_region][0]

    if not result:
        return('n')

    return(result)

def check_queries(sample, queries, output_dir, prediction_logfiles_folder, sequence_source):
    
    """
    Check whether a given region, for an individual has been predicted and logged.

    Parameters:
        sample: str 
            The name/id of an individual/sample
        queries: an iterable
            A batch of regions each in the genome in the form `chr_start_end`.
        output_dir: str (path)
            The folder where the predictions should have been logged. 
        prediction_logfiles_folder: a directory
            A directory within which the {sample}_log.csv file should have been logged. 
            If the file is not found, the regions are returned and will be logged after prediction.
        sequence_source: str ('personalized' or 'reference')
            where is the sequence sourced from? 
            This argument is needed to know what folders to look into within the output directory.
    
    Returns: dict
        'query': the query region if it has not been logged or predictions don't exist
        'logtype': whether it should be logged if it has not been logged i.e. 'y' or 'n'

    If predictions exist and the query has been logged, this function returns None.
    """
    import pandas as pd
    import os
    import numpy as np

    if isinstance(prediction_logfiles_folder, type(None)):
        output = [{'query':query, 'logtype':'y'} for query in queries]
        return(output)

    else:
        id_logfile = os.path.join(prediction_logfiles_folder, f'{sample}_log.csv')
        if os.path.isfile(id_logfile):
            #print(f'Logfile found for {sample}')
            try:
                id_logfile = pd.read_csv(id_logfile) 
                id_logfile = id_logfile.loc[id_logfile['sample'] == sample, : ]
                # check if the file is saved
                if sequence_source == 'personalized': # prediction must be present in two folders
                    queries_saved_h1 = [str(f'{output_dir}/{sample}/haplotype1/{query}_predictions.h5') for query in queries]
                    queries_saved_h2 = [str(f'{output_dir}/{sample}/haplotype2/{query}_predictions.h5') for query in queries]
                    queries_saved = np.array([os.path.isfile(q) for q in queries_saved_h1]) * np.array([os.path.isfile(q) for q in queries_saved_h2])
                elif sequence_source == 'reference':
                    queries_saved = [str(f'{output_dir}/{sample}/haplotype0/{query}_predictions.h5') for query in queries]
                    queries_saved = np.array([os.path.isfile(q) for q in queries_saved])
                # check if the file is in the logfile
                queries_logged = np.array([(query in id_logfile.region.values) for query in queries])
                # s
                if not queries_logged.shape[0] == queries_saved.shape[0]:
                    raise Exception("ERROR - Lengths of queries logged and saved conditions are not the same")
                queries_condition = queries_saved * queries_logged # true should not be predicted or logged; False should be
                queries_condition = queries_condition.tolist()
                # by default , log type is yes
                # id = queries.index('chr15_54740749_54740758')
                # print(queries_saved[id])
                # print(queries_logged[id])
                
                output = [{'query': queries[i], 'logtype':'y'} if qc is False else None for i, qc in enumerate(queries_condition)]

                # get indices that are nones
                none_ids = [i for i, q in enumerate(output) if q is None]

                # safely filter out the nones 
                #queries_saved = [q for i, q in enumerate(queries_saved) if i not in none_ids]
                queries_logged = [q for i, q in enumerate(queries_logged) if i not in none_ids]
                output = [q for i, q in enumerate(output) if i not in none_ids]
                #print(output)
                # assuming all predictions have not been logged hence logtype is y

                # change those regions where queries_logged == True to 'n'
                refined_output = [q_details.update({'logtype': 'n'}) if queries_logged[i] == True else q_details.update({'logtype': 'y'}) for i, q_details in enumerate(output)]

                return(output)

            except pd.errors.EmptyDataError:
                id_logfile = None
                output = [{'query':query, 'logtype':'y'} for query in queries]
                return(output)
        else:
            output = [{'query':query, 'logtype':'y'} for query in queries]
            return(output)
        
def check_predictions_and_logs(sample, predictions_folder, log_folder, interval_list_file, exclude_csv, sequence_source):

    import os
    import pandas as pd
    import numpy as np

    if not isinstance(sample, str):
        raise Exception(f'[ERROR] Sample argument should be a str of a valid sample name. You supplied a {type(sample).__name__} type')
    if not os.path.isdir(predictions_folder):
        raise Exception(f'[ERROR] Predictions folder does not exist. You supplied a {predictions_folder}')
    if not os.path.isdir(log_folder):
        raise Exception(f'[ERROR] Predictions log folder does not exist. You supplied a {log_folder}')
    if not isinstance(sequence_source, str):
        raise Exception(f'[ERROR] `sequence_source` argument should be a str of a valid sequence source (either `personalized` or `reference`). You supplied a {sequence_source}')

    if exclude_csv is None:
        exclude_these_regions = None
    elif os.path.isfile(exclude_csv):
        exclude_these_regions = pd.read_csv(exclude_csv)['region'].tolist()
    else:
        exclude_these_regions = None

    if os.path.isfile(interval_list_file):
        queries = pd.read_table(interval_list_file, sep=' ', header=None)[0].tolist()
        if exclude_these_regions is not None:
            queries = [q for q in queries if q not in exclude_these_regions]
    else:
        raise Exception(f'[ERROR] Interval list file {interval_list_file} does not exist.')

    # temporary solutions => remove nans
    queries = [q for q in queries if str(q) != 'nan']
    
    logged_queries = pd.read_csv(os.path.join(log_folder, f'{sample}_log.csv'))['region'].tolist()

    queries_logged = np.array([(query in logged_queries) for query in queries])
    
    if sequence_source == 'personalized': # prediction must be present in two folders
        queries_saved_h1 = [str(f'{predictions_folder}/{sample}/haplotype1/{query}_predictions.h5') for query in queries]
        queries_saved_h2 = [str(f'{predictions_folder}/{sample}/haplotype2/{query}_predictions.h5') for query in queries]
        queries_saved = np.array([os.path.isfile(q) for q in queries_saved_h1]) * np.array([os.path.isfile(q) for q in queries_saved_h2])
    elif sequence_source == 'reference':
        queries_saved = [str(f'{predictions_folder}/{sample}/haplotype0/{query}_predictions.h5') for query in queries]
        queries_saved = np.array([os.path.isfile(q) for q in queries_saved])

    if queries_logged.shape[0] != queries_saved.shape[0]:
        raise Exception("Lengths of queries logged and saved conditions are not the same")
    queries_condition = queries_saved * queries_logged # true should not be predicted or logged; False should be
    queries_condition = queries_condition.tolist()

    if all(queries_condition):
        message = f'SUCCESS - For {sample}, all predictions match all logged queries in the interval list files minus excluded regions, if any.'
        return({'logtype': 'info', 'logmessage':message, 'sample':sample})
    else:
        message = f'WARNING - For {sample}, either all predictions don\'t match all logged queries in the interval list files minus excluded regions or vice versa. This can happen if you have supplied a list of intervals but have chosen to predict on a subset. If this is the case, this behavior is normal. If you are unsure, please re-run the enformer prediction pipeline with the same parameters. You may supply a csv file of regions to exclude if available, but this should not matter.'
        return({'logtype': 'warning', 'logmessage':message, 'sample':sample})

def return_check_function(use_parsl, fxn=check_predictions_and_logs):
    '''
    Decorate or not the `run_batch_predictions` function based on whether `use_parsl` is true or false

    Returns: 
        function object
        The function if parsl is not used
        The parsl decorated function if parsl is meant to be used
    
    '''
    from parsl.app.app import python_app
    if use_parsl == True:
        return python_app(fxn)
    elif use_parsl == False:
        return fxn