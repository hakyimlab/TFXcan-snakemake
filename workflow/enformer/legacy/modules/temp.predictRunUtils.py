
# Author: Temi
# DAte: Wed Sept 25 2024
# Description: This script is used for ENFORMER inference
# The following functions are available:
# - run_enformer
# - return_prediction_function


def run_enformer(filtered_check_result, batch_directives, module_directives):

    def enformer_predict_on_batch(batch_regions, samples, logging_dictionary, batch_directives, module_directives):
        import gc
        print("INFO - Garbage collection thresholds:", gc.get_threshold())

        class AttrDict(dict):
            def __init__(self, *args, **kwargs):
                super(AttrDict, self).__init__(*args, **kwargs)
                self.__dict__ = self

        batch_directives = AttrDict(batch_directives)
        module_directives = AttrDict(module_directives)

        import tensorflow as tf
        import os, sys, numpy as np, time

        sys.path.append(module_directives.path_to_modules)

        try:
            import predictionUtils, sequencesUtils, collectUtils, checksUtils, saveUtils, loggerUtils
        except (ImportError, ModuleNotFoundError) as ie:
            print(f"ERROR - {type(ie).__name__} at {__name__}. Cannot locate `...` modules.")

        if batch_directives.debugging == True:
            batch_directives.bins_indices = None
            batch_directives.tracks_indices = None
            batch_directives.predictions_expected_shape = (896, 5313)
            batch_directives.aggregate_by_width = False

        # this could mean
        # - check_queries returned nothing because
            # - nothing should be returned - good ; and it should return none
            # - something is wrong with check_queries ; and I should fix that
        if (not batch_regions) or (batch_regions is None):
            raise Exception(f'INFO - There are no regions to predict on in this batch {batch_directives.batch_number}.')

        #print(f'GPU Memory at start of batch {batch_num} predict function is {loggerUtils.get_gpu_memory()}')

        if batch_directives.grow_memory == True:
            gpus = tf.config.experimental.list_physical_devices('GPU')
            if gpus:
                try:
                    # Currently, memory growth needs to be the same across GPUs
                    for gpu in gpus:
                        tf.config.experimental.set_memory_growth(gpu, True)
                except RuntimeError as e:
                    raise Exception(f'RUNTIME ERROR - Batch {batch_directives.batch_number} of {type(e)} in module {__name__}')
        try:
            #model = predictionUtils.get_model(model_path)
            #model = enformer_model # check global definitions
            enformer_model = predictionUtils.get_model(batch_directives.model_path) # get the model
            fasta_extractor = sequencesUtils.get_fastaExtractor(batch_directives.fasta_file) # get the fasta file

            # == how to aggregate the predictions ==
            if batch_directives.aggregate == True:
                aggregation_dict = {'by_width': batch_directives.aggregation_width, 'by_function': batch_directives.aggregation_function}
            elif batch_directives.aggregate == False:
                aggregation_dict = None
            elif batch_directives.aggregate == None:
                aggregation_dict = None
            else:
                raise Exception(f'ERROR - Unknown aggregation type: {batch_directives.aggregate}')

            #print('Fasta and model successfully loaded')
            logger_output = []
            
            # #dlist = {sample[0]: batch_regions}
            # print(dlist)
            if batch_directives.batch_save == True:
                output = {}

            for sample in samples: # input_region is e.g. chr1_10_20
                # filter the samples to predict on
                for input_region in batch_regions:
                    if batch_directives.debugging == True:
                        v_samples = samples
                    else:
                        if not logging_dictionary:
                            v_samples = samples
                        elif logging_dictionary:
                            dlist = {sample: [elem['query'] for elem in logging_dictionary[sample]] for sample in logging_dictionary.keys()}
                            v_samples = checksUtils.return_samples_to_predict_on(query=input_region, logging_list_per_sample=dlist)

                tic = time.perf_counter()

                samples_enformer_inputs = sequencesUtils.create_input_for_enformer(query_region=input_region, samples=v_samples, path_to_vcf=batch_directives.path_to_vcf, fasta_func=fasta_extractor, hap_type = 'both', resize_for_enformer=True, resize_length=None, write_log=batch_directives.write_log, sequence_source=batch_directives.sequence_source, reverse_complement=batch_directives.reverse_complement, error_file = batch_directives.invalid_queries)

                toc = time.perf_counter()

                retrieve_time = (toc - tic)/len(v_samples) #len(list(samples_enformer_inputs['sequence'].keys()))

                #print(f'Region {input_region} sequences successfully created')

                # check that all the samples are accounted for
                #print(sorted(list(samples_enformer_inputs['sequence'].keys())))
                if samples_enformer_inputs is None:
                    logger_output.append(2)
                    continue
                elif samples_enformer_inputs is not None:
                    if samples_enformer_inputs['metadata']['sequence_source'] == 'var':
                        if sorted(v_samples) != sorted(list(samples_enformer_inputs['sequence'].keys())):
                            missing_samples = [s for s in sorted(v_samples) if s not in sorted(list(samples_enformer_inputs['sequence'].keys()))]
                            # remove missing samples from v_samples
                            if missing_samples:
                                print(f"WARNING - Removing {len(missing_samples)} missing samples from the input samples list.")
                                v_samples = [s for s in v_samples if s not in missing_samples]
                            print(f"WARNING - Some samples cannot be found. But this job will continue")
                    #logging_info_list = [] # collect all logging information here
                    sample_list = []
                    for sample in v_samples:
                        if samples_enformer_inputs['metadata']['sequence_source'] == 'var':
                            tic = time.perf_counter()

                            unfiltered_sample_predictions = predictionUtils.enformer_predict_on_sequence(model=enformer_model, sample_input=samples_enformer_inputs['sequence'][sample])
                        elif samples_enformer_inputs['metadata']['sequence_source'] in ['ref', 'random']:
                            tic = time.perf_counter()

                            unfiltered_sample_predictions = predictionUtils.enformer_predict_on_sequence(model=enformer_model, sample_input=samples_enformer_inputs['sequence'])
                        
                        toc = time.perf_counter()
                        predict_time = toc - tic
                        
                        # check that the predictions have the appropriate shapes
                        sample_predictions = {}
                        for hap in unfiltered_sample_predictions.keys():

                            haplo_prediction_cur = np.squeeze(unfiltered_sample_predictions[hap], axis=0)
                            #print(haplo_prediction_cur)
                            sample_predictions[hap] = collectUtils.collect_bins_and_tracks(haplo_prediction_cur, batch_directives.bins_indices, batch_directives.tracks_indices)
                            #print(sample_predictions[hap])

                            # should you aggregate or not:
                            if batch_directives.aggregate == True:
                                if batch_directives.aggregate_by_width == True:
                                    binwidth = collectUtils.slice_bins_for_width(locus=input_region, bin_size=128, nbins=896, padding=batch_directives.aggregation_width)
                                    sample_predictions[hap] = sample_predictions[hap][binwidth[0]:binwidth[1], :]
                                    if batch_directives.aggregation_function == 'aggBySum':
                                        sample_predictions[hap] = sample_predictions[hap].sum(axis=0).reshape(batch_directives.predictions_expected_shape)
                                        
                                    elif batch_directives.aggregation_function == 'aggByMean':
                                        sample_predictions[hap] = sample_predictions[hap].mean(axis=0).reshape(batch_directives.predictions_expected_shape)
                                    elif batch_directives.aggregation_function == None:
                                        pass
                                    else:
                                        raise Exception(f'ERROR - Unknown aggregation function: {batch_directives.aggregation_function}')
                                elif batch_directives.aggregate_by_width == False:
                                    sample_predictions[hap] = sample_predictions[hap][batch_directives.bins_indices, batch_directives.tracks_indices]
                                    if batch_directives.aggregation_function == 'aggBySum':
                                        sample_predictions[hap] = sample_predictions[hap].sum(axis=0).reshape(batch_directives.predictions_expected_shape)
                                    elif batch_directives.aggregation_function == 'aggByMean':
                                        sample_predictions[hap] = sample_predictions[hap].mean(axis=0).reshape(batch_directives.predictions_expected_shape)
                                    elif batch_directives.aggregation_function == None:
                                        pass
                                    else:
                                        raise Exception(f'ERROR - Unknown aggregation function: {batch_directives.aggregation_function}')

                            # now ensure that the shape of the predictions is as expected
                            if sample_predictions[hap].shape != batch_directives.predictions_expected_shape:
                                #print(sample_predictions[hap])
                                raise Exception(f'ERROR - {sample}\'s {hap} predictions shape is {sample_predictions[hap].shape} and is not equal to expected shape {batch_directives.predictions_expected_shape}.')
                        if batch_directives.batch_save == True:
                            output[input_region] = dict()
                            output[input_region][sample] = sample_predictions
                            sample_list.append(sample)

            # save for each sample

            if batch_directives.batch_save == True:
                if batch_directives.debugging == False:
                    outfile = batch_directives.prediction_logfiles[sample]
                    combined_predictions = saveUtils.combine_predictions_from_dictionary(output, separator = ':', prefix = '')
                    split_predictions = saveUtils.split_predictions_into_metadata_and_array(combined_predictions, batch_directives.batch_number, separator = ':')
                    completion_status = saveUtils.batch_write_predictions_to_hdf5(split_predictions, output_file = os.path.join(f'{batch_directives.output_file_prefix}.h5'), compress = True, aggregation = aggregation_dict)
                    ll = saveUtils.write_completion_status_to_csv(output_file = outfile, completion_status = completion_status)
                elif batch_directives.debugging == True:
                    return(output)
            elif batch_directives.batch_save == False:
                return(logger_output)
                    #print(f'INFO - Predictions for {input_region} have been successfully made.')
                    
            if batch_directives.write_log['logtypes']['memory']:
                mem_use = loggerUtils.get_gpu_memory() # [123, 456]#
                msg_mem_log = f"[MEMORY] (GPU) at the end of batch {batch_directives.batch_number} prediction: free {mem_use[0]} mb, used {mem_use[1]} mb on " 
                if tf.config.list_physical_devices('GPU'):
                    MEMORY_LOG_FILE = os.path.join(batch_directives.write_log['logdir'], "memory_usage.log")
                    loggerUtils.write_logger(log_msg_type = 'memory', logfile = MEMORY_LOG_FILE, message = msg_mem_log)
            if batch_directives.write_log['logtypes']['cache']:
                msg_cac_log = f'[CACHE] (model) at batch {batch_directives.batch_number}: [{predictionUtils.get_model.cache_info()}]'
                CACHE_LOG_FILE = os.path.join(batch_directives.write_log['logdir'], 'cache_usage.log')
                loggerUtils.write_logger(log_msg_type = 'cache', logfile = CACHE_LOG_FILE, message = msg_cac_log)

        except (TypeError, AttributeError) as tfe:
            if batch_directives.write_log['logtypes']['error']:
                if tf.config.list_physical_devices('GPU'):
                    mem_use = loggerUtils.get_gpu_memory()
                    err_mem_log = f"[ERROR] GPU memory error of type {type(tfe).__name__} for batch {batch_directives.batch_number}): free {mem_use[0]} mb, used {mem_use[1]} mb on {loggerUtils.get_gpu_name()}"
                    MEMORY_ERROR_FILE = os.path.join(batch_directives.write_log['logdir'], 'error_details.log')
                    loggerUtils.write_logger(log_msg_type = 'error', logfile = MEMORY_ERROR_FILE, message = err_mem_log)
            else:
                raise Exception(f'[ERROR] of {type(tfe).__name__} at batch {batch_directives.batch_number}')
        
    import time
    # with open(f'{batch_directives}') as f:
    #     parameters = yaml.safe_load(f)
    #     batch_directives = datatypeUtils.DirectivesHolder(**parameters)

    # filter out nones
    # filtered_check_result = [r for r in check_results if r is not None]
    if not filtered_check_result: # i.e. if the list is empty
        return(1)
    else:
        pqueries = [v for k, v in filtered_check_result.items()]
        pqueries = [l for l in pqueries for l in l]
        pqueries = list(set([d['query'] for d in pqueries]))
        # print(f'{len(pqueries)}')
        #print(pqueries)

        if pqueries:
            samples = list(filtered_check_result.keys())
            tic = time.perf_counter()

            reg_prediction = enformer_predict_on_batch(batch_regions=pqueries, samples=samples, logging_dictionary=filtered_check_result, batch_directives = batch_directives, module_directives = module_directives)
            
            toc = time.perf_counter()
            print(f"INFO - Time to predict on batch {batch_directives['batch_number']} is {toc - tic}")
            return(reg_prediction) # returns 0 returned by enformer_predict
        else:
            return(1)

def return_prediction_function(use_parsl, fxn=run_enformer):
    '''
    Decorate or not the `run_batch_predictions` function based on whether `use_parsl` is true or false

    Returns: 
        function object
        The function if parsl is not used
        The parsl decorated function if parsl is meant to be used
    
    '''
    import parsl
    # from parsl.app.app import python_app
    if use_parsl == True:
        return parsl.app.app.python_app(fxn)
    elif use_parsl == False:
        return fxn
    




# if batch_directives.batch_save == True:
#                                 output[input_region] = dict()
#                                 output[input_region][sample] = sample_predictions
#                         elif batch_directives.batch_save == False: 
#                             if batch_directives.debugging == False:
#                                 # otherwise, you can save the predictions ; prediction will be reshaped to (17, 5313) here
#                                 sample_logging_info = saveUtils.save_haplotypes_h5_prediction(haplotype_predictions=sample_predictions, metadata=samples_enformer_inputs['metadata'], output_dir=batch_directives.output_directory, sample=sample, aggregate_by_width=batch_directives.aggregate_by_width)

#                                 # check logging info/dictionary for the sample and the region
#                                 logging_type = checksUtils.return_sample_logging_type(sample=sample, query_region=input_region, logging_dictonary=logging_dictionary)
#                                 #print(logging_type)

#                                 if logging_type == 'y':
#                                     if (sample_logging_info is not None) and (len(sample_logging_info) == 4):
#                                         predictions_log_file = os.path.join(batch_directives.prediction_logfiles_folder, f'{sample}_log.csv')
#                                         sample_logging_info.extend([predict_time, retrieve_time])
#                                         logger_output.append(loggerUtils.log_predictions(predictions_log_file=predictions_log_file, what_to_write=sample_logging_info))
#                                     #print(f'Sample {sample} {input_region} haplotypes predictions have been logged.')
#                                 elif logging_type == 'n':
#                                     logger_output.append(1)
#                                 continue
#                             elif batch_directives.batch_save == True:
#                                 continue