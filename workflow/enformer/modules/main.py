
# Author: Temi
# DAte: Wed Sept 25 2024
# Description: This script is used for ENFORMER inference
# The following functions are available:
# - main

import os, sys, yaml, re, parsl, time, dataclasses
import pandas as pd, numpy as np
import directives, datatypeUtils, collectUtils, predictionUtils, batchUtils, checksUtils, datatypeUtils, mergeUtils

# main 
def main(args, script_path):

    params_path = args.parameters

    if not os.path.isabs(params_path):
        params_path = os.path.abspath(params_path)
    
    print(f'INFO - Using parameters from {params_path}')

    # I need mainRun.py downstream and and it is quite tricky to use with parsl
    # p_two = os.path.join(script_path, 'modules', 'mainRun.py')

    with open(f'{params_path}') as f:
        parameters = yaml.safe_load(f)
        prediction_data_name = parameters['prediction_data_name']
        run_date = parameters['date'] #if parameters['date'] is not None else date.today().strftime("%Y-%m-%d")
        output_dir = os.path.join(os.path.abspath(parameters['output_dir']), f'{prediction_data_name}', f'predictions_{run_date}')

        interval_list_file = parameters['interval_list_file']
        predictions_log_dir = os.path.join(output_dir, parameters['write_log']['logdir'])
        predictions_dir = os.path.join(output_dir, 'predictions')

        if not os.path.isdir(predictions_dir):
            os.makedirs(predictions_dir)
        if not os.path.isdir(predictions_log_dir):
            os.makedirs(predictions_log_dir)

        job_log_dir = os.path.join(output_dir, parameters['write_log']['logdir'])
        n_regions = parameters["n_regions"]
        batch_regions = int(parameters['batch_regions'])
        use_parsl = parameters['use_parsl']
        parsl_parameters = parameters['parsl_parameters']
        sequence_source = parameters['sequence_source']
        exclude_regions = parameters["exclude_regions"]
        reverse_complement = parameters["reverse_complement"]
        metadata_dir = os.path.join(output_dir, 'metadata')
        if not os.path.isdir(metadata_dir):
            os.makedirs(metadata_dir)

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        if int(n_regions) == -1:
            n_regions = None
        elif int(n_regions) > 0:
            n_regions = (n_regions) if isinstance(n_regions, int) else None

        # personalized parameters 
        individuals = parameters['individuals'] if sequence_source == 'personalized' else None
        vcf_files_dict = parameters['vcf_files'] if sequence_source == 'personalized' else None
        invalid_queries = os.path.join(metadata_dir, 'invalid_queries.txt')

        if sequence_source == 'personalized':
             # use only the chromosomes that have been made available in the config file vcf params
            print(f'INFO - Sequence source is {sequence_source}. Using a reference genome + vcf files.')
            chromosomes = list(vcf_files_dict['files'].keys())

            batch_individuals = parameters["batch_individuals"]
            n_individuals = int(parameters['n_individuals'])
        # list of chromosomes (if the sequence source is reference)
        elif sequence_source == 'reference':
            print(f'INFO - Sequence source is {sequence_source}. Using a reference genome.')
            chromosomes = [f'chr{i}' for i in range(1, 23)]
            chromosomes.extend(['chrX'])

        if reverse_complement:
            print(f'INFO - Predicting on reverse complements too')

    # modify parsl parameters to add the working directory
    parsl_parameters['working_dir'] = output_dir

    if not os.path.isdir(job_log_dir):
        os.makedirs(job_log_dir)

    # set parsl directives
    directives.parsl_directives(use_parsl, parsl_parameters)
    
    # importing this module does not work; best to execute it here
    # predictionUtilities = os.path.join(script_path, 'modules', 'predictionRunUtils.py')
    # exec(open(predictionUtilities).read(), globals(), globals())

    exec(open(os.path.join(script_path, 'modules', 'predictionRunUtils.py')).read(), globals())
    exec(open(os.path.join(script_path, 'modules', 'datatypeUtils.py')).read(), globals())
    exec(open(os.path.join(script_path, 'modules', 'predictionUtils.py')).read(), globals())
    #print(os.path.isfile(os.path.join(script_path, 'modules', 'predictionRunUtils.py')))

    # decorate the prediction function with or without parsl
    prediction_fxn = return_prediction_function(use_parsl)

    # save the modules to pass to a parsl
    module_holder = dict(path_to_modules=os.path.join(script_path, 'modules'), predictionUtils=os.path.join(script_path, 'modules', 'predictionUtils'), checksUtils=os.path.abspath(checksUtils.__file__))

    # determine what individuals to predict on and all that
    if sequence_source == 'personalized':
        
        if isinstance(individuals, list):
            id_list = individuals
            pass
        elif isinstance(individuals, type('str')):
            if os.path.isfile(individuals):
                if n_individuals == -1:
                    id_list = pd.read_table(individuals, header=None)[0].tolist()[0:]
                elif n_individuals > 0:
                    id_list = pd.read_table(individuals, header=None)[0].tolist()[0:(n_individuals)]
            else:
                id_list = [individuals]
        print(f'INFO - Found {len(id_list)} individuals to predict on')

    elif sequence_source == 'reference':
        id_list = [prediction_data_name]
        print(f'INFO - Found one reference set named {id_list[0]} to predict on')
    elif sequence_source == 'random':
        id_list = [prediction_data_name]
        print(f'INFO - Prediction will be on a randomly generated set')

    # set log files to be put in a folder and touch the log files per sample
    prediction_logfiles_folder = predictions_log_dir
    if not os.path.isdir(prediction_logfiles_folder):
        os.makedirs(prediction_logfiles_folder)
    prediction_logfile = os.path.join(prediction_logfiles_folder, f"{prediction_data_name}.prediction_log.txt")
    if parameters['write_log']['overwrite'] == True:
        if os.path.isfile(prediction_logfile):
            os.remove(prediction_logfile)
        
    # list of intervals to be predicted on
    a = pd.read_table(interval_list_file, sep=' ', header=None).dropna(axis=0) #.drop_duplicates(subset=['region', 'sample', 'status', 'sequence_source'], keep='last')
    list_of_regions = a[0].tolist()[0:(n_regions)] # a list of queries
    print(f'INFO - Found {len(list_of_regions)} regions to be split into batches with at most {batch_regions} regions in each batch.')

    # filter the list of chromosomes to be compatible with the available regions
    chromosomes = list(set([r.split('_')[0] for r in list_of_regions]))
    #print(f'INFO - Chromosomes to predict on are: {chromosomes}')

    # should some regions be excluded?
    if exclude_regions == True:
        # seach for the invalid_regions.csv file
        exclude_file = os.path.join(job_log_dir, 'invalid_queries.csv')
        if os.path.isfile(exclude_file):
            exclude_these_regions = pd.read_csv(exclude_file)['region'].tolist()
            print(f'INFO - Found regions to be excluded from the input regions.')
            list_of_regions = [l for l in list_of_regions if l not in exclude_these_regions]  
            print(f'INFO - Updated number of regions to predict on is {len(list_of_regions)}')
        else:
            print(f'INFO - No regions to exclude yet. You either did not supply a file, this is the first run, or there are truly no regions to exclude')
            exclude_these_regions = None
    else:
        exclude_file = None
    
    # batch the samples too
    # if you have 1000 individuals, it may be too much
    if len(id_list) > 5:
        if batch_individuals is not None:
            if isinstance(batch_individuals, int):
                sample_batches = list(batchUtils.generate_batch_n_elems(id_list, n = batch_individuals)) # 5 samples in each batch
                print(f'INFO - Predictions will be done for every {batch_individuals} individuals.')
            else:
                raise Exception(f'ERROR - argument `batch_individuals` is not a str type. You supplied a {type(batch_individuals).__name__}')
        else:
            print(f'INFO - You have multiple individuals/samples and have not supplied how to batch them. For efficient use of resources, use the `batch_individuals` argument.')
    else:
        sample_batches = [id_list] # put the list in a list
        print(f'INFO - There seem to be just one sample i.e. {sample_batches}. No need to batch.')

    
    # ======== WHAT BINS AND TRACKS SHOULD BE SAVED? what is the expected shape of each prediction ========
    #bins_indices, tracks_indices = collectUtils.parse_bins_and_tracks(parameters['bins_to_save'], parameters['tracks_to_save'])
    # Check prediction size for correctness
    predictions_expected_shape = collectUtils.calculate_expected_shape(parameters["aggregation"])
    print(f'INFO - Expected shape of predictions is {predictions_expected_shape}')
    bins_indices, tracks_indices = parameters['aggregation']['bins_to_save'], parameters['aggregation']['tracks_to_save']
    by_width = True if isinstance(parameters['aggregation']['by_width'], int) else False
    width = parameters['aggregation']['by_width']
    if by_width:
        print(f"INFO - Aggregating by width, and padding with +/-{width}")
    by_function = 'aggByMean' if (parameters['aggregation']['by_function'] is None) or ("by_function" not in parameters['aggregation'].keys()) else parameters['aggregation']['by_function']

    # ==== to make this fast, pass multiple regions to one parsl app ======
    sample_app_futures = []
    count = 0
    output_files_list = []
    for sample_list in sample_batches:
        #print(f'INFO - Predicting on {len(sample_list)} samples')
        for chromosome in chromosomes:
            #print(f'INFO - Predicting on {chromosome}')
            #print(chromosome)
            chr_list_of_regions = [r for r in list_of_regions if r.startswith(f"{chromosome}_")]
            if sequence_source == 'personalized':
                if chromosome not in vcf_files_dict['files'].keys():
                    print(f'WARNING - {chromosome} VCF is not available.')
                    continue
                else:
                    chr_vcf_file = os.path.join(vcf_files_dict['folder'], vcf_files_dict['files'][chromosome])
            elif sequence_source == 'reference':
                chr_vcf_file = None

            if not chr_list_of_regions:
                print(f'WARNING - {chromosome} sites are not available.')
                continue

            # I want many regions to be put in a parsl app
            if len(chr_list_of_regions) > batch_regions:
                region_batches = batchUtils.generate_batch_n_elems(chr_list_of_regions, n=batch_regions) # batch_regions total batches
            else:
                region_batches = [chr_list_of_regions]
            
            for region_list in region_batches:
                #print(len(sample_list))
                #print(f'{len(region_list)} regions in {chromosome} for {len(sample_list)} samples')
                # sample_app_futures.append(prediction_fxn(batch_regions=list(region_list), samples=list(sample_list), path_to_vcf = chr_vcf_file, batch_num = count, script_path=script_path, output_dir=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, sequence_source=sequence_source, tmp_config_path=params_path, p_two=p_two))
                # 
                output_prefix = os.path.join(predictions_dir, f'{prediction_data_name}.{run_date}.batch_{count}')
                output_files_list.append(f'{output_prefix}.h5')
                batch_holder_dict = datatypeUtils.DirectivesHolder(use_parsl=use_parsl, batch_regions=list(region_list), 
                                                                 samples=list(sample_list), path_to_vcf=chr_vcf_file, batch_number=count, 
                                                                 script_path=script_path, output_directory=output_dir, prediction_logfiles_folder=prediction_logfiles_folder, sequence_source=sequence_source, debugging=False, grow_memory=True, write_log=parameters['write_log'], reverse_complement=reverse_complement, mainRunScript=None, tmp_config_path=params_path, 
                                                                 bins_indices = bins_indices, tracks_indices = tracks_indices, run_date=run_date, predictions_expected_shape=predictions_expected_shape, model_path=parameters['model_path'], fasta_file=parameters['fasta_file'],
                                                                 aggregate = parameters["aggregation"]["aggregate"], aggregate_by_width = by_width, aggregation_width = width, aggregation_function = by_function, output_file_prefix=output_prefix, prediction_logfile = prediction_logfile, invalid_queries=invalid_queries)
                batch_holder_dict = dataclasses.asdict(batch_holder_dict)
                queries, queries_directives = check_for_batches_to_run(batch_dictionary=batch_holder_dict, module_directives = module_holder)
                # print(queries)
                # print(queries_directives)
                sample_app_futures.append(prediction_fxn(queries, queries_directives, module_holder))
                count = count + 1

    if use_parsl == True:
        print(f'INFO - Executing {len(sample_app_futures)} parsl apps')
        # for q in sample_app_futures:
        #     print(q.result())
        exec_futures = [q.result() for q in sample_app_futures] 
        #print(sample_app_futures)
        #print(f'INFO - Finished predictions for all')
    elif use_parsl == False:
        print(f'INFO - Finished predictions for: {len(sample_app_futures)} batches')
    
    # == After predictions are complete, a yaml file will be written out to help with aggregation
    print(f'INFO - Writing `{prediction_data_name}.aggregation_config.yaml` file to {metadata_dir}')
    agg_dt = {'predictions_directory': f"{os.path.join(output_dir, 'predictions')}", 'predictions_logfile':prediction_logfile, 
              'expected_shape': predictions_expected_shape, 'prediction_data_name': prediction_data_name, 'sequence_source': sequence_source, 
              'run_date':run_date, 'individuals': None if sequence_source in ['reference', 'random'] else individuals, 
              'n_individuals': n_individuals if sequence_source == 'personalized' else None}

    mfile = os.path.join(metadata_dir, f'{prediction_data_name}.aggregation_config.yaml')
    with(open(mfile, mode='w')) as wj:
        yaml.dump(agg_dt, wj)







    # remove temporatry config file
    # print(f"INFO - Cleaning up: Removing temporary config file at {tmp_config_file}")
    # os.remove(tmp_config_file)


# # just so I don't have to deal with having too many resources, I can request a small amount of resource
    # check_fxn = return_check_function(use_parsl)
    # SUMMARY_FILE = os.path.join(job_log_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.summary')
    # summary_exec = []
    # for sample in id_list:
    #     if os.path.isfile(os.path.join(prediction_logfiles_folder, f"{sample}_log.csv")):
    #         summary_exec.append(check_fxn(sample=sample, predictions_folder=output_dir, log_folder=prediction_logfiles_folder, interval_list_file=interval_list_file, exclude_csv=exclude_file, sequence_source=sequence_source))

    # if use_parsl:
    #     summary_exec = [q.result() for q in summary_exec]
    #     parsl.clear() # end parsl

    # #summary_exec = list(set(summary_exec))
    # for i, qr in enumerate(summary_exec):
    #     loggerUtils.write_logger(log_msg_type=qr['logtype'], logfile=SUMMARY_FILE, message=qr['logmessage'])

    # # regex the summary file and save the failed ones e.t.c to csv
    # # --- there is a better way to do this but for now, this will do

    # warning_pattern = r"^\[WARNING.*For\s(\w+|\d+).*"
    # success_pattern = r"^\[INFO.*For\s(\w+|\d+).*"
    # with open(SUMMARY_FILE, 'r') as f:
    #     lines = list(set(f.readlines()))
    # # print(line)
    # warning_result = [re.search(warning_pattern, l).group(1) for l in lines if not re.search(warning_pattern, l) is None]
    # success_result = [re.search(success_pattern, l).group(1) for l in lines if not re.search(success_pattern, l) is None]
    # pd.DataFrame(list(set(warning_result))).to_csv(os.path.join(metadata_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.unsuccessful_predictions.csv'), index=False, header=False)
    # pd.DataFrame(list(set(success_result))).to_csv(os.path.join(metadata_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.successful_predictions.csv'), index=False, header=False)

    # collect the successfule predictions
    # successful_predictions = list(set([q['sample'] for q in summary_exec if q['logtype'] == 'INFO']))
    # unsuccessful_predictions = list(set([q['sample'] for q in summary_exec if q['logtype'] == 'WARNING']))
    # pd.DataFrame({'successful_predictions':successful_predictions}).to_csv(os.path.join(metadata_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.successful_predictions.csv'), index=False, header=False)
    # pd.DataFrame({'unsuccessful_predictions':unsuccessful_predictions}).to_csv(os.path.join(metadata_dir, f'{prediction_data_name}_{prediction_id}_{run_date}.unsuccessful_predictions.csv'), index=False, header=False)
