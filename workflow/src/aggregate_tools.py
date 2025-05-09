# Author: Temi
# Usage: python aggregate_tools.py, a module 
# Mofified: July 5: to aggregate single predictions
# Date: Mon Apr 17 2021


import pandas as pd
import numpy as np

def aggregate_enformer_predictions(each_id, log_data, predictions_path, prediction_id, save_dir, agg_types, sequence_source, batch_num=None):


    import h5py
    import numpy as np
    import os, sys
    import pandas as pd
    import multiprocessing
    import itertools

    # these are the bins/ positions
    # upstream = list(range(0, 7)) # or 0 to 7
    # center = 8
    # pre_center = 7
    # post_center = 9
    mean_center=[7,8,9]
    #global read_file

    # can aggregate by the mean of all bins, mean of the upstream and/or downstream alone, or just select the center
    def agg_by_mean(pred_tracks, use_bins=None):

        y = []
        X = []

        for k, v in pred_tracks.items():
            y.append(k) #if k.startswith('pos') else y.append(0)

            if isinstance(use_bins, type(None)):
                v = v.mean(axis=0)
            elif isinstance(use_bins, type([])):
                v = v[use_bins, :].mean(axis=0)
            v = np.expand_dims(v, axis=1).T
            X.append(v)

        y = np.expand_dims(np.array(y), axis=1)
        dt = np.hstack((y, np.vstack(X)))
        return dt

    def agg_by_center(pred_tracks, center):

        y = []
        X = []

        for k, v in pred_tracks.items():
            y.append(k)
            v = v[center, :]
            v = np.expand_dims(v, axis=1).T
            X.append(v)

        y = np.expand_dims(np.array(y), axis=1)
        dt = np.hstack((y, np.vstack(X)))

        return dt

    def agg_by_collect(pred_tracks, use_bins=None):
        y = []
        X = []
        for k, v in pred_tracks.items():
            y.append(k) #if k.startswith('pos') else y.append(0)
            if isinstance(use_bins, type(None)):
                v = v.reshape((1, 5313))
            X.append(v)
        y = np.expand_dims(np.array(y), axis=1)
        dt = np.hstack((y, np.vstack(X)))
        #dt = np.vstack(X)
        return dt
    
    def main(regions_list, predictions_path, haplotype): 
        global read_file

        def read_file(region, predictions_dir, each_haplotype):
            import os
            import h5py
            import numpy as np

            output = {}
            # print(f'INFO - Reading predictions for {region}')
            # print(f'INFO - Reading predictions for {each_haplotype}')
            # print(f'INFO - Reading predictions for {predictions_dir}')
            fle = os.path.join(predictions_dir, each_haplotype, f'{region}_predictions.h5')  #f'{dir}//{region}_predictions.h5'
            if os.path.isfile(fle):
                with h5py.File(fle, 'r') as f:
                    filekey = list(f.keys())[0]
                    output[region] = np.vstack(list(f[filekey]))
            else:
                print(f'ERROR - {region} predictions file does not exist.')
            return(output)
        iterargs = itertools.product(regions_list, [predictions_path], [haplotype])
        out = [read_file(*args) for args in iterargs]
        return(out)
        
        # with multiprocessing.Pool(4) as pool:
        #     #print(itertools.product(regions_list, [predictions_path], [haplotype_list]))
        #     out = pool.starmap(read_file, itertools.product(regions_list, [predictions_path], [haplotype]))
        #     return(out)

    if each_id in ['kawakami', 'cistrome']:
        pred_type = 'ref'
        haplotypes = ['haplotype0']
    elif each_id in ['random']:
        pred_type = 'random'
        haplotypes = ['haplotype0']
    elif each_id in ['personalized']:
        pred_type = 'var'
        haplotypes  = ['haplotype1', 'haplotype2']
    elif sequence_source == 'personalized':
        pred_type = 'var'
        haplotypes  = ['haplotype1', 'haplotype2']


    print(f'INFO - Seeing {multiprocessing.cpu_count()} CPUs')
    print(f'INFO - Starting to collect predictions for {each_id}')

    
    if __name__ == '__main__':
        
        regions_list = log_data.loc[log_data['sequence_source'] == pred_type, ].region.values.tolist()
        print(f'{regions_list[0:2]}...')
        
        pooled_dictionary = {}
        for haplotype in haplotypes:
            outputs_list = main(regions_list, predictions_path, haplotype)
            print(f'INFO - Reading predictions for {haplotype}')
            #outputs_list = read_file(regions_list, predictions_path, haplotype)
            # print(outputs_list)

            # with multiprocessing
            # pool = multiprocessing.Pool(4)
            # outputs_list = pool.starmap(read_file, itertools.product(regions_list, [predictions_path], [haplotype])) #pool.map(read_file, regions_list) # 'haplotype1

            # no multiprocessing
            # iterargs = itertools.product(regions_list, [predictions_path], [haplotype])
            # outputs_list = [read_file(*args) for args in iterargs]

            predictions =  {k: v for d in outputs_list for k, v in d.items()}
            pooled_dictionary[haplotype] = predictions
            print(f'INFO - Finished reading predictions for {haplotype}')

            #print({k: predictions[k] for k in predictions.keys()[:3]})

        if len(haplotypes) == 2:
            # check that the keys match
            match_condition = sorted(list(pooled_dictionary[haplotypes[0]].keys())) == sorted(list(pooled_dictionary[haplotypes[1]].keys()))
            if not match_condition:
                raise Exception(f'ERROR - Fatal: Haplotypes 1 and 2 regions are different. Haplotype1 length is {len(list(pooled_dictionary[haplotypes[0]].keys()))} and Haplotype2 length is {len(list(pooled_dictionary[haplotypes[1]].keys()))}')
            else:
                # one one of them
                ids = list(pooled_dictionary[haplotypes[0]].keys())
                summed_pooled_dictionary = {id: np.add(pooled_dictionary[haplotypes[0]][id], pooled_dictionary[haplotypes[1]][id]) for id in ids}

        elif len(haplotypes) == 1:
            summed_pooled_dictionary = pooled_dictionary[haplotypes[0]]
        print(f'INFO - Successfully read all files into pooled dictionary.')
        
        for agg_type in agg_types:
            try:
                data_dict = {}
                if agg_type == 'aggByMean': data_dict[agg_type] = agg_by_mean(summed_pooled_dictionary)
                if agg_type == 'aggByMeanCenter': data_dict[agg_type] = agg_by_mean(summed_pooled_dictionary, use_bins=mean_center)
                if agg_type == 'aggByCenter': data_dict[agg_type] = agg_by_center(summed_pooled_dictionary, center=center)
                if agg_type == 'aggByPreCenter': data_dict[agg_type] = agg_by_center(summed_pooled_dictionary, center=pre_center)
                if agg_type == 'aggByPostCenter': data_dict[agg_type] = agg_by_center(summed_pooled_dictionary, center=post_center)
                if agg_type == 'aggByUpstream': data_dict[agg_type] = agg_by_mean(summed_pooled_dictionary, use_bins=upstream)
                if agg_type == 'aggByDownstream': data_dict[agg_type] = agg_by_mean(summed_pooled_dictionary, use_bins=downstream)
                if agg_type == 'aggByUpstreamDownstream': data_dict[agg_type] = agg_by_mean(summed_pooled_dictionary, use_bins=upstream + downstream)
                if agg_type == 'aggByCollect': data_dict[agg_type] = agg_by_collect(summed_pooled_dictionary)
            
            except ValueError as ve:
                raise ValueError(f'ERROR - Problem with arrays for {each_id}, {agg_type}')

            #ty = pd.concat([pd.Series(list(predictions.keys())), pd.DataFrame(dt)], axis=1)
            ty = pd.DataFrame(data_dict[agg_type])
            print(f'INFO - Dimension of collected data is {ty.shape[0]} by {ty.shape[1]}')

            column_names = ['id']
            column_names.extend([f'f_{i}' for i in range(1, ty.shape[1])])

            #print(len(column_names))
            ty = ty.set_axis(column_names, axis=1)
            # print(ty.iloc[0:5, 0:5])

            print(f'INFO - Saving file to {save_dir}/{each_id}_{agg_type}_{prediction_id}.csv.gz')
            if batch_num is None:
                # writing to a gz file results in several issues
                outfile = f'{save_dir}/{each_id}_{agg_type}_{prediction_id}.csv.gz'
                ty.to_csv(path_or_buf=outfile, index=False, compression='gzip')
                print(f'INFO - Finished saving data for {each_id}')
            else:
                ty.to_csv(path_or_buf=f'{save_dir}/{each_id}_{agg_type}_{prediction_id}_batch_{batch_num}.csv.gz', index=False, compression='gzip')
                print(f'INFO - Finished saving data for {each_id} for batch {batch_num}')
            
    return(0)

def return_prediction_function(use_parsl, fxn=aggregate_enformer_predictions):
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

def generate_batch(lst, batch_n, len_lst = None):
    """
    Given a list, this function yields batches of an unspecified size but the number of batches is equal to `batch_n`
    E.g. generate_batch([0, 1, 2, 3, 4, 5, 6], batch_n=2) -> (0, 1, 2, 3), (4, 5, 6)
    
    Parameters:
        lst: list
        batch_n: int
            Number of batches to return
        len_lst: None or num (length of the input list)
    Yields
        `batch_n` batches of the list
    """
    import math
    # how many per batch
    if len_lst is not None:
        n_elems = math.ceil(len_lst/batch_n)
    else:
        n_elems = math.ceil(len(lst)/batch_n)
    for i in range(0, len(lst), n_elems):
        yield lst[i:(i + n_elems)]