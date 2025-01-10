

import numpy as np
import collectUtils


def combine_predictions_from_dictionary(dd, separator=':', prefix='') -> dict:
    # https://www.geeksforgeeks.org/python-convert-nested-dictionary-into-flattened-dictionary/
    """
    Input: dictionary of predictions: {region: {sample: {haplotype: 1 by 1 numpy matrix}}}
    output: dictionary of predictions: {metadata: pd.DataFrame, haplotype0: 3D numpy matrix}
    """
    res = {}
    for key, value in dd.items():
        if isinstance(value, dict):
            res.update(combine_predictions_from_dictionary(value, separator, prefix + key + separator))
        else:
            res[prefix + key] = value 
    res = {k: v for k,v in res.items() if v is not None}
    return(res)

def split_predictions_into_metadata_and_array(dd, batch_num, separator = ':') -> tuple:
    import numpy as np
    print(dd)
    if dd is None:
        raise Exception('ERROR - Predictions dictionary is empty')
    else:
        meta = [f'{d}{separator}batch_{batch_num}' for d in dd.keys()]
        data = np.stack([d for d in dd.values()], axis=0)

    return(meta, data)

def batch_write_predictions_to_hdf5(dd, output_file, compress = False, aggregation: dict=None) -> None:
    import os, h5py, numpy as np

    meta, data = dd

    # if aggregation is not None:
    #     if aggregation['by_width'] == False:
    #         if aggregation['by_function'] == 'aggBySum':
    #             data = data.sum(axis=2)
    #         elif aggregation['by_function'] == 'aggByMean':
    #             data = data.mean(axis=2)
    #         else:
    #             raise Exception(f'ERROR - Unknown aggregation function: {aggregation["by_function"]}')

    comp = 'gzip' if compress == True else None
    compopts = 9 if compress == True else None

    if os.path.isfile(output_file) == True:
        print(f'WARNING - {output_file} already exists; overwriting...') 
        os.remove(output_file)
    with h5py.File(output_file, 'w') as f:
        f.create_dataset('metadata', data=np.array(meta, dtype='S'), compression=comp, compression_opts=compopts)
        f.create_dataset('predictions', data=data, compression=comp, compression_opts=compopts)

    print(f'INFO - Predictions written to {output_file}')
    return(meta)


def write_completion_status_to_csv(output_file, completion_status) -> None:
    import os, csv

    open_mode = 'a' if os.path.isfile(output_file) else 'w'
    completion_status = [m.split(':') for m in completion_status]
    with open(output_file, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file, delimiter='\t')
        if open_mode == 'w':
            logwriter.writerow(['locus', 'sample', 'haplotype', 'batch']) # 
        logwriter.writerows(completion_status)

        running_log_file.flush()
        os.fsync(running_log_file)
    return(0)

def single_write_predictions_to_hdf5(haplotype_predictions, metadata, output_dir, sample, aggregate_by_width=True):

    import h5py, os, numpy, collectUtils

    region = metadata['region']
    if isinstance(aggregate_by_width, bool):
        if aggregate_by_width == True:
            bins_to_aggregate = collectUtils.slice_bins(region)
            print(f'INFO - aggregating by width: {bins_to_aggregate}')
        elif aggregate_by_width == False:
            bins_to_aggregate = [None, None]
            print(f'INFO - Not aggregating by width')
    elif isinstance(aggregate_by_width, int):
        bins_to_aggregate = collectUtils.slice_bins(region)
        bins_to_aggregate = [bins_to_aggregate[0] - aggregate_by_width, bins_to_aggregate[1] + aggregate_by_width]
        print(f'INFO - aggregating by width +/- {aggregate_by_width}: {bins_to_aggregate}')

    for key, values in haplotype_predictions.items():
        #print(f'[INFO] This is what is being saved {values.shape}')
        #print(values.shape)
        values = values[bins_to_aggregate[0]:bins_to_aggregate[1], : ].mean(axis=0)
        #print(values.shape)

        houtput = os.path.join(output_dir, sample, key)
        if not os.path.exists(houtput): os.makedirs(houtput, exist_ok=True)
        h5save = str(f'{houtput}/{region}_predictions.h5')
        with h5py.File(h5save, 'w') as hf:
            hf.create_dataset(region, data=values)

    output = [metadata['region'], sample, 'completed', metadata['sequence_source']]

    return(output)