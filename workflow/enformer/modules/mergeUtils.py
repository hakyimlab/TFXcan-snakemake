




def merge_hdf5_files(h5files: list, xshape: tuple, output_file: str, overwrite = True, compress = True):
    import h5py, os
    import numpy as np
    if os.path.isfile(output_file) == True and overwrite == True:
        #print(f'WARNING - {output_file} already exists; overwriting...') 
        os.remove(output_file)

    meta = []
    comp = 'gzip' if compress == True else None
    compopts = 9 if compress == True else None

    with h5py.File(output_file, mode='w') as h5fw:
        h5fw.create_dataset('predictions', dtype="f", shape=xshape, maxshape=(None, xshape[1], xshape[2]), compression=comp, compression_opts=compopts)
        #
        row0 = 0
        for h5name in h5files:
            if os.path.isfile(h5name) == False:
                print(f'ERROR - {h5name} does not exist; skipping...')
                continue
            else:
                #print(f'INFO - Processing {h5name}...')
                h5fr = h5py.File(h5name, 'r') 
                # dset1 = list(h5fr.keys())
                arr_metadata = h5fr['metadata'][:]
                arr_data = h5fr['predictions'][:]
                h5fr.close()
                meta.append(arr_metadata)
                if row0 == 0:
                    h5fw['predictions'][row0:arr_data.shape[0],:,:] = arr_data[:]
                else:
                    h5fw['predictions'][row0:row0+arr_data.shape[0],:,:] = arr_data[:]
                row0 += arr_data.shape[0]
        h5fw.create_dataset('metadata', data = np.array(np.concatenate(meta), dtype = 'S'), compression=comp, compression_opts=compopts)
    return(0)          