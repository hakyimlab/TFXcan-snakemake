# Author: Temi
# DAte: Wed Sept 25 2024
# Description: This script is used for ENFORMER inference
# The following functions are available:
# - slice_bins_for_width
# - collect_bins_and_tracks
# - calculate_expected_shape


import numpy as np

def slice_bins_for_width(locus, bin_size=128, nbins=896, padding=0):
    import math
    import numpy as np
    locus = locus.split('_')
    start = int(locus[1])
    end = int(locus[2])
    midn = math.ceil((start + end) / 2)
    nstart = midn - ((bin_size*nbins) / 2) # (128*896) / 2
    sstart = (start - nstart)
    send = sstart + ((end - start) - 1)
    bins = list(range(0, nbins*bin_size, bin_size))
    out = np.digitize([sstart, send], bins=bins).tolist()

    if((end - start) <= 128):
        pass
    elif(((end - start) > 128) or (end - start) % bin_size > 0):
        out[1] = out[1] + 1

    # if [448], make it [448, 449]
    # if [448, 448] make it [448, 449]
    
    if len(out) == 1:
        out.append(out[0] + 1)
    elif len(out) == 2:
        if out[0] == out[1]:
            out[1] = out[1] + 1

    out[0] = out[0] - padding
    out[1] = out[1] + padding
    return(out)

def collect_bins_and_tracks(predictions, bins_indices, tracks_indices):

    # print(predictions.shape)
    # print(bins_indices)
    # print(tracks_indices)

    # print(type(predictions))
    # print(type(bins_indices))
    # print(type(tracks_indices))

    if bins_indices is not None and tracks_indices is not None:
        return(predictions[bins_indices, :][:, tracks_indices])
    else:
        if bins_indices is None:
            if tracks_indices is None:
                return(predictions[:, :])
            else:
                return(predictions[:, tracks_indices])
        else:
            if tracks_indices is None:
                return(predictions[bins_indices, :])
    
    # if bins_indices is None:
    #     if tracks_indices is None:
    #         return(predictions[:, :])
    #     else:
    #         return(predictions[:, tracks_indices])
    # else:
    #     if tracks_indices is None:
    #         return(predictions[bins_indices, :])
    #     else:
    #         print(predictions.ix_(bins_indices, tracks_indices))
    #         return(predictions.ix_(bins_indices, tracks_indices))
        
def calculate_expected_shape(params):
    if params['aggregate'] == True:
        #params["by_function"] = 'aggByMean' if (params["by_function"] is None) or ("by_function" not in params.keys()) else params["by_function"]
        if params["by_width"] == True:
            if isinstance(params["tracks_to_save"], list):
                return((1, len(params['tracks_to_save'])))
            elif params["tracks_to_save"] is None:
                return((1, 5313))
        elif (params["by_width"] is None) or params["by_width"] == False:
            if isinstance(params['tracks_to_save'], list):
                return((1, len(params['tracks_to_save'])))
            elif params['tracks_to_save'] is None:
                return((1, 5313))
        else:
            return((1, 5313))
    elif params['aggregate'] == False:
        if params['tracks_to_save'] is None:
            if params['bins_to_save'] is None:
                return((896, 5313))
        elif isinstance(params['bins_to_save'], list):
            if isinstance(params['tracks_to_save'], list):
                return((len(params['bins_to_save']), len(params['tracks_to_save'])))
            elif params['tracks_to_save'] is None:
                return((len(params['bins_to_save']), 5313))
    


# import numpy as np
# x = np.array(((1,2,3),(4,5,6)))

# b_inds = "0-1"
# t_inds = "2-2"

# bins_indices,tracks_indices = parse_bins_and_tracks(b_inds,t_inds)
# out = collect_bins_and_tracks(x,bins_indices,tracks_indices)
# print(out)




# def parse_bins_and_tracks(unparsed_bins, unparsed_tracks):
#     def split_range(range_str):
#         if "-" in range_str:
#             spl_range_str = [int(x) for x in range_str.split("-")]
#             return [x for x in range(spl_range_str[0],spl_range_str[1]+1)]
#         else:
#             return [int(range_str)]
    
#     if unparsed_bins == -1:
#         bins_indices = None
#     else:
#         bins_indices = []
#         bins_split = unparsed_bins.split(",")
#         for b in bins_split:
#             bins_indices += split_range(b)
#         bins_indices.sort()
    
#     if unparsed_tracks == -1:
#         tracks_indices = None
#     else:
#         tracks_indices = []
#         tracks_split = unparsed_tracks.split(",")
#         for t in tracks_split:
#             tracks_indices += split_range(t)
#         tracks_indices.sort()

#         # with open(mapping_path,"w") as mp:
#         #     mp.write("stored_track_index\tmodel_track_index\n")
#         #     for i,x in enumerate(tracks_indices):
#         #         mp.write("\t".join([str(i),str(x)])+"\n")

#     return bins_indices,tracks_indices