# Author: Temi
# DAte: Wed Sept 25 2024
# Description: This script is used for ENFORMER inference
# The following functions are available:
# - batch
# - generate_n_batches
# - generate_batch_n_elems

def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]

def generate_n_batches(lst, batch_n, len_lst = None):
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

def generate_batch_n_elems(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]