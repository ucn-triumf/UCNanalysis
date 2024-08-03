# read UCN data or a set of UCN data
# Derek Fujimoto
# Aug 2024

import glob
from multiprocessing import cpu_count, Pool
from .ucndata import ucndata
from tqdm import tqdm
import numpy as np
from functools import partial

def read(path, header_only=False):
    """Read out single or multiple UCN run files from ROOT

    Args:
        path (str|list): path to file, may include wildcards, may be a list of paths which may include wildcards
        header_only (bool): if true, read only the header

    Returns:
        dict: keyed by run number, contains ucndata objects
    """

    # normalize input
    if type(path) is str:
        path = [path]

    # expand wildcards
    pathlist = []
    for p in path:
        pathlist.extend(glob.glob(p))

    # read out the data
    nproc = max(cpu_count()-1, 1)
    with Pool(nproc) as pool:
        fn = partial(ucndata, header_only=header_only)
        iterable = tqdm(pool.imap_unordered(fn, pathlist),
                        leave=False,
                        desc='Reading')

        data = np.fromiter(iterable, dtype=object)

    # sort result
    run_numbers = [d.run_number for d in data]
    idx = np.argsort(run_numbers)

    return data[idx]
