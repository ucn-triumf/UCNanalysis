# Merge a list of ucndata objects into a single object
# Derek Fujimoto
# Aug 2024

from .udata import udata
from rootloader import tfile, ttree, th1, th2
import numpy as np
import pandas as pd

def merge(ucnlist):
    """Merge a list of ucndata objects into a single object

    Args:
        ucnlist (list): iterable of ucndata objects

    Returns:
        ucndata: single object with all data inside of it
    """

    # sort by run number, assume run numbers are in chronological order
    ucnlist = np.array(ucnlist)
    run_num = [d.run_number for d in ucnlist]
    idx = np.argsort(run_num)
    ucnlist = ucnlist[idx]

    # initialize output object
    ucnmerged = ucndata(None)

    # set header items as arrays
    for key in ucnlist[0].__dict__.keys():

        # skip some things
        if key in ('tfile', 'cycle'):
            continue

        # make and set arrays
        values = np.array([getattr(data, key) for data in ucnlist])
        setattr(ucnmerged, key, values)

    # set cycle
    ucnmerged.cycle = None

    obj_names = [] # list of all contained objects
    for dat in ucnlist:

        # convert tfiles to dataframes
        dat.tfile.to_dataframe()

        # get names of objects
        obj_names.extend(list(dat.tfile.keys()))

    # get unique object names
    obj_names = np.unique(obj_names)

    # set tfile
    ucnmerged.tfile = tfile(None)

    # merge
    for name in obj_names:

        name = str(name)
        lst = [dat.tfile[name] for dat in ucnlist if name in dat.tfile.keys()]
        df = pd.concat(lst, axis='index')

        first_attrs = lst[0].attrs # first item for its metadata

        # merge tree
        if first_attrs['type'] is ttree:
            ucnmerged.tfile[name] = df.copy()
            ucnmerged.tfile[name].attrs['type'] = ttree

        # merge th1
        elif first_attrs['type'] is th1:

            # error pre-treatment
            df[first_attrs['ylabel'] + " error"] **= 2

            # sum
            values = df.groupby(first_attrs['xlabel']).sum()

            # errors post treatment
            values[first_attrs['ylabel'] + " error"] **= 0.5

            # reset index
            values.reset_index(inplace=True)

            # set
            ucnmerged.tfile[name] = values.copy()

            # common attrs
            ucnmerged.tfile[name].attrs['type'] = th1
            for key in ('name', 'title', 'xlabel', 'ylabel', 'base_class', 'nbins'):
                ucnmerged.tfile[name].attrs[key] = first_attrs[key]

            # summed attrs
            ucnmerged.tfile[name].attrs['sum'] = sum([i.attrs['sum'] for i in lst])
            ucnmerged.tfile[name].attrs['entries'] = int(sum([i.attrs['entries'] for i in lst]))

        # merge th2
        elif first_attrs['type'] is th2:

            # error pre-treatment
            df[first_attrs['zlabel'] + " error"] **= 2

            # sum
            values = df.groupby([first_attrs['xlabel'], first_attrs['ylabel']]).sum()

            # errors post treatment
            values[first_attrs['zlabel'] + " error"] **= 0.5

            # set
            ucnmerged.tfile[name] = values.copy()

            # common attrs
            ucnmerged.tfile[name].attrs['type'] = th2
            for key in ('name', 'title', 'xlabel', 'ylabel', 'zlabel',
                        'base_class', 'nbinsx', 'nbinsy'):
                ucnmerged.tfile[name].attrs[key] = first_attrs[key]

            # summed attrs
            ucnmerged.tfile[name].attrs['sum'] = sum([i.attrs['sum'] for i in lst])
            ucnmerged.tfile[name].attrs['entries'] = int(sum([i.attrs['entries'] for i in lst]))

    # back to original objects
    ucnmerged.tfile.from_dataframe()

    return ucnmerged