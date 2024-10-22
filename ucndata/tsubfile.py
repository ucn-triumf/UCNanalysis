# Fetch a subrange of tfile object
# Derek Fujimoto
# Oct 2024

from .exceptions import *
import numpy as np
import pandas as pd
from rootloader import tfile

class tsubfile(tfile):
    """Wrapper for tfile which restricts access to values only within given times

    Args:
        tfileobj (tfile): object to wrap
        start (int): starting epoch time
        stop (int): stopping epoch time
    """

    def __init__(self, tfileobj, start, stop):

        for key, value in tfileobj.items():
            self[key] = value

        self._start = start
        self._stop = stop

    def __getitem__(self, key):

        # get the data
        val = super().__getitem__(key)

        # convert to dataframe
        is_dataframe = type(val) is pd.DataFrame
        if not is_dataframe:
            try:
                val = val.to_dataframe()
            except AttributeError:
                return val

        # get sub range
        try:
            index_name = val.index.name.lower()
        except AttributeError:
            pass
        else:
            if 'time' in index_name:
                start = np.min(val.index[val.index >= self._start]) # not sure why needed
                stop = np.max(val.index[val.index <= self._stop]) # it should work without it
                val = val.loc[start:stop]

        # convert back
        if not is_dataframe:
            return val.attrs['type'](val)
        else:
            return val

    def __getattr__(self, name):

        if name in self.keys():
            return self[name]
        else:
            return self.__getattribute__(name)

