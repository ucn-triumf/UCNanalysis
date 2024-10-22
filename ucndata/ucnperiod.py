# Open and analyze UCN data for a one cycle
# Derek Fujimoto
# Oct 2024

"""
    TODO List, things which haven't been ported from WS code

    * get temperature
    * get vapour pressure
    * data checks for periods
    * check that period durations match between detector frontends
"""

from .exceptions import *
from . import settings
from .ucnbase import ucnbase
from .tsubfile import tsubfile

import numpy as np
import os

class ucnperiod(ucnbase):
    """Stores the data from a single UCN period from a single cycle

    Args:
        ucycle (ucncycle): object to pull period from
        period (int): period number
    """

    def __init__(self, ucycle, period):

        # get start and stop time
        if period > 0:      start = int(ucycle.cycle_param.period_end_times[period-1])
        else:               start = int(ucycle.cycle_start)

        stop = int(ucycle.cycle_param.period_end_times[period])

        # copy data
        for key, value in ucycle.__dict__.items():
            if key == 'tfile':
                setattr(self, key, tsubfile(value, start, stop))
            elif hasattr(value, 'copy'):
                setattr(self, key, value.copy())
            else:
                setattr(self, key, value)

        # trim cycle parameters
        self.cycle_param.period_durations_s = self.cycle_param.period_durations_s[period]
        self.cycle_param.period_end_times = self.cycle_param.period_end_times[period]

        self.period = period
        self.period_start = start
        self.period_stop = stop

    def __repr__(self):
        klist = [d for d in self.__dict__.keys() if d[0] != '_']
        if klist:

            # sort without caps
            klist.sort(key=lambda x: x.lower())

            # get number of columns based on terminal size
            maxsize = max((len(k) for k in klist)) + 2
            terminal_width = os.get_terminal_size().columns
            ncolumns = int(np.floor(terminal_width / maxsize))
            ncolumns = min(ncolumns, len(klist))

            # split into chunks
            needed_len = int(np.ceil(len(klist) / ncolumns)*ncolumns) - len(klist)
            klist = np.concatenate((klist, np.full(needed_len, '')))
            klist = np.array_split(klist, ncolumns)

            # print
            cyc_str = f' (cycle {self.cycle}, period {self.period})'
            s = f'run {self.run_number}{cyc_str}:\n'
            for key in zip(*klist):
                s += '  '
                s += ''.join(['{0: <{1}}'.format(k, maxsize) for k in key])
                s += '\n'
            return s
        else:
            return self.__class__.__name__ + "()"

    def get_counts(self, detector, bkgd=None, dbkgd=None, norm=None, dnorm=None):
        """Get sum of ucn hits

        Args:
            detector (str): one of the keys to settings.DET_NAMES
            bkgd (float|None): background counts
            dbkgd(float|None): error in background counts
            norm (float|None): normalize to this value
            dnorm (float|None): error in normalization

        Returns:
            tuple: (count, error) number of hits
        """
        hit_tree = self.get_hits(detector)
        counts = len(hit_tree.index)
        dcounts = np.sqrt(counts) # error assumed poissonian

        # subtract background, but no less than 0 counts
        if bkgd is not None:

            # check if iterable, else fetch only for this period
            try:                iter(bkgd)
            except TypeError:   b = bkgd
            else:               b = bkgd[self.period]
            counts = max(counts-b, 0)

            # error correction
            if dbkgd is not None:

                # check if iterable, else fetch only for this period
                try:                iter(dbkgd)
                except TypeError:   db = dbkgd
                else:               db = dbkgd[self.period]
                dcounts = (dcounts**2 + db**2)**0.5

        # normalize with error corretion
        if dnorm is not None:

            # check if iterable, else fetch only for this period
            try:
                iter(dnorm)
            except TypeError:
                dn = dnorm
                n  = norm
            else:
                dn = dnorm[self.period]
                n  = norm[self.period]

            # normalize
            dcounts = counts*((dcounts/counts)**2 + (dn/n)**2)**0.5
            counts /= n

        # normalize without error correction
        elif norm is not None:

            try:                iter(norm)
            except TypeError:   n = norm
            else:               n = norm[self.period]
            counts /= n

        return (counts, dcounts)

    def get_hits(self, detector):
        """Get times of ucn hits

        Args:
            detector (str): one of the keys to settings.DET_NAMES

        Returns:
            pd.DataFrame: hits tree as a dataframe, only the values when a hit is registered
        """

        # get the tree
        hit_tree = super().get_hits(detector)

        ## filter pileup for period data

        # get thresholds
        dt = settings.DATA_CHECK_THRESH['pileup_within_first_s']
        count_thresh = settings.DATA_CHECK_THRESH['pileup_cnt_per_ms']

        # make histogram
        t = hit_tree.index.values
        counts, edges = np.histogram(t, bins=int(1/0.001),
                                        range=(min(t), min(t)+dt))

        # delete bad count ranges
        ncounts_total = len(t)
        for i, count in enumerate(counts):
            if count > count_thresh:
                hit_tree = hit_tree.loc[(hit_tree.index < edges[i]) | (hit_tree.index > edges[i+1])]
        ncounts_removed = ncounts_total - len(hit_tree.index)

        if ncounts_removed > 0:
            print(f'Removed {ncounts_removed} pileup counts ({int(ncounts_removed/ncounts_total*100):d}%) from run{self.run_number} (cycle {self.cycle}, period {self.period})')

        return hit_tree

    def get_rate(self, detector, bkgd=None, dbkgd=None, norm=None, dnorm=None):
        """Get sum of ucn hits per unit time of period

        Args:
            detector (str): one of the keys to settings.DET_NAMES
            bkgd (float|None): background counts
            dbkgd(float|None): error in background counts
            norm (float|None): normalize to this value
            dnorm (float|None): error in normalization

        Returns:
            tuple: (count rate, error)
        """
        # get counts
        counts, dcounts = self.get_counts(detector=detector,
                                          bkgd=bkgd,
                                          dbkgd=dbkgd,
                                          norm=norm,
                                          dnorm=dnorm,
                                          )

        # get rate
        duration = self.cycle_param.period_durations_s

        counts /= duration
        dcounts /= duration

        return (counts, dcounts)

