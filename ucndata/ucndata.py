# Open and analyze UCN data
# Derek Fujimoto
# June 2024

from rootloader import tfile, attrdict, ttree, tdirectory
from .exceptions import MissingDataException
import ROOT
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class ucndata(object):
    """UCN run data. Cleans data and performs analysis

    Args:
        cycle (int|None): indicates cycle number, none if full file
        cycle_start (float): epoch time cycle start time (only if single cycle)
        cycle_stop (float) : epoch time cycle stop time (only if single cycle)
        filename (str): path to file to open
        header_only (bool): if true, read only the header

    Attributes:
        tfile (tfile): stores tfile raw readback
    """

    def __init__(self, filename, header_only=False):

        if filename is None:
            return

        # read
        if header_only:
            fid = ROOT.TFile(filename, 'READ')
            head = ttree(fid.Get('header'))
            fid.Close()
            self._head = head # needed to keep value in memory
        else:
            self.tfile = tfile(filename, empty_ok=False, quiet=True)
            head = self.tfile['header']

        # fix header values in tfile
        for key, value in self.tfile.header.items():
            self.tfile.header[key] = str(value[0])

        # reformat header and move to top level
        for k, val in head.items():
            setattr(self, k.replace(' ', '_').lower(), val)
        self.run_number = int(self.run_number)

        # reformat tfile branch names to remove spaces
        for key in tuple(self.tfile.keys()):
            if ' ' in key:
                self.tfile[key.replace(' ', '_')] = self.tfile[key]
                del self.tfile[key]

        # stop
        if header_only: return

        # set detector default
        self._names = attrdict()
        self.set_li6()

        # cycle number
        self.cycle = None

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
            cyc_str = '' if self.cycle is None else f' (cycle {self.cycle})'
            s = f'run {self.run_number}{cyc_str}:\n'
            for key in zip(*klist):
                s += '  '
                s += ''.join(['{0: <{1}}'.format(k, maxsize) for k in key])
                s += '\n'
            return s
        else:
            return self.__class__.__name__ + "()"

    def check_data(self):
        """Run some checks to determine if the data is ok.

        Checks:
            Do the following trees exist and have entries?
                BeamlineEpics
                UCN2Epics
                SequencerTree
                LNDDetectorTree
            Are there nonzero counts in UCNHits?
        """

        # check data trees
        keys = tuple(self.tfile.keys())
        for tree in ('BeamlineEpics', 'SequencerTree', 'LNDDetectorTree'):#, 'UCN2Epics'):

            # does tree exist?
            if tree not in keys:
                raise MissingDataException(f'Missing ttree "{tree}"')

            # does tree have entries?
            if self.tfile[tree].entries == 0:
                raise MissingDataException(f'Zero entries found in "{tree}" ttree')

        # check for nonzero counts
        if not self.tfile[self._names.hits].tIsUCN.any():
            raise MissingDataException(f'No UCN hits in "{self._names.hits}" ttree')

    def copy(self):
        """Return a copy of this objet"""
        copy = ucndata(None)

        for key, value in self.__dict__.items():
            if hasattr(value, 'copy'):
                setattr(copy, key, value.copy())
            else:
                setattr(copy, key, value)
        return copy

    def get_cycle(self, cycle):
        """Return a copy of this object, but trees are trimmed to only one cycle.

        Note that this process converts all objects to dataframes

        Args:
            cycle (int): cycle number

        Returns:
            ucndata: a copy of this object but with data from only one cycle.
        """

        # make copy
        copy = self.copy()

        # get cycles to keep
        cycles = copy.get_cycles_times()
        start = int(cycles.loc[cycle, 'start'])
        stop = int(cycles.loc[cycle, 'stop'])

        # trim the trees
        for key, value in copy.tfile.items():
            if key == 'header': continue

            # trim ttree
            if type(value) is ttree:
                value = value.to_dataframe()
                if value.index.name != '':
                    idx = (value.index < stop) & (value.index > start)
                    value = value.loc[idx]
                copy.tfile[key] = ttree(value)

            # trim dataframe
            elif type(value) is pd.DataFrame:
                if value.index.name != '':
                    idx = (value.index < stop) & (value.index > start)
                    copy.tfile[key] = value.loc[idx].copy()

        copy.cycle = cycle
        copy.cycle_start = start
        copy.cycle_stop = stop
        return copy

    def get_cycles_times(self):
        """Get start and end times of each cycle from the sequencer

        Returns:
            pd.DataFrame: with columns "start", "stop", and "duration (s)". Values are in epoch time. Indexed by cycle id
        """

        # get data
        df = self.tfile.SequencerTree
        if type(df) == ttree:
            df = df.to_dataframe()
        df = df.diff()

        # get start and end times
        times = {'start': df.index[df.inCycle == 1],
                 'stop': df.index[df.inCycle == -1]}

        # convert to dataframe
        times = pd.DataFrame(times)
        times['duration (s)'] = times.stop - times.start
        times.index.name = 'cycle'

        return times

    def get_hits_histogram(self, bin_ms=100):
        """Get histogram of UCNHits ttree times

        Args:
            bin_ms (int): histogram bin size in milliseconds

        Returns:
            tuple: (bin_centers, histogram counts)
        """

        # index_col = 'tUnixTimePrecise'

        # get cycle start and stop times from ucn counts histogram

        # get data
        df = self.tfile[self._names.hits].to_dataframe()
        index_col = df.index.name
        df.reset_index(inplace=True)

        # purge bad timestamps
        df = df.loc[df[index_col] > 15e8]

        # combine timestamps which are identical
        df = df.groupby(index_col).sum()

        # get timesteps for which there is an ucn
        times = df.index[df.tIsUCN.values.astype(bool)].values
        times = np.sort(times)

        # get histogram bin edges
        bins = np.arange(times.min(), times.max()+bin_ms/1000, bin_ms/1000)
        bins -= bin_ms/1000/2

        # histogram
        hist, bins = np.histogram(times, bins=bins)
        bin_centers = (bins[1:] + bins[:-1])/2

        return (bin_centers, hist)

    def set_he3(self):
        """Set name patterns to match He3 detector"""
        self._names['hits'] =        'UCNHits_He3'
        self._names['charge'] =      'He3_Charge'
        self._names['rate'] =        'He3_Rate'
        self._names['transitions'] = 'RunTransitions_He3'
        self._names['hitsseq'] =     'hitsinsequence_he3'
        self._names['hitsseqcumul'] ='hitsinsequencecumul_he3'

    def set_li6(self):
        """Set name patterns to match li6 detector"""
        self._names['hits'] =        'UCNHits_Li-6'
        self._names['charge'] =      'Li6_Charge'
        self._names['rate'] =        'Li6_Rate'
        self._names['transitions'] = 'RunTransitions_Li-6'
        self._names['hitsseq'] =     'hitsinsequence_li6'
        self._names['hitsseqcumul'] ='hitsinsequencecumul_li6'

    def to_dataframe(self):
        """Convert self.tfile contents to pd.DataFrame"""

        for key, value in self.tfile.items():

            # non-header items
            if key not in ('header',):

                # tdirectory is converted inplace
                if type(value) is tdirectory:
                    value.to_dataframe()

                # everything else produces a copy
                else:
                    self.tfile[key] = value.to_dataframe()

            # convert header
            else:
                self.tfile[key] = pd.Series(value)