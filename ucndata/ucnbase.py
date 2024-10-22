# Base class for ucnrun, ucncycle, and ucnperiod
# Derek Fujimoto
# Oct 2024

from rootloader import ttree
from .exceptions import *
import ucndata.constants as const
import numpy as np
import pandas as pd

class ucnbase(object):
    """UCN run data. Cleans data and performs analysis

    Args:
        run (int|str): if int, generate filename with settings.datadir
            elif str then run is the path to the file
        header_only (bool): if true, read only the header

    Attributes:
        comment (str): comment input by users
        cycle (int|none): cycle number, none if no cycle selected
        cycle_param (attrdict): cycle parameters from sequencer settings
        experiment_number (str): experiment number input by users
        month (int): month of run start
        run_number (int): run number
        run_title (str): run title input by users
        shifter (str): experimenters on shift at time of run
        start_time (str): start time of the run
        stop_time (str): stop time of the run
        supercycle (int|none): supercycle number, none if no cycle selected
        tfile (tfile): stores tfile raw readback
        year (int): year of run start

    Notes:
        Can access attributes of tfile directly from top-level object
        Need to define the values in ucndata.settings if you want non-default
        behaviour
        Object is indexed as [cycle, period] for easy access to sub time frames
    """

    # detector names
    DET_NAMES = {'He3':{'hits':         'UCNHits_He3',
                        'charge':       'He3_Charge',
                        'rate':         'He3_Rate',
                        'transitions':  'RunTransitions_He3',
                        'hitsseq':      'hitsinsequence_he3',
                        'hitsseqcumul': 'hitsinsequencecumul_he3',
                        },
                 'Li6':{'hits':         'UCNHits_Li-6',
                        'charge':       'Li6_Charge',
                        'rate':         'Li6_Rate',
                        'transitions':  'RunTransitions_Li-6',
                        'hitsseq':      'hitsinsequence_li6',
                        'hitsseqcumul': 'hitsinsequencecumul_li6',
                        },
                }

    # needed slow control trees
    SLOW_TREES = ('BeamlineEpics', 'SequencerTree', 'LNDDetectorTree')

    # data thresholds for checking data
    DATA_CHECK_THRESH = {'beam_min_current': 0.1, # uA
                         'beam_max_current_std': 0.02, # uA
                         'max_bkgd_count_rate': 4, # fractional increase over DET_BKGD values
                         'min_total_counts': 100, # number of counts total
                         'pileup_cnt_per_ms': 3, # if larger than this, then pileup and delete
                         'pileup_within_first_s': 1, # time frame for pileup in each period
                         }

    # default detector backgrounds - from 2019
    DET_BKGD = {'Li6':     1.578,
                'Li6_err': 0.009,
                'He3':     0.0349,
                'He3_err': 0.0023}

    def __iter__(self):
        # setup iteration
        self._iter_current = 0
        return self

    def __next__(self):
        # permit iteration over object like it was a list

        if self._iter_current < self.cycle_param.ncycles:

            # skip elements that should be filtered
            if self.cycle_param.filter is not None:

                # skip
                while not self.cycle_param.filter[self._iter_current]:
                    self._iter_current += 1

                    # end condition
                    if self._iter_current >= self.cycle_param.ncycles:
                        raise StopIteration()

            # iterate
            cyc = self[self._iter_current]
            self._iter_current += 1
            return cyc

        # end of iteration
        else:
            raise StopIteration()

    def _get_beam_duration(self, on=True):
        # Get beam on/off durations

        # get needed info
        cycle_times = self.cycle_param.cycle_times

        try:
            beam = self.tfile.BeamlineEpics
        except AttributeError:
            raise MissingDataError("No saved ttree named BeamlineEpics")

        # setup storage
        beam_dur = []
        epics_val = 'B1V_KSM_RDBEAMON_VAL1' if on else 'B1V_KSM_RDBEAMOFF_VAL1'

        # get as dataframe
        if isinstance(beam, ttree):
            beam = beam.to_dataframe()

        # get durations closest to cycle start time
        for start in cycle_times.start:

            start_times = abs(beam.index.values - start)
            durations = beam[epics_val].values*const.beam_bucket_duration_s
            idx = np.argmin(start_times)
            beam_dur.append(durations[idx])

        out = pd.Series(beam_dur, index=cycle_times.index)
        out.index.name = cycle_times.index.name

        if len(out) == 1:
            return float(out.values[0])
        return out

    def apply(self, fn_handle):
        """Apply function to each cycle

        Args:
            fn_handle (function handle): function to be applied to each cycle

        Returns:
            np.ndarray: output of the function
        """
        return [fn_handle(c) for c in self.cycles()]

    def copy(self):
        """Return a copy of this objet"""
        copy = ucnrun(None)

        for key, value in self.__dict__.items():
            if hasattr(value, 'copy'):
                setattr(copy, key, value.copy())
            else:
                setattr(copy, key, value)
        return copy

    def get_hits(self, detector):
        """Get times of ucn hits

        Args:
            detector (str): one of the keys to self.DET_NAMES

        Returns:
            pd.DataFrame: hits tree as a dataframe, only the values when a hit is registered
        """

        # get the tree
        hit_tree = self.tfile[self.DET_NAMES[detector]['hits']] # maybe should be a copy?
        if type(hit_tree) is not pd.DataFrame:
            hit_tree = hit_tree.to_dataframe()

        # get times only when a hit is registered
        hit_tree = hit_tree.loc[hit_tree.tIsUCN.astype(bool)]

        return hit_tree

    def get_hits_histogram(self, detector, bin_ms=100):
        """Get histogram of UCNHits ttree times

        Args:
            detector (str): Li6|He3
            bin_ms (int): histogram bin size in milliseconds

        Returns:
            tuple: (bin_centers, histogram counts)
        """

        # get data
        df = self.tfile[self.DET_NAMES[detector]['hits']]

        if not isinstance(df, pd.DataFrame):
            df = df.to_dataframe()

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

    def from_dataframe(self):
        """Convert self.tfile contents to rootfile struture types"""
        self.tfile.from_dataframe()

    def to_dataframe(self):
        """Convert self.tfile contents to pd.DataFrame"""
        self.tfile.to_dataframe()

    # quick access properties
    @property
    def beam_current_uA(self):

        if type(self.tfile.BeamlineEpics) is pd.DataFrame:
            df = self.tfile.BeamlineEpics
        else:
            df = self.tfile.BeamlineEpics.to_dataframe()

        # PREDCUR is the predicted current in beamline 1U.
        # PREDCUR is calculated by using the beamline 1V extraction foil current
        # (the current as it leaves the cyclotron) and multiplid by the fraction
        # of beam that is going to the 1U beamline (as opposed to 1A beamline).
        # So if the extraction foil current is 100uA and we are kicking 1 bucket
        # out of 10 buckets to 1U, then PREDCUR will be 10uA
        predcur = df.B1V_KSM_PREDCUR

        # BONPRD is a bool, which indicates if there is beam down 1U
        bonprd = df.B1V_KSM_BONPRD

        # current in the 1U beamline
        return predcur*bonprd

    @property
    def beam_on_s(self): return self._get_beam_duration(on=True)

    @property
    def beam_off_s(self): return self._get_beam_duration(on=False)
