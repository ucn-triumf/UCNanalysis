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
from .applylist import applylist
from .ucnbase import ucnbase
from .ucnperiod import ucnperiod
from .tsubfile import tsubfile

import numpy as np
import os


class ucncycle(ucnbase):
    """View for the data from a single UCN cycle

    Notes:
        Any changes to the data frame affects all periods (for the time steps
        contained in that period) and the containing run

    Args:
        urun (ucnrun): object to pull cycle from
        cycle (int): cycle number
    """

    def __init__(self, urun, cycle):

        # get cycles to keep
        cycles = urun.cycle_param.cycle_times
        start = int(cycles.loc[cycle, 'start'])
        stop = int(cycles.loc[cycle, 'stop'])
        supercycle = int(cycles.loc[cycle, 'supercycle'])

        # copy data
        for key, value in urun.__dict__.items():
            if key == 'tfile':
                setattr(self, key, tsubfile(value, start, stop))
            elif hasattr(value, 'copy'):
                setattr(self, key, value.copy())
            else:
                setattr(self, key, value)

        # trim cycle parameters
        self.cycle_param.period_durations_s = self.cycle_param.period_durations_s[cycle]
        self.cycle_param.period_end_times = self.cycle_param.period_end_times[cycle]

        self.cycle = cycle
        self.supercycle = supercycle
        self.cycle_start = start
        self.cycle_stop = stop

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

    def __getitem__(self, key):
        # get cycle or period based on slicing indexes

        # get a single key
        if isinstance(key, (np.integer, int)):
            if key > self.cycle_param.nperiods:
                raise IndexError(f'Run {self.run_number}, cycle {self.cycle}: Index larger than number of periods ({self.cycle_param.nperiods})')

            return self.get_period(key)

        # slice on cycles
        if isinstance(key, slice):
            period = self.get_period()[:self.cycle_param.nperiods]
            return applylist(period[key])

        raise IndexError('Cycles indexable only as a 1-dimensional object')

    def check_data(self, period_production=None, period_count=None, period_background=None,
                   raise_error=False, quiet=False):
        """Run some checks to determine if the data is ok.

        Args:
            period_production (int): index of period where the beam should be stable. Enables checks of beam stability
            period_count (int): index of period where we count ucn. Enables checks of data quantity
            period_background (int): index of period where we do not count ucn. Enables checks of background
            raise_error (bool): if true, raise an error if check fails, else return false
            quiet (bool): if true don't print or raise exception

        Returns:
            bool: true if check passes, else false.

        Checks:
            Do the following trees exist and have entries?
                BeamlineEpics
                UCN2Epics
                SequencerTree
                LNDDetectorTree
            Are there nonzero counts in UCNHits?
        """
        # setup error message
        msg = f'Run {self.run_number}, cycle {self.cycle}:'

        # setup raise or warn
        if quiet:
            def warn(error, message):
                return False
        elif raise_error:
            def warn(error, message):
                raise error(message)
        else:
            def warn(error, message):
                print(message)
                return False

        ## overall data quality checks ---------------------------------------

        # beam data exists
        if len(self.beam_current_uA) == 0:
            return warn(BeamError, f'{msg} No beam data saved')

        # total duration
        if self.cycle_stop - self.cycle_start <= 0:
            return warn(DataError, f'{msg} Cycle duration nonsensical: {self.cycle_stop - self.cycle_start} s')

        # skip cycle because cycle sequence is longer than irradiation interval
        # TODO

        # valve states
        if not self.cycle_param.valve_states.any().any():
            return warn(ValveError, f'{msg} No valves operated')

        # has counts
        if not any([self.tfile[self.DET_NAMES[det]['hits']].tIsUCN.sum() > 1 for det in self.DET_NAMES.keys()]):
            return warn(DataError, f'{msg} No counts detected')

        ## production period checks ------------------------------------------
        if period_production is not None:
            period = self.get_period(period_production)
            beam_current = period.beam_current_uA

            # beam data exists during production
            if len(beam_current) == 0:
                return warn(BeamError, f'{msg} No beam data during production period')

            # beam dropped too low
            if beam_current.min() < self.DATA_CHECK_THRESH['beam_min_current']:
                return warn(BeamError, f'{msg} Beam current dropped to {beam_current.min()} uA')

            # beam current unstable
            if beam_current.std() > self.DATA_CHECK_THRESH['beam_max_current_std']:
                return warn(BeamError, f'{msg} Beam current fluctuated by {beam_current.std()} uA')

        ## background period checks ------------------------------------------
        if period_background is not None:

            period = self.get_period(period_background)

            for det in self.DET_NAMES.keys():
                counts, _ = period.get_counts(det)

                # background count rate too high
                rate = counts / period.cycle_param.period_durations_s
                if rate / self.DET_BKGD[det] > self.DATA_CHECK_THRESH['max_bkgd_count_rate']:
                    raise warn(DataError, f'{msg} Background count rate in {det} detector is {rate / self.DET_BKGD[det]} times larger than expected ({self.DET_BKGD[det]} counts/s)')

                # background counts missing
                if counts == 0:
                    raise warn(DataError, f'{msg} No detected background counts in {det} detector')

        ## count period checks -----------------------------------------------
        if period_count is not None:

            period = self.get_period(period_count)
            for det in self.DET_NAMES.keys():
                counts = period.get_counts(det)[0]
                if counts > self.DATA_CHECK_THRESH['min_total_counts']:
                    raise warn(DataError, f'{msg} Too few counts in {det} detector during counting period ({counts} counts)')


        return True

    # TODO: FIX THIS
    def get_counts(self, detector, period=None, bkgd=None, norm=None):
        """Get counts for each period
        Args:
            detector (str): one of the keys to self.DET_NAMES
            period (None|int):  if None get for entire cycle
                                elif < 0 get for each period
                                elif >=0 get for that period
            bkgd (tuple|None): if not None subtract this as the background (value, error)
                               bkgd.shape = (2, nperiods) if period < 0 else (2, 1)
            norm (tuple|None): if not None normalize to this value (value, error)
                                norm.shape = (2, nperiods) if period < 0 else (2, 1)

        Returns:
            np.ndarray: number of hits for each period and error
        """

        # check input
        nperiods = self.cycle_param.nperiods
        if period is not None and period > nperiods:
            raise RuntimeError(f"Run {self.run_number}, cycle {self.cycle}: Period index must be less than {self.cycle_param.nperiods}")

        # get ucn hits
        hit_tree = self.get_hits(detector)

        if period is None:
            counts = len(hit_tree.index)
        else:

            # make histogram of counts
            edges = np.concatenate(([self.cycle_start], self.cycle_param.period_end_times))
            counts, _ = np.histogram(hit_tree.index, bins=edges)

            # trim to nperiods
            if period < 0:
                counts = counts[:nperiods]
                edges = edges[:nperiods+1]

            # select single period
            else:
                counts = counts[period]
                edges = edges[period:period+2]

        # error assumed poissonian
        dcounts = np.sqrt(counts)

        # subtract background
        if bkgd is not None:
            if isinstance(counts (int, np.int64)):  zero = 0
            else:                                   zero = np.zeros(len(counts))
            counts = np.max(counts-bkgd[0], zero)
            dcounts = (dcounts**2 + bkgd[1]**2)**0.5

        # normalize
        if norm is not None:
            dcounts = counts*((dcounts/counts)**2 + (norm[1]/norm[0])**2)**0.5
            counts /= norm[0]

        return (counts, dcounts)

    def get_period(self, period=None):
        """Return a copy of this object, but trees are trimmed to only one period.

        Notes:
            This process converts all objects to dataframes
            Must be called for a single cycle only

        Args:
            period (int): period number, if None, get all periods
            cycle (int|None) if cycle not specified then specify a cycle

        Returns:
            run:
                if period > 0: a copy of this object but with data from only one period.
                if period < 0 | None: a list of copies of this object for all periods for a single cycle
        """

        # get all periods
        if period is None or period < 0:
            nperiods = self.cycle_param.nperiods
            return applylist(map(self.get_period, range(nperiods)))
        else:
            return ucnperiod(self, period)

    #TODO: FIX THIS
    def get_rate(self, detector, bkgd=True, norm=False):
        """Get count rate for each period
        Args:
            detector (str): one of the keys to self.DET_NAMES
            bkgd (tuple|None): if not None subtract this as the background (value, error)
            norm (tuple|None): if not None normalize to this value (value, error)


        Returns:
            np.ndarray: count rate each period and error
        """
        rate = [p.get_rate(detector, bkgd, norm) for p in self.periods()]
        return np.array(rate)
