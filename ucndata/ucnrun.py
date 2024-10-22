# Open and analyze UCN data for a whole run
# Derek Fujimoto
# June 2024

"""
    TODO List, things which haven't been ported from WS code

    * get temperature
    * get vapour pressure
    * data checks for periods
    * check that period durations match between detector frontends
"""

from rootloader import tfile, ttree, attrdict
from .exceptions import *
from .applylist import applylist
from .ucnbase import ucnbase
from .ucncycle import ucncycle
from . import settings
import ROOT
import numpy as np
import pandas as pd
import itertools, warnings, os

ROOT.gROOT.SetBatch(1)

class ucnrun(ucnbase):
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

    def __init__(self, run, header_only=False):

        # check if copying
        if run is None:
            return

        # make filename from defaults
        elif type(run) is int:
            filename = os.path.join(settings.datadir, f'ucn_run_{run:0>8d}.root')

        # fetch from specified path
        elif type(run) is str:
            filename = run

        # read
        if header_only:
            fid = ROOT.TFile(filename, 'READ')
            head = ttree(fid.Get('header'))
            fid.Close()
            head = {k:str(val[0]) for k, val in head.items()}
            self._head = head # needed to keep value in memory

        else:
            self.tfile = tfile(filename, empty_ok=False, quiet=True,
                               key_filter=settings.keyfilter)
            head = self.tfile['header']

            # fix header values in tfile
            for key, value in self.tfile.header.items():
                self.tfile.header[key] = str(value[0])

        # reformat header and move to top level
        for k, val in head.items():
            setattr(self, k.replace(' ', '_').lower(), val)

        if type(self.run_number) is pd.Series:
            self.run_number = int(self.run_number[0])
        else:
            self.run_number = int(self.run_number)

        # set cycle parameters
        self._get_cycle_param()

        # set other header items
        date = pd.to_datetime(self.start_time)
        self.year = date.year
        self.month = date.month

        # stop
        if header_only:
            return

        # reformat tfile branch names to remove spaces
        for key in tuple(self.tfile.keys()):
            if ' ' in key:
                self.tfile[key.replace(' ', '_')] = self.tfile[key]
                del self.tfile[key]

        # get cycle times
        self.set_cycle_times()

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
            s = f'run {self.run_number}:\n'
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
            if key > self.cycle_param.ncycles:
                raise IndexError(f'Run {self.run_number}: Index larger than number of cycles ({self.cycle_param.ncycles})')
            return self.get_cycle(key)

        # slice on cycles
        if isinstance(key, slice):
            cycles = self.get_cycle()[:self.cycle_param.ncycles]

            # no filter
            if self.cycle_param.filter is None or all(self.cycle_param.filter):
                cyc = cycles[key]

            # yes filter
            else:

                # fetch the filter and slice in the same way as the return value
                cfilter = self.cycle_param.filter[key]

                # fetch cycles and slice, then apply filter
                cyc = np.array(cycles[key])
                cyc = cyc[cfilter]

            return applylist(cyc)

        # slice on periods
        if isinstance(key, tuple):
            cycles = self[key[0]]
            if isinstance(cycles, (np.ndarray, applylist, list)):
                return applylist([c[key[1]] for c in cycles])
            else:
                return cycles[key[1]]

        raise IndexError(f'Run {self.run_number} given an unknown index type ({type(key)})')

    def _get_cycle_param(self):
        # set self.cycle_param dict

        # reformat cycle param tree
        paramtree = self.tfile.CycleParamTree
        self.cycle_param = {'nperiods': paramtree.nPeriods[0],
                            'nsupercyc': paramtree.nSuperCyc[0],
                            'enable': bool(paramtree.enable[0]),
                            'inf_cyc_enable': bool(paramtree.infCyclesEnable[0]),
                            }
        self.cycle_param = attrdict(self.cycle_param)

        # setup cycle paramtree array outputs from transition trees
        for detector in settings.DET_NAMES.values():
            if detector['transitions'] in self.tfile.keys():
                tree = self.tfile[detector['transitions']]
                break
        tree = tree.to_dataframe()

        # cycle and supercycle indices
        self.cycle_param['cycle'] = tree['cycleIndex'].astype(int)
        self.cycle_param['supercycle'] = tree['superCycleIndex'].astype(int)

        # convert the array in each cell into a dataframe
        def item_to_df(x):
            s = pd.DataFrame(np.array(x).copy(), index=np.arange(len(x))+1)
            return s

        # valve states -------------------------------------------------------
        df = tree[[col for col in tree.columns if 'valveState' in col]]
        col_map = {col:int(col.replace("valveStatePeriod", "")) for col in df.columns}
        df = df.rename(columns=col_map)
        df.columns.name = 'period'
        df.index.name = 'cycle_idx'

        # valve states should not change across cycles
        df2 = df.loc[0]
        df2 = pd.concat([item_to_df(df2[period]) for period in df2.index], axis='columns')

        # rename columns and index
        df2.columns = np.arange(len(df2.columns))
        df2.columns.name = 'period'
        df2.index.name = 'valve'
        self.cycle_param['valve_states'] = df2.transpose()

        # period end times ---------------------------------------------------
        df = tree[[col for col in tree.columns if 'cyclePeriod' in col]]
        col_map = {col:int(col.replace("cyclePeriod", "").replace("EndTime", "")) for col in df.columns}

        # rename columns and index
        df = df.rename(columns=col_map)
        df.columns.name = 'period'
        df.index.name = 'cycle'
        self.cycle_param['period_end_times'] = df.transpose()

        # period durations ---------------------------------------------------
        cycle_start = tree.cycleStartTime
        df_diff = df.diff(axis='columns')
        df_diff[0] = df[0] - cycle_start

        self.cycle_param['period_durations_s'] = df_diff.transpose()

        # number of cycles
        self.cycle_param['ncycles'] = len(df.index)

        # cycle filter
        self.cycle_param['filter'] = None

    def check_data(self, raise_error=False):
        """Run some checks to determine if the data is ok.

        Args:
            raise_error (bool): if true, raise an error if check fails, else return false

        Returns:
            bool: true if check passes, else false.

        Checks:
            Do the settings.SLOW_TREES exist and have entries?
            Are there nonzero counts in UCNHits?
        """

        # check some necessary data trees
        for tree in settings.SLOW_TREES:

            msg = None

            # does tree exist?
            if tree not in self.tfile.keys():
                msg = f'Missing ttree "{tree}" in run {self.run_number}'

            # does tree have entries?
            elif len(self.tfile[tree]) == 0:
                msg = f'Zero entries found in "{tree}" ttree in run {self.run_number}'

            # raise error or return
            if msg is not None:
                if raise_error:
                    raise MissingDataError(msg)
                else:
                    print(msg)
                    return False

        for name, det in settings.DET_NAMES.items():

            # check for nonzero counts
            if not self.tfile[det['hits']].tIsUCN.any():
                msg = f'No UCN hits in "{name}" ttree in run {self.run_number}'

            # check if sequencer was enabled but no run transitions
            elif any(self.tfile.SequencerTree.sequencerEnabled):

                if len(self.tfile[det['transitions']]) == 0:
                    msg = 'No cycles found in run {self.run_number}, although sequencer was active'

            # raise error or return
            if msg is not None:
                if raise_error:
                    raise MissingDataError(msg)
                else:
                    print(msg)
                    return False

        return True

    def gen_cycle_filter(self, period_production=None, period_count=None,
                         period_background=None, quiet=False):
        """Generate filter array for cycles. Use with self.set_cycle_filter to filter cycles.

        Args:
            period_production (int): index of period where the beam should be stable. Enables checks of beam stability
            period_count (int): index of period where we count ucn. Enables checks of data quantity
            period_background (int): index of period where we do not count ucn. Enables checks of background
            quiet (bool): if true don't print or raise exception

        Returns:
            np.array: of bool, true if keep cycle, false if discard

        Notes:
            calls ucncycle.check_data on each cycle
        """

        cycles = self.get_cycle()
        cfilter = [c.check_data(period_background=period_background,
                                period_count=period_count,
                                period_production=period_production,
                                quiet=quiet,
                                raise_error=False) for c in cycles]
        return np.array(cfilter)

    def get_cycle(self, cycle=None):
        """Return a copy of this object, but trees are trimmed to only one cycle.

        Note that this process converts all objects to dataframes

        Args:
            cycle (int): cycle number, if None, get all cycles

        Returns:
            ucncycle:
                if cycle > 0:  ucncycle object
                if cycle < 0 | None: a list ucncycle objects for all cycles
        """

        if cycle is None or cycle < 0:
            ncycles = len(self.cycle_param.cycle_times.index)
            return applylist(map(self.get_cycle, range(ncycles)))
        else:
            return ucncycle(self, cycle)

    def set_cycle_filter(self, cfilter=None):
        """Set filter for which cycles to fetch when slicing or iterating

        Notes:
            Filter is ONLY applied when fetching cycles as a slice or as an iterator. ucnrun.get_cycle() always returns unfiltered cycles.

        Examples where the filter is applied:
            * run[:]
            * run[3:10]
            * run[:3]
            * for c in run: print(c)

        Examples where the filter is not applied:
            * run[2]
            * run.get_cycle()
            * run.get_cycle(2)

        Args:
            cfilter (None|iterable): list of bool, True if keep cycle, False if reject.
                if None then same as if all True

        Returns:
            None: sets self.cycle_param.filter
        """

        # check input
        cfilter = np.array(cfilter).astype(bool)

        if len(cfilter) != self.cycle_param.ncycles:
            raise RuntimeError(f'Run {self.run_number}: Length of cycle filter ({len(cfilter)}) does not match expected number of cycles ({self.cycle_param.ncycles})')

        # set
        self.cycle_param.filter = cfilter

    def set_cycle_times(self, mode='matched'):
        """Get start and end times of each cycle from the sequencer and save
        into self.cycle_param.cycle_times

        Run this if you want to change how cycle start times are calculated

        Args:
            mode (str): matched|sequencer|he3|li6
                if matched: look for identical timestamps in RunTransitions from detectors
                if sequencer: look for inCycle timestamps in SequencerTree
                if he3: use He3 detector cycle start times
                if li6: use Li6 detector cycle start times

        Notes:
            - If run ends before sequencer stop is called, a stop is set to final timestamp.
            - If the sequencer is disabled mid-run, a stop is set when disable ocurrs.
            - If sequencer is not enabled, then make the entire run one cycle
            - For matched mode,
                - set run stops as start of next transition
                - set offset as start_He3 - start_Li6
                - set start/stop/duration based on start_He3
            - If the object reflects a single cycle, return from cycle_start, cycle_stop

        Returns:
            pd.DataFrame: with columns "start", "stop", "offset" and "duration (s)". Values are in epoch time. Indexed by cycle id. Offset is the difference in detector start times: he3_start-li6_start
        """

        # check if single cycle
        if hasattr(self, 'cycle'):
            return pd.DataFrame({'start':[self.cycle_start],
                                 'stop':[self.cycle_stop],
                                 'duration (s)': [self.cycle_stop-self.cycle_start],
                                 'offset (s)': [0.0],
                                 'supercycle': [self.supercycle]},
                                 index=[self.cycle])

        # get data
        df = self.tfile.SequencerTree
        if type(df) == ttree:
            df = df.to_dataframe()

        # get run end time from control trees - used in matched and detector cycles times
        run_stop = -np.inf
        for treename in settings.SLOW_TREES:
            try:
                idx = self.tfile[treename].to_dataframe().index
            except AttributeError as err:
                if type(self.tfile[treename]) is pd.DataFrame:
                    idx = self.tfile[treename].index
                else:
                    raise err from None
            run_stop = max((idx.max(), run_stop))

        ## if no sequencer, make the whole run a single cycle
        if not any(df.sequencerEnabled):

            times = {'start': np.inf,
                     'stop': -np.inf,
                     'supercycle': 0}

            # use timestamps from slow control trees to determine timestamps
            for treename in settings.SLOW_TREES:
                idx = self.tfile[treename].to_dataframe().index
                times['start'] = min((idx.min(), times['start']))
                times['stop']  = max((idx.max(), times['stop']))

        ## get matched timesteps from He3 and Li6 RunTransitions
        elif mode in 'matched':
            hestart = self.tfile[settings.DET_NAMES['He3']['transitions']].cycleStartTime
            listart = self.tfile[settings.DET_NAMES['Li6']['transitions']].cycleStartTime
            scycle = self.tfile[settings.DET_NAMES['He3']['transitions']].superCycleIndex

            # drop duplicate timestamps
            hestart = hestart.drop_duplicates()
            listart = listart.drop_duplicates()

            # get all possible pairs and sort by time difference
            pairs = sorted(itertools.product(hestart, listart),
                            key = lambda t: abs(t[0] - t[1]))

            # save output
            matchedhe3 = []
            matchedli6 = []

            # go through all possible pairs
            for pair in pairs:

                # if none of the two start times are already in the matched list, add the pair to the matched list
                if pair[0] not in matchedhe3 and pair[1] not in matchedli6:
                    offset = pair[0] - pair[1]

                    # discard if time difference is too large
                    if abs(offset) > 20:
                        warnings.warn(f'He3 cycle start time ({pair[0]}) too distant from Li6 start ({pair[1]}) in run {self.run_number}', CycleWarning)
                    else:
                        matchedhe3.append(pair[0])
                        matchedli6.append(pair[1])

            matchedhe3 = np.sort(matchedhe3)
            matchedli6 = np.sort(matchedli6)

            # warnings for unmatched cycles
            unmatched = [t not in matchedhe3 for t in hestart]
            if any(unmatched):
                warnings.warn(f'Found no match to He3 cycles at {unmatched} in run {self.run_number}', CycleWarning)

            unmatched = [t not in matchedhe3 for t in listart]
            if any(unmatched):
                warnings.warn(f'Found no match to Li6 cycles at {unmatched} in run {self.run_number}', CycleWarning)

            # setup output
            times = {'start': matchedhe3,
                     'duration (s)': np.concatenate((np.diff(matchedhe3), [run_stop])),
                     'offset (s)': matchedhe3-matchedli6}
            times['stop'] = times['start'] + times['duration (s)']
            times['supercycle'] = scycle

        ## get timestamps from sequencer
        elif mode in 'sequencer':

            # if sequencer is not enabled cause a stop transition
            df.inCycle *= df.sequencerEnabled

            # start counting only after first start flag
            df = df.loc[df.loc[df.cycleStarted > 0].index[0]:]

            # get start and end times
            df = df.diff()
            times = {'start': df.index[df.inCycle == 1],
                    'stop': df.index[df.inCycle == -1],
                    'supercycle': 0}

            # check lengths
            if len(times['start']) > len(times['stop']):
                times['stop'].append(df.index[-1])

        ## detector start times
        elif mode in 'he3':

            start = self.tfile[settings.DET_NAMES['He3']['transitions']].cycleStartTime
            start = start.drop_duplicates()

            # setup output
            times = {'start': start,
                     'duration (s)': np.concatenate((np.diff(start), [run_stop]))
                    }
            times['stop'] = times['start'] + times['duration (s)']
            times['supercycle'] = self.tfile[settings.DET_NAMES['He3']['transitions']].superCycleIndex

        ## detector start times
        elif mode in 'li6':

            start = self.tfile[settings.DET_NAMES['Li6']['transitions']].cycleStartTime
            start = start.drop_duplicates()

            # setup output
            times = {'start': start,
                     'duration (s)': np.concatenate((np.diff(start), [run_stop])),
                    }
            times['stop'] = times['start'] + times['duration (s)']
            times['supercycle'] = self.tfile[settings.DET_NAMES['Li6']['transitions']].superCycleIndex

        # convert to dataframe
        times = pd.DataFrame(times)
        times['duration (s)'] = times.stop - times.start
        times.index.name = 'cycle'

        # save
        self.cycle_param['cycle_times'] = times
        return times
