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
import ucndata.settings as default_settings
import ucndata.constants as const
import ROOT
import numpy as np
import pandas as pd
import itertools, warnings, os

ROOT.gROOT.SetBatch(1)

class ucnrun(object):
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

    def __init__(self, run, header_only=False):

        # check if copying
        if run is None:
            return

        # make filename from defaults
        elif type(run) is int:
            try:
                _dirname = datadir
            except NameError:
                _dirname = default_settings.datadir

            filename = os.path.join(_dirname, f'ucn_run_{run:0>8d}.root')

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

            try:
                _keyfilter = keyfilter
            except NameError:
                _keyfilter = default_settings.keyfilter

            self.tfile = tfile(filename, empty_ok=False, quiet=True,
                               key_filter=_keyfilter)
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
        for detector in self.DET_NAMES.values():
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

    def apply(self, fn_handle):
        """Apply function to each cycle

        Args:
            fn_handle (function handle): function to be applied to each cycle

        Returns:
            np.ndarray: output of the function
        """
        return [fn_handle(c) for c in self.cycles()]

    def check_data(self, raise_error=False):
        """Run some checks to determine if the data is ok.

        Args:
            raise_error (bool): if true, raise an error if check fails, else return false

        Returns:
            bool: true if check passes, else false.

        Checks:
            Do the self.SLOW_TREES exist and have entries?
            Are there nonzero counts in UCNHits?
        """

        # check some necessary data trees
        for tree in self.SLOW_TREES:

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

        for name, det in self.DET_NAMES.items():

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

    def copy(self):
        """Return a copy of this objet"""
        copy = ucnrun(None)

        for key, value in self.__dict__.items():
            if hasattr(value, 'copy'):
                setattr(copy, key, value.copy())
            else:
                setattr(copy, key, value)
        return copy

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

        # filter pileup for period data
        if type(self) is ucnperiod:

            # get thresholds
            dt = self.DATA_CHECK_THRESH['pileup_within_first_s']
            count_thresh = self.DATA_CHECK_THRESH['pileup_cnt_per_ms']

            # make histogram
            t = hit_tree.index.values
            counts, edges = np.histogram(t,
                                         bins=int(1/0.001),
                                         range=(min(t),
                                                min(t)+dt))

            # delete bad count ranges
            ncounts_total = len(t)
            for i, count in enumerate(counts):
                if count > count_thresh:
                    hit_tree= hit_tree.loc[(hit_tree.index < edges[i]) & (hit_tree.index > edges[i+1])]
            ncounts_removed = ncounts_total - len(hit_tree.index)

            if ncounts_removed > 0:
                print(f'Removed {ncounts_removed} pileup counts ({int(ncounts_removed/ncounts_total*100):d}%) from run{self.run_number} (cycle {self.cycle}, period {self.period})')

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
        for treename in self.SLOW_TREES:
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
            for treename in self.SLOW_TREES:
                idx = self.tfile[treename].to_dataframe().index
                times['start'] = min((idx.min(), times['start']))
                times['stop']  = max((idx.max(), times['stop']))

        ## get matched timesteps from He3 and Li6 RunTransitions
        elif mode in 'matched':
            hestart = self.tfile[self.DET_NAMES['He3']['transitions']].cycleStartTime
            listart = self.tfile[self.DET_NAMES['Li6']['transitions']].cycleStartTime
            scycle = self.tfile[self.DET_NAMES['He3']['transitions']].superCycleIndex

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

            start = self.tfile[self.DET_NAMES['He3']['transitions']].cycleStartTime
            start = start.drop_duplicates()

            # setup output
            times = {'start': start,
                     'duration (s)': np.concatenate((np.diff(start), [run_stop]))
                    }
            times['stop'] = times['start'] + times['duration (s)']
            times['supercycle'] = self.tfile[self.DET_NAMES['He3']['transitions']].superCycleIndex

        ## detector start times
        elif mode in 'li6':

            start = self.tfile[self.DET_NAMES['Li6']['transitions']].cycleStartTime
            start = start.drop_duplicates()

            # setup output
            times = {'start': start,
                     'duration (s)': np.concatenate((np.diff(start), [run_stop])),
                    }
            times['stop'] = times['start'] + times['duration (s)']
            times['supercycle'] = self.tfile[self.DET_NAMES['Li6']['transitions']].superCycleIndex

        # convert to dataframe
        times = pd.DataFrame(times)
        times['duration (s)'] = times.stop - times.start
        times.index.name = 'cycle'

        # save
        self.cycle_param['cycle_times'] = times
        return times

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

class ucncycle(ucnrun):
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

class ucnperiod(ucncycle):
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
            detector (str): one of the keys to self.DET_NAMES
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

    # TODO: FIX THIS
    def get_rate(self, detector, bkgd=True, norm=False):
        """Get sum of ucn hits per unit time of period

        Args:
            detector (str): one of the keys to self.DET_NAMES
            bkgd (tuple|None): if not None subtract this as the background (value, error)
            norm (tuple|None): if not None normalize to this value (value, error)

        Returns:
            float: count rate
        """
        # get counts
        counts, dcounts = self.get_counts(detector=detector,
                                          bkgd=bkgd,
                                          norm=norm)

        # get rate
        duration = self.cycle_param.period_durations_s

        counts /= duration
        dcounts /= duration

        return (counts, dcounts)

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

