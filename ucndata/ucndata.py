# Open and analyze UCN data for a whole run
# Derek Fujimoto
# June 2024

"""
    TODO List, things which haven't been ported from WS code

    * skip cycles that have total duration of 0s
    * skip cycle because cycle sequence is longer than irradiation interval
    * skip cycles with no detector counts
    * get temperature
    * get vapour pressure
    * data checks for periods
"""

from rootloader import tfile, ttree, attrdict
from .exceptions import *
import ucndata.settings as default_settings
import ucndata.constants as const
import ROOT
import numpy as np
import pandas as pd
import itertools, warnings, os
from tqdm import tqdm

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

    def _get_beam_duration(self, on=True):
        # Get beam on/off durations

        # get needed info
        cycle_times = self.get_cycle_times()

        try:
            beam = self.tfile.BeamlineEpics
        except AttributeError:
            raise MissingDataError("No saved ttree named BeamlineEpics")

        # setup storage
        beam_dur = []
        epics_val = 'B1V_KSM_RDBEAMON_VAL1' if on else 'B1V_KSM_RDBEAMOFF_VAL1'

        # get durations closest to cycle start time
        for start in cycle_times.start:

            start_times = abs(beam.timestamp - start)
            durations = getattr(beam, epics_val)*const.beam_bucket_duration_s
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

    def _trim(self, start, stop):
        # trim trees such that timestamps are between start and stop

        for key, value in self.tfile.items():
            if key == 'header': continue

            # trim ttree
            if type(value) is ttree:
                value = value.to_dataframe()
                if value.index.name is not None:
                    idx = (value.index < stop) & (value.index > start)
                    value = value.loc[idx]
                self.tfile[key] = ttree(value)

            # trim dataframe
            elif type(value) is pd.DataFrame:
                if value.index.name is not None:
                    idx = (value.index < stop) & (value.index > start)
                    self.tfile[key] = value.loc[idx].copy()

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
            elif self.tfile[tree].entries == 0:
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
                if self.tfile[det['transitions']].entries == 0:
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

    def cycles(self):
        """Cycles generator, calls get_cycle"""
        for i in range(self.cycle_param.ncycles):
            yield self.get_cycle(i)

    def get_cycle(self, cycle=None, **cycle_times_args):
        """Return a copy of this object, but trees are trimmed to only one cycle.

        Note that this process converts all objects to dataframes

        Args:
            cycle (int): cycle number, if None, get all cycles
            cycle_times_args: passed to get_cycle_times

        Returns:
            ucncycle:
                if cycle > 0:  ucncycle object
                if cycle < 0 | None: a list ucncycle objects for all cycles
        """

        if cycle is None or cycle < 0:
            ncycles = len(self.get_cycle_times(**cycle_times_args).index)
            return list(map(self.get_cycle, tqdm(range(ncycles),
                                                 total=ncycles,
                                                 leave=False,
                                                 desc='Fetch all cycles')
                            )
                        )
        else:
            return ucncycle(self, cycle, **cycle_times_args)

    def get_cycle_times(self, mode='matched'):
        """Get start and end times of each cycle from the sequencer

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
            idx = self.tfile[treename].to_dataframe().index
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

        return times

    def get_hits(self, detector):
        """Get times of ucn hits

        Args:
            detector (str): one of the keys to self.DET_NAMES

        Returns:
            pd.DataFrame: hits tree as a dataframe, only the values when a hit is registered
        """

        # get the tree
        hit_tree = self.tfile[self.DET_NAMES[detector]['hits']]
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

        # index_col = 'tUnixTimePrecise'

        # get cycle start and stop times from ucn counts histogram

        # get data
        df = self.tfile[self.DET_NAMES[detector]['hits']].to_dataframe()
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

    @property
    def souce_temperature_k(self):
        raise NotImplementedError()

    @property
    def source_pressure_kpa(self):
        raise NotImplementedError()

class ucncycle(ucnrun):
    """Stores the data from a single UCN cycle

    Args:
        urun (ucnrun): object to pull cycle from
        cycle (int): cycle number
        cycle_times_args: passed to urun.get_cycle_times
    """

    def __init__(self, urun, cycle, **cycle_times_args):

        # copy data
        for key, value in urun.__dict__.items():
            if hasattr(value, 'copy'):
                setattr(self, key, value.copy())
            else:
                setattr(self, key, value)

        # get cycles to keep
        cycles = urun.get_cycle_times(**cycle_times_args)
        start = int(cycles.loc[cycle, 'start'])
        stop = int(cycles.loc[cycle, 'stop'])
        supercycle = int(cycles.loc[cycle, 'supercycle'])

        # trim the trees
        self._trim(start, stop)

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

    def check_data(self, raise_error=False):
        """Run some checks to determine if the data is ok.

        Args:
            raise_error (bool): if true, raise an error if check fails, else return false

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

        # run full data checks
        super().check_data(raise_error=raise_error)

        # setup error message
        msg = None
        run_msg = f'Run {self.run}, cycle {self.cycle} failure:'

        # beam data exists
        beam_current = self.beam_current_uA
        if len(beam_current) == 0:
            msg = f'{run_msg} No beam data saved'
            err = BeamError

        # beam current too low
        elif beam_current.min() < self.DATA_CHECK_THRESH['beam_min_current']:
            msg = f'{run_msg} Beam current dropped to {beam_current.min()} uA'
            err = BeamError

        # beam current unstable
        elif beam_current.std() > self.DATA_CHECK_THRESH['beam_max_current_std']:
            msg = f'{run_msg} Beam current fluctuated by {beam_current.std()} uA'
            err = BeamError

        # valve states
        elif not self.cycle_param.valve_states.any().any():
            msg = f'{run_msg} No valves operated'
            err = ValveError

        # has counts
        elif not any([self.tfile[self.DET_NAMES['hits']].tIsUCN.sum() > 1 for det in self.DET_NAMES.keys()]):
            msg = f'{run_msg} No counts detected'
            err = DataError

        # raise error or return value
        if msg is not None:
            if raise_error:
                raise err(msg)
            else:
                print(msg)
                return False

        return True

    def cycles(self, *args, **kwargs):
        """Cannot get cycle from current cycle"""
        raise RuntimeError('Object already reflets a single cycle')

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
        if period > nperiods:
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
            if type(counts) in (int, np.int64): zero = 0
            else:                   zero = np.zeros(len(counts))
            counts = np.max(counts-bkgd[0], zero)
            dcounts = (dcounts**2 + bkgd[1]**2)**0.5

        # normalize
        if norm is not None:
            dcounts = counts*((dcounts/counts)**2 + (norm[1]/norm[0])**2)**0.5
            counts /= norm[0]

        return (counts, dcounts)

    def get_cycle(self, *args, **kwargs):
        """Cannot get cycle from current cycle"""
        raise RuntimeError('Object already reflets a single cycle')

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
            return list(map(self.get_period, tqdm(range(nperiods),
                                                 total=nperiods,
                                                 leave=False,
                                                 desc='Fetch all periods')
                            )
                        )
        else:
            return ucnperiod(self, period)

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

    def periods(self):
        """Periods generator, calls get_period"""
        for i in range(self.cycle_param.nperiods):
            yield self.get_period(i)

class ucnperiod(ucncycle):
    """Stores the data from a single UCN period from a single cycle

    Args:
        ucycle (ucncycle): object to pull period from
        period (int): period number
    """

    def __init__(self, ucycle, period):

        # copy data
        for key, value in ucycle.__dict__.items():
            if hasattr(value, 'copy'):
                setattr(self, key, value.copy())
            else:
                setattr(self, key, value)

        # get start and stop time
        if period > 0:      start = int(self.cycle_param.period_end_times[period-1])
        else:               start = int(self.cycle_start)

        stop = int(self.cycle_param.period_end_times[period])

        # trim the trees
        self._trim(start, stop)

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

    def get_bkgd(self, detector):
        """Get default background counts

        Args:
            detector (str): one of the keys to self.DET_NAMES

        Returns:
            np.ndarray: (background, error)
        """

         # get background counts
        brate = self.DET_BKGD[detector]
        brate_err = self.DET_BKGD[detector+'_err']

        bcount = brate * self.cycle_param.period_durations_s
        bcount_err = brate_err * self.cycle_param.period_durations_s

        return np.array((bcount, bcount_err))

    def get_counts(self, detector, bkgd=None, norm=None):
        """Get sum of ucn hits

        Args:
            detector (str): one of the keys to self.DET_NAMES
            bkgd (tuple|None): if not None subtract this as the background (value, error)
            norm (tuple|None): if not None normalize to this value (value, error)

        Returns:
            tuple: (count, error) number of hits
        """
        hit_tree = self.get_hits(detector)
        counts = len(hit_tree.index)
        dcounts = np.sqrt(counts) # error assumed poissonian

        # subtract background
        if bkgd is not None:
            counts = max(counts-bkgd[0], 0)
            dcounts = (dcounts**2 + bkgd[1]**2)**2

        # normalize
        if norm is not None:
            dcounts = counts*((dcounts/counts)**2 + (norm[1]/norm[0])**2)**0.5
            counts /= norm[0]

        return (counts, dcounts)

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
