# Open and analyze UCN data
# Derek Fujimoto
# June 2024

"""
    TODO List, things which haven't been ported from WS code

    * Filter pileup
    * skip cycles that have total duration of 0s
    * skip cycle because cycle sequence is longer than irradiation interval
    * get supercycle index
    * skip cycles with no detector counts
    * get temperature
    * get vapour pressure
    * period extraction and data checks
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

class udata(object):
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
    DATA_CHECK_THRESH = {'beam_min_current':        0.1,    # uA
                         'beam_max_current_std':    0.02,   # uA
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

        # reformat cycle param tree
        paramtree = self.tfile.CycleParamTree
        self.cycle_param = {'ncycles':  paramtree.nCycles[0],
                            'nperiods': paramtree.nPeriods[0],
                            'nsupercyc': paramtree.nSuperCyc[0],
                            'enable': bool(paramtree.enable[0]),
                            'inf_cyc_enable': bool(paramtree.infCyclesEnable[0]),
                            }
        self.cycle_param = attrdict(self.cycle_param)

        # setup cycle paramtree array outputs
        df_param = paramtree.to_dataframe()
        df_param.set_index('periodNumber', inplace=True)
        df_param.index.name = 'period'

        # get valve states
        df = df_param[[col for col in df_param.columns if 'Valve' in col]]

        col_map = {col:int(col.replace("Valve", "").replace("State", "")) for col in df.columns}
        df = df.rename(columns=col_map)
        df.columns.name = 'valve'
        df = df.astype(bool)
        self.cycle_param['valve_states'] = df

        # get period durations
        df = df_param[[col for col in df_param.columns if 'periodDurationInCycle' in col]]
        col_map = {col:int(col.replace("periodDurationInCycle", "")) for col in df.columns}
        df = df.rename(columns=col_map)
        df.columns.name = 'cycle'
        df = df[sorted(df.columns)]
        self.cycle_param['period_durations'] = df

        # set other header items
        self.cycle = None
        self.supercycle = None
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
            cyc_str = '' if self.cycle is None else f' (cycle {self.cycle})'
            s = f'run {self.run_number}{cyc_str}:\n'
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

        # individual cycle checks
        if self.cycle is not None:

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
            # elif self.tfile.???.UCN_UGD_IV1_STATON.max() < 1:
            #     msg = f'{run_msg} IV1 never opened'
            #     err = ValveError

            # raise error or return value
            if msg is not None:
                if raise_error:
                    raise err(msg)
                else:
                    print(msg)
                    return False

        return True

    def copy(self):
        """Return a copy of this objet"""
        copy = udata(None)

        for key, value in self.__dict__.items():
            if hasattr(value, 'copy'):
                setattr(copy, key, value.copy())
            else:
                setattr(copy, key, value)
        return copy

    def correct_bkgd(self, detector, rate=None, rate_err=None):
        """Subtract background and normalize to the average beam current

        Args:
            detector (str): He3|Li6
            rate (float): background rate. If None, use self.DET_BKGD
            rate_err (float): error in background rate. If None, use self.DET_BKGD
        """

        # get rates
        if rate is None: rate = self.DET_BKGD[detector]
        if rate_err is None: rate_err = self.DET_BKGD[detector+'_err']

        # TODO: needs finishing. See UCN.py: SubtractBackgroundAndNormalize

    def get_cycle(self, cycle=None, nproc=-1, **cycle_times_args):
        """Return a copy of this object, but trees are trimmed to only one cycle.

        Note that this process converts all objects to dataframes

        Args:
            cycle (int): cycle number, if None, get all cycles
            nproc (int): number of processors to use
            cycle_times_args: passed to get_cycle_times

        Returns:
            udata:
                if cycle > 0: a copy of this object but with data from only one cycle.
                if cycle < 0 | None: a list of copies of this object for all cycles
        """

        # get all cycles
        if cycle is None or cycle < 0:
            ncycles = len(self.get_cycle_times(**cycle_times_args).index)
            return list(map(self.get_cycle, tqdm(range(ncycles),
                                                 total=ncycles,
                                                 leave=False,
                                                 desc='Fetch all cycles')
                            )
                        )

        # make copy
        copy = self.copy()

        # get cycles to keep
        cycles = copy.get_cycle_times(**cycle_times_args)
        start = int(cycles.loc[cycle, 'start'])
        stop = int(cycles.loc[cycle, 'stop'])
        supercycle = int(cycles.loc[cycle, 'supercycle'])

        # trim the trees
        for key, value in copy.tfile.items():
            if key == 'header': continue

            # trim ttree
            if type(value) is ttree:
                value = value.to_dataframe()
                if value.index.name is not None:
                    idx = (value.index < stop) & (value.index > start)
                    value = value.loc[idx]
                copy.tfile[key] = ttree(value)

            # trim dataframe
            elif type(value) is pd.DataFrame:
                if value.index.name is not None:
                    idx = (value.index < stop) & (value.index > start)
                    copy.tfile[key] = value.loc[idx].copy()

        # trim cycle parameters
        copy.cycle_param.period_durations = copy.cycle_param.period_durations[cycle]

        copy.cycle = cycle
        copy.supercycle = supercycle
        copy.cycle_start = start
        copy.cycle_stop = stop
        return copy

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
        if self.cycle is not None:
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
            times['supercycle'] = self.tfile[self.DET_NAMES['He3']['transitions']].superCycleIndex

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


