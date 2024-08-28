# Open and analyze UCN data
# Derek Fujimoto
# June 2024

"""
    TODO List, things which haven't been ported from WS code

    * Filter pileup
    * skip cycles that have total duration of 0s
    * skip cycle because cycle sequence is longer than irradiation interval
    * get supercycle index
    * skip cycles with no beam data
    * skip cycles with fluctuating beam current
    * skip cycles with low beam current
    * skip cycles with no detector counts
    * skip cycles no IV1 open during measurement
    * get temperature
    * get vapour pressure
"""

from rootloader import tfile, ttree
from .exceptions import *
import ucndata.constants as const
import ROOT
import numpy as np
import pandas as pd
import itertools, warnings, os

ROOT.gROOT.SetBatch(1)

class ucndata(object):
    """UCN run data. Cleans data and performs analysis

    Args:
        cycle (int|None): indicates cycle number, none if full file
        cycle_start (float): epoch time cycle start time (only if single cycle)
        cycle_stop (float) : epoch time cycle stop time (only if single cycle)
        filename (str): path to file to open
        header_only (bool): if true, read only the header

    Attributes:
        comment (str): comment input by users
        cycle (int|none): cycle number, none if no cycle selected
        experiment_number (str): experiment number input by users
        month (int): month of run start
        run_number (int): run number
        run_title (str): run title input by users
        shifter (str): experimenters on shift at time of run
        start_time (str): start time of the run
        stop_time (str): stop time of the run
        tfile (tfile): stores tfile raw readback
        year (int): year of run start

    Notes:
        Can access attributes of tfile directly from top-level object
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

    def __init__(self, filename, header_only=False):

        if filename is None:
            return

        # read
        if header_only:
            fid = ROOT.TFile(filename, 'READ')
            head = ttree(fid.Get('header'))
            fid.Close()
            head = {k:str(val[0]) for k, val in head.items()}
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

        if type(self.run_number) is pd.Series:
            self.run_number = int(self.run_number[0])
        else:
            self.run_number = int(self.run_number)

        # set other header items
        self.cycle = None
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

            # unsure why 0.00088801 is needed...
            start_times = abs(beam.timestamp - start)
            durations = getattr(beam, epics_val)*const.beam_bucket_duration_s
            idx = np.argmin(start_times)
            beam_dur.append(durations[idx])

        out = pd.Series(beam_dur, index=cycle_times.index)
        out.index.name = cycle_times.index.name
        return out

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

        # check some necessary data trees
        for tree in self.SLOW_TREES:

            # does tree exist?
            if tree not in self.tfile.keys():
                raise MissingDataError(f'Missing ttree "{tree}" in run {self.run_number}')

            # does tree have entries?
            if self.tfile[tree].entries == 0:
                raise MissingDataError(f'Zero entries found in "{tree}" ttree in run {self.run_number}')

        for name, det in self.DET_NAMES.items():

            # check for nonzero counts
            if not self.tfile[det['hits']].tIsUCN.any():
                raise MissingDataError(f'No UCN hits in "{name}" ttree in run {self.run_number}')

            # check if sequencer was enabled but no run transitions
            if any(self.tfile.SequencerTree.sequencerEnabled):
                if self.tfile[det['transitions']].entries == 0:
                    raise MissingDataError('No cycles found in run {self.run_number}, although sequencer was active')

    def copy(self):
        """Return a copy of this objet"""
        copy = ucndata(None)

        for key, value in self.__dict__.items():
            if hasattr(value, 'copy'):
                setattr(copy, key, value.copy())
            else:
                setattr(copy, key, value)
        return copy

    def get_cycle(self, cycle=None, **cycle_times_args):
        """Return a copy of this object, but trees are trimmed to only one cycle.

        Note that this process converts all objects to dataframes

        Args:
            cycle (int): cycle number, if None, get all cycles
            cycle_times_args: passed to get_cycle_times

        Returns:
            ucndata:
                if cycle > 0: a copy of this object but with data from only one cycle.
                if cycle < 0: a list of copies of this object for all cycles
        """

        # get all cycles
        if cycle is None:
            ncycles = len(self.get_cycle_times(**cycle_times_args).index)
            return [self.get_cycle(c) for c in range(ncycles)]

        # make copy
        copy = self.copy()

        # get cycles to keep
        cycles = copy.get_cycle_times(**cycle_times_args)
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

    def get_cycle_times(self, mode='matched'):
        """Get start and end times of each cycle from the sequencer

        Args:
            mode (str): matched|sequencer
                if matched: look for identical timestamps in RunTransitions from detectors
                if sequencer: look for inCycle timestamps in SequencerTree

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
            pd.DataFrame: with columns "start", "stop", and "duration (s)". Values are in epoch time. Indexed by cycle id
        """

        # check if single cycle
        if self.cycle is not None:
            return pd.DataFrame({'start':[self.cycle_start],
                                 'stop':[self.cycle_stop],
                                 'duration (s)': [self.cycle_stop-self.cycle_start],
                                 'offset (s)': [0.0]},
                                 index=[self.cycle])

        # get data
        df = self.tfile.SequencerTree
        if type(df) == ttree:
            df = df.to_dataframe()

        ## if no sequencer, make the whole run a single cycle
        if not any(df.sequencerEnabled):

            times = {'start': np.inf,
                     'stop': -np.inf}

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

            # get run end time from control trees
            run_stop = -np.inf
            for treename in self.SLOW_TREES:
                idx = self.tfile[treename].to_dataframe().index
                run_stop = max((idx.max(), run_stop))

            # setup output
            times = {'start': matchedhe3,
                     'duration (s)': np.concatenate((np.diff(matchedhe3), [run_stop])),
                     'offset (s)': matchedhe3-matchedli6}
            times['stop'] = times['start'] + times['duration (s)']

        ## get timestamps from sequencer
        elif mode in 'sequencer':

            # if sequencer is not enabled cause a stop transition
            df.inCycle *= df.sequencerEnabled

            # start counting only after first start flag
            df = df.loc[df.loc[df.cycleStarted > 0].index[0]:]

            # get start and end times
            df = df.diff()
            times = {'start': df.index[df.inCycle == 1],
                    'stop': df.index[df.inCycle == -1]}

            # check lengths
            if len(times['start']) > len(times['stop']):
                times['stop'].append(df.index[-1])

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