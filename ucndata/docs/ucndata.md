# Ucndata

[Ucndata Index](./README.md#ucndata-index) / Ucndata

> Auto-generated documentation for [ucndata](../ucndata.py) module.

- [Ucndata](#ucndata)
  - [tsubfile](#tsubfile)
  - [ucncycle](#ucncycle)
    - [ucncycle().check_data](#ucncycle()check_data)
    - [ucncycle().get_counts](#ucncycle()get_counts)
    - [ucncycle().get_period](#ucncycle()get_period)
    - [ucncycle().get_rate](#ucncycle()get_rate)
  - [ucnperiod](#ucnperiod)
    - [ucnperiod().get_counts](#ucnperiod()get_counts)
    - [ucnperiod().get_rate](#ucnperiod()get_rate)
  - [ucnrun](#ucnrun)
    - [ucnrun().apply](#ucnrun()apply)
    - [ucnrun().beam_current_uA](#ucnrun()beam_current_ua)
    - [ucnrun().beam_off_s](#ucnrun()beam_off_s)
    - [ucnrun().beam_on_s](#ucnrun()beam_on_s)
    - [ucnrun().check_data](#ucnrun()check_data)
    - [ucnrun().copy](#ucnrun()copy)
    - [ucnrun().from_dataframe](#ucnrun()from_dataframe)
    - [ucnrun().get_cycle](#ucnrun()get_cycle)
    - [ucnrun().get_hits](#ucnrun()get_hits)
    - [ucnrun().get_hits_histogram](#ucnrun()get_hits_histogram)
    - [ucnrun().set_cycle_filter](#ucnrun()set_cycle_filter)
    - [ucnrun().set_cycle_times](#ucnrun()set_cycle_times)
    - [ucnrun().to_dataframe](#ucnrun()to_dataframe)

## tsubfile

[Show source in ucndata.py:1209](../ucndata.py#L1209)

Wrapper for tfile which restricts access to values only within given times

#### Arguments

- `tfileobj` *tfile* - object to wrap
- `start` *int* - starting epoch time
- `stop` *int* - stopping epoch time

#### Signature

```python
class tsubfile(tfile):
    def __init__(self, tfileobj, start, stop): ...
```



## ucncycle

[Show source in ucndata.py:773](../ucndata.py#L773)

View for the data from a single UCN cycle

#### Notes

Any changes to the data frame affects all periods (for the time steps
contained in that period) and the containing run

#### Arguments

- `urun` *ucnrun* - object to pull cycle from
- `cycle` *int* - cycle number

#### Signature

```python
class ucncycle(ucnrun):
    def __init__(self, urun, cycle): ...
```

#### See also

- [ucnrun](#ucnrun)

### ucncycle().check_data

[Show source in ucndata.py:857](../ucndata.py#L857)

Run some checks to determine if the data is ok.

#### Arguments

- `period_production` *int* - index of period where the beam should be stable. Enables checks of beam stability
- `period_count` *int* - index of period where we count ucn. Enables checks of data quantity
- `period_background` *int* - index of period where we do not count ucn. Enables checks of background
- `raise_error` *bool* - if true, raise an error if check fails, else return false
- `quiet` *bool* - if true don't print or raise exception

#### Returns

- `bool` - true if check passes, else false.

Checks:
    Do the following trees exist and have entries?
        BeamlineEpics
        UCN2Epics
        SequencerTree
        LNDDetectorTree
    Are there nonzero counts in UCNHits?

#### Signature

```python
def check_data(
    self,
    period_production=None,
    period_count=None,
    period_background=None,
    raise_error=False,
    quiet=False,
): ...
```

### ucncycle().get_counts

[Show source in ucndata.py:962](../ucndata.py#L962)

Get counts for each period

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES
- `period` *None|int* - if None get for entire cycle
                    elif < 0 get for each period
                    elif >=0 get for that period
- `bkgd` *tuple|None* - if not None subtract this as the background (value, error)
                   bkgd.shape = (2, nperiods) if period < 0 else (2, 1)
- `norm` *tuple|None* - if not None normalize to this value (value, error)
                    norm.shape = (2, nperiods) if period < 0 else (2, 1)

#### Returns

- `np.ndarray` - number of hits for each period and error

#### Signature

```python
def get_counts(self, detector, period=None, bkgd=None, norm=None): ...
```

### ucncycle().get_period

[Show source in ucndata.py:1021](../ucndata.py#L1021)

Return a copy of this object, but trees are trimmed to only one period.

#### Notes

This process converts all objects to dataframes
Must be called for a single cycle only

#### Arguments

- `period` *int* - period number, if None, get all periods
cycle (int|None) if cycle not specified then specify a cycle

#### Returns

run:
    if period > 0: a copy of this object but with data from only one period.
    if period < 0 | None: a list of copies of this object for all periods for a single cycle

#### Signature

```python
def get_period(self, period=None): ...
```

### ucncycle().get_rate

[Show source in ucndata.py:1046](../ucndata.py#L1046)

Get count rate for each period

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES
- `bkgd` *tuple|None* - if not None subtract this as the background (value, error)
- `norm` *tuple|None* - if not None normalize to this value (value, error)

#### Returns

- `np.ndarray` - count rate each period and error

#### Signature

```python
def get_rate(self, detector, bkgd=True, norm=False): ...
```



## ucnperiod

[Show source in ucndata.py:1060](../ucndata.py#L1060)

Stores the data from a single UCN period from a single cycle

#### Arguments

- `ucycle` *ucncycle* - object to pull period from
- `period` *int* - period number

#### Signature

```python
class ucnperiod(ucncycle):
    def __init__(self, ucycle, period): ...
```

#### See also

- [ucncycle](#ucncycle)

### ucnperiod().get_counts

[Show source in ucndata.py:1122](../ucndata.py#L1122)

Get sum of ucn hits

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES
- `bkgd` *float|None* - background counts
- `dbkgd(float|None)` - error in background counts
- `norm` *float|None* - normalize to this value
- `dnorm` *float|None* - error in normalization

#### Returns

- `tuple` - (count, error) number of hits

#### Signature

```python
def get_counts(self, detector, bkgd=None, dbkgd=None, norm=None, dnorm=None): ...
```

### ucnperiod().get_rate

[Show source in ucndata.py:1185](../ucndata.py#L1185)

Get sum of ucn hits per unit time of period

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES
- `bkgd` *tuple|None* - if not None subtract this as the background (value, error)
- `norm` *tuple|None* - if not None normalize to this value (value, error)

#### Returns

- `float` - count rate

#### Signature

```python
def get_rate(self, detector, bkgd=True, norm=False): ...
```



## ucnrun

[Show source in ucndata.py:26](../ucndata.py#L26)

#### Attributes

- `DET_NAMES` - detector names: {'He3': {'hits': 'UCNHits_He3', 'charge': 'He3_Charge', 'rate': 'He3_Rate', 'transitions': 'RunTransitions_He3', 'hitsseq': 'hitsinsequence_he3', 'hitsseqcumul': 'hitsinsequencecumul_he3'}, 'Li6': {'hits': 'UCNHits_Li-6', 'charge': 'Li6_Charge', 'rate': 'Li6_Rate', 'transitions': 'RunTransitions_Li-6', 'hitsseq': 'hitsinsequence_li6', 'hitsseqcumul': 'hitsinsequencecumul_li6'}}

- `SLOW_TREES` - needed slow control trees: ('BeamlineEpics', 'SequencerTree', 'LNDDetectorTree')

- `DATA_CHECK_THRESH` - data thresholds for checking data: {'beam_min_current': 0.1, 'beam_max_current_std': 0.02, 'max_bkgd_count_rate': 4, 'min_total_counts': 100, 'pileup_cnt_per_ms': 3, 'pileup_within_first_s': 1}

- `DET_BKGD` - default detector backgrounds - from 2019: {'Li6': 1.578, 'Li6_err': 0.009, 'He3': 0.0349, 'He3_err': 0.0023}


UCN run data. Cleans data and performs analysis

#### Arguments

- `run` *int|str* - if int, generate filename with settings.datadir
    elif str then run is the path to the file
- `header_only` *bool* - if true, read only the header

#### Attributes

- `comment` *str* - comment input by users
- `cycle` *int|none* - cycle number, none if no cycle selected
- `cycle_param` *attrdict* - cycle parameters from sequencer settings
- `experiment_number` *str* - experiment number input by users
- `month` *int* - month of run start
- `run_number` *int* - run number
- `run_title` *str* - run title input by users
- `shifter` *str* - experimenters on shift at time of run
- `start_time` *str* - start time of the run
- `stop_time` *str* - stop time of the run
- `supercycle` *int|none* - supercycle number, none if no cycle selected
- `tfile` *tfile* - stores tfile raw readback
- `year` *int* - year of run start

#### Notes

Can access attributes of tfile directly from top-level object
Need to define the values in ucndata.settings if you want non-default
behaviour
Object is indexed as [cycle, period] for easy access to sub time frames

#### Signature

```python
class ucnrun(object):
    def __init__(self, run, header_only=False): ...
```

### ucnrun().apply

[Show source in ucndata.py:362](../ucndata.py#L362)

Apply function to each cycle

#### Arguments

fn_handle (function handle): function to be applied to each cycle

#### Returns

- `np.ndarray` - output of the function

#### Signature

```python
def apply(self, fn_handle): ...
```

### ucnrun().beam_current_uA

[Show source in ucndata.py:745](../ucndata.py#L745)

#### Signature

```python
@property
def beam_current_uA(self): ...
```

### ucnrun().beam_off_s

[Show source in ucndata.py:770](../ucndata.py#L770)

#### Signature

```python
@property
def beam_off_s(self): ...
```

### ucnrun().beam_on_s

[Show source in ucndata.py:767](../ucndata.py#L767)

#### Signature

```python
@property
def beam_on_s(self): ...
```

### ucnrun().check_data

[Show source in ucndata.py:373](../ucndata.py#L373)

Run some checks to determine if the data is ok.

#### Arguments

- `raise_error` *bool* - if true, raise an error if check fails, else return false

#### Returns

- `bool` - true if check passes, else false.

Checks:
    Do the self.SLOW_TREES exist and have entries?
    Are there nonzero counts in UCNHits?

#### Signature

```python
def check_data(self, raise_error=False): ...
```

### ucnrun().copy

[Show source in ucndata.py:430](../ucndata.py#L430)

Return a copy of this objet

#### Signature

```python
def copy(self): ...
```

### ucnrun().from_dataframe

[Show source in ucndata.py:545](../ucndata.py#L545)

Convert self.tfile contents to rootfile struture types

#### Signature

```python
def from_dataframe(self): ...
```

### ucnrun().get_cycle

[Show source in ucndata.py:441](../ucndata.py#L441)

Return a copy of this object, but trees are trimmed to only one cycle.

Note that this process converts all objects to dataframes

#### Arguments

- `cycle` *int* - cycle number, if None, get all cycles

#### Returns

ucncycle:
    if cycle > 0:  ucncycle object
    if cycle < 0 | None: a list ucncycle objects for all cycles

#### Signature

```python
def get_cycle(self, cycle=None): ...
```

### ucnrun().get_hits

[Show source in ucndata.py:461](../ucndata.py#L461)

Get times of ucn hits

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES

#### Returns

- `pd.DataFrame` - hits tree as a dataframe, only the values when a hit is registered

#### Signature

```python
def get_hits(self, detector): ...
```

### ucnrun().get_hits_histogram

[Show source in ucndata.py:505](../ucndata.py#L505)

Get histogram of UCNHits ttree times

#### Arguments

- `detector` *str* - Li6|He3
- `bin_ms` *int* - histogram bin size in milliseconds

#### Returns

- `tuple` - (bin_centers, histogram counts)

#### Signature

```python
def get_hits_histogram(self, detector, bin_ms=100): ...
```

### ucnrun().set_cycle_filter

[Show source in ucndata.py:549](../ucndata.py#L549)

Set filter for which cycles to fetch when slicing or iterating

#### Arguments

- `cfilter` *None|iterable* - list of bool, True if keep cycle, False if reject.
    if None then same as if all True

#### Returns

- `None` - sets self.cycle_param.filter

#### Signature

```python
def set_cycle_filter(self, cfilter=None): ...
```

### ucnrun().set_cycle_times

[Show source in ucndata.py:569](../ucndata.py#L569)

Get start and end times of each cycle from the sequencer and save
into self.cycle_param.cycle_times

Run this if you want to change how cycle start times are calculated

#### Arguments

- `mode` *str* - matched|sequencer|he3|li6
    - `if` *matched* - look for identical timestamps in RunTransitions from detectors
    - `if` *sequencer* - look for inCycle timestamps in SequencerTree
    - `if` *he3* - use He3 detector cycle start times
    - `if` *li6* - use Li6 detector cycle start times

#### Notes

- If run ends before sequencer stop is called, a stop is set to final timestamp.
- If the sequencer is disabled mid-run, a stop is set when disable ocurrs.
- If sequencer is not enabled, then make the entire run one cycle
- For matched mode,
    - set run stops as start of next transition
    - set offset as start_He3 - start_Li6
    - set start/stop/duration based on start_He3
- If the object reflects a single cycle, return from cycle_start, cycle_stop

#### Returns

- `pd.DataFrame` - with columns "start", "stop", "offset" and "duration (s)". Values are in epoch time. Indexed by cycle id. Offset is the difference in detector start times: he3_start-li6_start

#### Signature

```python
def set_cycle_times(self, mode="matched"): ...
```

### ucnrun().to_dataframe

[Show source in ucndata.py:740](../ucndata.py#L740)

Convert self.tfile contents to pd.DataFrame

#### Signature

```python
def to_dataframe(self): ...
```