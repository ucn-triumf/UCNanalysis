# Ucndata

[Ucndata Index](./README.md#ucndata-index) / Ucndata

> Auto-generated documentation for [ucndata](../ucndata.py) module.

- [Ucndata](#ucndata)
  - [ucncycle](#ucncycle)
    - [ucncycle().check_data](#ucncycle()check_data)
    - [ucncycle().cycles](#ucncycle()cycles)
    - [ucncycle().get_counts](#ucncycle()get_counts)
    - [ucncycle().get_cycle](#ucncycle()get_cycle)
    - [ucncycle().get_period](#ucncycle()get_period)
    - [ucncycle().get_rate](#ucncycle()get_rate)
    - [ucncycle().periods](#ucncycle()periods)
  - [ucnperiod](#ucnperiod)
    - [ucnperiod().get_counts](#ucnperiod()get_counts)
    - [ucnperiod().get_period](#ucnperiod()get_period)
    - [ucnperiod().get_rate](#ucnperiod()get_rate)
    - [ucnperiod().periods](#ucnperiod()periods)
  - [ucnrun](#ucnrun)
    - [ucnrun().beam_current_uA](#ucnrun()beam_current_ua)
    - [ucnrun().beam_off_s](#ucnrun()beam_off_s)
    - [ucnrun().beam_on_s](#ucnrun()beam_on_s)
    - [ucnrun().check_data](#ucnrun()check_data)
    - [ucnrun().copy](#ucnrun()copy)
    - [ucnrun().cycles](#ucnrun()cycles)
    - [ucnrun().from_dataframe](#ucnrun()from_dataframe)
    - [ucnrun().get_cycle](#ucnrun()get_cycle)
    - [ucnrun().get_cycle_times](#ucnrun()get_cycle_times)
    - [ucnrun().get_hits](#ucnrun()get_hits)
    - [ucnrun().get_hits_histogram](#ucnrun()get_hits_histogram)
    - [ucnrun().souce_temperature_k](#ucnrun()souce_temperature_k)
    - [ucnrun().source_pressure_kpa](#ucnrun()source_pressure_kpa)
    - [ucnrun().to_dataframe](#ucnrun()to_dataframe)

## ucncycle

[Show source in ucndata.py:689](../ucndata.py#L689)

Stores the data from a single UCN cycle

#### Arguments

- `urun` *ucnrun* - object to pull cycle from
- `cycle` *int* - cycle number
- `cycle_times_args` - passed to urun.get_cycle_times

#### Signature

```python
class ucncycle(ucnrun):
    def __init__(self, urun, cycle, **cycle_times_args): ...
```

#### See also

- [ucnrun](#ucnrun)

### ucncycle().check_data

[Show source in ucndata.py:754](../ucndata.py#L754)

Run some checks to determine if the data is ok.

#### Arguments

- `raise_error` *bool* - if true, raise an error if check fails, else return false

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
def check_data(self, raise_error=False): ...
```

### ucncycle().cycles

[Show source in ucndata.py:815](../ucndata.py#L815)

Cannot get cycle from current cycle

#### Signature

```python
def cycles(self, *args, **kwargs): ...
```

### ucncycle().get_counts

[Show source in ucndata.py:819](../ucndata.py#L819)

Get counts for each period

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES
- `subtr_bkgd` *bool* - if true subtract background
- `norm_dur` *bool* - if true normalize to beam current

#### Returns

- `np.ndarray` - number of hits for each period and error

#### Signature

```python
def get_counts(self, detector, subtr_bkgd=True, norm_beam=False): ...
```

### ucncycle().get_cycle

[Show source in ucndata.py:832](../ucndata.py#L832)

Cannot get cycle from current cycle

#### Signature

```python
def get_cycle(self, *args, **kwargs): ...
```

### ucncycle().get_period

[Show source in ucndata.py:836](../ucndata.py#L836)

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

[Show source in ucndata.py:865](../ucndata.py#L865)

Get count rate for each period

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES
- `subtr_bkgd` *bool* - if true subtract background
- `norm_dur` *bool* - if true normalize to beam current

#### Returns

- `np.ndarray` - count rate each period and error

#### Signature

```python
def get_rate(self, detector, subtr_bkgd=True, norm_beam=False): ...
```

### ucncycle().periods

[Show source in ucndata.py:878](../ucndata.py#L878)

Periods generator, calls get_period

#### Signature

```python
def periods(self): ...
```



## ucnperiod

[Show source in ucndata.py:883](../ucndata.py#L883)

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

[Show source in ucndata.py:946](../ucndata.py#L946)

Get sum of ucn hits

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES
- `subtr_bkgd` *bool* - if true subtract background
- `norm_dur` *bool* - if true normalize to beam current

#### Returns

- `float` - number of hits

#### Signature

```python
def get_counts(self, detector, subtr_bkgd=True, norm_beam=False): ...
```

### ucnperiod().get_period

[Show source in ucndata.py:1010](../ucndata.py#L1010)

Cannot get period from current period

#### Signature

```python
def get_period(self, *args, **kwargs): ...
```

### ucnperiod().get_rate

[Show source in ucndata.py:986](../ucndata.py#L986)

Get sum of ucn hits per unit time of period

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES
- `subtr_bkgd` *bool* - if true subtract background
- `norm_dur` *bool* - if true normalize to beam current

#### Returns

- `float` - count rate

#### Signature

```python
def get_rate(self, detector, subtr_bkgd=True, norm_beam=False): ...
```

### ucnperiod().periods

[Show source in ucndata.py:1014](../ucndata.py#L1014)

Cannot get period from current period

#### Signature

```python
def periods(self): ...
```



## ucnrun

[Show source in ucndata.py:28](../ucndata.py#L28)

#### Attributes

- `DET_NAMES` - detector names: {'He3': {'hits': 'UCNHits_He3', 'charge': 'He3_Charge', 'rate': 'He3_Rate', 'transitions': 'RunTransitions_He3', 'hitsseq': 'hitsinsequence_he3', 'hitsseqcumul': 'hitsinsequencecumul_he3'}, 'Li6': {'hits': 'UCNHits_Li-6', 'charge': 'Li6_Charge', 'rate': 'Li6_Rate', 'transitions': 'RunTransitions_Li-6', 'hitsseq': 'hitsinsequence_li6', 'hitsseqcumul': 'hitsinsequencecumul_li6'}}

- `SLOW_TREES` - needed slow control trees: ('BeamlineEpics', 'SequencerTree', 'LNDDetectorTree')

- `DATA_CHECK_THRESH` - data thresholds for checking data: {'beam_min_current': 0.1, 'beam_max_current_std': 0.02, 'pileup_cnt_per_ms': 3, 'pileup_within_first_s': 1}

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

#### Signature

```python
class ucnrun(object):
    def __init__(self, run, header_only=False): ...
```

### ucnrun().beam_current_uA

[Show source in ucndata.py:653](../ucndata.py#L653)

#### Signature

```python
@property
def beam_current_uA(self): ...
```

### ucnrun().beam_off_s

[Show source in ucndata.py:678](../ucndata.py#L678)

#### Signature

```python
@property
def beam_off_s(self): ...
```

### ucnrun().beam_on_s

[Show source in ucndata.py:675](../ucndata.py#L675)

#### Signature

```python
@property
def beam_on_s(self): ...
```

### ucnrun().check_data

[Show source in ucndata.py:302](../ucndata.py#L302)

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

[Show source in ucndata.py:358](../ucndata.py#L358)

Return a copy of this objet

#### Signature

```python
def copy(self): ...
```

### ucnrun().cycles

[Show source in ucndata.py:369](../ucndata.py#L369)

Cycles generator, calls get_cycle

#### Signature

```python
def cycles(self): ...
```

### ucnrun().from_dataframe

[Show source in ucndata.py:644](../ucndata.py#L644)

Convert self.tfile contents to rootfile struture types

#### Signature

```python
def from_dataframe(self): ...
```

### ucnrun().get_cycle

[Show source in ucndata.py:374](../ucndata.py#L374)

Return a copy of this object, but trees are trimmed to only one cycle.

Note that this process converts all objects to dataframes

#### Arguments

- `cycle` *int* - cycle number, if None, get all cycles
- `cycle_times_args` - passed to get_cycle_times

#### Returns

ucncycle:
    if cycle > 0:  ucncycle object
    if cycle < 0 | None: a list ucncycle objects for all cycles

#### Signature

```python
def get_cycle(self, cycle=None, **cycle_times_args): ...
```

### ucnrun().get_cycle_times

[Show source in ucndata.py:400](../ucndata.py#L400)

Get start and end times of each cycle from the sequencer

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
def get_cycle_times(self, mode="matched"): ...
```

### ucnrun().get_hits

[Show source in ucndata.py:560](../ucndata.py#L560)

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

[Show source in ucndata.py:604](../ucndata.py#L604)

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

### ucnrun().souce_temperature_k

[Show source in ucndata.py:681](../ucndata.py#L681)

#### Signature

```python
@property
def souce_temperature_k(self): ...
```

### ucnrun().source_pressure_kpa

[Show source in ucndata.py:685](../ucndata.py#L685)

#### Signature

```python
@property
def source_pressure_kpa(self): ...
```

### ucnrun().to_dataframe

[Show source in ucndata.py:648](../ucndata.py#L648)

Convert self.tfile contents to pd.DataFrame

#### Signature

```python
def to_dataframe(self): ...
```