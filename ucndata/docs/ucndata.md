# ucndata

[Ucndata Index](./README.md#ucndata-index) / ucndata

> Auto-generated documentation for [ucndata](../ucndata.py) module.

- [ucndata](#ucndata)
  - [ucndata](#ucndata-1)
    - [ucndata().beam_current_uA](#ucndata()beam_current_ua)
    - [ucndata().beam_off_s](#ucndata()beam_off_s)
    - [ucndata().beam_on_s](#ucndata()beam_on_s)
    - [ucndata().check_data](#ucndata()check_data)
    - [ucndata().copy](#ucndata()copy)
    - [ucndata().from_dataframe](#ucndata()from_dataframe)
    - [ucndata().get_cycle](#ucndata()get_cycle)
    - [ucndata().get_cycle_times](#ucndata()get_cycle_times)
    - [ucndata().get_hits_histogram](#ucndata()get_hits_histogram)
    - [ucndata().souce_temperature_k](#ucndata()souce_temperature_k)
    - [ucndata().source_pressure_kpa](#ucndata()source_pressure_kpa)
    - [ucndata().to_dataframe](#ucndata()to_dataframe)

## ucndata

[Show source in ucndata.py:31](../ucndata.py#L31)

#### Attributes

- `DET_NAMES` - detector names: {'He3': {'hits': 'UCNHits_He3', 'charge': 'He3_Charge', 'rate': 'He3_Rate', 'transitions': 'RunTransitions_He3', 'hitsseq': 'hitsinsequence_he3', 'hitsseqcumul': 'hitsinsequencecumul_he3'}, 'Li6': {'hits': 'UCNHits_Li-6', 'charge': 'Li6_Charge', 'rate': 'Li6_Rate', 'transitions': 'RunTransitions_Li-6', 'hitsseq': 'hitsinsequence_li6', 'hitsseqcumul': 'hitsinsequencecumul_li6'}}

- `SLOW_TREES` - needed slow control trees: ('BeamlineEpics', 'SequencerTree', 'LNDDetectorTree')


UCN run data. Cleans data and performs analysis

#### Arguments

- `cycle` *int|None* - indicates cycle number, none if full file
- `cycle_start` *float* - epoch time cycle start time (only if single cycle)
cycle_stop (float) : epoch time cycle stop time (only if single cycle)
- `filename` *str* - path to file to open
- `header_only` *bool* - if true, read only the header

#### Attributes

- `comment` *str* - comment input by users
- `cycle` *int|none* - cycle number, none if no cycle selected
- `experiment_number` *str* - experiment number input by users
- `month` *int* - month of run start
- `run_number` *int* - run number
- `run_title` *str* - run title input by users
- `shifter` *str* - experimenters on shift at time of run
- `start_time` *str* - start time of the run
- `stop_time` *str* - stop time of the run
- `tfile` *tfile* - stores tfile raw readback
- `year` *int* - year of run start

#### Notes

Can access attributes of tfile directly from top-level object

#### Signature

```python
class ucndata(object):
    def __init__(self, filename, header_only=False): ...
```

### ucndata().beam_current_uA

[Show source in ucndata.py:453](../ucndata.py#L453)

#### Signature

```python
@property
def beam_current_uA(self): ...
```

### ucndata().beam_off_s

[Show source in ucndata.py:478](../ucndata.py#L478)

#### Signature

```python
@property
def beam_off_s(self): ...
```

### ucndata().beam_on_s

[Show source in ucndata.py:475](../ucndata.py#L475)

#### Signature

```python
@property
def beam_on_s(self): ...
```

### ucndata().check_data

[Show source in ucndata.py:181](../ucndata.py#L181)

Run some checks to determine if the data is ok.

Checks:
    Do the following trees exist and have entries?
        BeamlineEpics
        UCN2Epics
        SequencerTree
        LNDDetectorTree
    Are there nonzero counts in UCNHits?

#### Signature

```python
def check_data(self): ...
```

### ucndata().copy

[Show source in ucndata.py:215](../ucndata.py#L215)

Return a copy of this objet

#### Signature

```python
def copy(self): ...
```

### ucndata().from_dataframe

[Show source in ucndata.py:444](../ucndata.py#L444)

Convert self.tfile contents to rootfile struture types

#### Signature

```python
def from_dataframe(self): ...
```

### ucndata().get_cycle

[Show source in ucndata.py:226](../ucndata.py#L226)

Return a copy of this object, but trees are trimmed to only one cycle.

Note that this process converts all objects to dataframes

#### Arguments

- `cycle` *int* - cycle number, if None, get all cycles
- `cycle_times_args` - passed to get_cycle_times

#### Returns

ucndata:
    if cycle > 0: a copy of this object but with data from only one cycle.
    if cycle < 0: a list of copies of this object for all cycles

#### Signature

```python
def get_cycle(self, cycle=None, **cycle_times_args): ...
```

### ucndata().get_cycle_times

[Show source in ucndata.py:277](../ucndata.py#L277)

Get start and end times of each cycle from the sequencer

#### Arguments

- `mode` *str* - matched|sequencer
    - `if` *matched* - look for identical timestamps in RunTransitions from detectors
    - `if` *sequencer* - look for inCycle timestamps in SequencerTree

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

- `pd.DataFrame` - with columns "start", "stop", and "duration (s)". Values are in epoch time. Indexed by cycle id

#### Signature

```python
def get_cycle_times(self, mode="matched"): ...
```

### ucndata().get_hits_histogram

[Show source in ucndata.py:404](../ucndata.py#L404)

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

### ucndata().souce_temperature_k

[Show source in ucndata.py:481](../ucndata.py#L481)

#### Signature

```python
@property
def souce_temperature_k(self): ...
```

### ucndata().source_pressure_kpa

[Show source in ucndata.py:485](../ucndata.py#L485)

#### Signature

```python
@property
def source_pressure_kpa(self): ...
```

### ucndata().to_dataframe

[Show source in ucndata.py:448](../ucndata.py#L448)

Convert self.tfile contents to pd.DataFrame

#### Signature

```python
def to_dataframe(self): ...
```