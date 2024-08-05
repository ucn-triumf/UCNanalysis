# ucndata

[Ucndata Index](./README.md#ucndata-index) / ucndata

> Auto-generated documentation for [ucndata](../ucndata.py) module.

- [ucndata](#ucndata)
  - [ucndata](#ucndata-1)
    - [ucndata().check_data](#ucndata()check_data)
    - [ucndata().copy](#ucndata()copy)
    - [ucndata().from_dataframe](#ucndata()from_dataframe)
    - [ucndata().get_cycle](#ucndata()get_cycle)
    - [ucndata().get_cycles_times](#ucndata()get_cycles_times)
    - [ucndata().get_hits_histogram](#ucndata()get_hits_histogram)
    - [ucndata().set_he3](#ucndata()set_he3)
    - [ucndata().set_li6](#ucndata()set_li6)
    - [ucndata().to_dataframe](#ucndata()to_dataframe)

## ucndata

[Show source in ucndata.py:14](../ucndata.py#L14)

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

#### Signature

```python
class ucndata(object):
    def __init__(self, filename, header_only=False): ...
```

### ucndata().check_data

[Show source in ucndata.py:117](../ucndata.py#L117)

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

[Show source in ucndata.py:145](../ucndata.py#L145)

Return a copy of this objet

#### Signature

```python
def copy(self): ...
```

### ucndata().from_dataframe

[Show source in ucndata.py:284](../ucndata.py#L284)

Convert self.tfile contents to rootfile struture types

#### Signature

```python
def from_dataframe(self): ...
```

### ucndata().get_cycle

[Show source in ucndata.py:156](../ucndata.py#L156)

Return a copy of this object, but trees are trimmed to only one cycle.

Note that this process converts all objects to dataframes

#### Arguments

- `cycle` *int* - cycle number

#### Returns

- [ucndata](#ucndata) - a copy of this object but with data from only one cycle.

#### Signature

```python
def get_cycle(self, cycle): ...
```

### ucndata().get_cycles_times

[Show source in ucndata.py:199](../ucndata.py#L199)

Get start and end times of each cycle from the sequencer

#### Returns

- `pd.DataFrame` - with columns "start", "stop", and "duration (s)". Values are in epoch time. Indexed by cycle id

#### Signature

```python
def get_cycles_times(self): ...
```

### ucndata().get_hits_histogram

[Show source in ucndata.py:223](../ucndata.py#L223)

Get histogram of UCNHits ttree times

#### Arguments

- `bin_ms` *int* - histogram bin size in milliseconds

#### Returns

- `tuple` - (bin_centers, histogram counts)

#### Signature

```python
def get_hits_histogram(self, bin_ms=100): ...
```

### ucndata().set_he3

[Show source in ucndata.py:262](../ucndata.py#L262)

Set name patterns to match He3 detector

#### Signature

```python
def set_he3(self): ...
```

### ucndata().set_li6

[Show source in ucndata.py:271](../ucndata.py#L271)

Set name patterns to match li6 detector

#### Signature

```python
def set_li6(self): ...
```

### ucndata().to_dataframe

[Show source in ucndata.py:280](../ucndata.py#L280)

Convert self.tfile contents to pd.DataFrame

#### Signature

```python
def to_dataframe(self): ...
```