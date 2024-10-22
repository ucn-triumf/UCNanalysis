# ucnbase

[Ucndata Index](./README.md#ucndata-index) / ucnbase

> Auto-generated documentation for [ucnbase](../ucnbase.py) module.

- [ucnbase](#ucnbase)
  - [ucnbase](#ucnbase-1)
    - [ucnbase.apply](#ucnbaseapply)
    - [ucnbase.beam_current_uA](#ucnbasebeam_current_ua)
    - [ucnbase.beam_off_s](#ucnbasebeam_off_s)
    - [ucnbase.beam_on_s](#ucnbasebeam_on_s)
    - [ucnbase.copy](#ucnbasecopy)
    - [ucnbase.from_dataframe](#ucnbasefrom_dataframe)
    - [ucnbase.get_hits](#ucnbaseget_hits)
    - [ucnbase.get_hits_histogram](#ucnbaseget_hits_histogram)
    - [ucnbase.to_dataframe](#ucnbaseto_dataframe)

## ucnbase

[Show source in ucnbase.py:11](../ucnbase.py#L11)

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
class ucnbase(object): ...
```

### ucnbase.apply

[Show source in ucnbase.py:140](../ucnbase.py#L140)

Apply function to each cycle

#### Arguments

fn_handle (function handle): function to be applied to each cycle

#### Returns

- `np.ndarray` - output of the function

#### Signature

```python
def apply(self, fn_handle): ...
```

### ucnbase.beam_current_uA

[Show source in ucnbase.py:231](../ucnbase.py#L231)

#### Signature

```python
@property
def beam_current_uA(self): ...
```

### ucnbase.beam_off_s

[Show source in ucnbase.py:256](../ucnbase.py#L256)

#### Signature

```python
@property
def beam_off_s(self): ...
```

### ucnbase.beam_on_s

[Show source in ucnbase.py:253](../ucnbase.py#L253)

#### Signature

```python
@property
def beam_on_s(self): ...
```

### ucnbase.copy

[Show source in ucnbase.py:151](../ucnbase.py#L151)

Return a copy of this objet

#### Signature

```python
def copy(self): ...
```

### ucnbase.from_dataframe

[Show source in ucnbase.py:222](../ucnbase.py#L222)

Convert self.tfile contents to rootfile struture types

#### Signature

```python
def from_dataframe(self): ...
```

### ucnbase.get_hits

[Show source in ucnbase.py:162](../ucnbase.py#L162)

Get times of ucn hits

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES

#### Returns

- `pd.DataFrame` - hits tree as a dataframe, only the values when a hit is registered

#### Signature

```python
def get_hits(self, detector): ...
```

### ucnbase.get_hits_histogram

[Show source in ucnbase.py:182](../ucnbase.py#L182)

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

### ucnbase.to_dataframe

[Show source in ucnbase.py:226](../ucnbase.py#L226)

Convert self.tfile contents to pd.DataFrame

#### Signature

```python
def to_dataframe(self): ...
```