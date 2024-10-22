# ucncycle

[Ucndata Index](./README.md#ucndata-index) / ucncycle

> Auto-generated documentation for [ucncycle](../ucncycle.py) module.

- [ucncycle](#ucncycle)
  - [ucncycle](#ucncycle-1)
    - [ucncycle.check_data](#ucncyclecheck_data)
    - [ucncycle.get_counts](#ucncycleget_counts)
    - [ucncycle.get_period](#ucncycleget_period)
    - [ucncycle.get_rate](#ucncycleget_rate)

## ucncycle

[Show source in ucncycle.py:24](../ucncycle.py#L24)

View for the data from a single UCN cycle

#### Notes

Any changes to the data frame affects all periods (for the time steps
contained in that period) and the containing run

#### Arguments

- `urun` *ucnrun* - object to pull cycle from
- `cycle` *int* - cycle number

#### Signature

```python
class ucncycle(ucnbase):
    def __init__(self, urun, cycle): ...
```

#### See also

- [ucnbase](./ucnbase.md#ucnbase)

### ucncycle.check_data

[Show source in ucncycle.py:108](../ucncycle.py#L108)

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

### ucncycle.get_counts

[Show source in ucncycle.py:212](../ucncycle.py#L212)

Get counts for a/each period

#### Arguments

- `detector` *str* - one of the keys to settings.DET_NAMES
- `period` *None|int* - if None get for entire cycle
                    elif < 0 get for each period
                    elif >=0 get for that period
- `bkgd` *float|None* - background counts
- `dbkgd(float|None)` - error in background counts
- `norm` *float|None* - normalize to this value
- `dnorm` *float|None* - error in normalization

#### Returns

- `tuple` - number of hits for each period and error

#### Signature

```python
def get_counts(
    self, detector, period=None, bkgd=None, dbkgd=None, norm=None, dnorm=None
): ...
```

### ucncycle.get_period

[Show source in ucncycle.py:276](../ucncycle.py#L276)

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

### ucncycle.get_rate

[Show source in ucncycle.py:300](../ucncycle.py#L300)

Get count rate for each period

#### Arguments

- `detector` *str* - one of the keys to settings.DET_NAMES
- `bkgd` *float|None* - background counts
- `dbkgd(float|None)` - error in background counts
- `norm` *float|None* - normalize to this value
- `dnorm` *float|None* - error in normalization

#### Returns

- `np.ndarray` - count rate each period and error

#### Signature

```python
def get_rate(self, detector, bkgd=None, dbkgd=None, norm=None, dnorm=None): ...
```