# ucnperiod

[Ucndata Index](./README.md#ucndata-index) / ucnperiod

> Auto-generated documentation for [ucnperiod](../ucnperiod.py) module.

- [ucnperiod](#ucnperiod)
  - [ucnperiod](#ucnperiod-1)
    - [ucnperiod.get_counts](#ucnperiodget_counts)
    - [ucnperiod.get_hits](#ucnperiodget_hits)
    - [ucnperiod.get_rate](#ucnperiodget_rate)

## ucnperiod

[Show source in ucnperiod.py:21](../ucnperiod.py#L21)

Stores the data from a single UCN period from a single cycle

#### Arguments

- `ucycle` *ucncycle* - object to pull period from
- `period` *int* - period number

#### Signature

```python
class ucnperiod(ucnbase):
    def __init__(self, ucycle, period): ...
```

#### See also

- [ucnbase](./ucnbase.md#ucnbase)

### ucnperiod.get_counts

[Show source in ucnperiod.py:83](../ucnperiod.py#L83)

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

### ucnperiod.get_hits

[Show source in ucnperiod.py:145](../ucnperiod.py#L145)

Get times of ucn hits

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES

#### Returns

- `pd.DataFrame` - hits tree as a dataframe, only the values when a hit is registered

#### Signature

```python
def get_hits(self, detector): ...
```

### ucnperiod.get_rate

[Show source in ucnperiod.py:184](../ucnperiod.py#L184)

Get sum of ucn hits per unit time of period

#### Arguments

- `detector` *str* - one of the keys to self.DET_NAMES
- `bkgd` *float|None* - background counts
- `dbkgd(float|None)` - error in background counts
- `norm` *float|None* - normalize to this value
- `dnorm` *float|None* - error in normalization

#### Returns

- `tuple` - (count rate, error)

#### Signature

```python
def get_rate(self, detector, bkgd=None, dbkgd=None, norm=None, dnorm=None): ...
```