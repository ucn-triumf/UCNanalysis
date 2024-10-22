# Loading Runs

[**Back to Index**](index.md)\
[**Next Page: Accessing Cycles and Periods**](cycandperiods.md)

---

## Specification

There are two main ways to specify which run to load:

1. Pass a filename
2. Pass a run number

The latter is convenient and useful for making readable analysis scripts:

```python
from ucndata import ucnrun

# we want to go from this:
run = ucnrun('/data3/ucn/root_files/ucn_run_00002050.root')

# to this:
run = ucnrun(2050)
```

The way to do this is to set the default path in the `ucnrun.settings` file:

```python
from ucndata import ucnrun, settings

# redefine settings.datadir
settings.datadir = '/data3/ucn/root_files'

# now this works
run = ucnrun(2050)
```

## Efficient Loading

We don't need to load all the data in the root file. This can slow the loading process by quite a bit. The `settings` file also defines the function `keyfilter` which tells the rootloader which objects to read into memory. The default is as follows:

```python
def keyfilter(name):
    """Don't load all the data in each file, only that which is needed"""

    name = name.replace(' ', '_').lower()

    # reject some keys based on partial matches
    reject_keys = ('v1725', 'v1720', 'v792', 'tv1725', 'charge', 'edge_diff',
                   'pulse_widths', 'iv2', 'iv3')
    for key in reject_keys:
        if key in name:
            return False

    return True
```

But one can define this in the same way as with the data directory:

```python
from ucndata import settings

# this function loads all the data in the file
settings.keyfilter = lambda x: True
```

Note that by default empty trees and histograms are not loaded into memory.

If one needs to read many runs it is recommended that one uses the [read] function:

```python
from ucndata import read

runlist = read([2050, 2051, 2052])
```

---

[**Back to Index**](index.md)\
[**Next Page: Accessing Cycles and Periods**](cycandperiods.md)

[tfile]: #tfile
[DataFrame]: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html
[ttree]:https://github.com/ucn-triumf/rootloader/blob/main/docs/rootloader/ttree.md
[attrdict]:https://github.com/ucn-triumf/rootloader/blob/main/docs/rootloader/attrdict.md
[rootloader]: https://github.com/ucn-triumf/rootloader
[ucnrun]: ../docs/ucndata.md#ucnrun
[ucncycle]: ../docs/ucndata.md#ucncycle
[ucnperiod]: ../docs/ucndata.md#ucnperiod
[applylist]: ../docs/applylist.md
[read]: ../docs/read.md
[merge]: ../docs/merge.md
