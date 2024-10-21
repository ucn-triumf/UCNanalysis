# Filtering

In the case of imperfect data we often want to not load every cycle present in the file. Of course one can use the [slicing](../docs/cycandperiods.md#slicing-and-indexing) feature of the [ucnrun] or [applylist] objects to achieve this, however we offer an alternative.

The `ucnrun.cycle_param` dictionary has a `filter` item which sets the cycles filter. Effectively, this is a list of length `ucnrun.cycle_param.ncycles` containing boolean values which indicate whether a cycle should be kept or not (True to keep the cycle).

This filter should be set using the [ucnrun.set_cycle_filter()](../docs/ucndata.md#set_cycle_filter) method. A default filter can be generated with the [ucnrun.gen_cycle_filter()](../docs/ucndata.md#gen_cycle_filter) method.

## When is filtering applied

Filtering is applied **ONLY** when slicing as described in the [cycles and periods page](cycandperiods.md#slicing-and-indexing) or when iterating through the run in the following way:

```python
for cycle in run:
    # do something...
```

It does **NOT** apply when attempting to access a single cycle or using the [`get_cycle()`](ucndata.md#ucnrunget_cycle) function.

## How it works

When slicing, both the filter list and the list of cycles are sliced, then the reduced filter is applied to the reduced list of cycles. Some examples:

```python
# setup filter to permit fetching of all cycles
run.set_cycle_filter(np.full(run.cycle_param.ncycles, True))

# print the cycle numbers fetched
print(run[:].cycle)
# output: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

# don't allow cycle 1 to be fetched
run.cycle_param.filter[1] = False

# print the cycle numbers fetched
print(run[:].cycle)
# output: [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

# this can result in some odd behaviour
run[1:3] # this fetches only cycle 2, since cycle 1 is filtered
run[:4] # this returns a list that is of length 3, not 4
```

## Why??

This lets a user rerun an analysis with a new cycle filter with minimal changes to their code.





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
