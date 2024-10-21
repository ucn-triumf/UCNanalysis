# Using ApplyLists

An [applylist] is an object defined in the ucndata package. It inherits from the python [list](https://docs.python.org/3/tutorial/datastructures.html) object. However we often want to get an attribute or call a function from the contents of lists of objects. For example, if you want to get the beam current for period 0 of each cycle, and the counts in period 2, you would need to do something like:

```python
from ucndata import ucnrun
run = ucnrun('ucn_run_00001846.root')

beam_current_mean = [c[0].beam_current_uA.mean() for c in run]
beam_current_std =  [c[0].beam_current_uA.std() for c in run]
counts = [c[2].get_counts(detector='Li6') for c in run]
```

Fetching more parameters would result in many list comprehensions. Since the slicing operation on [ucnrun] and [ucncycle] objects returns an [applylist] we can instead do the following:

```python
beam_current_mean = run[:, 0].beam_current_uA.mean()
beam_current_std  = run[:, 0].beam_current_uA.std()
counts = run[:, 2].get_counts(detector='Li6')
```

Breaking it down:

1. `run[:, 0]` returns an [applylist] of [ucnperiod]s corresponding to period 0 of each cycle.
2. Since `beam_current_uA` is not an attribute of the [applylist] itself, the list instead fetches this attribute from each of the contained items and returns a new [applylist] with the result (a pandas [DataFrame]).
3. Since `mean()` is also not an attribute of the [applylist] it instead tries to call this method on each of the contained items and returns a new [applylist] with the result (a float).

The [applylist] object works similarly to a [numpy array](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html) in that it enables slicing, index-wise comparison and arithmetic, transpose, and type conversions; but also has an `apply()` function borrowed from pandas [DataFrame]s which applies a function to each element.

Note that this works on nested [applylist]s. Therefore the following how to fetch the mean beam current for each period in each cycle:

```python
run[:, :].beam_current_uA.mean()
```


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
