# Cycles and Periods

[**Back to Index**](index.md)\
[**Next Page: Filtering Cycles**](filter.md)

---

Each run which uses the sequencer makes use of many cycles, each of which is composed of up to 10 periods. Each period has a set of valves open or closed for a fixed amount of time. Each cycle can have different period settings, but each supercycle reflects a repetition of the full cycle set. It is therefore prudent to access the data for each cycles and period in the run in an easy manner.

## Table of Contents

* [Basics](#basics-and-introduction)
* [Slicing and indexing](#slicing-and-indexing)
* [Determining cycle start and end times](#cycle-start-and-end-times)

## Basics and Introduction

The overarching concept is to access the entire contents of the data set in a new object each time a cycle or period is accessed. These new objects, [ucncycle] and [ucnperiod] respectively, each keep track of their own timings and identity but do not copy the contents of [tfile]. Rather, when the user attempts to access the contents of [tfile] the [ucncycle] or [ucnperiod] object instead fetches only part of the contents as needed. Thus, no large copying of data is needed, improving run times significantly. This means, however, that modifying data in [ucncycle] or [ucnperiod] objects may modify the data elsewhere, since this is a shared property. Otherwise, each of the [ucncycle] or [ucnperiod] objects behave for the most part the same as the containing ucnrun object.

A simple example of usage:

```python
In [0]: from ucndata import ucnrun
In [1]: run = ucnrun('ucn_run_00001846.root')

In [2]: run.get_cycle(0)
Out[2]:
run 1846 (cycle 0):
  comment            cycle_start        month              shifters           supercycle
  cycle              cycle_stop         run_number         start_time         tfile
  cycle_param        experiment_number  run_title          stop_time          year
```

As you can see the cycle knows it's cycle 0, and will tell the user as such. It also gains the attributes

* `cycle_start`: the start time of the cycle in epoch time
* `cycle_stop`: the stop time of the cycle in epoch time
* `cycle`: the cycle id index

However its contents are now restricted to the time frame associated with that cycle:

```python
In [3]: run.get_cycle(0).tfile.BeamlineEpics
Out[3]:
            B1UT_CM01_RDCOND  B1UT_CM02_RDCOND  ...  B1V_KSM_RDMODE_VAL1  B1_FOIL_ADJCUR
timestamp                                       ...
1572461640          0.000000          0.000000  ...                  0.0       40.876598
1572461645          0.009375          0.015625  ...                  0.0       40.876598
1572461650          0.000000          0.003125  ...                  0.0       40.876598
1572461655          0.000000          0.000000  ...                  0.0       40.876598
...

[58 rows x 49 columns]
```

Similarly, once a cycle is fetched, one can then access the periods within with [`ucncycle.get_period()`](../docs/ucndata.md#ucncycle). One can fetch all the cycles/periods by passing no parameter (or None) to the function.

## Slicing and Indexing

Since accessing the cycle and period views of the [ucnrun] are such a common thing to need in an analysis, the [ucnrun] object can be indexed as if it were a 2-dimensional array. The indexing follows the scheme of `[cycle, period]` and employs slicing. Thus,

```python
# the following statement
run.get_cycle(0)

# is equivalent to
run[0]
```

Similarly, fetching the entire cycle list is easily reduced:

```python
# the following statement
run.get_cycle()

# is equivalent to
run[:]
```

While this doesn't yet seem to be too beneficial, the true advantage comes when we want to start fetching periods:

```python
# the following statement
run.get_cycle(0).get_period(0)

# is equivalent to
run[0, 0]
```

And more so when we want to get a list of all the periods and cycles:

```python
# the following statements
list_of_cycles = []
for i in range(run.cycle_param.ncycles):
    list_of_periods = []
    for j in range(run.cycle_param.nperiods):
        list_of_periods.append(run.get_cycle(i).get_period(j))
    list_of_cycles.append(list_of_periods)

# is equivalent to
list_of_cycles = run[:, :]
```

All the usual slicing rules and syntax apply. Some examples:

```python
run[0, :]       # fetch all periods from cycle 0
run[:, 0]       # fetch period 0 from all cycles
run[2:5]        # fetch cycles 2, 3, and 4
run[2:5, 0]     # fetch period 0 from cycles 2, 3, and 4
run[2:5, :2]    # fetch periods 0 and 1 from cycles 2, 3, and 4
```

You can also treat the [ucnrun] object as an iterator for notational simplicity:

```python
for cycle in run:
    print(cycle)
```

## Cycle Start and End Times

Cycle start and end times can be calculated in a few different ways. When the [ucnrun] object is created, during its construction it calls [ucnrun.set_cycle_times()](../docs/ucndata.md#ucnrunset_cycle_times) which in turn calculates the start and end times of each cycle. Fetching cycles uses the output of this function (saved to `ucnrun.cycle_param.cycle_times`) to determine what time range is associated with each cycle. To change how each cycle is determined, call [ucnrun.set_cycle_times()](../docs/ucndata.md#ucnrunset_cycle_times) with different parameters. The strategies are as follows:

* **matched (default)**: look at He3 and Li6 detector transitions and find pairs which are the closest in time to each other. Start times are then set to the He3 transition state and the difference is saved as the offset. Raises a warning if unmatched pairs exist.
* **sequencer**: Look at the `inCycle` flag of the `SequencerTree` and determine times from this parameter
* **he3**: Uses the He3 detector transitions only (`RunTransitions_He3`)
* **li6**: Uses the Li6 detector transitions only (`RunTransitions_Li-6`)

For runs without a sequencer, this function returns the run start and stop times as the cycle timing.

---

[**Back to Index**](index.md)\
[**Next Page: Filtering Cycles**](filter.md)

[tfile]: #tfile
[DataFrame]: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html
[ttree]:https://github.com/ucn-triumf/rootloader/blob/main/docs/rootloader/ttree.md
[attrdict]:https://github.com/ucn-triumf/rootloader/blob/main/docs/rootloader/attrdict.md
[rootloader]: https://github.com/ucn-triumf/rootloader
[ucnrun]: ../docs/ucnrun.md
[ucncycle]: ../docs/ucncycle.md
[ucnperiod]: ../docs/ucnperiod.md
[applylist]: ../docs/applylist.md
[read]: ../docs/read.md
[merge]: ../docs/merge.md
