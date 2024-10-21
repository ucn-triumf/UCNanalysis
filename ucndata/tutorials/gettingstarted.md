# Getting Started

The `ucndata` package is based around the [ucnrun] object. This object reads a root file into memory and organizes its contents for easy scripting of future analyses. These root files can be generated with the `midas2root` program, which is a part of the [ucn_detector_analyzer](https://github.com/ucn-triumf/ucn_detector_analyzer/tree/2024) code set.

### Table of Contents

* [`ucnrun.cycle_param`](#cycle_param)
* [`ucnrun.tfile`](#tfile)
* [Accessing cycles and periods](#cycles-and-periods)

## Basics

Here is a minimum working example, assuming we have the file `ucn_run_00001846.root` in our current working directory.

```python
In [1]: from ucndata import ucnrun
In [2]: run = ucnrun('ucn_run_00001846.root')
```

This prompts a loading bar showing the progress of the file load. File inspection is easy: in the interpreter, simply type out the loaded run and it will print its contents to screen:

```python
In [3]: run
Out[3]:
run 1846:
  comment            month              shifters           tfile
  cycle_param        run_number         start_time         year
  experiment_number  run_title          stop_time
```

Most of these variables are simple strings or integers but there are two important attributes you should pay attention to: `cycle_param` and `tfile`.

## `cycle_param`

This object is an [attrdict] defined in the [rootloader] package. This is simply a dictionary whose contents can be accessed either in the normal way, or as an attribute. This allows for tab-completion in the interpreter.

Let's inspect the contents of this attribute:

```python
In [4]: run.cycle_param
Out[4]: attrdict: {'nperiods', 'nsupercyc', 'enable', 'inf_cyc_enable', 'cycle', 'supercycle', 'valve_states', 'period_end_times', 'period_durations_s', 'ncycles', 'filter', 'cycle_times'}
```

These are the various settings and properties of each cycle and period.

* `nperiods`: Number of periods in each cycle
* `nsupercyc`: Number of supercycles contained in the run
* `enable`: Enable status of the sequencer
* `inf_cyc_enable`: Enable status of infinite cycles
* `cycle`: Cycle ID numbers
* `supercycle`: Supercycle ID numbers
* `valve_states`: Valve states in each period and cycle
* `period_end_times`: End time of each period in each cycle in epoch time
* `period_durations_s`: Duration in sections of each period in each cycle
* `ncycles`: Number of total cycles contained in the run
* `filter`: A list indicating how we should filter cycles. More on that in [filters](filters.md)
* `cycle_time`: The start and end times of each cycle

## `tfile`

This is a [tfile](https://github.com/ucn-triumf/rootloader/blob/main/docs/rootloader/tfile.md) object from the [rootloader] package. It is effectively an [attrdict] with some embellishment for reading root files. This is the object which contains all your data read from the file. Its contents can be either [ttree] objects (again, based on the [attrdict] object) or a pandas [DataFrame]. These can be converted relatively easily. First, inspecting with its contents as [ttree]s:

```python
In [5]: run.tfile
Out[5]:
contents:
    BeamlineEpics            RunTransitions_Li-6      hitsinsequence_he3
    CycleParamTree           SequencerTree            hitsinsequence_li6
    LNDDetectorTree          UCNHits_He3              hitsinsequencecumul_he3
    RunTransitions_He3       UCNHits_Li-6             hitsinsequencecumul_li6
    RunTransitions_He3Det2   header

In [6]: run.tfile.BeamlineEpics
Out[6]:
ttree branches:
    B1UT_CM01_RDCOND        B1U_HARP2_RDUPDATE      B1U_TNIM2_10MINAVG      B1U_YCB1_RDCUR
    B1UT_CM02_RDCOND        B1U_IV0_STATON          B1U_TNIM2_10MINTRIP     B1V_KICK_RDHICUR
    B1UT_LM50_RDLVL         B1U_IV2_STATON          B1U_TNIM2_10SECAVG      B1V_KSM_BONPRD
    B1UT_PT01_RDPRESS       B1U_PNG0_RDVAC          B1U_TNIM2_10SECTRIP     B1V_KSM_INSEQ
    B1UT_PT02_RDPRESS       B1U_PNG2_RDVAC          B1U_TNIM2_1SECAVG       B1V_KSM_PREDCUR
    B1UT_PT50_RDPRESS       B1U_Q1_STATON           B1U_TNIM2_1SECTRIP      B1V_KSM_RDBEAMOFF_VAL1
    B1U_B0_RDCUR            B1U_Q1_VT_RDCUR         B1U_TNIM2_5MINAVG       B1V_KSM_RDBEAMON_VAL1
    B1U_B0_STATON           B1U_Q2_RDCUR            B1U_TNIM2_RAW           B1V_KSM_RDFRCTN_VAL1
    B1U_COL2DOWN_RDTEMP     B1U_Q2_STATON           B1U_WTEMP1_RDTEMP       B1V_KSM_RDMODE_VAL1
    B1U_COL2LEFT_RDTEMP     B1U_SEPT_RDCUR          B1U_WTEMP2_RDTEMP       B1_FOIL_ADJCUR
    B1U_COL2RIGHT_RDTEMP    B1U_SEPT_STATON         B1U_XCB1_RDCUR          timestamp
    B1U_COL2UP_RDTEMP       B1U_TGTTEMP1_RDTEMP     B1U_YCB0_RDCUR
    B1U_HARP0_RDUPDATE      B1U_TGTTEMP2_RDTEMP     B1U_YCB0_STATON
```

Then converting to [DataFrame]s in-place:

```python
In [7]: run.to_dataframe()

In [8]: run.tfile.BeamlineEpics
Out[8]:
            B1UT_CM01_RDCOND  B1UT_CM02_RDCOND  ...  B1V_KSM_RDMODE_VAL1  B1_FOIL_ADJCUR
timestamp                                       ...
1572460997          0.018750           0.00000  ...                  0.0        0.000000
1572461002          0.000000           0.01875  ...                  0.0        2.151400
1572461007          0.021875           0.01250  ...                  0.0        2.151400
1572461012          0.012500           0.00000  ...                  0.0        2.151400
1572461017          0.000000           0.00000  ...                  0.0        2.151400
...                      ...               ...  ...                  ...             ...
1572466463          0.000000           0.01250  ...                  0.0       38.294899
1572466468          0.000000           0.00000  ...                  0.0       38.294899
1572466473          0.018750           0.00000  ...                  0.0       37.864700
1572466478          0.034375           0.00000  ...                  0.0       37.864700
1572466479          0.000000           0.01250  ...                  0.0       38.294899

[1093 rows x 49 columns]
```

We see that the [ucnrun] object knows to set the `timestamp` as the index. Thus, one can access values directly from their timestamp:

```python
In [9]: run.tfile.BeamlineEpics.loc[1572460997:1572461017]
Out[9]:
            B1UT_CM01_RDCOND  B1UT_CM02_RDCOND  ...  B1V_KSM_RDMODE_VAL1  B1_FOIL_ADJCUR
timestamp                                       ...
1572460997          0.018750           0.00000  ...                  0.0          0.0000
1572461002          0.000000           0.01875  ...                  0.0          2.1514
1572461007          0.021875           0.01250  ...                  0.0          2.1514
1572461012          0.012500           0.00000  ...                  0.0          2.1514
1572461017          0.000000           0.00000  ...                  0.0          2.1514

[5 rows x 49 columns]
```

In general if each [ttree] structure is simple enough to convert to a [DataFrame] it is recommended that one does so.


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
