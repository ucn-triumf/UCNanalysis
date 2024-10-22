# Settings

[Ucndata Index](./README.md#ucndata-index) / Settings

> Auto-generated documentation for [settings](../settings.py) module.

#### Attributes

- `datadir` - path to the directory which contains the root files: '/data3/ucn/root_files'

- `DET_NAMES` - detector tree names: {'He3': {'hits': 'UCNHits_He3', 'charge': 'He3_Charge', 'rate': 'He3_Rate', 'transitions': 'RunTransitions_He3', 'hitsseq': 'hitsinsequence_he3', 'hitsseqcumul': 'hitsinsequencecumul_he3'}, 'Li6': {'hits': 'UCNHits_Li-6', 'charge': 'Li6_Charge', 'rate': 'Li6_Rate', 'transitions': 'RunTransitions_Li-6', 'hitsseq': 'hitsinsequence_li6', 'hitsseqcumul': 'hitsinsequencecumul_li6'}}

- `SLOW_TREES` - needed slow control trees: for checking data quality: ('BeamlineEpics', 'SequencerTree', 'LNDDetectorTree')

- `DATA_CHECK_THRESH` - data thresholds for checking data: {'beam_min_current': 0.1, 'beam_max_current_std': 0.02, 'max_bkgd_count_rate': 4, 'min_total_counts': 100, 'pileup_cnt_per_ms': 3, 'pileup_within_first_s': 1}

- `DET_BKGD` - default detector backgrounds - from 2019: {'Li6': 1.578, 'Li6_err': 0.009, 'He3': 0.0349, 'He3_err': 0.0023}


- [Settings](#settings)
  - [keyfilter](#keyfilter)

## keyfilter

[Show source in settings.py:43](../settings.py#L43)

Don't load all the data in each file, only that which is needed

#### Signature

```python
def keyfilter(name): ...
```