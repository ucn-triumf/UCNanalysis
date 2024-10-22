# Overall settings. Import this, or override its values at every execution
# Derek Fujimoto
# Oct 2024

# path to the directory which contains the root files
datadir = "/data3/ucn/root_files"

# detector tree names
DET_NAMES = {'He3':{'hits':         'UCNHits_He3',
                    'charge':       'He3_Charge',
                    'rate':         'He3_Rate',
                    'transitions':  'RunTransitions_He3',
                    'hitsseq':      'hitsinsequence_he3',
                    'hitsseqcumul': 'hitsinsequencecumul_he3',
                    },
                'Li6':{'hits':         'UCNHits_Li-6',
                    'charge':       'Li6_Charge',
                    'rate':         'Li6_Rate',
                    'transitions':  'RunTransitions_Li-6',
                    'hitsseq':      'hitsinsequence_li6',
                    'hitsseqcumul': 'hitsinsequencecumul_li6',
                    },
            }

# needed slow control trees: for checking data quality
SLOW_TREES = ('BeamlineEpics', 'SequencerTree', 'LNDDetectorTree')

# data thresholds for checking data
DATA_CHECK_THRESH = {'beam_min_current': 0.1, # uA
                     'beam_max_current_std': 0.02, # uA
                     'max_bkgd_count_rate': 4, # fractional increase over DET_BKGD values
                     'min_total_counts': 100, # number of counts total
                     'pileup_cnt_per_ms': 3, # if larger than this, then pileup and delete
                     'pileup_within_first_s': 1, # time frame for pileup in each period
                    }

# default detector backgrounds - from 2019
DET_BKGD = {'Li6':     1.578,
            'Li6_err': 0.009,
            'He3':     0.0349,
            'He3_err': 0.0023}

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