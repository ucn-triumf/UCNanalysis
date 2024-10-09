# Overall settings. Import this, or override its values at every execution
# Derek Fujimoto
# Oct 2024

# path to the directory which contains the root files
datadir = "root_data"

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