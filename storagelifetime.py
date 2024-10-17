# Get storage lifetime
# Derek Fujimoto
# Oct 2024

# TODO: Fix beam current output (why all zero?)
# TODO: Fix background fetching
# TODO: use memoryviews for cycles and periods
# TODO: improve fetching of value from run level

from ucndata import ucnrun, settings, read, merge
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# settings
settings.datadir = 'root_files'
production_period = 0
storageperiod = 1
countperiod = 2
backgroundperiod = 1
detector = 'Li6'

run_numbers = [1846] #, 1873, 1875, 1892, 1902, 1905, 1908, 1918, 1933, 1940, 1960, 1971, 1982, 1988, 1995, 2008, 2009, 2017, 2022, 2038, 2039]

# get all runs
runs = read(run_numbers)

# get all cycles
cycles = np.concatenate([r.get_cycle() for r in runs])

# get beam current and means
beam_currents = [c.get_period(production_period).beam_current_uA for c in cycles]
dbeam_currents = np.array([b.std() for b in beam_currents])
beam_currents = np.array([b.mean() for b in beam_currents])

# get storage durations and associated counts
storage_duration = [c.cycle_param.period_durations_s[storageperiod] for c in cycles]
counts_bkgd = np.array([c.get_counts(detector, period=backgroundperiod)for c in cycles])
dcounts_bkgd = counts_bkgd[:,1]
counts_bkgd = counts_bkgd[:,0]
counts2 = [c.get_counts(detector, period=countperiod,
                                               bkgd=(counts_bkgd[i], dcounts_bkgd[i]),
                                               norm=(beam_currents[i],dbeam_currents[i])
                                               ) for i, c in enumerate(cycles)]

storage_duration = np.array(storage_duration)
counts = np.array(counts2)

# sort
idx = np.argsort(storage_duration)
storage_duration = storage_duration[idx]
counts = counts[idx]

idx = counts[:,0] > 0
storage_duration = storage_duration[idx]
counts = counts[idx]

# fit
fn = lambda x, p0, tau: p0*np.exp(-x/tau)
par, cov = curve_fit(fn, storage_duration, counts[:,0], sigma=counts[:,1], absolute_sigma=True)
std = np.diag(cov)**0.5

# draw
plt.errorbar(storage_duration, counts[:, 0], counts[:, 1], fmt='.')
plt.plot(storage_duration, fn(storage_duration, *par))
plt.yscale('log')


# debugging
r = runs[0]