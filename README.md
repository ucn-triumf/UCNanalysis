# UCNanalysis2018

This a collection of Python scripts to analyze data from the 2018 run.

## extractcycles.py:

Takes ROOT files (generated from midas runs using midas2root) as command-line parameters. Extracts the relevant data for all cycles in each run and stores it in ucn_output.root.
The output file contains a folder "DetectorRates" containing histograms of detector rates from all cycles and a tree "cycledata" with one entry for each cycle with the following structure:

* start (integer, timestamp of cycle start time [s])
* runnumber (integer)
* cyclenumber (integer)
* beamonduration (float, [s])
* beamoffduration (float, [s])
* timingoffset (float, offset between Li6 and He3 timestamps in the raw data for debugging [s])
* periods (branch containing data for each period in cycle)
  * periodsLen (integer, length of arrays in this branch for this cycle, equal to number of periods (11))
  * durations (array of 11 floats, duration of each period [s])
  * countsLi6 (array of 11 integers, total counts in Li6 detector during each period)
  * countsHe3 (array of 11 integers, total counts in He3 detector during each period)
  * valveXstate (array of 11 integers, setting for up to 8 valves in each period [0,1])
* Li6 (branch containing Li6 hits)
  * Li6Len (integer, length of arrays in this branch for this cycle)
  * hits (array of Li6Len floats, timestamps of Li6 hits during this cycle)
* He3 (branch containing He3 hits, similar to Li6 branch)
* Source (branch containing Epics data from source)
  * SourceLen (integer, length of arrays in this branch for this cycle)
  * timestamp (array of SourceLen integers, timestamps when the Epics variables were sampled)
  * UCN_EPICS_VARIABLE (array of SourceLen floats, contains data for Epics variables, e.g. IV1 state, temperatures, vapor pressures)
* Beamline (branch containing Epics data from beamline, similar to Source branch)
* SCM (branch containing data from SCM, similar to Source branch)
  * SCMVoltages3 (array of SCMLen floats, voltage drop of SCM current flowing through a 250uOhm shunt resistor [V])
* LND (branch containing data from LND detector, similar to Source branch)

Run with
`python extractcycles.py /data/ucn/root_files_20190114/ucn_tree_0000*.root`

The latest output, generated from files in /data/ucn/root_files_20190114 on daq01, is found at https://ucn.triumf.ca/ucn-source/ucnanalysis2018/ucn_output_20190123.root

## transmission.py

This script takes transmission experiments with two periods per cycle (irradiation + counting) that were performed with a monitor detector available during irradiation.
It takes the counts in the He3 detector during irradiation and the counts in the Li6 detector during counting.
It subtracts a fixed background rate in the Li6 detector of 2.16 +/- 0.01 per second.
Then it determines the ratio of background-corrected Li6 counts to 3He counts and prints the weighted average over all cycles and saves a pdf file showing the ratio and average for all cycles.
The He3 detector is assumed to be background-free.

It also plots the background in the Li6 detector during the last ten seconds of all cycles in each run.

Run with
`python transmission.py ucn_output_20190123.root`

The latest output, generated from ucn_output_20190123.root, is found in the folder transmission.

## storagelifetime.py

This script takes storage-lifetime experiments with three periods per cycle (irradiation, storage, counting) that were performed without a monitor detector available during irradiation.
It takes the counts in either the He3 or the Li6 detector (whichever was used for counting) and determines the storage lifetime in two ways:
1. Subtract the background rate, determined during the storage period and averaged over all cycles, and divide the background-corrected detector counts by the average beam current for ach cycle. Plot against duration of the storage period. A single-exponential fit determines the storage lifetime.
2. Plot the uncorrected and unnormalized counts against duration of the storage period and fit a single-exponential with background.

Results are plotted into pdf files.
It also plots the standard storagelifetimes over time (TCN18-015), storage lifetime vs temperature/pressure (TCN18-300), storage lifetime while spoiling the source (TCN18-170), and a histogram of backgrounds from all analyzed runs.

Run with
`python storagelifetime.py ucn_output_20190123.root`

The latest output, generated from ucn_out_20190123.root, is found in the folder storagelifetime.

## storagelifetime_with_monitor.py

This script takes storage-lifetime experiments with three periods per cycle (irradiation, storage, counting) with the He3 detector used as a monitor detector during irradiation.
The He3 detector is assumed to be background-free.

The storage lifetime is determined in three ways:
1. Take the counts in the Li6 detector during counting and subtract the background rate determined during the storage period and averaged over all cycles. Divide by the counts in the He3 detector during irradiation and plot the ratio against the duration of the storage period. Single-exponential fits (including or excluding the measurement with 0s storage time) or double exponential fits fit determine the storage lifetimes.
2. Plot the uncorrected and unnormalized Li6 counts against the duration of the storage period. Fit a single exponential with background.
3. Fit the raw rate in the He3 detector during storage with a single exponential (pinhole method). This only makes sense if the monitor detector is actually connected to the storage volume (e.g. it is useless for storage between IV2 and IV3).

The results are printed to pdf files.
It also plots the background rate in the Li6 detector during the storage time from all analyzed runs.

Run with
`python storagelifetime_with_monitor.py ucn_output_20190123.root`

The latest output, generated from ucn_output_20190123.root, is found in the folder storagelifetime_with_monitor.

## time_of_flight.py

This script takes the same transmission experiments as transmission.py and plots the rate in the Li6 detector normalized to the counts in the He3 detector during irradiation, averaged over all cycles of each experiment.
The result is basically a time-of-flight spectrum. All spectra are printed to pdfs. It also can divide time-of-flight spectra, e.g. to normalize a transmission spectrum to a reference experiment.

Run with
`python time_of_flight.py ucn_output_20190123.root`

The latest output, generated from ucn_output_20190123.root, is found in the folder time_of_flight.

## pyROOT crash course

```python
import ROOT # import ROOT library

f = ROOT.TFile('ucn_output.root') # open ROOT file
cycles = f.Get('cycledata') # get tree

for cycle in cycles: # loop over entries in tree
  start = cycle.start # access leaves of tree
  Li6hits = getattr(cycle, 'Li6/hits') # syntax cycle.Li6.hits to access sub-branches unfortunately not (yet?) supported by pyROOT.
  rate = ROOT.TH1I('rate', 'rate', 180, 0, 180) # create histogram and fill with hits
  for h in Li6hits:
    rate.Fill(h)
  c = ROOT.TCanvas('c', 'c') # create canvas
  rate.Draw() # draw histogram
  c.Print('plot.pdf') # print to pdf
```
