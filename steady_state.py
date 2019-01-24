import ROOT
import sys
import numpy
import inspect

f = ROOT.TFile(sys.argv[1])

f.cycledata.SetBranchStatus('Li6',1)
print(f.cycledata)

for cycle in f.cycledata:
  print(cycle.runnumber)
  hits = [h for h in getattr(cycle, 'Li6/hits')]
  Ttime = [t for t in getattr(cycle, 'Source/timestamp')]
  hist = numpy.histogram(hits, Ttime)
  rate = hist[0]/numpy.diff(Ttime)
