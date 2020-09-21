import sys
import ROOT

# Collect all valve-timing histograms from root files passed as command-line parameters
# Adds up all timing histograms and prints them to pdfs

hists = {}
for filename in sys.argv[1:]:
  print(filename)
  f = ROOT.TFile(filename, 'READ')
  for valve in ['IV2', 'IV3']:
    for diff in ['DriveClosedDiff', 'DriveOpenDiff', 'NoDriveClosedDiff', 'NoDriveOpenDiff']:
      name = valve + diff
      hist = f.Get(name)
      if name not in hists:
        hists[name] = hist
        hists[name].SetDirectory(0)
      else:
        hists[name].Add(hist)

c = ROOT.TCanvas('c', 'c')
for name in hists:
  hists[name].Draw()
  c.Print(name + '.pdf')
