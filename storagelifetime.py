import ROOT
import sys
import numpy
import math
import datetime
import itertools
import scipy.optimize

# calculate 4He vapor pressure from temperature
def HeVaporPressure(T):
  # from Clement, Logan, Gaffney, Phys. Rev. 100, 743
  # https://doi.org/10.1103/PhysRev.100.743
  I = 4.6202
  A = 6.399
  B = 2.541
  C = 0.00612
  D = 0.5197
  a = 7.
  b = 14.14
  lnP = I - A/T + B*math.log(T) + C/2*T**2 - D*(a*b/(b**2 + 1) - 1./T)*math.atan(a*T - b) - a*D/2/(b**2 + 1)*math.log(T**2/(1 + (a*T - b)**2))
  return math.exp(lnP)

# use scipy solver to invert vapor pressure formula to calculate temperature from vapor pressure
def HeTemperature(P):
  return scipy.optimize.brentq(lambda T: HeVaporPressure(T) - P, 0.5, 4)

# analyze storage time from list of runs
def StorageLifetime(infile, exp, runs, excycles):
  print('\nAnalyzing {0} in run {2}, excluding cycles {1}'.format(exp, excycles, run))
  starts = []
  li6counts = []
  he3counts = []
  countdurations = []
  backgroundcounts = []
  backgrounddurations = []
  storagedurations = []
  beamcurrent = []
  temperature = []
  vaporpress = []
  beamaverage = []

  # loop over cycledata
  tcycles = infile.Get('cycledata')
  for cycle in tcycles:
    d = cycle.durations
    det = cycle.countsLi6
    
    # use only cycles that are part of listed runs and not in exclusion list
    if cycle.runnumber in runs and (cycle.runnumber not in excycles or cycle.cyclenumber not in excycles[cycle.runnumber]):

      # filter useless runs
      if cycle.countsLi6[countperiod] == 0 and cycle.countsHe3[countperiod] == 0:
        print('SKIPPING cycle {0} in run {1} because it contains no detector counts'.format(cycle.cyclenumber, cycle.runnumber))
        continue
      if cycle.countsLi6[countperiod] > 0 and cycle.countsLi6[10] == 0: # indicates that run was stopped during counting
        print('SKIPPING cycle {0} in run {1} because not all periods contain Li6detector counts'.format(cycle.cyclenumber, cycle.runnumber))
        continue
      beam = [cur for cur in cycle.B1V_KSM_PREDCUR]
      if min(beam) < 0.1:
        print('SKIPPING cycle {0} in run {1} because beam current dropped to {2}uA!'.format(cycle.cyclenumber, cycle.runnumber, min(beam)))
        continue
      if numpy.std(beam) > 0.02:
        print('SKIPPING cycle {0} in run {1} because beam current fluctuated by {2}uA!'.format(cycle.cyclenumber, cycle.runnumber, numpy.std(beam)))
        continue
      if max(cycle.UCN_UGD_IV1_STATON) < 1:
        print('SKIPPING cycle {0} in run {1} because IV1 never opened!'.format(cycle.cyclenumber, cycle.runnumber))
        continue

      # collect data in lists
      starts.append(cycle.start)
      beamcurrent = beamcurrent + beam
      beamaverage.append([numpy.average(beam), numpy.std(beam)])
      temperature.append([min(min(cycle.UCN_ISO_TS11_RDTEMP), min(cycle.UCN_ISO_TS12_RDTEMP), min(cycle.UCN_ISO_TS14_RDTEMP)),\
                          max(max(cycle.UCN_ISO_TS11_RDTEMP), max(cycle.UCN_ISO_TS12_RDTEMP), max(cycle.UCN_ISO_TS14_RDTEMP))])
      if max(cycle.UCN_ISO_PG9L_RDPRESS) >= 2.:   
        vaporpress = vaporpress + [p for p in cycle.UCN_ISO_PG9H_RDPRESS]
      else:
        vaporpress = vaporpress + [p for p in cycle.UCN_ISO_PG9L_RDPRESS]
      li6counts.append(cycle.countsLi6[countperiod])
      he3counts.append(cycle.countsHe3[countperiod])
      countdurations.append(d[countperiod])
      backgroundcounts.append(cycle.countsLi6[backgroundperiod])
      backgrounddurations.append(d[backgroundperiod])
      storagedurations.append(d[storageperiod])

  # report start time
  print(datetime.datetime.fromtimestamp(starts[0]))

  # calculate background rate averaged over all cycles
  brate = sum(backgroundcounts)/sum(backgrounddurations)
#  print(backgroundcounts)
  brateerr = math.sqrt(sum(backgroundcounts))/sum(backgrounddurations)
  print('Detector background rate: {0} +/- {1} 1/s'.format(brate, brateerr))
 
  # report range of beam current
  print('Beam current from {0} to {1} uA'.format(min(beamcurrent), max(beamcurrent)))

  # report range of temperature and vapor pressure
  mintemp = min([t[0] for t in temperature])
  maxtemp = max([t[1] for t in temperature])
  print('Temperatures from {0} to {1} K'.format(mintemp, maxtemp))
  print('Vapor pressure from {0} to {1} torr'.format(min(vaporpress), max(vaporpress)))

  x = storagedurations
  xerr = [0. for _ in storagedurations]
  # subtract background from UCN counts
  bgsub = [c - brate*cd for c, cd in zip(li6counts, countdurations)]
  bgsuberr = [math.sqrt(c + brateerr**2*cd) for c,cd in zip(li6counts, countdurations)]

  # normalize to beam current
  y = [bgs/cur[0] for bgs, cur in zip(bgsub, beamaverage)]
  yerr = [math.sqrt((bgserr/cur[0])**2 + (cur[1]*bgs/cur[0]**2)**2) for bgserr, bgs, cur in zip(bgsuberr, bgsub, beamaverage)]

  # plot normalized Li6 counts vs storage time
  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph.SetTitle(exp + ' (Li6 detector, single exponential fit, with background subtracted, normalized to beam current)')
  graph.GetXaxis().SetTitle('Storage time (s)')
  graph.GetYaxis().SetTitle('UCN count (#muA^{-1})')
#  graph.SetMarkerStyle(20)
  # do single exponential fit
  f = ROOT.TF1('f','[0]*exp(-x/[1])')
  f.SetParameters(1000, 15)
  f.SetParName(1, '#tau')
  graph.Fit(f, 'Q')
  li6tau = [f.GetParameter(1), f.GetParError(1)]
  canvas = ROOT.TCanvas('li6', 'li6')
  #canvas.SetLogy()
  graph.Draw('AP')
  canvas.Print('{0}.pdf('.format(exp))
  print('{0} +/- {1} (Li6 detector, single exponential fit, with background subtracted, normalized to beam current)'.format(li6tau[0], li6tau[1]))
  
  # plot uncorrected UCN counts
  y = [float(c) for c in li6counts]
  yerr = [math.sqrt(c) for c in li6counts]
  graph1 = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph1.SetTitle(exp + ' (Li6 detector, single exponential fit + background, unnormalized)')
  graph1.GetXaxis().SetTitle('Storage time (s)')
  graph1.GetYaxis().SetTitle('UCN count')
  # do single exponential fit with background
  f1 = ROOT.TF1('f1', '[0]*exp(-x/[1]) + [2]')
  f1.SetParameters(1000, 15, 200)
  f1.SetParName(1, '#tau')
  f1.SetParName(2, 'Background')
  graph1.Fit(f1, 'Q')
  canvas1 = ROOT.TCanvas('li6_2', 'li6_2')
  graph1.Draw('AP')
  canvas1.Print('{0}.pdf'.format(exp))
  print('{0} +/- {1} (single exponential fit with {2} +/- {3} background, unnormalized)'.format(f1.GetParameter(1), f1.GetParError(1), f1.GetParameter(2), f1.GetParError(2)))

  # plot beam-normalized He3 counts vs storage time
  y = [c/cur[0] for c, cur in zip(he3counts, beamaverage)]
  yerr = [math.sqrt(c/cur[0]**2 + (cur[1]*c/cur[0]**2)**2) for c, cur in zip(he3counts, beamaverage)]
  graph2 = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph2.SetTitle(exp + ' (He3 detector, single exponential fit, normalized to beam current)')
  graph2.GetXaxis().SetTitle('Storage time (s)')
  graph2.GetYaxis().SetTitle('UCN count (#muA^{-1})')
  # do single exponential fit
  f.SetParameters(1000, 15)
  graph2.Fit(f, 'Q')
  he3tau = [f.GetParameter(1), f.GetParError(1)]
  canvas2 = ROOT.TCanvas('he3', 'he3')
  graph2.Draw('AP')
  canvas2.Print('{0}.pdf)'.format(exp))
  print('{0} +/- {1} (He3 detector, single exponential fit, normalized to beam current)'.format(he3tau[0], he3tau[1]))

  # return result from primary detector
  if max(li6counts) > max(he3counts):
    return starts[0], li6tau[0], li6tau[1], min(vaporpress), max(vaporpress), mintemp, maxtemp, brate, brateerr
  else:
    return starts[0], he3tau[0], he3tau[1], min(vaporpress), max(vaporpress), mintemp, maxtemp, brate, brateerr
    



runs = {}
excycles = {}

# list runs for different source-storage experiments
runs['TCN18-015'] = [[869], [870], [872], [895], [900], [907], [911], [923], [931], [941], [953], [968], [975], [984], [988], [998], [1011], [1019], [1030], [1049], [1053], [1088], [1123], [1137], [1151], [1167], [1179], [1193]] #Daily storage lifetimes
excycles['TCN18-015']= {968: [0], # not sure what happened here, only background in Li6. IV2/3 did not open?
                        975: [0,1,2,3,4]} # IG5/6 were turned on #872: [4,17], 900: [10], 953: [24], 975: [0,1,2,3,23,29], 1123: [8], 1179: [8]}

runs['TCN18-300'] = [[1153], [1154], [1155], [1156], [1158], [1159], [1160], [1161], [1167]]
excycles['TCN18-300'] = {1161: [0]} # huge background. IG5 turned on? #1153: [0,9], 1155: [0], 1156: [2,5,6,14], 1158: [7], 1159: [5], 1160: [0], 1161: [0,1,2,3,15]}

runs['TCN18-170'] = [[1194], [1195, 1196, 1197], [1198], [1199], [1200], [1201], [1202], [1203], [1204], [1205]]
excycles['TCN18-170'] = {}#1195: [5,6], 1196: [3], 1198: [8], 1202: [8], 1203: [8], 1204: [8]}
injectedpress =      [1.04,  1.995,               4.236,  7.953, 16.284, 32.425,  62.84,    315,    628,   1365]

countperiod = 2
monitorperiod = 0
backgroundperiod = 1
storageperiod = 1
ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1

infile = ROOT.TFile(sys.argv[1])
bg = ROOT.TH1D('bg','bg',100,0,5)

taus = []
# loop over daily storage lifetime measurements
for run in runs['TCN18-015']:
  taus.append(StorageLifetime(infile, 'TCN18-015_{0}'.format(run[0]), run, excycles['TCN18-015']))

# plot storage lifetime vs time
ctaus = ROOT.TCanvas('taus', 'taus')
grtaus = ROOT.TGraphErrors(len(taus), numpy.array([tau[0] for tau in taus]), numpy.array([tau[1] for tau in taus]), numpy.array([0. for _ in taus]), numpy.array([tau[2] for tau in taus]))
grtaus.GetXaxis().SetTimeDisplay(1)
grtaus.GetXaxis().SetNdivisions(10, 10, 0)
grtaus.GetXaxis().SetTitle('Date')
grtaus.GetYaxis().SetTitle('Storage lifetime (s)')
grtaus.GetYaxis().SetRangeUser(0,40)
grtaus.Draw('AP')
ctaus.Print('dailytau.pdf')
# add background rates to histogram
for tau in taus:
  bg.Fill(tau[7])

# loop over storage lifetime measurements
taus = []
for run in runs['TCN18-300']:
  taus.append(StorageLifetime(infile, 'TCN18-300_{0}'.format(run[0]), run, excycles['TCN18-300']))

# plot storage lifetime vs temperature
ofile = ROOT.TFile('tauvsvaporpress.root', 'RECREATE')
ctemp = ROOT.TCanvas('temp', 'temp')
grtemp = ROOT.TGraphErrors(len(taus), numpy.array([(tau[5] + tau[6])/2 for tau in taus]), numpy.array([tau[1] for tau in taus]), numpy.array([(tau[6] - tau[5])/2 for tau in taus]), numpy.array([tau[2] for tau in taus]))
grtemp.GetXaxis().SetTitle('Temperature (K)')
grtemp.GetYaxis().SetTitle('Storage lifetime (s)')
grtemp.SetLineColor(ROOT.kRed)
grtemp.Draw('AP')
ctemp.Print('tauvstemp.pdf')
ofile.cd()
grtemp.Write()

# plot storage lifetime vs vapor pressure
cpress = ROOT.TCanvas('press', 'press')
#grpress = ROOT.TGraphErrors(len(taus), numpy.array([(HeTemperature(tau[3]) + HeTemperature(tau[4]))/2 for tau in taus]), numpy.array([tau[1] for tau in taus]), numpy.array([(HeTemperature(tau[4]) - HeTemperature(tau[3]))/2 for tau in taus]), numpy.array([tau[2] for tau in taus]))
grpress = ROOT.TGraphErrors(len(taus), numpy.array([(HeTemperature(tau[3]) + HeTemperature(tau[4]))/2 for tau in taus]),\
                                       numpy.array([tau[1] for tau in taus]),\
                                       numpy.array([(HeTemperature(tau[4]) - HeTemperature(tau[3]))/2 for tau in taus]),\
                                       numpy.array([tau[2] for tau in taus]))
grpress.GetXaxis().SetTitle('Temperature derived from vapor pressure (K)')
grpress.GetYaxis().SetTitle('Storage lifetime (s)')
grpress.SetLineColor(ROOT.kBlue)
grpress.Draw('AP')
#cpress.SetLogx()
cpress.Print('tauvsvaporpress.pdf')
ofile.cd()
grpress.Write()
ofile.Close()
# add background rates to histogram
for tau in taus:
  bg.Fill(tau[7])

# loop over storage lifetime measurements
taus = []
for run, P in zip(runs['TCN18-170'], injectedpress):
  taus.append(StorageLifetime(infile, 'TCN18-170_{0}'.format(run[0]), run, excycles['TCN18-170']))

# plot storage lifetime vs source spoilage
cspoil = ROOT.TCanvas('spoil', 'spoil')
grspoil = ROOT.TGraphErrors(len(taus), numpy.cumsum([p*0.1 for p in injectedpress]), numpy.array([tau[1] for tau in taus]), numpy.array([0. for _ in taus]), numpy.array([tau[2] for tau in taus]))
grspoil.GetXaxis().SetTitle('Spoiling gas injected (torr L)')
grspoil.GetYaxis().SetTitle('Storage Lifetime (s)')
grspoil.Draw('AP')
cspoil.SetLogx()
cspoil.Print('tauvsspoilage.pdf')
# add background rates to histogram
for tau in taus:
  bg.Fill(tau[7])

# plot histogram of background counts
cbg = ROOT.TCanvas('bg','bg')
bg.Draw()
cbg.Print('bg.pdf')
