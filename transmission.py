import ROOT
import sys
import numpy
import math


def Transmission(infile, exp, runs):
  print('\nAnalyzing ' + exp)
  cyclenumbers = []
  counts = []
  countdurations = []
  monitorcounts = []
  backgroundcounts = []
  backgrounddurations = []
  beamcurrent = []
  temperature = []
  first = True

  tcycles = infile.Get('cycledata')
  # loop over cycles in cycledata
  for cycle in tcycles:
    d = cycle.durations
    li6 = cycle.countsLi6
    he3 = cycle.countsHe3
    if cycle.runnumber not in runs:
      continue # only include runs in list for experiment

    # filter useless runs
    if he3[monitorperiod] < 1000:
      print('SKIPPING cycle {0} in run{1} because there are less than 1000 monitor counts ({2})'.format(cycle.cyclenumber, cycle.runnumber, he3[monitorperiod]))
      continue
    if min(cycle.B1V_KSM_PREDCUR) < 0.1:
      print('SKIPPING cycle {0} in run {1} because beam current dropped below 0.1uA ({2}uA)'.format(cycle.cyclenumber, cycle.runnumber, min(cycle.B1V_KSM_PREDCUR)))
      continue
    if numpy.std(cycle.B1V_KSM_PREDCUR) > 0.02:
      print('SKIPPING cycle {0} in run {1} because beam current fluctuated by {2}'.format(cycle.cyclenumber, cycle.runnumber, numpy.std(cycle.B1V_KSM_PREDCUR)))
      continue
    if li6[10] == 0:
      print('SKIPPING cycle {0} in run {1} because not all periods contain Li6 data'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    if li6[countperiod] < 10*d[countperiod]:
      print('SKIPPING cycle {0} in run {1} because Li6 seems to see only background ({2}/s)'.format(cycle.cyclenumber, cycle.runnumber, li6[countperiod]/d[countperiod]))
      continue
    if li6[backgroundperiod] > 10*d[backgroundperiod]:
      print('SKIPPING cycle {0} in run {1} because Li6 sees high background ({2}/s)'.format(cycle.cyclenumber, cycle.runnumber, li6[backgroundperiod]/d[backgroundperiod]))
      continue
    
    # collect relevant data in lists
    cyclenumbers.append(float(cycle.cyclenumber))
    beamcurrent.append([numpy.mean(cycle.B1V_KSM_PREDCUR), numpy.std(cycle.B1V_KSM_PREDCUR)])
    temperature.append([min(min(cycle.UCN_ISO_TS11_RDTEMP), min(cycle.UCN_ISO_TS12_RDTEMP), min(cycle.UCN_ISO_TS14_RDTEMP)),\
                        max(max(cycle.UCN_ISO_TS11_RDTEMP), max(cycle.UCN_ISO_TS12_RDTEMP), max(cycle.UCN_ISO_TS14_RDTEMP))]) 
    counts.append(li6[countperiod])
    countdurations.append(d[countperiod])
    monitorcounts.append(he3[monitorperiod])
    backgroundcounts.append(li6[backgroundperiod])
    backgrounddurations.append(d[backgroundperiod])

  # calculate background rate averaged over all cycles in the experiment
  brate = sum(backgroundcounts)/sum(backgrounddurations)
  brateerr = math.sqrt(sum(backgroundcounts))/sum(backgrounddurations)
  print 'Li6 background rate: {0} +/- {1} 1/s'.format(brate, brateerr)

  # report average monitor counts
  monitoravg = numpy.average(monitorcounts, None, [1./m for m in monitorcounts], True)
  print 'Monitor counts: {0} +/- {1}'.format(monitoravg[0], 1./math.sqrt(monitoravg[1]))

  # report average beam current
  beamaverage = numpy.average([bc[0] for bc in beamcurrent], None, [1./bc[1]**2 if bc[1] > 0 else 0 for bc in beamcurrent], True)
  print 'Beam current: {0} +/- {1} uA'.format(beamaverage[0], 1./math.sqrt(beamaverage[1]))

  # report He-II temperature range
  print 'Temperatures from {0} to {1} K'.format(min([t[0] for t in temperature]), max([t[1] for t in temperature]))

  # subtract background from UCN counts
  bgsub = [c - brate*cd for c, cd in zip(counts, countdurations)]
  bgsuberr = [math.sqrt(c + (brateerr*cd)**2) for c,cd in zip(counts, countdurations)]

  # normalize to monitor counts
  y = [bgs/m for bgs, m in zip(bgsub, monitorcounts)]
  yerr = [math.sqrt((bgserr/m)**2 + m*(bgs/m**2)**2) for bgserr, m, bgs in zip(bgsuberr, monitorcounts, bgsub)]

  # plot ratio of background-corrected counts to monitor counts
  canvas = ROOT.TCanvas('c', 'c')
  graph = ROOT.TGraphErrors(len(cyclenumbers), numpy.array(cyclenumbers), numpy.array(y), numpy.array([0. for _ in cyclenumbers]), numpy.array(yerr))
  graph.SetTitle(exp)
  graph.GetXaxis().SetTitle('Cycle')
  graph.GetYaxis().SetTitle('UCN-count-to-monitor ratio')
  graph.SetMarkerStyle(20)
  f = graph.Fit('pol0', 'QS')
  #canvas.SetLogy()
  graph.Draw('AP')
  canvas.Print('{0}.pdf'.format(exp))
  print('{0} +/- {1}'.format(f.GetParams()[0], f.GetErrors()[0]))

  # return average ratio and background rate, including uncertainties
  return f.GetParams()[0], f.GetErrors()[0], brate, brateerr


### Main program starts here ###

# list of runs belonging to each experiment
runs = {}
runs['TCN18-029'] = [934]
runs['TCN18-031'] = [938]
runs['TCN18-035'] = [944]
runs['TCN18-043'] = [954]
runs['TCN18-045'] = [964]
runs['TCN18-080 (production up to IV1)'] = [973]
runs['TCN18-053'] = [985]
runs['TCN18-085 (MV open)'] = [990]
runs['TCN18-085'] = [993]
runs['TCN18-090'] = [1000]
runs['TCN18-290'] = [1009]
runs['TCN18-060'] = [1012,1013]

runs['TCN18-065_0A'] = [1054]
runs['TCN18-065_25A'] = [1055]
runs['TCN18-065_200A'] = [1056]
runs['TCN18-065_100A'] = [1057]
runs['TCN18-065_150A'] = [1058]
runs['TCN18-065_175A'] = [1059]
runs['TCN18-065_125A'] = [1064]
runs['TCN18-065_75A'] = [1065]
runs['TCN18-065_50A'] = [1066]

runs['TCN18-265_0A'] = [1081]
runs['TCN18-265_200A'] = [1082]
runs['TCN18-265_100A'] = [1083]
runs['TCN18-265_150A'] = [1084]
runs['TCN18-265_50A'] = [1085]
runs['TCN18-265_75A'] = [1086]
runs['TCN18-265_25A'] = [1087]

runs['TCN18-115'] = [1125]
runs['TCN18-245'] = [1129]
runs['TCN18-480'] = [1131]
runs['TCN18-480 (production up to IV1)'] = [1132, 1133]
runs['TCN18-057'] = [1141]
runs['TCN18-302'] = [1165]
runs['TCN18-240'] = [1176]
runs['TCN18-215'] = [1181]
runs['TCN18-380'] = [1188]
runs['TCN18-310'] = [1192]

countperiod = 1
monitorperiod = 0
backgroundperiod = 0
ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
#ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1

infile = ROOT.TFile(sys.argv[1])

h = ROOT.TH1D('bg', 'bg', 100, 0, 30)
# loop over experiments and analyze transmission in runs
for tcn in runs:
  result = Transmission(infile, tcn, runs[tcn])
  h.Fill(result[2]) # fill histogram with the determined background rate
c = ROOT.TCanvas('c', 'c')
h.Draw()
c.Print('bg.pdf')
