import ROOT
import sys
import numpy
import math
import itertools

# analyze storage time from list of runs
def StorageLifetime(infile, exp, runs):
  print('\nAnalyzing ' + exp)
  counts = []
  countdurations = []
  monitorcounts = []
  backgroundcounts = []
  backgrounddurations = []
  storagedurations = []
  beamcurrent = []
  temperature = []
  monitortau = []
  monitortauerr = []
  first = True

  # loop over cycles
  tcycles = infile.Get('cycledata')
  for cycle in tcycles:
    d = cycle.durations
    li6 = cycle.countsLi6
    he3 = cycle.countsHe3
    if cycle.runnumber not in runs:
      continue # only use runs in list

    # filter useless cycles
    if he3[monitorperiod] < 1000:
      print('SKIPPING cycle {0} in run {1} because He3 saw less than 1000 monitor counts ({2})'.format(cycle.cyclenumber, cycle.runnumber, he3[monitorperiod]))
      continue
    if min(cycle.B1V_KSM_PREDCUR) < 0.1:
      print('SKIPPING cycle {0} in run {1} because beam dropped below 0.1uA ({2}uA)'.format(cycle.cyclenumber, cycle.runnumber, min(cycle.B1V_KSM_PREDCUR)))
      continue
    if numpy.std(cycle.B1V_KSM_PREDCUR) > 0.02:
      print('SKIPPING cycle {0} in run {1} because beam fluctuated by {2}uA'.format(cycle.cyclenumber, cycle.runnumber, numpy.std(cycle.B1V_KSM_PREDCUR)))
      continue
    if li6[10] == 0:
      print('SKIPPING cycle {0} in run {1} because Li6 does not contain data in all periods'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    if d[backgroundperiod] > 0 and li6[backgroundperiod]/d[backgroundperiod] > 10:
      print('SKIPPING cycle {0} in run {1} because of high Li6 background ({2}/s)'.format(cycle.cyclenumber, cycle.runnumber, li6[backgroundperiod]/d[backgroundperiod]))
      continue

    # store data in lists
    beamcurrent.append([numpy.mean(cycle.B1V_KSM_PREDCUR), numpy.std(cycle.B1V_KSM_PREDCUR)])
    temperature.append([min(min(cycle.UCN_ISO_TS11_RDTEMP), min(cycle.UCN_ISO_TS12_RDTEMP), min(cycle.UCN_ISO_TS14_RDTEMP)),\
                        max(max(cycle.UCN_ISO_TS11_RDTEMP), max(cycle.UCN_ISO_TS12_RDTEMP), max(cycle.UCN_ISO_TS14_RDTEMP))]) 
    counts.append(li6[countperiod])
    countdurations.append(d[countperiod])
    monitorcounts.append(he3[monitorperiod])
    backgroundcounts.append(li6[backgroundperiod])
    backgrounddurations.append(d[backgroundperiod])
    storagedurations.append(d[storageperiod])

    # fit single exponential to rate in He3 detector
    he3rate = ROOT.TH1I('He3_{0}_{1}'.format(cycle.runnumber, cycle.cyclenumber), 'He3 detector rate', int(math.ceil(sum(d))), 0., sum(d))
    for h in getattr(cycle, 'He3/hits'):
      he3rate.Fill(h)
    fitstart = sum(itertools.islice(d, storageperiod)) + 5
    fitend = sum(itertools.islice(d, storageperiod + 1))
    if fitend > fitstart + 10 and he3rate.Integral(he3rate.FindBin(fitstart), he3rate.FindBin(fitend)) > 0:
      rf = ROOT.TF1('rf', '[0]*exp(-x/[1])')
      rf.SetParameters(10, 10)
      he3rate.Fit(rf, 'Q', '', fitstart, fitend)
      monitortau.append(rf.GetParameter(1))
      monitortauerr.append(rf.GetParError(1))
      che3 = ROOT.TCanvas('{0}_{1}'.format(cycle.runnumber, cycle.cyclenumber),'che3')
      he3rate.Draw()
      # print fitted He3 rate to pdf
      che3.Print('{0}.pdf{1}'.format(exp, '(' if first else ''))
      first = False

  # print average storage lifetime from He3 fits to pdf
  if len(monitortau) > 0:
    che3 = ROOT.TCanvas('che3','che3')
    tauavg = numpy.average(monitortau, None, [1./dt**2 for dt in monitortauerr], True)
    label = ROOT.TLatex(0.1, 0.5, '#tau = {0} +/- {1}'.format(tauavg[0], 1./tauavg[1]**2))
    label.Draw()
    che3.Print('{0}.pdf'.format(exp))
    print('{0} +/- {1} (single exponential fit to rate in monitor detector during storage period)'.format(tauavg[0], 1./tauavg[1]**2))
  
  # calculate background rate averaged over all cycles of experiment
  brate = sum(backgroundcounts)/sum(backgrounddurations)
  brateerr = numpy.std(backgroundcounts)/sum(backgrounddurations)
  print 'Li6 background rate: {0} +/- {1} 1/s'.format(brate, brateerr)

  # report average monitor counts, average beam current, range of He-II temperature
  monitoravg = numpy.average(monitorcounts, None, [1./m for m in monitorcounts], True)
  print 'Monitor counts: {0} +/- {1}'.format(monitoravg[0], 1./math.sqrt(monitoravg[1]))
  beamaverage = numpy.average([bc[0] for bc in beamcurrent], None, [1./bc[1]**2 if bc[1] > 0 else 0 for bc in beamcurrent], True)
  print 'Beam current: {0} +/- {1} uA'.format(beamaverage[0], 1./math.sqrt(beamaverage[1]))
  print 'Temperatures from {0} to {1} K'.format(min([t[0] for t in temperature]), max([t[1] for t in temperature]))

  x = storagedurations
  xerr = [0. for _ in storagedurations]
  # subtract background from UCN counts
  bgsub = [c - brate*cd for c, cd in zip(counts, countdurations)]
  bgsuberr = [math.sqrt(c + brateerr**2) for c in counts]
  # normalize to monitor counts
  y = [bgs/m for bgs, m in zip(bgsub, monitorcounts)]
  yerr = [math.sqrt((bgserr/m)**2 + m*(bgs/m**2)**2) for bgserr, m, bgs in zip(bgsuberr, monitorcounts, bgsub)]
  # plot normalized, background corrected counts vs storage time
  canvas = ROOT.TCanvas('c', 'c')
  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph.SetTitle(exp + ' (single exponential fit, with background subtracted, normalized to monitor detector)')
  graph.GetXaxis().SetTitle('Storage time (s)')
  graph.GetYaxis().SetTitle('UCN-count-to-monitor ratio')
#  graph.SetMarkerStyle(20)
  # do single exponential fit
  f = ROOT.TF1('f','[0]*exp(-x/[1])')
  f.SetParameters(1, 3)
  f.SetParName(1, '#tau')
  graph.Fit(f, 'Q')
  #canvas.SetLogy()
  graph.Draw('AP')
  canvas.Print('{0}.pdf'.format(exp))
  print('{0} +/- {1} (single exponential fit, with background subtracted, normalized to monitor detector)'.format(f.GetParameter(1), f.GetParError(1)))

  # do double exponential fit
  f2 = ROOT.TF1('f2', '[0]*exp(-x/[1]) + [2]*exp(-x/[3])')
  f2.SetParameters(1, 10, 0.1, 50)
  graph.SetTitle(exp + ' (double exponential fit, with background subtracted, normalized to monitor detector)')
  graph.Fit(f2, 'Q')
  graph.Draw('AP')
  canvas.Print('{0}.pdf'.format(exp))
  print('{0} +/- {1}, {2} +/- {3} (double exponential fit, with background subtracted, normalized to monitor detector)'.format(f2.GetParameter(1), f2.GetParError(1), f2.GetParameter(3), f2.GetParError(3)))
  
  # plot uncorrected UCN counts
  y = [float(c) for c in counts]
  yerr = [math.sqrt(c) for c in counts]
  canvas1 = ROOT.TCanvas(exp, exp)
  graph1 = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph1.SetTitle(exp + ' (single exponential fit + background, unnormalized)')
  graph1.GetXaxis().SetTitle('Storage time (s)')
  graph1.GetYaxis().SetTitle('UCN count')
  # do single exponential fit with background
  f1 = ROOT.TF1('f1', '[0]*exp(-x/[1]) + [2]')
  f1.SetParameters(10000, 10, 30)
  f1.SetParName(1, '#tau')
  f1.SetParName(2, 'Background')
  graph1.Fit(f1, 'Q')
  graph1.Draw('AP')
  canvas1.Print('{0}.pdf)'.format(exp))
  print('{0} +/- {1} (single exponential fit with {2} +/- {3} background, unnormalized)'.format(f1.GetParameter(1), f1.GetParError(1), f1.GetParameter(2), f1.GetParError(2)))

  # return storage lifetime from single exponential fit, background rate, background from exponential fit with background
  return f.GetParameter(1), f.GetParError(1), brate, f1.GetParameter(2)/numpy.mean(countdurations)


# list runs for each experiment
runs = {}
runs['TCN18-025'] = [932] # TCN18-025
runs['TCN18-026'] = [933] # TCN18-026
runs['TCN18-032'] = [939] # TCN18-032
runs['TCN18-033'] = [940] # TCN18-033
runs['TCN18-036'] = [949, 955, 956, 957, 958] #[947, 948, 949, 955, 956, 957, 958] # TCN18-036
runs['TCN18-040'] = [950] # TCN18-040
runs['TCN18-041'] = [951] # TCN18-041
runs['TCN18-042'] = [952] # TCN18-042
runs['TCN18-046'] = [961, 962, 967, 969] # TCN18-046
runs['TCN18-050'] = [966] # TCN18-050
runs['TCN18-051'] = [965] # TCN18-051
runs['TCN18-052'] = [970] # TCN18-052
runs['TCN18-081'] = [976, 977] # TCN18-081
runs['TCN18-082'] = [974] # TCN18-082
runs['TCN18-055'] = [983]
runs['TCN18-054'] = [986]
runs['TCN18-087 (MV open)'] = [989]
runs['TCN18-086 (MV open)'] = [991]
runs['TCN18-087'] = [992]
runs['TCN18-092'] = [999]
runs['TCN18-091'] = [1001, 1002, 1003, 1004]
runs['TCN18-291'] = [1008]
runs['TCN18-292'] = [1010]
runs['TCN18-061'] = [1014, 1015, 1016, 1018]
runs['TCN18-062'] = [1017]
runs['TCN18-066_0A'] = [1067]
runs['TCN18-066_200A'] = [1068]
runs['TCN18-066_150A'] = [1069]
runs['TCN18-066_100A'] = [1070]
runs['TCN18-066_50A'] = [1071]
runs['TCN18-068_0A'] = [1072, 1073, 1074]
runs['TCN18-068_200A'] = [1075]
runs['TCN18-068_150A'] = [1076]
runs['TCN18-068_50A'] = [1077]
runs['TCN18-068_100A'] = [1078]
runs['TCN18-268_0A'] = [1089]
runs['TCN18-268_200A'] = [1090]
runs['TCN18-268_100A'] = [1091]
runs['TCN18-268_150A'] = [1092]
runs['TCN18-268_50A'] = [1093]
runs['TCN18-266_0A'] = [1094]
runs['TCN18-266_200A'] = [1095]
runs['TCN18-266_100A'] = [1096]
runs['TCN18-266_50A'] = [1097]
runs['TCN18-266_150A'] = [1098]
runs['TCN18-125'] = [1118, 1119]
runs['TCN18-126'] = [1122]
runs['TCN18-127'] = [1124]
runs['TCN18-116'] = [1126]
runs['TCN18-117'] = [1127]
runs['TCN18-481'] = [1136]
runs['TCN18-482'] = [1134,1135]
runs['TCN18-058'] = [1142]
runs['TCN18-059'] = [1143]
runs['TCN18-216'] = [1182, 1183, 1184]
runs['TCN18-217'] = [1185, 1186, 1187]
runs['TCN18-381'] = [1189, 1190]
runs['TCN18-382'] = [1191]

countperiod = 2
monitorperiod = 0
backgroundperiod = 1
storageperiod = 1
ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1

infile = ROOT.TFile(sys.argv[1])

result = {}
bg = ROOT.TH1D('bg','background', 100, 0, 10)
bg2 = ROOT.TH1D('bg2','background', 100, 0, 20)

# loop over experiments
for tcn in runs:
  result[tcn] = StorageLifetime(infile, tcn, runs[tcn])
  # fill histograms with determined background rates
  bg.Fill(result[tcn][2])
  bg2.Fill(result[tcn][3])

# draw background histograms
c = ROOT.TCanvas('bg','bg')
bg.Draw()
c.Print('background.pdf(')
c2 = ROOT.TCanvas('bg','bg')
bg2.Draw()
c2.Print('background.pdf)')

