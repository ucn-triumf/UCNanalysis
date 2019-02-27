import ROOT
import sys
import numpy
import math
import itertools
import UCN

def ReadCycles(infile, experiments):
  countperiod = 2
  monitorperiod = 0
  backgroundperiod = 1
  storageperiod = 1

  for ex in experiments:
    ex['start'] = []
    ex['beamcurrent'] = []
    ex['li6counts'] = []
    ex['countduration'] = []
    ex['monitorcounts'] = []
    ex['monitorduration'] = []
    ex['li6background'] = []
    ex['backgroundduration'] = []
    ex['li6irradiation'] = []
    ex['irradiationduration'] = []
    ex['storageduration'] = []
    ex['beamcurrent'] = []
    ex['mintemperature'] = []
    ex['maxtemperature'] = []
    ex['minvaporpressure'] = []
    ex['maxvaporpressure'] = []
    ex['SCMcurrent'] = []
    ex['he3rate'] = []

  for cycle in infile.cycledata:
    run = cycle.runnumber
  
    if not any(run in ex['runs'] for ex in experiments): # if there is no experiment using this cycle
      continue
    
    Li6 = cycle.countsLi6
    He3 = cycle.countsHe3
    d = cycle.durations
   
    # filter useless cycles
    if min(cycle.B1V_KSM_PREDCUR) < 0.1:
      print('SKIPPING cycle {0} in run {1} because beam dropped below 0.1uA ({2}uA)'.format(cycle.cyclenumber, cycle.runnumber, min(cycle.B1V_KSM_PREDCUR)))
      continue
    if numpy.std(cycle.B1V_KSM_PREDCUR) > 0.02:
      print('SKIPPING cycle {0} in run {1} because beam fluctuated by {2}uA'.format(cycle.cyclenumber, cycle.runnumber, numpy.std(cycle.B1V_KSM_PREDCUR)))
      continue
    if Li6[10] == 0:
      print('SKIPPING cycle {0} in run {1} because Li6 does not contain data in all periods'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    if He3[monitorperiod] < 1000:
      print('SKIPPING cycle {0} in run {1} because He3 saw less than 1000 monitor counts ({2})'.format(cycle.cyclenumber, cycle.runnumber, He3[monitorperiod]))
      continue
    if d[backgroundperiod] > 0 and Li6[backgroundperiod]/d[backgroundperiod] > 10:
      print('SKIPPING cycle {0} in run {1} because of high Li6 background ({2}/s)'.format(cycle.cyclenumber, cycle.runnumber, Li6[backgroundperiod]/d[backgroundperiod]))
      continue
    if any([1e-7 < ig5 < 1e-2 for ig5 in cycle.UCN_EXP_IG5_RDVAC]) or any([1e-7 < ig6 < 1e-2 for ig6 in cycle.UCN_EXP_IG6_RDVAC]):
      print ('SKIPPING cycle {0} in run {1} because IG5 or IG6 were on!'.format(cycle.cyclenumber, cycle.runnumber))
      continue

    if (cycle.valve0state[0] != 1 or cycle.valve0state[1] != 0 or cycle.valve0state[2] != 0 or # IV1 should be open during irradiation, closed after
       cycle.valve1state[0] != 1 or # IV2 should be open during irradiation, state after depends on storage mode (between IV1+IV3 or between IV2+IV3)
       cycle.valve2state[0] != 0 or cycle.valve2state[1] != 0 or cycle.valve2state[2] != 1): # IV3 should be closed during irradiation and storage, open during counting
      print('Abnormal valve configuration in cycle {0} of run {1}'.format(cycle.cyclenumber, cycle.runnumber))
    
    for ex in experiments:
      if run not in ex['runs']:
        continue

      ex['start'].append(cycle.start)
      ex['beamcurrent'].append([cur for cur in cycle.B1V_KSM_PREDCUR])
      ex['mintemperature'].append(min([min(getattr(cycle,'UCN_ISO_{0}_RDTEMP'.format(TS))) for TS in ['TS11', 'TS12', 'TS14']]))
      ex['maxtemperature'].append(max([max(getattr(cycle,'UCN_ISO_{0}_RDTEMP'.format(TS))) for TS in ['TS11', 'TS12', 'TS14']]))
      if max(cycle.UCN_ISO_PG9L_RDPRESS) >= 2.:
        ex['minvaporpressure'].append(min(cycle.UCN_ISO_PG9H_RDPRESS))
        ex['maxvaporpressure'].append(max(cycle.UCN_ISO_PG9H_RDPRESS))
      else:
        ex['minvaporpressure'].append(min(cycle.UCN_ISO_PG9L_RDPRESS))
        ex['maxvaporpressure'].append(max(cycle.UCN_ISO_PG9L_RDPRESS))
      ex['SCMcurrent'].append([v/250e-6 for v in cycle.SCMVoltages3])
      ex['li6counts'].append(Li6[countperiod])
      ex['countduration'].append(d[countperiod])
      ex['monitorcounts'].append(He3[monitorperiod])
      ex['monitorduration'].append(d[monitorperiod])
      ex['li6background'].append(Li6[backgroundperiod])
      ex['backgroundduration'].append(d[backgroundperiod])
      ex['storageduration'].append(d[storageperiod])
      ex['irradiationduration'].append(d[0])
      ex['li6irradiation'].append(Li6[0])

      he3rate = ROOT.TH1I('He3_{0}_{1}'.format(cycle.runnumber, cycle.cyclenumber), 'He3 detector rate', int(math.floor(sum(d))), 0., math.floor(sum(d)))
      for h in getattr(cycle, 'He3/hits'):
        he3rate.Fill(h)
      he3rate.GetXaxis().SetTitle('Time (s)')
      he3rate.GetYaxis().SetTitle('He3 rate (s^{-1})')
      he3rate.SetDirectory(0)
      ex['he3rate'].append(he3rate)
	  
	  
# analyze storage time from list of runs
def StorageLifetime(ex):
  print('\nAnalyzing TCN' + ex['TCN'])
  
  if len(ex['start']) == 0:
    print('Found no cycles with run numbers {0}!'.format(ex['runs']))
    return

  ex['li6backgroundrate'], ex['li6backgroundrateerr'] = UCN.BackgroundRate(ex['li6background'], ex['backgroundduration'])
  print('Li6 detector background rate: {0} +/- {1} 1/s'.format(ex['li6backgroundrate'], ex['li6backgroundrateerr']))
  beam = [numpy.mean(cur) for cur in ex['beamcurrent']], [numpy.std(cur) for cur in ex['beamcurrent']]
  ex['li6irradiationrate'], ex['li6irradiationrateerr'] = UCN.SubtractBackgroundAndNormalizeRate(ex['li6irradiation'], ex['irradiationduration'], 'li6', beam[0], beam[1])

  # report average monitor counts, range of beam current, range of He-II temperature
  monitoravg = numpy.average(ex['monitorcounts'], None, [1./m for m in ex['monitorcounts']], True)
  print 'Monitor counts: {0} +/- {1}'.format(monitoravg[0], 1./math.sqrt(monitoravg[1]))
  print('Beam current from {0} to {1} uA'.format(min(min(c) for c in ex['beamcurrent']), max(max(c) for c in ex['beamcurrent'])))
  print 'Temperatures from {0} to {1} K'.format(min(ex['mintemperature']), max(ex['maxtemperature']))

  x = ex['storageduration']
  xerr = [0. for _ in ex['storageduration']]  

  y, yerr = UCN.SubtractBackgroundAndNormalize(ex['li6counts'], ex['countduration'], 'li6', ex['monitorcounts'], [math.sqrt(m) for m in ex['monitorcounts']])
  
  # plot normalized, background corrected counts vs storage time
  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph.SetTitle('TCN{0} (single exponential fit, with background subtracted, normalized to monitor detector)'.format(ex['TCN']))
  graph.GetXaxis().SetTitle('Storage time (s)')
  graph.GetYaxis().SetTitle('UCN-count-to-monitor ratio')
#  graph.SetMarkerStyle(20)

  # do single exponential fit
  f = graph.Fit(UCN.SingleExpo(), 'SQB', '', 0., 1000.)

  canvas = ROOT.TCanvas('c', 'c')
  canvas.SetLogy()
  graph.Draw('AP')
  pdf = 'TCN{0}.pdf'.format(ex['TCN'])
  canvas.Print(pdf + '(')
  ex['tau'] = f.GetParams()[1]
  ex['tauerr'] = f.GetErrors()[1]
  print('{0} +/- {1} (single exponential fit, with background subtracted, normalized to monitor detector)'.format(ex['tau'], ex['tauerr']))
  
  #do single exponential fit with data point at 0 excluded
  f = graph.Fit(UCN.SingleExpo(), 'SQB', '', 1., 1000.)
  graph.SetTitle('TCN{0} (single exponential fit, with background subtracted, normalized to monitor detector, 0s excluded)'.format(ex['TCN']))
  graph.Draw('AP')
  canvas.Print(pdf)

  # do double exponential fit
  graph.SetTitle('TCN{0} (double exponential fit, with background subtracted, normalized to monitor detector)'.format(ex['TCN']))
  f = graph.Fit(UCN.DoubleExpo(), 'SQB')
  graph.Draw('AP')
  canvas.Print(pdf)
  print('{0} +/- {1}, {2} +/- {3} (double exponential fit, with background subtracted, normalized to monitor detector)'.format(f.GetParams()[1], f.GetErrors()[1], f.GetParams()[3], f.GetErrors()[3]))
  
  # plot uncorrected UCN counts
  y = [float(c) for c in ex['li6counts']]
  yerr = [math.sqrt(c) for c in ex['li6counts']]
  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph.SetTitle('TCN{0} (single exponential fit + background, unnormalized)'.format(ex['TCN']))
  graph.GetXaxis().SetTitle('Storage time (s)')
  graph.GetYaxis().SetTitle('UCN count')
  # do single exponential fit with background
  f = graph.Fit(UCN.SingleExpoWithBackground(), 'SQB')
  graph.Draw('AP')
  canvas.Print(pdf)
  print('{0} +/- {1} (single exponential fit with {2} +/- {3} background, unnormalized)'.format(f.GetParams()[1], f.GetErrors()[1], f.GetParams()[2], f.GetErrors()[2]))

  mtau = []
  mtauerr = []
  for he3rate, m, s in zip(ex['he3rate'], ex['monitorduration'], ex['storageduration']):
    fitstart = m + 5
    fitend = m + s
    if fitend > fitstart + 10 and he3rate.Integral(he3rate.FindBin(fitstart), he3rate.FindBin(fitend)) > 0:
      f = he3rate.Fit(UCN.SingleExpo(), 'SQB', '', fitstart, fitend)
      he3rate.Draw()
      canvas.Print(pdf) # print fitted He3 rate to pdf
      mtau.append(f.GetParams()[1])
      mtauerr.append(f.GetErrors()[1])
	  
  # print average storage lifetime from He3 fits to pdf
  canvas = ROOT.TCanvas('c','c')
  if len(mtau) > 0:
    tauavg = numpy.average(mtau, None, [1./dt**2 for dt in mtauerr], True)
    label = ROOT.TLatex(0.1, 0.6, 'Average lifetime from He3 data:')
    label.Draw()
    label2 = ROOT.TLatex(0.1, 0.5, '#tau = {0} +/- {1} s'.format(tauavg[0], 1./tauavg[1]**2))
    label2.Draw()
    print('{0} +/- {1} (single exponential fit to rate in monitor detector during storage period)'.format(tauavg[0], 1./tauavg[1]**2))
  canvas.Print(pdf + ')')
  

ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1

# list runs for each experiment
experiments = [{'TCN': '18-025 (source-IV2, no elbow)', 'runs': [ 932]},
               {'TCN': '18-026 (IV1-IV2, no elbow)', 'runs': [ 933]},
               {'TCN': '18-032 (source-IV2, with elbow)', 'runs': [ 939]},
               {'TCN': '18-033 (IV1-IV2, with elbow)', 'runs': [ 940]},
               {'TCN': '18-036 (IV2-IV3, non-O-ring side)', 'runs': [ 949, 955, 956, 957, 958]}, #[947, 948, 949, 955, 956, 957, 958] # TCN18-036
               {'TCN': '18-040 (source-IV3, O-ring downstream)', 'runs': [ 950]},
               {'TCN': '18-041 (IV1-IV3, O-ring downstream)', 'runs': [ 951]},
               {'TCN': '18-042 (source-IV2, O-ring upstream)', 'runs': [ 952]},
               {'TCN': '18-046 (IV2-IV3, O-ring side)', 'runs': [ 961, 962, 967, 969]},
               {'TCN': '18-050 (source-IV3, O-ring upstream)', 'runs': [ 966]},
               {'TCN': '18-051 (IV1-IV3, O-ring upstream)', 'runs': [ 965]},
               {'TCN': '18-052 (IV1-IV2, O-ring downstream)', 'runs': [ 970]},
               {'TCN': '18-081 (UGD22+2, IV2-IV3)', 'runs': [ 976, 977]},
               {'TCN': '18-082 (UGD22+2, IV1-IV3)', 'runs': [ 974]},
               {'TCN': '18-055 (burst disk+UGD2, IV1-IV3)', 'runs': [ 983]},
               {'TCN': '18-054 (burst disk+UGD2, IV2-IV3)', 'runs': [ 986]},
               {'TCN': '18-087 (UGD22+19, IV1-IV3, MV open)', 'runs': [ 989]},
               {'TCN': '18-086 (UGD22+19, IV2-IV3, MV open)', 'runs': [ 991]},
               {'TCN': '18-087 (UGD22+19, IV1-IV3)', 'runs': [ 992]},
               {'TCN': '18-092 (UGD22+UGA11+UGG3+UGA5, IV1-IV3)', 'runs': [ 999]},
               {'TCN': '18-091 (UGD22+UGA11+UGG3+UGA5, IV2-IV3)', 'runs': [1001, 1002, 1003, 1004]},
               {'TCN': '18-291 (UGD22+UGA5+UGG3+UGA6, IV2-IV3)', 'runs': [1008]},
               {'TCN': '18-292 (UGD22+UGA5+UGG3+UGA6, IV1-IV3)', 'runs': [1010]},
               {'TCN': '18-061 (UGD10+17+11, IV2-IV3)', 'runs': [1014, 1015, 1016, 1018]},
               {'TCN': '18-062 (UGD10+17+11, IV1-IV3)', 'runs': [1017]},
			   
               {'TCN': '18-066_0A', 'runs': [1067]},
               {'TCN': '18-066_200A', 'runs': [1068]},
               {'TCN': '18-066_150A', 'runs': [1069]},
               {'TCN': '18-066_100A', 'runs': [1070]},
               {'TCN': '18-066_50A', 'runs': [1071]},
			   
               {'TCN': '18-068_0A', 'runs': [1072, 1073, 1074]},
               {'TCN': '18-068_200A', 'runs': [1075]},
               {'TCN': '18-068_150A', 'runs': [1076]},
               {'TCN': '18-068_50A', 'runs': [1077]},
               {'TCN': '18-068_100A', 'runs': [1078]},
			   
               {'TCN': '18-268_0A', 'runs': [1089]},
               {'TCN': '18-268_200A', 'runs': [1090]},
               {'TCN': '18-268_100A', 'runs': [1091]},
               {'TCN': '18-268_150A', 'runs': [1092]},
               {'TCN': '18-268_50A', 'runs': [1093]},
			   
               {'TCN': '18-266_0A', 'runs': [1094]},
               {'TCN': '18-266_200A', 'runs': [1095]},
               {'TCN': '18-266_100A', 'runs': [1096]},
               {'TCN': '18-266_50A', 'runs': [1097]},
               {'TCN': '18-266_150A', 'runs': [1098]},
			   
               {'TCN': '18-125 (NiP bottle, unbaked)', 'runs': [1118, 1119]},
               {'TCN': '18-126 (NiP bottle, baked 100C)', 'runs': [1122]},
               {'TCN': '18-127 (NiP bottle, baked 150C)', 'runs': [1124]},
               {'TCN': '18-116 (UGD22+2+Ti, IV2-IV3)', 'runs': [1126]},
               {'TCN': '18-117 (UGD22+2+Ti, IV1-IV3)', 'runs': [1127]},
               {'TCN': '18-481 (UGD22+2, IV2-IV3)', 'runs': [1136]},
               {'TCN': '18-482 (UGD22+2, IV1-IV3)', 'runs': [1134, 1135]},
               {'TCN': '18-058 (spider+UGD2, IV2-IV3)', 'runs': [1142]},
               {'TCN': '18-059 (spider+UGD2, IV1-IV3)', 'runs': [1143]},
               {'TCN': '18-216 (UGD22+20+Ti, high pos, IV2-IV3)', 'runs': [1182, 1183, 1184]},
               {'TCN': '18-217 (UGD22+20+Ti, high pos, IV1-IV3)', 'runs': [1185, 1186, 1187]},
               {'TCN': '18-381 (UGD22+20, high pos, IV2-IV3)', 'runs': [1189, 1190]},
               {'TCN': '18-382 (UGD22+20, high pos, IV1-IV3)', 'runs': [1191]}
			  ]

ReadCycles(ROOT.TFile(sys.argv[1]), experiments)
			  
# loop over experiments
for ex in experiments:
  StorageLifetime(ex)

UCN.PrintBackground(experiments)

canvas = ROOT.TCanvas('c','c')
for tcn in ['18-066', '18-068', '18-268', '18-266']:
  SCMex = [ex for ex in experiments if ex['TCN'].startswith(tcn)]
  gr = ROOT.TGraphErrors()
  for ex in SCMex:
    SCMcurrent = numpy.concatenate(ex['SCMcurrent'])
    i = gr.GetN()
    gr.SetPoint(i, numpy.mean(SCMcurrent), ex['tau'])
    gr.SetPointError(i, numpy.std(SCMcurrent)/math.sqrt(len(SCMcurrent)), ex['tauerr'])
  gr.SetTitle('TCN' + tcn)
  gr.GetXaxis().SetTitle('SCM current (A)')
  gr.GetYaxis().SetTitle('Storage lifetime (s)')
  gr.Draw('AP')
  canvas.Print('TCN{0}.pdf'.format(tcn))
