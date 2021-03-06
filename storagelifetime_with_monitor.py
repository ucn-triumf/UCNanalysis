import ROOT
import sys
import numpy
import math
import itertools
import UCN

# Read cycles from file and map them to TCN numbers as given in 'experiments' array of dictionaries
# Fill array of dictionaries with data
def ReadCycles(infile, experiments):
  countperiod = 2
  monitorperiod = 0
  backgroundperiod = 1
  storageperiod = 1

  for ex in experiments:
    # initialize dictionary for each experiment
    ex['start'] = []                # list cycle-start times
    ex['cyclenumber'] = []          # list of cycle numbers
    ex['beamcurrent'] = []          # list of arrays of beam-current samples
    ex['li6counts'] = []            # list of Li6 counts during countperiod
    ex['countduration'] = []        # list of countperiod durations
    ex['monitorcounts'] = []        # list of He3 counts during monitorperiod
    ex['monitorduration'] = []      # list of monitorperiod durations
    ex['monitorcounts2'] = []       # list of He3 counts during time window at end of monitorperiod
    ex['monitorduration2'] = []     # list of time windows for monitorcounts2
    ex['li6background'] = []        # list of Li6 counts during background period
    ex['backgroundduration'] = []   # list of backgroundperiod durations
    ex['li6irradiation'] = []       # list of Li6 counts during irradiationperiod
    ex['irradiationduration'] = []  # list of irradiation durations
    ex['storageduration'] = []      # list of storage durations
    ex['mintemperature'] = []       # list of minimum He-II temperature during cycles (from temperature sensors)
    ex['maxtemperature'] = []       # list of maximum He-II temperature during cycles (from temperature sensors)
    ex['minvaporpressure'] = []     # list of minimum He-II temperature during cycles (from tvapor pressure)
    ex['maxvaporpressure'] = []     # list of maximum He-II temperature during cycles (from tvapor pressure)
    ex['SCMcurrent'] = []           # list of arrays of SCM current samples
    ex['li6rate'] = []              # list of Li6 count-rate histograms
    ex['he3rate'] = []              # list of He3 count-rate histograms
    ex['channels'] = ROOT.TH1D('TCN{0}_ch'.format(ex['TCN']), 'TCN{0};Channel;Count'.format(ex['TCN']), 10, 0, 10) # histogram of Li6 channels
    ex['channels'].SetDirectory(0)
    ex['IV1IV2'] = []               # list of boolean values that indicate if storage was between IV2 and IV3 (i.e. not connected to He3 detector)

  for cycle in infile.cycledata: # loop over all cycles
    run = cycle.runnumber
  
    if not any(run in ex['runs'] for ex in experiments): # skip if there is no experiment using this cycle
      continue
    
    Li6 = cycle.countsLi6
    He3 = cycle.countsHe3
    d = cycle.durations
   
    # filter useless cycles
    beam = [cur*bon for cur, bon, t in zip(cycle.B1V_KSM_PREDCUR, cycle.B1V_KSM_BONPRD, getattr(cycle, 'Beamline/timestamp')) if 1 < t - cycle.start < 59]
    if min(beam) < 0.1:
      print('SKIPPING cycle {0} in run {1} because beam dropped below 0.1uA ({2}uA)'.format(cycle.cyclenumber, cycle.runnumber, min(beam)))
      continue
    if numpy.std(beam) > 0.02:
      print('SKIPPING cycle {0} in run {1} because beam fluctuated by {2}uA'.format(cycle.cyclenumber, cycle.runnumber, numpy.std(beam)))
      continue
    if Li6[10] == 0:
      print('SKIPPING cycle {0} in run {1} because Li6 does not contain data in all periods'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    if max(cycle.UCN_UGD_IV1_STATON) < 1:
      print('SKIPPING cycle {0} in run {1} because IV1 never opened!'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    if He3[monitorperiod] < 1000:
      print('SKIPPING cycle {0} in run {1} because He3 saw less than 1000 monitor counts ({2})'.format(cycle.cyclenumber, cycle.runnumber, He3[monitorperiod]))
      continue
    if d[backgroundperiod] > 0 and Li6[backgroundperiod]/d[backgroundperiod] > 10:
      print('SKIPPING cycle {0} in run {1} because of high Li6 background ({2}/s)'.format(cycle.cyclenumber, cycle.runnumber, Li6[backgroundperiod]/d[backgroundperiod]))
      continue
    if (cycle.runnumber == 956 and cycle.cyclenumber == 56) or (cycle.runnumber == 1003 and cycle.cyclenumber == 41) or (cycle.runnumber == 1003 and cycle.cyclenumber == 45):
      print('SKIPPING cycle {0} in run {1}!'.format(cycle.cyclenumber, cycle.runnumber))
      continue
#    chist, bin_edges = numpy.histogram(getattr(cycle, 'Li6/channel'), 10, (0., 10.), True)
#    if sum(chist) > 1000 and not (chist[0] > chist[1] and chist[1] < chist[2] and chist[2] > chist[3] and chist[3] < chist[4] and chist[4] > chist[5] and chist[5] < chist[6] and chist[7] == 0 and chist[8] < chist[9]):
#      print('SKIPPING cycle {0} in run {1} because Li6 channel distribution looks weird!'.format(cycle.cyclenumber, cycle.runnumber))
#      continue
      
#    if (#any([1e-7 < ig5 < 1e-2 for ig5 in cycle.UCN_EXP_IG5_RDVAC]) or 
#        any([1e-7 < ig6 < 1e-2 for ig6 in cycle.UCN_EXP_IG6_RDVAC])):
#      print ('SKIPPING cycle {0} in run {1} because IG6 was on!'.format(cycle.cyclenumber, cycle.runnumber))
#      continue

    if (cycle.valve0state[0] != 1 or cycle.valve0state[1] != 0 or cycle.valve0state[2] != 0 or # IV1 should be open during irradiation, closed after
       cycle.valve1state[0] != 1 or # IV2 should be open during irradiation, state after depends on storage mode (between IV1+IV3 or between IV2+IV3)
       cycle.valve2state[0] != 0 or cycle.valve2state[1] != 0 or cycle.valve2state[2] != 1): # IV3 should be closed during irradiation and storage, open during counting
      print('Abnormal valve configuration in cycle {0} of run {1}'.format(cycle.cyclenumber, cycle.runnumber))
    
    for ex in experiments:
      if run not in ex['runs']: # find experiment this cycle belongs to
        continue

      # add data to experiment
      ex['start'].append(cycle.start)
      ex['cyclenumber'].append(float(cycle.cyclenumber))
      ex['beamcurrent'].append(beam)
      ex['mintemperature'].append(min([min(getattr(cycle,'UCN_ISO_{0}_RDTEMP'.format(TS))) for TS in ['TS11', 'TS12', 'TS14']]))
      ex['maxtemperature'].append(max([max(getattr(cycle,'UCN_ISO_{0}_RDTEMP'.format(TS))) for TS in ['TS11', 'TS12', 'TS14']]))
      if max(cycle.UCN_ISO_PG9L_RDPRESS) >= 2.:
        ex['minvaporpressure'].append(min(cycle.UCN_ISO_PG9H_RDPRESS))
        ex['maxvaporpressure'].append(max(cycle.UCN_ISO_PG9H_RDPRESS))
      else:
        ex['minvaporpressure'].append(min(cycle.UCN_ISO_PG9L_RDPRESS))
        ex['maxvaporpressure'].append(max(cycle.UCN_ISO_PG9L_RDPRESS))
      ex['SCMcurrent'].append([v/250e-6 for v in cycle.SCMVoltages3]) # calculate SCM current from voltage drop over 250uOhm shunt resistor
      ex['li6counts'].append(Li6[countperiod])
      ex['countduration'].append(d[countperiod])
      ex['monitorcounts'].append(He3[monitorperiod])
      ex['monitorduration'].append(d[monitorperiod])
      he3window = (-10., 0.)
      ex['monitorcounts2'].append(len([t for t in getattr(cycle, 'He3/hits') if d[monitorperiod] + he3window[0] < t < d[monitorperiod] + he3window[1]]))
      ex['monitorduration2'].append(he3window[1] - he3window[0])
      ex['li6background'].append(Li6[backgroundperiod])
      ex['backgroundduration'].append(d[backgroundperiod])
      ex['storageduration'].append(d[storageperiod])
      ex['irradiationduration'].append(d[0])
      ex['li6irradiation'].append(Li6[0])

      li6rate = ROOT.TH1I('Li6_{0}_{1}'.format(cycle.runnumber, cycle.cyclenumber), 'Li6 detector rate', int(math.floor(sum(d))), 0., math.floor(sum(d)))
      for h in getattr(cycle, 'Li6/hits'): # fill Li6 count-rate histogram
        li6rate.Fill(h)
      li6rate.GetXaxis().SetTitle('Time (s)')
      li6rate.GetYaxis().SetTitle('Li6 rate (s^{-1})')
      li6rate.SetDirectory(0)
      ex['li6rate'].append(li6rate)

      if cycle.valve1state[storageperiod] == 0 and cycle.valve0state[storageperiod] == 0:
        ex['IV1IV2'].append(True)
      he3rate = ROOT.TH1I('He3_{0}_{1}'.format(cycle.runnumber, cycle.cyclenumber), 'He3 detector rate', int(math.floor(sum(d))), 0., math.floor(sum(d)))
      for h in getattr(cycle, 'He3/hits'): # fill He3 count-rate histogram
        he3rate.Fill(h)
      he3rate.GetXaxis().SetTitle('Time (s)')
      he3rate.GetYaxis().SetTitle('He3 rate (s^{-1})')
      he3rate.SetDirectory(0)
      ex['he3rate'].append(he3rate)

      for c in getattr(cycle, 'Li6/channel'):
        ex['channels'].Fill(c)

  print('Read {0} cycles\n'.format(len(numpy.concatenate([ex['start'] for ex in experiments]))))
	  
	  
# analyze storage time from list of runs
def StorageLifetime(ex):
  print('\nAnalyzing TCN' + ex['TCN'])
  
  if len(ex['start']) == 0:
    print('Found no cycles with run numbers {0}!'.format(ex['runs']))
    return

  canvas = ROOT.TCanvas('c', 'c')
  pdf = 'TCN{0}.pdf'.format(ex['TCN'])

  # subtract background from Li6 counts and normalize to monitor counts
  y, yerr = UCN.SubtractBackgroundAndNormalize(ex['li6counts'], ex['countduration'], 'li6', ex['monitorcounts2'], [math.sqrt(m) for m in ex['monitorcounts2']])

  x = ex['storageduration']
  xerr = [0. for _ in ex['storageduration']]  

  # plot normalized, background corrected counts vs storage time
  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph.SetTitle('TCN{0} (single exponential fit, with background subtracted, normalized to monitor detector)'.format(ex['TCN']))
  graph.GetXaxis().SetTitle('Storage time (s)')
  graph.GetYaxis().SetTitle('UCN-count-to-monitor ratio')
  graph.SetMarkerStyle(20)

  canvas.SetLogy()

  #do single exponential fit with data point at 0 excluded
  f = graph.Fit(UCN.SingleExpo(), 'SQB', '', 1., 1000.)
  graph.SetTitle('TCN{0} (single exponential fit, with background subtracted, normalized to monitor detector, 0s excluded)'.format(ex['TCN']))
  graph.Draw('AP')
  canvas.Print(pdf + '(')
  ex['tau'] = f.GetParams()[1]
  ex['tauerr'] = f.GetErrors()[1]*max(math.sqrt(f.Chi2()/f.Ndf()), 1.0)

  # do double exponential fit
  graph.SetTitle('TCN{0} (double exponential fit, with background subtracted, normalized to monitor detector, 0s excluded)'.format(ex['TCN']))
  f = graph.Fit(UCN.DoubleExpo(), 'SQB', '', 1., 1000.)
  graph.Draw('AP')
  canvas.Print(pdf)
  print('{0} +/- {1}, {2} +/- {3} (double exponential fit, with background subtracted, normalized to monitor detector, 0s excluded)'.format(f.GetParams()[1], f.GetErrors()[1], f.GetParams()[3], f.GetErrors()[3]))
  
  # plot uncorrected UCN counts
#  y = [float(c) for c in ex['li6counts']]
#  yerr = [math.sqrt(c) for c in ex['li6counts']]
#  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
#  graph.SetTitle('TCN{0} (single exponential fit + background, unnormalized)'.format(ex['TCN']))
#  graph.GetXaxis().SetTitle('Storage time (s)')
#  graph.GetYaxis().SetTitle('UCN count')
  # do single exponential fit with background
#  f = graph.Fit(UCN.SingleExpoWithBackground(), 'SQB')
#  graph.Draw('AP')
#  canvas.Print(pdf)
#  print('{0} +/- {1} (single exponential fit with {2} +/- {3} background, unnormalized)'.format(f.GetParams()[1], f.GetErrors()[1], f.GetParams()[2], f.GetErrors()[2]))

  # draw plot of temperature during each cycle
  UCN.PrintTemperatureVsCycle(ex, pdf)

  mtau = []
  mtauerr = []
  # fit single exponential to He3 count-rate histogram during storage period (pinhole method)
  for he3rate, m, s, c in zip(ex['he3rate'], ex['monitorduration'], ex['storageduration'], ex['countduration']):
    fitstart = m + 5
    fitend = m + s
#    if not ex['pinhole']:
#      fitend = fitend + c
    if fitend > fitstart + 10 and he3rate.Integral(he3rate.FindBin(fitstart), he3rate.FindBin(fitend)) > 0:
      f = he3rate.Fit(UCN.SingleExpo(), 'SQB', '', fitstart, fitend)
      he3rate.SetTitle('TCN{0} (He3 rate)'.format(ex['TCN']))
      he3rate.Draw()
      canvas.Print(pdf) # print fitted He3 rate to pdf
      mtau.append(f.GetParams()[1])
      mtauerr.append(f.GetErrors()[1]*max(math.sqrt(f.Chi2()/f.Ndf()), 1.0))
	  
  # print average storage lifetime from He3 fits to pdf
  canvas = ROOT.TCanvas('c','c')
  if len(mtau) > 0:
    he3tau = ROOT.TGraphErrors(len(mtau), numpy.array(ex['cyclenumber']), numpy.array(mtau), numpy.array([0. for _ in mtau]), numpy.array(mtauerr))
    fit = he3tau.Fit('pol0', 'SQ')
    ex['pinholetau'] = fit.Parameter(0)
    ex['pinholetauerr'] = fit.ParError(0)*max(math.sqrt(f.Chi2()/f.Ndf()), 1.0)
    he3tau.SetMarkerStyle(20)
    he3tau.GetXaxis().SetTitle('Cycle')
    he3tau.GetYaxis().SetTitle('Pinhole storage lifetime (s)')
    he3tau.SetTitle('')
    he3tau.Draw('AP')
    print('{0} +/- {1} (single exponential fit to rate in monitor detector during storage period)'.format(fit.Parameter(0), fit.ParError(0)))
  else:
    ex['pinholetau'] = 0.
    ex['pinholetauerr'] = 0.
  canvas.Print(pdf)

  # draw plot of Li6 background rate during each cycle
  ex['li6backgroundrate'], ex['li6backgroundrateerr'] = UCN.PrintBackgroundVsCycle(ex, pdf, 'li6')
  print('Li6 detector background rate: {0} +/- {1} 1/s'.format(ex['li6backgroundrate'], ex['li6backgroundrateerr']))
  beam = [numpy.mean(cur) for cur in ex['beamcurrent']], [numpy.std(cur) for cur in ex['beamcurrent']]
  # subtract background from Li6 counts during irradiation and normalize to beam current, draw plot for each cycle
  ex['li6irradiationrate'], ex['li6irradiationrateerr'] = UCN.SubtractBackgroundAndNormalizeRate(ex['li6irradiation'], ex['irradiationduration'], 'li6', beam[0], beam[1])
  UCN.PrintIrradiationBackgroundVsCycle(ex, pdf, 'li6')

  # report average monitor counts, range of beam current, range of He-II temperature
  monitoravg = numpy.average(ex['monitorcounts2'], None, [1./m for m in ex['monitorcounts2']], True)
  print('Monitor counts: {0} +/- {1}'.format(monitoravg[0], 1./math.sqrt(monitoravg[1])))
  print('Beam current from {0} to {1} uA'.format(min(min(c) for c in ex['beamcurrent']), max(max(c) for c in ex['beamcurrent'])))
  print('Temperatures from {0} to {1} K'.format(min(ex['mintemperature']), max(ex['maxtemperature'])))

  # draw Li6 channel histogram
  ex['channels'].Draw()
  canvas.Print(pdf + ')')
  

# show additional fit statistics
ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptFit(1111)
# run ROOT silently
ROOT.gROOT.SetBatch(1)
# suppress stupid ROOT warnings
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

# read all data from file
ReadCycles(ROOT.TFile(sys.argv[1]), experiments)
			  
# loop over experiments
for ex in experiments:
  StorageLifetime(ex)

# draw plot of average background during each experiment
UCN.PrintBackground(experiments)

# plot storage lifetime vs SCM current for all SCM measurements
canvas = ROOT.TCanvas('c','c')
for tcn in ['18-066', '18-068', '18-268', '18-266']:
  SCMex = [ex for ex in experiments if ex['TCN'].startswith(tcn)]
  mg = ROOT.TMultiGraph()
  mg.SetTitle('TCN' + tcn)
  mg.GetXaxis().SetTitle('SCM current (A)')
  mg.GetYaxis().SetTitle('Storage lifetime (s)')
  for measurement in ['tau', 'pinholetau']: # plot both results from storage and pinhole measurement
    gr = ROOT.TGraphErrors()
    if measurement == 'pinholetau':
      if any(ex['IV1IV2']): # skip pinhole measurement if storage was between IV2 and IV3
        continue
      gr.SetLineColor(ROOT.kRed)
    for ex in SCMex:
      SCMcurrent = numpy.concatenate(ex['SCMcurrent'])
      i = gr.GetN()
      gr.SetPoint(i, numpy.mean(SCMcurrent), ex[measurement])
      gr.SetPointError(i, numpy.std(SCMcurrent)/math.sqrt(len(SCMcurrent)), ex[measurement+'err'])
    mg.Add(gr)
  mg.Draw('AP')
  canvas.Print('TCN{0}.pdf'.format(tcn))

# plot all storage lifetimes determined from pinhole measurements between IV1 and IV2
exps = [ex for ex in experiments if any(ex['IV1IV2'])]
pinhole = ROOT.TGraphErrors(len(exps),
                            numpy.array([float(min(ex['runs'])) for ex in exps]),
                            numpy.array([ex['pinholetau'] for ex in exps]),
                            numpy.array([0. for _ in exps]),
                            numpy.array([ex['pinholetauerr'] for ex in exps]))
pinhole.SetTitle(';Run;Storage lifetime between IV1 and IV2 (s)')
pinhole.Fit('pol0','Q','',960,1150)
pinhole.GetYaxis().SetRangeUser(15,40)
pinhole.Draw('AP')
canvas.Print('pinhole.pdf')

