import ROOT
import sys
import numpy
import math
import UCN


def ReadCycles(infile, experiments):
  countperiod = 1
  monitorperiod = 0
  backgroundperiod = 10
  tofrange = (60, 120)
  binspersec = 10

  for ex in experiments:
    ex['start'] = []
    ex['cyclenumber'] = []
    ex['beamcurrent'] = []
    ex['li6counts'] = []
    ex['countduration'] = []
    ex['monitorcounts'] = []
    ex['monitorduration'] = []
    ex['li6background'] = []
    ex['he3background'] = []
    ex['li6irradiation'] = []
    ex['he3irradiation'] = []
    ex['irradiationduration'] = []
    ex['backgroundduration'] = []
    ex['beamcurrent'] = []
    ex['mintemperature'] = []
    ex['maxtemperature'] = []
    ex['minvaporpressure'] = []
    ex['maxvaporpressure'] = []
    ex['SCMcurrent'] = []

    ex['tofspectrum'] = ROOT.TH1D('TCN{0}'.format(ex['TCN']), 'TCN{0}'.format(ex['TCN']), (tofrange[1] - tofrange[0])*binspersec, tofrange[0], tofrange[1])
    ex['tofspectrum'].GetXaxis().SetTitle('Time in cycle (s)')
    ex['tofspectrum'].GetYaxis().SetTitle('UCN-rate-to-monitor ratio (per {0} s)'.format(1./binspersec))
    ex['tofspectrum'].SetBit(ROOT.TH1.kIsAverage)
    ex['tofspectrum'].SetDirectory(0)

  for cycle in infile.cycledata:
    run = cycle.runnumber
  
    if not any(run in ex['runs'] for ex in experiments): # if there is no experiment using this cycle
      continue
    
    Li6 = cycle.countsLi6
    He3 = cycle.countsHe3
    d = cycle.durations
   
    # filter useless runs
    if min(cycle.B1V_KSM_PREDCUR) < 0.1:
      print('SKIPPING cycle {0} in run {1} because beam current dropped below 0.1uA ({2}uA)'.format(cycle.cyclenumber, cycle.runnumber, min(cycle.B1V_KSM_PREDCUR)))
      continue
    if numpy.std(cycle.B1V_KSM_PREDCUR) > 0.02:
      print('SKIPPING cycle {0} in run {1} because beam current fluctuated by {2}'.format(cycle.cyclenumber, cycle.runnumber, numpy.std(cycle.B1V_KSM_PREDCUR)))
      continue
    if Li6[10] == 0:
      print('SKIPPING cycle {0} in run {1} because not all periods contain Li6 data'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    if He3[monitorperiod] < 1000:
      print('SKIPPING cycle {0} in run{1} because there are less than 1000 monitor counts ({2})'.format(cycle.cyclenumber, cycle.runnumber, He3[monitorperiod]))
      continue
    if any([1e-7 < ig5 < 1e-2 for ig5 in cycle.UCN_EXP_IG5_RDVAC]) or any([1e-7 < ig6 < 1e-2 for ig6 in cycle.UCN_EXP_IG6_RDVAC]):
      print ('SKIPPING cycle {0} in run {1} because IG5 or IG6 were on!'.format(cycle.cyclenumber, cycle.runnumber))
      continue  
    if Li6[countperiod] < 10*d[countperiod]:
      print('SKIPPING cycle {0} in run {1} because Li6 seems to see only background ({2}/s)'.format(cycle.cyclenumber, cycle.runnumber, Li6[countperiod]/d[countperiod]))
      continue
    if Li6[backgroundperiod] > 10*d[backgroundperiod]:
      print('SKIPPING cycle {0} in run {1} because Li6 sees high background ({2}/s)'.format(cycle.cyclenumber, cycle.runnumber, Li6[backgroundperiod]/d[backgroundperiod]))
      continue

    if cycle.valve0state[0] != 1 or cycle.valve0state[1] != 1 or cycle.valve1state[0] != 0 or cycle.valve1state[1] != 1:
      print('Abnormal valve configuration in cycle {0} of run {1}'.format(cycle.cyclenumber, cycle.runnumber))

    for ex in experiments:
      if run not in ex['runs']:
        continue

      ex['start'].append(cycle.start)
      ex['cyclenumber'].append(float(cycle.cyclenumber))
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
      ex['li6irradiation'].append(Li6[0])
      ex['irradiationduration'].append(d[0])

      li6hits, bins = numpy.histogram([hit for hit in getattr(cycle, 'Li6/hits')], (tofrange[1] - tofrange[0])*binspersec, (tofrange[0], tofrange[1]))
      norm, normerr = UCN.SubtractBackgroundAndNormalize(li6hits, [1./binspersec for _ in li6hits], 'li6', [He3[monitorperiod] for _ in li6hits], [math.sqrt(He3[monitorperiod]) for _ in li6hits])
      tof = ROOT.TH1D('tofspectrum', 'tofspectrum', (tofrange[1] - tofrange[0])*binspersec, tofrange[0], tofrange[1])
      for n, ne, binlo, binhi in zip(norm, normerr, bins[:-1], bins[1:]):
        b = tof.FindBin((binlo + binhi)/2.)
        tof.SetBinContent(b, n)
        tof.SetBinError(b, ne)
      tof.SetBit(ROOT.TH1.kIsAverage)
      ex['tofspectrum'].Add(tof)

	  
def Transmission(ex):
  print('\nAnalyzing TCN{0}'.format(ex['TCN']))

  if len(ex['start']) == 0:
    print('Found no cycles with run numbers {0}!'.format(ex['runs']))
    return

  ex['li6backgroundrate'], ex['li6backgroundrateerr'] = UCN.BackgroundRate(ex['li6background'], ex['backgroundduration'])
  ex['li6irradiationrate'], ex['li6irradiationrateerr'] = UCN.SubtractBackgroundAndNormalizeRate(ex['li6irradiation'], ex['irradiationduration'], 'li6', \
                                                                                                 [numpy.mean(cur) for cur in ex['beamcurrent']], [numpy.std(cur) for cur in ex['beamcurrent']])
  print 'Li6 background rate: {0} +/- {1} 1/s'.format(ex['li6backgroundrate'], ex['li6backgroundrateerr'])

  # report average monitor counts
  monitoravg = numpy.average(ex['monitorcounts'], None, [1./m for m in ex['monitorcounts']], True)
  print 'Monitor counts: {0} +/- {1}'.format(monitoravg[0], 1./math.sqrt(monitoravg[1]))

  # report range of beam current
  print('Beam current from {0} to {1} uA'.format(min(min(c) for c in ex['beamcurrent']), max(max(c) for c in ex['beamcurrent'])))

  # report He-II temperature range
  print 'Temperatures from {0} to {1} K'.format(min(ex['mintemperature']), max(ex['maxtemperature']))

  y, yerr = UCN.SubtractBackgroundAndNormalize(ex['li6counts'], ex['countduration'], 'li6', ex['monitorcounts'], [math.sqrt(m) for m in ex['monitorcounts']])

  # plot ratio of background-corrected counts to monitor counts
  canvas = ROOT.TCanvas('c', 'c')
  graph = ROOT.TGraphErrors(len(y), numpy.array(ex['cyclenumber']), numpy.array(y), numpy.array([0. for _ in ex['cyclenumber']]), numpy.array(yerr))
  graph.SetTitle('TCN' + ex['TCN'])
  graph.GetXaxis().SetTitle('Cycle')
  graph.GetYaxis().SetTitle('UCN-count-to-monitor ratio')
  graph.SetMarkerStyle(20)
  f = graph.Fit('pol0', 'QS')
  #canvas.SetLogy()
  graph.Draw('AP')
  canvas.Print('TCN{0}.pdf('.format(ex['TCN']))
  print('{0} +/- {1}'.format(f.GetParams()[0], f.GetErrors()[0]))
  
  ex['transmission'] = f.GetParams()[0]
  ex['transmissionerr'] = f.GetErrors()[0]

  ex['tofspectrum'].SetMinimum(0.)
  ex['tofspectrum'].Draw('')
  canvas.Print('TCN{0}.pdf)'.format(ex['TCN']))


# normalize one time-of-flight spectrum to another and print to pdf
def Normalize(experiments, transtcn, reftcn):
  trans = next((ex for ex in experiments if ex['TCN'].startswith(transtcn)), None) # find experiment with given TCN number
  ref = next((ex for ex in experiments if ex['TCN'].startswith(reftcn)), None) # find reference experiment with given TCN number
  if not trans or not ref:
    return 0., 0.

  transmission = trans['transmission']/ref['transmission']
  transmissionerr = math.sqrt((trans['transmissionerr']/trans['transmission'])**2 + (ref['transmissionerr']/ref['transmission'])**2)*transmission

  tofspec = trans['tofspectrum'].Clone() # make copy of tof spectrum
  tofspec.Divide(ref['tofspectrum']) # normalize to reference spectrum
  tofspec.GetYaxis().SetRangeUser(0, 1.5)
  tofspec.SetTitle('TCN{0} normalized to TCN{1}: {2} +/- {3}'.format(transtcn, reftcn, transmission, transmissionerr))
  tofspec.GetYaxis().SetTitle('Transmission')
  c = ROOT.TCanvas('c', 'c')
  tofspec.Draw('')
  c.Print('TCN{0}_TCN{1}.pdf'.format(transtcn, reftcn)) # print to pdf

  return transmission, transmissionerr


### Main program starts here ###

ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
#ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1

# list of runs belonging to each experiment
experiments = [{'TCN': '18-029 (IV2+UGD17, no elbow)', 'runs': [ 934]},
               {'TCN': '18-031 (IV2+UGD17+elbow)', 'runs': [ 938]},
               {'TCN': '18-035 (IV2+UGD22+IV3, O-rings out)', 'runs': [ 944]},
               {'TCN': '18-043 (IV2+UGD22+IV3, O-rings out, IV2 open during production)', 'runs': [ 954]},
               {'TCN': '18-045 (IV2+UGD22+IV3, O-rings in)', 'runs': [ 964]},
               {'TCN': '18-080 (UGD22+2, IV1 closed after production)', 'runs': [ 973]},
               {'TCN': '18-053 (burst disk+UGD2)', 'runs': [ 985]},
               {'TCN': '18-085 (UGD22+19, MV open)', 'runs': [ 990]},
               {'TCN': '18-085 (UGD22+19)', 'runs': [ 993]},
               {'TCN': '18-090 (UGD22+UGA11+UGG3+UGA5)', 'runs': [1000]},
               {'TCN': '18-290 (UGD22+UGA5+UGG3+UGA6)', 'runs': [1009]},
               {'TCN': '18-060 (UGD10+17+11)', 'runs': [1013]},
			   
               {'TCN': '18-065_0A', 'runs': [1054]},
               {'TCN': '18-065_25A', 'runs': [1055]},
               {'TCN': '18-065_200A', 'runs': [1056]},
               {'TCN': '18-065_100A', 'runs': [1057]},
               {'TCN': '18-065_150A', 'runs': [1058]},
               {'TCN': '18-065_175A', 'runs': [1059]},
               {'TCN': '18-065_125A', 'runs': [1064]},
               {'TCN': '18-065_75A', 'runs': [1065]},
               {'TCN': '18-065_50A', 'runs': [1066]},
			   
               {'TCN': '18-265_0A', 'runs': [1081]},
               {'TCN': '18-265_200A', 'runs': [1082]},
               {'TCN': '18-265_100A', 'runs': [1083]},
               {'TCN': '18-265_150A', 'runs': [1084]},
               {'TCN': '18-265_50A', 'runs': [1085]},
               {'TCN': '18-265_75A', 'runs': [1086]},
               {'TCN': '18-265_25A', 'runs': [1087]},
			   
               {'TCN': '18-115 (UGD22+2+Ti)', 'runs': [1125]},
               {'TCN': '18-245 (IV2+UGD22+IV3, O-rings in)', 'runs': [1129]},
               {'TCN': '18-480 (UGD22+2)', 'runs': [1131]},
               {'TCN': '18-480 (UGD22+2, IV1 closed after production)', 'runs': [1132, 1133]},
               {'TCN': '18-057 (spider+UGD2)', 'runs': [1141]},
               {'TCN': '18-302 (high position: IV2+elbow+UGD10+18)', 'runs': [1165]},
               {'TCN': '18-240 (high position: IV2+elbow+UGD10+Al+UGD18)', 'runs': [1176]},
               {'TCN': '18-215 (high position: UGD22+20+Ti)', 'runs': [1181]},
               {'TCN': '18-380 (high position: UGD22+20)', 'runs': [1188]},
               {'TCN': '18-310 (high position: UGD22+20, smooth elbow)', 'runs': [1192]}]

ReadCycles(ROOT.TFile(sys.argv[1]), experiments)

# loop over experiments and analyze transmission in runs
for ex in experiments:
  Transmission(ex)

UCN.PrintBackground(experiments, 'li6')

Normalize(experiments, '18-215', '18-380') # Ti foil in high position
Normalize(experiments, '18-115', '18-480') # Ti foil in low position
Normalize(experiments, '18-240', '18-302') # Al foil above detector
Normalize(experiments, '18-310', '18-380') # smooth elbow
Normalize(experiments, '18-057', '18-480') # spider
Normalize(experiments, '18-053', '18-480') # burst disk
Normalize(experiments, '18-085', '18-480') # NiP-coated guide
Normalize(experiments, '18-090', '18-480') # NiMo-coated glass guide
Normalize(experiments, '18-290', '18-480') # NiMo-coated glass guide
Normalize(experiments, '18-065_0A', '18-060') # warm bore without foil
Normalize(experiments, '18-265_0A', '18-060') # warm bore with foil

canvas = ROOT.TCanvas('c','c')
for tcn in ['18-065', '18-265']: # normalize all the SCM measurements to zero current and plot transmission vs. SCMcurrent
  gr = ROOT.TGraphErrors()
  for cur in ['0', '25', '50', '75', '100', '125', '150', '175', '200']:
    fulltcn = '{0}_{1}A'.format(tcn, cur)
    ex = next((ex for ex in experiments if ex['TCN'].startswith(fulltcn)), None)
    if not ex:
      continue
    t, te = Normalize(experiments, fulltcn, tcn+'_0A')
    SCMcurrent = numpy.concatenate(ex['SCMcurrent'])
    i = gr.GetN()
    gr.SetPoint(i, numpy.mean(SCMcurrent), t)
    gr.SetPointError(i, numpy.std(SCMcurrent)/math.sqrt(len(SCMcurrent)), te)
  gr.SetTitle('TCN' + tcn)
  gr.GetXaxis().SetTitle('SCM current (A)')
  gr.GetYaxis().SetTitle('Transmission')
  canvas = ROOT.TCanvas('c','c')
  gr.Draw('AP')
  canvas.Print('TCN{0}.pdf'.format(tcn))
