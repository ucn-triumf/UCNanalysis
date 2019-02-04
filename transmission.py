import ROOT
import sys
import numpy
import math
import UCN


def ReadCycles(infile, experiments):
  countperiod = 1
  monitorperiod = 0
  backgroundperiod = 10

  for ex in experiments:
    ex['start'] = []
    ex['cyclenumber'] = []
    ex['beamcurrent'] = []
    ex['li6counts'] = []
    ex['countduration'] = []
    ex['monitorcounts'] = []
    ex['monitorduration'] = []
    ex['li6background'] = []
    ex['backgroundduration'] = []
    ex['beamcurrent'] = []
    ex['mintemperature'] = []
    ex['maxtemperature'] = []
    ex['minvaporpressure'] = []
    ex['maxvaporpressure'] = []
    ex['SCMcurrent'] = []

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
    if Li6[countperiod] < 10*d[countperiod]:
      print('SKIPPING cycle {0} in run {1} because Li6 seems to see only background ({2}/s)'.format(cycle.cyclenumber, cycle.runnumber, Li6[countperiod]/d[countperiod]))
      continue
    if Li6[backgroundperiod] > 10*d[backgroundperiod]:
      print('SKIPPING cycle {0} in run {1} because Li6 sees high background ({2}/s)'.format(cycle.cyclenumber, cycle.runnumber, Li6[backgroundperiod]/d[backgroundperiod]))
      continue
    
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

	  
def Transmission(ex):
  print('\nAnalyzing TCN{0}'.format(ex['TCN']))

  ex['li6backgroundrate'] = sum(ex['li6background'])/sum(ex['backgroundduration'])
  ex['li6backgroundrateerr'] = math.sqrt(sum(ex['li6background']))/sum(ex['backgroundduration'])
  print 'Li6 background rate: {0} +/- {1} 1/s'.format(ex['li6backgroundrate'], ex['li6backgroundrateerr'])

  # report average monitor counts
  monitoravg = numpy.average(ex['monitorcounts'], None, [1./m for m in ex['monitorcounts']], True)
  print 'Monitor counts: {0} +/- {1}'.format(monitoravg[0], 1./math.sqrt(monitoravg[1]))

  # report range of beam current
  print('Beam current from {0} to {1} uA'.format(min(min(c) for c in ex['beamcurrent']), max(max(c) for c in ex['beamcurrent'])))

  # report He-II temperature range
  print 'Temperatures from {0} to {1} K'.format(min(ex['mintemperature']), max(ex['maxtemperature']))

  y, yerr = UCN.SubtractBackgroundAndNormalizeToMonitor(ex['li6counts'], ex['countduration'], 2.16, 0.01, ex['monitorcounts'])

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
  canvas.Print('TCN{0}.pdf'.format(ex['TCN']))
  print('{0} +/- {1}'.format(f.GetParams()[0], f.GetErrors()[0]))
  
  ex['transmission'] = f.GetParams()[0]
  ex['transmissionerr'] = f.GetErrors()[0]


### Main program starts here ###

ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
#ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1

# list of runs belonging to each experiment
experiments = [{'TCN': '18-029', 'runs': [ 934]},
               {'TCN': '18-031', 'runs': [ 938]},
               {'TCN': '18-035', 'runs': [ 944]},
               {'TCN': '18-043', 'runs': [ 954]},
               {'TCN': '18-045', 'runs': [ 964]},
               {'TCN': '18-080 (production up to IV1)', 'runs': [ 973]},
               {'TCN': '18-053', 'runs': [ 985]},
               {'TCN': '18-085 (MV open)', 'runs': [ 990]},
               {'TCN': '18-085', 'runs': [ 993]},
               {'TCN': '18-090', 'runs': [1000]},
               {'TCN': '18-290', 'runs': [1009]},
               {'TCN': '18-060', 'runs': [1012, 1013]},
			   
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
			   
               {'TCN': '18-115', 'runs': [1125]},
               {'TCN': '18-245', 'runs': [1129]},
               {'TCN': '18-480', 'runs': [1131]},
               {'TCN': '18-480 (production up to IV1)', 'runs': [1132, 1133]},
               {'TCN': '18-057', 'runs': [1141]},
               {'TCN': '18-302', 'runs': [1165]},
               {'TCN': '18-240', 'runs': [1176]},
               {'TCN': '18-215', 'runs': [1181]},
               {'TCN': '18-380', 'runs': [1188]},
               {'TCN': '18-310', 'runs': [1192]}]

ReadCycles(ROOT.TFile(sys.argv[1]), experiments)

# loop over experiments and analyze transmission in runs
for ex in experiments:
  Transmission(ex)

UCN.PrintBackground(experiments)

canvas = ROOT.TCanvas('c','c')
for tcn in ['18-065', '18-265']:
  SCMex = [ex for ex in experiments if ex['TCN'].startswith(tcn)]
  gr = ROOT.TGraphErrors()
  for ex in SCMex:
    SCMcurrent = numpy.concatenate(ex['SCMcurrent'])
    i = gr.GetN()
    gr.SetPoint(i, numpy.mean(SCMcurrent), ex['transmission'])
    gr.SetPointError(i, numpy.std(SCMcurrent)/math.sqrt(len(SCMcurrent)), ex['transmissionerr'])
  gr.SetTitle('TCN' + tcn)
  gr.GetXaxis().SetTitle('SCM current (A)')
  gr.GetYaxis().SetTitle('Transmission')
  gr.Draw('AP')
  canvas.Print('TCN{0}.pdf'.format(tcn))
