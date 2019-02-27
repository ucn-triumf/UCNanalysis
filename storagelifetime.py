import ROOT
import sys
import numpy
import math
import datetime
import itertools
import scipy.optimize
import UCN

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



def ReadCycles(infile, experiments):
  countperiod = 2
  monitorperiod = 0
  backgroundperiod = 1
  storageperiod = 1

  for ex in experiments:
    ex['start'] = []
    ex['beamcurrent'] = []
    ex['li6counts'] = []
    ex['he3counts'] = []
    ex['countduration'] = []
    ex['li6background'] = []
    ex['he3background'] = []
    ex['backgroundduration'] = []
    ex['li6irradiation'] = []
    ex['he3irradiation'] = []
    ex['irradiationduration'] = []
    ex['storageduration'] = []
    ex['beamcurrent'] = []
    ex['mintemperature'] = []
    ex['maxtemperature'] = []
    ex['minvaporpressure'] = []
    ex['maxvaporpressure'] = []

  for cycle in infile.cycledata:
    run = cycle.runnumber
  
    if not any(run in ex['runs'] for ex in experiments): # if there is no experiment using this cycle
      continue
    Li6 = cycle.countsLi6
    He3 = cycle.countsHe3
    d = cycle.durations
   
    # filter useless cycles
    if Li6[countperiod] == 0 and He3[countperiod] == 0:
      print('SKIPPING cycle {0} in run {1} because it contains no detector counts'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    if Li6[countperiod] > 0 and Li6[10] == 0: # indicates that run was stopped during counting
      print('SKIPPING cycle {0} in run {1} because not all periods contain Li6detector counts'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    beam = [cur for cur, bon in zip(cycle.B1V_KSM_PREDCUR, cycle.B1V_KSM_BONPRD)]
    if min(beam) < 0.1:
      print('SKIPPING cycle {0} in run {1} because beam current dropped to {2}uA!'.format(cycle.cyclenumber, cycle.runnumber, min(beam)))
      continue
    if numpy.std(beam) > 0.02:
      print('SKIPPING cycle {0} in run {1} because beam current fluctuated by {2}uA!'.format(cycle.cyclenumber, cycle.runnumber, numpy.std(beam)))
      continue
    if max(cycle.UCN_UGD_IV1_STATON) < 1:
      print('SKIPPING cycle {0} in run {1} because IV1 never opened!'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    if any([1e-7 < ig5 < 1e-2 for ig5 in cycle.UCN_EXP_IG5_RDVAC]) or any([1e-7 < ig6 < 1e-2 for ig6 in cycle.UCN_EXP_IG6_RDVAC]):
      print ('SKIPPING cycle {0} in run {1} because IG5 or IG6 were on!'.format(cycle.cyclenumber, cycle.runnumber))
      continue

    # IV1 should be closed during irradiation and storage, all valves should be open during counting
    if (cycle.valve0state[0] != 0 or cycle.valve0state[1] != 0 or cycle.valve0state[2] != 1 or cycle.valve1state[2] != 1 or cycle.valve2state[2] != 1):
      print('Abnormal valve configuration in cycle {0} of run {1}'.format(cycle.cyclenumber, cycle.runnumber))
    
    
    for ex in experiments:
      if run not in ex['runs']:
        continue

      ex['start'].append(cycle.start)
      ex['beamcurrent'].append(beam)
      ex['mintemperature'].append(min([min(getattr(cycle,'UCN_ISO_{0}_RDTEMP'.format(TS))) for TS in ['TS11', 'TS12', 'TS14']]))
      ex['maxtemperature'].append(max([max(getattr(cycle,'UCN_ISO_{0}_RDTEMP'.format(TS))) for TS in ['TS11', 'TS12', 'TS14']]))
      if max(cycle.UCN_ISO_PG9L_RDPRESS) >= 2.:
        ex['minvaporpressure'].append(min(cycle.UCN_ISO_PG9H_RDPRESS))
        ex['maxvaporpressure'].append(max(cycle.UCN_ISO_PG9H_RDPRESS))
      else:
        ex['minvaporpressure'].append(min(cycle.UCN_ISO_PG9L_RDPRESS))
        ex['maxvaporpressure'].append(max(cycle.UCN_ISO_PG9L_RDPRESS))
      ex['li6counts'].append(Li6[countperiod])
      ex['he3counts'].append(He3[countperiod])
      ex['countduration'].append(d[countperiod])
      ex['li6background'].append(Li6[backgroundperiod])
      ex['he3background'].append(He3[backgroundperiod])
      ex['backgroundduration'].append(d[backgroundperiod])
      ex['storageduration'].append(d[storageperiod])
      ex['li6irradiation'].append(Li6[0])
      ex['he3irradiation'].append(He3[0])
      ex['irradiationduration'].append(d[0])


# analyze storage time from experiment
def StorageLifetime(ex):
  print('\nAnalyzing TCN{0} in runs {1}'.format(ex['TCN'], ex['runs']))

  if len(ex['start']) == 0:
    print('Found no cycles with run numbers {0}!'.format(ex['runs']))
    return

  # report start time
  print(datetime.datetime.fromtimestamp(min(ex['start'])))

  # report range of beam current
  print('Beam current from {0} to {1} uA'.format(min(min(c) for c in ex['beamcurrent']), max(max(c) for c in ex['beamcurrent'])))

  # report range of temperature and vapor pressure
  print('Temperatures from {0} to {1} K'.format(min(ex['mintemperature']), max(ex['maxtemperature'])))
  print('Vapor pressure from {0} to {1} torr'.format(min(ex['minvaporpressure']), max(ex['maxvaporpressure'])))

  beam = [numpy.mean(cur) for cur in ex['beamcurrent']], [numpy.std(cur) for cur in ex['beamcurrent']]
  for det in ['li6', 'he3']:
    ex[det + 'backgroundrate'], ex[det + 'backgroundrateerr'] = UCN.BackgroundRate(ex[det + 'background'], ex['backgroundduration'])
    ex[det + 'irradiationrate'], ex[det + 'irradiationrateerr'] = UCN.SubtractBackgroundAndNormalizeRate(ex[det + 'irradiation'], ex['irradiationduration'], det, beam[0], beam[1])
    print(det + ' detector background rate: {0} +/- {1} 1/s'.format(ex[det + 'backgroundrate'], ex[det + 'backgroundrateerr']))

  x = ex['storageduration']
  xerr = [0. for _ in x]
  # subtract background from UCN counts
  y, yerr = UCN.SubtractBackgroundAndNormalize(ex['li6counts'], ex['countduration'], 'li6', beam[0], beam[1])
 
  # plot normalized Li6 counts vs storage time
  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph.SetTitle('TCN{0} (Li6 detector, single exponential fit, background subtracted, normalized to beam current)'.format(ex['TCN']))
  graph.GetXaxis().SetTitle('Storage time (s)')
  graph.GetYaxis().SetTitle('UCN count (#muA^{-1})')
#  graph.SetMarkerStyle(20)
  # do single exponential fit
  f = graph.Fit(UCN.SingleExpo(), 'SQB')
  li6tau = [f.GetParams()[1], f.GetErrors()[1]]
  canvas = ROOT.TCanvas('c', 'c')
  canvas.SetLogy()
  graph.Draw('AP')
  pdf = 'TCN{0}_{1}.pdf'.format(ex['TCN'], ex['runs'][0])
  canvas.Print(pdf + '(')
  print('{0} +/- {1} (Li6 detector, single exponential fit, background subtracted, normalized to beam current)'.format(li6tau[0], li6tau[1]))
  
  # plot uncorrected UCN counts
  y = [float(c) for c in ex['li6counts']]
  yerr = [math.sqrt(c) for c in ex['li6counts']]
  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph.SetTitle('TCN{0} (Li6 detector, single exponential fit + background, unnormalized)'.format(ex['TCN']))
  graph.GetXaxis().SetTitle('Storage time (s)')
  graph.GetYaxis().SetTitle('UCN count')
  # do single exponential fit with background
  f = graph.Fit(UCN.SingleExpoWithBackground(), 'SQB')
  graph.Draw('AP')
  canvas.Print(pdf)
  print('{0} +/- {1} (single exponential fit with {2} +/- {3} background, unnormalized)'.format(
           f.GetParams()[1], f.GetErrors()[1], f.GetParams()[2], f.GetErrors()[2])
        )

  # plot beam-normalized He3 counts vs storage time
  y, yerr = UCN.SubtractBackgroundAndNormalize(ex['he3counts'], ex['countduration'], 'he3', beam[0], beam[1])
  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph.SetTitle('TCN{0} (He3 detector, single exponential fit, background subtracted, normalized to beam current)'.format(ex['TCN']))
  graph.GetXaxis().SetTitle('Storage time (s)')
  graph.GetYaxis().SetTitle('UCN count (#muA^{-1})')
  # do single exponential fit
  f = graph.Fit(UCN.SingleExpo(), 'SQB')
  he3tau = [f.GetParams()[1], f.GetErrors()[1]]
  graph.Draw('AP')
  canvas.Print(pdf + ')')
  print('{0} +/- {1} (He3 detector, single exponential fit, background subtracted, normalized to beam current)'.format(he3tau[0], he3tau[1]))

  # return result from primary detector
  if max(ex['li6counts']) > max(ex['he3counts']):
    ex['tau'] = li6tau[0]
    ex['tauerr'] = li6tau[1]
  else:
    ex['tau'] = he3tau[0]
    ex['tauerr'] = he3tau[1]
    


ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1


# list runs for different source-storage experiments
#runs['TCN18-015'] = [[869], [870], [872], [895], [900], [907], [911], [923], [931], [941], [953], [968], [975], [984], [988], [998], [1011], [1019], [1030], [1049], [1053], [1088], [1123], [1137], [1151], [1167], [1179], [1193]] #Daily storage lifetimes
#runs['TCN18-300'] = [[1153], [1154], [1155], [1156], [1158], [1159], [1160], [1161], [1167]]
#runs['TCN18-170'] = [[1194], [1195, 1196, 1197], [1198], [1199], [1200], [1201], [1202], [1203], [1204], [1205]]
experiments = [{'TCN': '18-015', 'runs': [ 869]},
               {'TCN': '18-015', 'runs': [ 870]},
               {'TCN': '18-015', 'runs': [ 872]},
               {'TCN': '18-015', 'runs': [ 895]},
               {'TCN': '18-015', 'runs': [ 900]},
               {'TCN': '18-015', 'runs': [ 907]},
               {'TCN': '18-015', 'runs': [ 911]},
               {'TCN': '18-015', 'runs': [ 923]},
               {'TCN': '18-015', 'runs': [ 931]},
               {'TCN': '18-015', 'runs': [ 941]},
               {'TCN': '18-015', 'runs': [ 953]},
               {'TCN': '18-015', 'runs': [ 968]},
               {'TCN': '18-015', 'runs': [ 975]},
               {'TCN': '18-015', 'runs': [ 984]},
               {'TCN': '18-015', 'runs': [ 988]},
               {'TCN': '18-015', 'runs': [ 998]},
               {'TCN': '18-015', 'runs': [1011]},
               {'TCN': '18-015', 'runs': [1019]},
               {'TCN': '18-015', 'runs': [1030]},
               {'TCN': '18-015', 'runs': [1049]},
               {'TCN': '18-015', 'runs': [1053]},
               {'TCN': '18-015', 'runs': [1088]},
               {'TCN': '18-015', 'runs': [1123]},
               {'TCN': '18-015', 'runs': [1137]},
               {'TCN': '18-015', 'runs': [1151]},
               {'TCN': '18-015', 'runs': [1167]},
               {'TCN': '18-015', 'runs': [1179]},
               {'TCN': '18-015', 'runs': [1193]},

               {'TCN': '18-300', 'runs': [1153]},
               {'TCN': '18-300', 'runs': [1154]},
               {'TCN': '18-300', 'runs': [1155]},
               {'TCN': '18-300', 'runs': [1156]},
               {'TCN': '18-300', 'runs': [1158]},
               {'TCN': '18-300', 'runs': [1159]},
               {'TCN': '18-300', 'runs': [1160]},
               {'TCN': '18-300', 'runs': [1161]},
               {'TCN': '18-300', 'runs': [1167]},

               {'TCN': '18-170', 'runs': [1194]},
               {'TCN': '18-170', 'runs': [1195, 1196, 1197]},
               {'TCN': '18-170', 'runs': [1198]},
               {'TCN': '18-170', 'runs': [1199]},
               {'TCN': '18-170', 'runs': [1200]},
               {'TCN': '18-170', 'runs': [1201]},
               {'TCN': '18-170', 'runs': [1202]},
               {'TCN': '18-170', 'runs': [1203]},
               {'TCN': '18-170', 'runs': [1204]},
               {'TCN': '18-170', 'runs': [1205]}
              ]

injectedpress = [1.04, 1.995, 4.236, 7.953, 16.284, 32.425, 62.84, 315, 628, 1365]

ReadCycles(ROOT.TFile(sys.argv[1]), experiments)
for ex in experiments:
  StorageLifetime(ex)

UCN.PrintBackground(experiments, 'li6', 930, 1206)
UCN.PrintBackground(experiments, 'he3')
canvas = ROOT.TCanvas('c','c')

# plot storage lifetime vs time
dailytaus = [ex for ex in experiments if ex['TCN'] == '18-015' and len(ex['start']) > 0]
grtaus = ROOT.TGraphErrors(len(dailytaus),
                           numpy.array([min(ex['start']) for ex in dailytaus]),
                           numpy.array([ex['tau']        for ex in dailytaus]),
                           numpy.array([0.               for _  in dailytaus]),
                           numpy.array([ex['tauerr']     for ex in dailytaus]))
grtaus.GetXaxis().SetTimeDisplay(1)
grtaus.GetXaxis().SetNdivisions(10, 10, 0)
grtaus.GetXaxis().SetTitle('Date')
grtaus.GetYaxis().SetTitle('Storage lifetime (s)')
grtaus.GetYaxis().SetRangeUser(0,40)
grtaus.Draw('AP')
canvas.Print('dailytau.pdf')

# plot storage lifetime vs temperature
tauvstemp = [ex for ex in experiments if ex['TCN'] == '18-300']
grtemp = ROOT.TGraphErrors(len(tauvstemp),
                           numpy.array([(max(ex['maxtemperature']) + min(ex['mintemperature']))/2 for ex in tauvstemp]),
                           numpy.array([ex['tau']                                                 for ex in tauvstemp]),
                           numpy.array([(max(ex['maxtemperature']) - min(ex['mintemperature']))/2 for ex in tauvstemp]),
                           numpy.array([ex['tauerr']                                              for ex in tauvstemp]))
grtemp.GetXaxis().SetTitle('Temperature (K)')
grtemp.GetYaxis().SetTitle('Storage lifetime (s)')
grtemp.SetLineColor(ROOT.kRed)
grtemp.Draw('AP')
canvas.Print('tauvstemp.pdf')

# plot storage lifetime vs vapor pressure
mintemp = [HeTemperature(min(ex['minvaporpressure']))   for ex in tauvstemp]
maxtemp = [HeTemperature(max(ex['maxvaporpressure']))   for ex in tauvstemp]
grpress = ROOT.TGraphErrors(len(tauvstemp),
                            numpy.array([(maxT + minT)/2 for maxT, minT in zip(maxtemp, mintemp)]),
                            numpy.array([ex['tau']       for ex in tauvstemp]),
                            numpy.array([(maxT - minT)/2 for maxT, minT in zip(maxtemp, mintemp)]),
                            numpy.array([ex['tauerr']    for ex in tauvstemp]))
grpress.GetXaxis().SetTitle('Temperature derived from vapor pressure (K)')
grpress.GetYaxis().SetTitle('Storage lifetime (s)')
grpress.SetLineColor(ROOT.kBlue)
grpress.Draw('AP')
canvas.Print('tauvsvaporpress.pdf')

# plot storage lifetime vs source spoilage
tauvsspoil = [ex for ex in experiments if ex['TCN'] == '18-170']
grspoil = ROOT.TGraphErrors(len(tauvsspoil),
                            numpy.cumsum([p*0.1       for p  in injectedpress]),
                            numpy.array([ex['tau']    for ex in tauvsspoil]),
                            numpy.array([0.           for _  in tauvsspoil]),
                            numpy.array([ex['tauerr'] for ex in tauvsspoil]))
canvas.SetLogx()
grspoil.GetXaxis().SetTitle('Spoiling gas injected (torr L)')
grspoil.GetYaxis().SetTitle('Storage Lifetime (s)')
grspoil.Draw('AP')
canvas.Print('tauvsspoilage.pdf')

