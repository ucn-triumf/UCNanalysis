import ROOT
import sys
import numpy
import math
import time
import datetime
import itertools
import UCN

def UCNyield(x, p):
  scale = p[0] if x[0] >= 0 else p[5]
  T = UCN.HeTemperature(abs(x[0]))
  tau1 = 1./(1./p[1] + p[2]*(1. - p[3]*T)*T**7)
  tau2 = 1./(1. + p[4]*T**7)
#  print(tau1, tau2)
  if tau1 <= 0 or tau2 <= 0:
    return 0.
  return scale*tau1*tau2*(1. - math.exp(-x[2]*p[6]/tau1))*math.exp(-x[1]/tau1)


def ReadCycles(infile, experiments):
  countperiod = 2
  monitorperiod = 0
  backgroundperiod = 1
  storageperiod = 1

  for ex in experiments:
    ex['start'] = []
    ex['cyclenumber'] = []
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
    ex['meanvaporpressure'] = []
    ex['vaporpressurestd'] = []
    ex['channels'] = ROOT.TH1D('TCN{0}_ch'.format(ex['TCN']), 'TCN{0};Channel;Count'.format(ex['TCN']), 10, 0, 10)
    ex['channels'].SetDirectory(0)

  for cycle in infile.cycledata:
    run = cycle.runnumber
  
    if not any(run in ex['runs'] for ex in experiments): # if there is no experiment using this cycle
      continue
    Li6 = cycle.countsLi6
    He3 = cycle.countsHe3
    d = cycle.durations
   
    # filter useless cycles
    beam = [cur*bon for cur, bon, t in zip(cycle.B1V_KSM_PREDCUR, cycle.B1V_KSM_BONPRD, getattr(cycle, 'Beamline/timestamp')) if 1 < t - cycle.start < 59]
    if min(beam) < 0.1:
      print('SKIPPING cycle {0} in run {1} because beam current dropped to {2}uA!'.format(cycle.cyclenumber, cycle.runnumber, min(beam)))
      continue
    if numpy.std(beam) > 0.02:
      print('SKIPPING cycle {0} in run {1} because beam current fluctuated by {2}uA!'.format(cycle.cyclenumber, cycle.runnumber, numpy.std(beam)))
      continue
    if Li6[countperiod] > 0 and Li6[10] == 0: # indicates that run was stopped during counting
      print('SKIPPING cycle {0} in run {1} because not all periods contain Li6detector counts'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    if max(cycle.UCN_UGD_IV1_STATON) < 1:
      print('SKIPPING cycle {0} in run {1} because IV1 never opened!'.format(cycle.cyclenumber, cycle.runnumber))
      continue
#    if Li6[countperiod] > 80000:
#      print('SKIPPING cycle {0} in run {1} because Li6 saw too many counts ({2})'.format(cycle.cyclenumber, cycle.runnumber, Li6[countperiod]))
#      continue
#    if Li6[countperiod] > 0 and (any([1e-7 < ig5 < 1e-2 for ig5 in cycle.UCN_EXP_IG5_RDVAC])
#                             or any([1e-7 < ig6 < 1e-2 for ig6 in cycle.UCN_EXP_IG6_RDVAC])):
#      print ('SKIPPING cycle {0} in run {1} because IG5 and/or IG6 were on!'.format(cycle.cyclenumber, cycle.runnumber))
#      continue
    if Li6[countperiod] == 0 and He3[countperiod] == 0:
      print('SKIPPING cycle {0} in run {1} because it contains no detector counts'.format(cycle.cyclenumber, cycle.runnumber))
      continue

    # IV1 should be closed during irradiation and storage, all valves should be open during counting
    if (cycle.valve0state[0] != 0 or cycle.valve0state[1] != 0 or cycle.valve0state[2] != 1):
      print('Abnormal valve configuration in cycle {0} of run {1}'.format(cycle.cyclenumber, cycle.runnumber))
    
    
    for ex in experiments:
      if run not in ex['runs']:
        continue

      ex['start'].append(cycle.start)
      ex['cyclenumber'].append(float(cycle.cyclenumber))
      ex['beamcurrent'].append(beam)
      ex['mintemperature'].append(min([min(getattr(cycle,'UCN_ISO_{0}_RDTEMP'.format(TS))) for TS in ['TS11', 'TS12', 'TS14']]))
      ex['maxtemperature'].append(max([max(getattr(cycle,'UCN_ISO_{0}_RDTEMP'.format(TS))) for TS in ['TS11', 'TS12', 'TS14']]))
      vp = [PG9L if PG9L < 2. else PG9H + 0.18 for PG9L, PG9H, t in zip(cycle.UCN_ISO_PG9L_RDPRESS, cycle.UCN_ISO_PG9H_RDPRESS, getattr(cycle, 'Source/timestamp')) if t - cycle.start <= d[0]]
      ex['minvaporpressure'].append(min(vp))
      ex['maxvaporpressure'].append(max(vp))
      ex['meanvaporpressure'].append(numpy.mean(vp))
      ex['vaporpressurestd'].append(numpy.std(vp))
      for c in getattr(cycle, 'Li6/channel'):
        ex['channels'].Fill(c)
#      if False:(any([1e-7 < ig5 < 1e-2 for ig5 in cycle.UCN_EXP_IG5_RDVAC])
          #or
          #any([1e-7 < ig6 < 1e-2 for ig6 in cycle.UCN_EXP_IG6_RDVAC])
          #):
#      if ex['channels'].GetBinContent(ex['channels'].FindBin(3)) > ex['channels'].GetEntries()/10.:
#        print ('EXCLUDING Li6 data from cycle {0} in run {1} because channel distribution looks weird!'.format(cycle.cyclenumber, cycle.runnumber))
#        ex['li6counts'].append(0)
#        ex['li6background'].append(0)
#        ex['li6irradiation'].append(0)
#      else:
      ex['li6counts'].append(Li6[countperiod])
      ex['li6background'].append(Li6[backgroundperiod])
      ex['li6irradiation'].append(Li6[0])
      ex['he3counts'].append(He3[countperiod])
      ex['countduration'].append(d[countperiod])
      ex['he3background'].append(He3[backgroundperiod])
      ex['backgroundduration'].append(d[backgroundperiod])
      ex['storageduration'].append(d[storageperiod])
      ex['he3irradiation'].append(He3[0])
      ex['irradiationduration'].append(d[0])

  print('Read {0} cycles\n'.format(len(numpy.concatenate([ex['start'] for ex in experiments]))))


def DoCombinedFit(experiments, pdf, PreviousFit = None):
  if sum(sum(ex['li6counts']) + sum(ex['he3counts']) for ex in experiments) < 1000:
    for ex in experiments:
      ex['tau_wall'] = 0.
      ex['tau_wallerr'] = 0.
    return

  data = ROOT.Fit.BinData(2*sum([len(ex['start']) for ex in experiments]), 3, ROOT.Fit.BinData.kCoordError)
  for ex in experiments:
    for P, dP, ts, ti, li6, dli6, he3, dhe3 in zip(ex['meanvaporpressure'], ex['vaporpressurestd'], ex['storageduration'], ex['irradiationduration'],
                                                    ex['li6counts_normalized'], ex['li6counts_normalized_err'], ex['he3counts_normalized'], ex['he3counts_normalized_err']):
      if li6 > 0:
        data.Add(numpy.array([P, ts, ti]), li6, numpy.array([dP, 0., 0.]), dli6)
      if he3 > 0:
        data.Add(numpy.array([-P, ts, ti]), he3, numpy.array([dP, 0., 0.]), dhe3)


  f3 = ROOT.TF3('f3', UCNyield, -20., 20., 0., 180., 0., 120., 7)
  yieldfunc = ROOT.Math.WrappedMultiTF1(f3)
  fitter = ROOT.Fit.Fitter()
  fitter.SetFunction(yieldfunc, False)
  for i, setting in enumerate(zip([2800, 40, 0.3, 0.1, 0.1, 28, 5.], [10000, 100, 1, 1, 1, 10000, 30.])):
    fitter.Config().ParSettings(i).SetValue(setting[0])
    fitter.Config().ParSettings(i).SetLimits(0., setting[1])
  if PreviousFit:
    for i in [2, 3, 4, 6]:
      fitter.Config().ParSettings(i).SetValue(PreviousFit.Parameter(i))
      fitter.Config().ParSettings(i).Fix()
  fitter.Config().SetMinimizer('Minuit2')
  fitter.Fit(data)
  r = ROOT.TFitResult(fitter.Result())
#  r.Print()
  for ex in experiments:
    ex['tau_wall'] = r.Parameter(1)
    ex['tau_wallerr'] = r.Error(1)*max(math.sqrt(r.Chi2()/r.Ndf()), 1.)

  canvas = ROOT.TCanvas('cfit', 'cfit')
  canvas.SetLogy()
  n = data.Size()
  ddata = numpy.zeros(n)
  r.GetConfidenceIntervals(data, ddata, 0.683)
  grpress = {}
  for i in range(n):
    x = data.Coords(i)
    if x[2] == 30.:
      if x[1] not in grpress:
        grpress[x[1]] = (ROOT.TGraphErrors(), ROOT.TGraphErrors())
      gr = grpress[x[1]][0]
      grfit = grpress[x[1]][1]
      N = gr.GetN()
      gr.Set(N + 1)
      gr.SetPoint(N, abs(x[0]), data.Value(i))
      gr.SetPointError(N, data.CoordErrors(i)[0], data.Error(i))
      grfit.Set(N + 1)
      grfit.SetPoint(N, abs(x[0]), yieldfunc(x))
      grfit.SetPointError(N, 0., ddata[i])
    else:
      if -60. not in grpress:
        grpress[-60.] = (ROOT.TGraphErrors(), ROOT.TGraphErrors())
      gr = grpress[-60.][0]
      grfit = grpress[-60.][1]
      N = gr.GetN()
      gr.Set(N + 1)
      gr.SetPoint(N, x[1], data.Value(i))
      gr.SetPointError(N, data.CoordErrors(i)[1], data.Error(i))
      grfit.Set(N + 1)
      grfit.SetPoint(N, x[1], yieldfunc(x))
      grfit.SetPointError(N, 0., ddata[i])
    
  for i, ts in enumerate(grpress):
    gr = grpress[ts][0]
    grfit = grpress[ts][1]
    if ts != -60.:
      gr.SetTitle('{0}s irradiation, {1}s storage'.format(30., ts))
#      gr.SetMinimum(1e-3)
#      gr.SetMaximum(50000.)
      gr.GetXaxis().SetTitle('Vapor pressure (torr)')
    else:
      gr.SetTitle('{0}s irradiation'.format(60.))
      gr.GetXaxis().SetTitle('Storage time (s)')
    gr.GetYaxis().SetTitle('UCN count (#muA^{-1})')
    gr.Draw('AP')
    grfit.SetLineColor(ROOT.kRed)
    grfit.SetMarkerColor(ROOT.kRed)
    grfit.Draw('SAMEP')
    l = ROOT.TLatex()
    l.SetTextSize(0.03)
    for p in range(7):
      l.DrawLatexNDC(0.15, 0.21 + 0.03*(6 - p), 'p{0} = {1:.4g} #pm {2:.2g}'.format(p, r.Parameter(p), r.Error(p)))
    l.DrawLatexNDC(0.15, 0.18, '#chi^{{2}}/ndf = {0:.3g}'.format(r.Chi2()/r.Ndf()))
    l.DrawLatexNDC(0.15, 0.15, 'N(T, t_{s}, t_{i}) = p0|p5 #tau #frac{1}{1 + p4 T^{7}} (1 - exp(-#frac{t_{i} p6}{#tau})) exp(-#frac{t_{s}}{#tau})')
    if len(grpress) > 1 and i == 0:
      canvas.Print(pdf + '(')
    elif len(grpress) > 1 and i == len(grpress) - 1:
      canvas.Print(pdf + ')')
    else:
      canvas.Print(pdf)

  return r


# analyze storage time from experiment
def StorageLifetime(ex, FitResult = None):
  print('\nAnalyzing TCN{0} in runs {1}'.format(ex['TCN'], ex['runs']))

  if len(ex['start']) == 0:
    print('Found no cycles with run numbers {0}!'.format(ex['runs']))
    return

  # report start time
  print(datetime.datetime.fromtimestamp(min(ex['start'])))

  # report range of beam current
  print('Beam current from {0:.3} to {1:.3} uA'.format(min(min(c) for c in ex['beamcurrent']), max(max(c) for c in ex['beamcurrent'])))

  # report range of temperature and vapor pressure
  print('Temperatures from {0:.3} to {1:.3} K'.format(min(ex['mintemperature']), max(ex['maxtemperature'])))
  print('Vapor pressure from {0:.3} to {1:.3} torr'.format(min(ex['minvaporpressure']), max(ex['maxvaporpressure'])))

  beam = [numpy.mean(cur) for cur in ex['beamcurrent']], [numpy.std(cur) for cur in ex['beamcurrent']]
  canvas = ROOT.TCanvas('c1', 'c1')
  canvas.SetLogy()
  for det in ['li6', 'he3']:
    x = numpy.array(ex['storageduration'])
    xerr = numpy.array([0. for _ in x])
    # subtract background from UCN counts
    ex[det + 'counts_normalized'], ex[det + 'counts_normalized_err'] = UCN.SubtractBackgroundAndNormalize(ex[det + 'counts'], ex['countduration'], det, beam[0], beam[1])
    y = numpy.array(ex[det + 'counts_normalized'])
    yerr = numpy.array(ex[det + 'counts_normalized_err'])
 
    # plot normalized Li6 counts vs storage time
    graph = ROOT.TGraphErrors(len(x), x, y, xerr, yerr)
    graph.SetTitle('TCN{0} ({1} detector, single exponential fit, background subtracted, normalized to beam current)'.format(ex['TCN'], det))
    graph.GetXaxis().SetTitle('Storage time (s)')
    graph.GetYaxis().SetTitle('UCN count (#muA^{-1})')
#    graph.SetMarkerStyle(20)
    # do single exponential fit
    f = graph.Fit(UCN.SingleExpo(), 'SQB')
    if sum(y) > 0 and f.Error(1) < 5.:
      ex[det + 'tau'] = f.Parameter(1)
      ex[det + 'tauerr'] = f.Error(1)*max(math.sqrt(f.Chi2()/f.Ndf()), 1.)
    else:
      print('SKIPPING lifetime measurement from {0} detector in run(s) {1} because there were no counts detected or the error is larger than 5 s.'.format(det, ex['runs']))
      ex[det + 'tau'] = 0.
      ex[det + 'tauerr'] = 0.
    print('{0:.4} +/- {1:.2} ({2} detector, single exponential fit, background subtracted, normalized to beam current)'.format(ex[det + 'tau'], ex[det + 'tauerr'], det))
    graph.Draw('AP')
    pdf = 'TCN{0}_{1}.pdf'.format(ex['TCN'], ex['runs'][0])
    if det == 'li6':
      canvas.Print(pdf + '(')
    else:
      canvas.Print(pdf)

  if FitResult:
    r = DoCombinedFit([ex], pdf, FitResult)
    print('{0:.4} +/- {1:.2} (wall-storage lifetime from combined fit)'.format(ex['tau_wall'], ex['tau_wallerr']))

  canvas.SetLogy(0)

  UCN.PrintTemperatureVsCycle(ex, pdf)
  for det in ['li6', 'he3']:
    ex[det + 'backgroundrate'], ex[det + 'backgroundrateerr'] = UCN.PrintBackgroundVsCycle(ex, pdf, det)
    ex[det + 'irradiationrate'], ex[det + 'irradiationrateerr'] = UCN.SubtractBackgroundAndNormalizeRate(ex[det + 'irradiation'], ex['irradiationduration'], det, beam[0], beam[1])
    UCN.PrintIrradiationBackgroundVsCycle(ex, pdf, det)
    print(det + ' detector background rate: {0:.4} +/- {1:.2} 1/s'.format(ex[det + 'backgroundrate'], ex[det + 'backgroundrateerr']))

  ex['channels'].Draw()
  canvas.Print(pdf + ')')

  # return result from primary detector
  if max(ex['li6counts']) > max(ex['he3counts']):
    ex['tau'] = ex['li6tau']
    ex['tauerr'] = ex['li6tauerr']
  else:
    ex['tau'] = ex['he3tau']
    ex['tauerr'] = ex['he3tauerr']





ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1


# list runs for different source-storage experiments
experiments = [{'TCN': '19-010', 'runs': [1846]},
               {'TCN': '19-020', 'runs': [1873, 1875]},
               {'TCN': '19-190', 'runs': [1892]},
               {'TCN': '19-190', 'runs': [1902, 1905]},
               {'TCN': '19-191', 'runs': [1908]},
               {'TCN': '19-192', 'runs': [1918]},
               {'TCN': '19-250', 'runs': [1933]},
               {'TCN': '19-260', 'runs': [1940]}
              ]


#IsopureFillTimes = [time.mktime(dt.timetuple()) for dt in [datetime.datetime(2018, 11, 20, 16, 54), 
#                                                           datetime.datetime(2018, 11, 20, 22, 14),
#                                                           datetime.datetime(2018, 11, 24, 11, 27),
#                                                           datetime.datetime(2018, 11, 24, 18, 3)]]

ReadCycles(ROOT.TFile(sys.argv[1]), experiments)

for ex in experiments:
  StorageLifetime(ex)

UCN.PrintBackground(experiments, 'li6', 930, 1206)
UCN.PrintBackground(experiments, 'he3')

def TauVsTime(experiments, parameter, variable, timeformat, color):
  exs = [ex for ex in experiments if ex[variable] > 0.]
  x = numpy.array([float(min(ex[parameter])) for ex in exs])
  y = numpy.array([ex[variable] for ex in exs])
  yerr = numpy.array([ex[variable + 'err'] for ex in exs])
  grtaus = ROOT.TGraphErrors(len(exs), x, y, numpy.array([0. for _  in exs]), yerr)
  if timeformat:
    grtaus.GetXaxis().SetTimeDisplay(1)
    grtaus.GetXaxis().SetTimeFormat('%m-%d%F2019-01-01 00:00:00')
    grtaus.GetXaxis().SetNdivisions(10, 10, 0)
    grtaus.GetXaxis().SetTitle('Date')
  grtaus.SetTitle('Source-storage lifetime')
  grtaus.SetLineColor(color)
  grtaus.SetMarkerColor(color)
  grtaus.GetXaxis().SetTitle(parameter)
  grtaus.GetYaxis().SetTitle('Storage lifetime (s)')
  grtaus.GetYaxis().SetRangeUser(0,50)
  return grtaus


# plot storage lifetime vs time
canvas = ROOT.TCanvas('c','c')
l = ROOT.TLatex()
l.SetTextSize(0.03)
dailytaus = [ex for ex in experiments if len(ex['start']) > 0]

mg = ROOT.TMultiGraph('mg', ';Run;Lifetime (s)')
mg.Add(TauVsTime(dailytaus, 'runs', 'li6tau', False, ROOT.kBlack))
mg.Add(TauVsTime(dailytaus, 'runs', 'he3tau', False, ROOT.kRed))
mg.Draw('AP')
canvas.Print('dailytau.pdf(')

mg = ROOT.TMultiGraph('mg', ';Date;Lifetime (s)')
mg.Add(TauVsTime(dailytaus, 'start', 'li6tau', False, ROOT.kBlack))
mg.Add(TauVsTime(dailytaus, 'start', 'he3tau', False, ROOT.kRed))
mg.GetXaxis().SetTimeDisplay(1)
mg.GetXaxis().SetTimeFormat('%m-%d%F2019-01-01 00:00:00')
mg.GetXaxis().SetNdivisions(10, 10, 0)
mg.Draw('AP')
#box = ROOT.TBox(IsopureFillTimes[0], mg.GetHistogram().GetMinimum(), IsopureFillTimes[1], mg.GetHistogram().GetMaximum())
#box.Draw('SAME')
#box2 = ROOT.TBox(IsopureFillTimes[2], mg.GetHistogram().GetMinimum(), IsopureFillTimes[3], mg.GetHistogram().GetMaximum())
#box2.Draw('SAME')
canvas.Print('dailytau.pdf')

minvp = [min(ex['minvaporpressure']) for ex in dailytaus]
maxvp = [max(ex['maxvaporpressure']) for ex in dailytaus]
grtaus = ROOT.TGraphErrors(len(dailytaus),
                           numpy.array([float(min(ex['start'])) for ex in dailytaus]),
                           numpy.array([(ma + mi)/2. for ma, mi in zip(maxvp, minvp)]),
                           numpy.array([0.              for _  in dailytaus]),
                           numpy.array([(ma - mi)/math.sqrt(12) for ma, mi in zip(maxvp, minvp)]))
grtaus.GetXaxis().SetTitle('Date')
grtaus.GetXaxis().SetTimeDisplay(1)
grtaus.GetXaxis().SetTimeFormat('%m-%d%F2019-01-01 00:00:00')
grtaus.GetXaxis().SetNdivisions(10, 10, 0)
grtaus.GetYaxis().SetTitle('Vapor pressure (torr)')
grtaus.Draw('AP')
canvas.Print('dailytau.pdf)')

def TauVsTemp(experiments, parameter, variable, color = ROOT.kBlack, convert = False):
  exs = [ex for ex in experiments if ex[variable] > 0.]
  if not convert:
    x = numpy.array([(max(ex['max' + parameter]) + min(ex['min' + parameter]))/2 for ex in exs])
    xerr = numpy.array([(max(ex['max' + parameter]) - min(ex['min' + parameter]))/2 for ex in exs])
  else:
    x = numpy.array([(UCN.HeTemperature(max(ex['max' + parameter])) + UCN.HeTemperature(min(ex['min' + parameter])))/2 for ex in exs])
    xerr = numpy.array([(UCN.HeTemperature(max(ex['max' + parameter])) - UCN.HeTemperature(min(ex['min' + parameter])))/2 for ex in exs])
  gr = ROOT.TGraphErrors(len(exs), x, numpy.array([ex[variable] for ex in exs]), xerr, numpy.array([ex[variable + 'err'] for ex in exs]))
  gr.SetTitle('')

  if parameter == 'temperature':
    gr.GetXaxis().SetTitle('Temperature (K)')
  elif parameter == 'vaporpressure' and not convert:
    gr.GetXaxis().SetTitle('Vapor pressure (Torr)')
  elif parameter == 'vaporpressure' and convert:
    gr.GetXaxis().SetTitle('Temperature (K)')
  else:
    assert(True)
  gr.GetYaxis().SetTitle('Storage lifetime (s)')
  gr.SetLineColor(color)
  gr.SetMarkerColor(color)
  return gr


def FormatTaxis(Taxis):
  Taxis.SetTitle('Temperature (K)')
  Taxis.SetLabelFont(42)
  Taxis.SetLabelSize(0.035)
#  Taxis.SetLabelOffset(0.045)
  Taxis.SetTitleSize(0.035)
  Taxis.SetTitleFont(42)
#  Taxis.SetTitleOffset(1)

mg = ROOT.TMultiGraph('mg', ';Temperature (K);Storage Lifetime (s)')
mg.Add(TauVsTemp(dailytaus, 'temperature', 'tau', ROOT.kBlack)) # initial geometry
mg.Draw('AP')
canvas.Print('tauvstemp.pdf(')

canvas.SetLogx()
mg = ROOT.TMultiGraph('mg', ';Vapor pressure (Torr);Storage lifetime (s)')
mg.Add(TauVsTemp(dailytaus, 'vaporpressure', 'tau', ROOT.kBlack))
mg.Draw('AP')
fHeTemperature = ROOT.TF1('HeTemperature', lambda x: UCN.HeTemperature(x[0]), UCN.HeTemperature(mg.GetXaxis().GetXmin()), UCN.HeTemperature(mg.GetXaxis().GetXmax()))
Taxis = ROOT.TGaxis(mg.GetXaxis().GetXmin(), mg.GetHistogram().GetMaximum(), mg.GetXaxis().GetXmax(), mg.GetHistogram().GetMaximum(), 'HeTemperature', 510, '-')
FormatTaxis(Taxis)
Taxis.Draw()
canvas.Print('tauvstemp.pdf)')

