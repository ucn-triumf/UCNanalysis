import ROOT
import math
import numpy
import scipy.optimize

DetectorBackground = {'li6': (1.5, 0.1), 'he3': (0.03, 0.003)}


# calculate 4He vapor pressure from temperature
def HeVaporPressure(T):
  # from Clement, Logan, Gaffney, Phys. Rev. 100, 743
  # https://doi.org/10.1103/PhysRev.100.743
  if not 0.66 <= T <= 5.2:
    raise Exception('Tried to evaluate vapor pressure at T = {0}. Formula only valid between 0.66 and 5.2K!'.format(T))
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
  if P == 0.:
    return 0.
  return scipy.optimize.brentq(lambda T: HeVaporPressure(T) - P, 0.66, 5.2)


def SingleExpo():
  SingleExpo = ROOT.TF1('SingleExpo', '[0]*exp(-x/[1])')
  SingleExpo.SetParameters(10, 30)
  SingleExpo.SetParName(1, '#tau')
  SingleExpo.SetParLimits(0, 0, 1e6)
  SingleExpo.SetParLimits(1, 0, 1000)
  return SingleExpo

def SingleExpoWithBackground():
  SingleExpoWithBackground = ROOT.TF1('SingleExpoWithBackground', '[0]*exp(-x/[1]) + [2]')
  SingleExpoWithBackground.SetParameters(1000, 15, 200)
  SingleExpoWithBackground.SetParName(1, '#tau')
  SingleExpoWithBackground.SetParName(2, 'Background')
  SingleExpoWithBackground.SetParLimits(0, 0, 1e6)
  SingleExpoWithBackground.SetParLimits(1, 0, 1000)
  SingleExpoWithBackground.SetParLimits(2, 0, 1e6)
  return SingleExpoWithBackground

def DoubleExpo():
  DoubleExpo = ROOT.TF1('DoubleExpo', '[0]*exp(-x/[1]) + [2]*exp(-x/[3])')
  for i, param in enumerate([('N_{1}', 1, 1e6), ('#tau_{1}', 10, 1e6), ('N_{2}', 0.1, 1e6), ('#tau_{2}', 50, 1e6)]):
    DoubleExpo.SetParName(i, param[0])
    DoubleExpo.SetParameter(i, param[1])
    DoubleExpo.SetParLimits(i, 0, param[2])
  return DoubleExpo


def SubtractBackgroundAndNormalize(counts, countdurations, detector, normalization, normalizationerr, backgroundrate = None, backgroundrateerr = None):
  if not backgroundrate:
    backgroundrate = DetectorBackground[detector][0]
  if not backgroundrateerr:
    backgroundrateerr = DetectorBackground[detector][1]
  bgsub = [c - backgroundrate*cd if c > 0 else 0. for c, cd in zip(counts, countdurations)]
  bgsuberr = [math.sqrt(c + backgroundrateerr**2*cd**2) if c > 0 else 0. for c,cd in zip(counts, countdurations)]
 
  norm = [bgs/m for bgs, m in zip(bgsub, normalization)]
  normerr = [math.sqrt((bgserr/m)**2 + (dm*bgs/m**2)**2) for bgserr, bgs, m, dm in zip(bgsuberr, bgsub, normalization, normalizationerr)]

  return norm, normerr


def SubtractBackgroundAndNormalizeRate(counts, countdurations, detector, normalization, normalizationerr, backgroundrate = None, backgroundrateerr = None):
  if not backgroundrate:
    backgroundrate = DetectorBackground[detector][0]
  if not backgroundrateerr:
    backgroundrateerr = DetectorBackground[detector][1]
  norm, normerr = SubtractBackgroundAndNormalize(counts, countdurations, detector, normalization, normalizationerr, backgroundrate, backgroundrateerr)
  return [n/d for n, d in zip(norm, countdurations)], [ne/d for ne, d in zip(normerr, countdurations)]


def PrintBackground(experiments, detector = 'li6', fitmin = 0, fitmax = 0):
  canvas = ROOT.TCanvas('c', 'c')
  bgexps = [ex for ex in experiments if detector + 'backgroundrate' in ex and ex[detector + 'backgroundrateerr'] > 0]
  if len(bgexps) > 0:
    bg = ROOT.TGraphErrors(len(bgexps), numpy.array([float(min(ex['runs'])) for ex in bgexps]), 
                                             numpy.array([ex[detector + 'backgroundrate'] for ex in bgexps]),
                                             numpy.array([0. for _ in bgexps]),
                                             numpy.array([ex[detector + 'backgroundrateerr'] for ex in bgexps]))
    bg.SetTitle(detector + ' background')
    bg.GetXaxis().SetTitle('Run')
    bg.GetYaxis().SetTitle('Background rate (s^{-1})')
#    bg.SetMarkerColor(ROOT.kRed)
    bg.SetMarkerStyle(20)
    bg.Fit('pol0', 'Q')
    bg.Draw('AP')

#  lowbackground = [ex for ex in bgexps if ex[detector + 'backgroundrate'] < 2.5]
#  if len(lowbackground) > 0:
#    lowbg = ROOT.TGraphErrors(len(lowbackground), numpy.array([float(min(ex['runs'])) for ex in lowbackground]), 
#                                                  numpy.array([ex[detector + 'backgroundrate'] for ex in lowbackground]),
#                                                  numpy.array([0. for _ in lowbackground]),
#                                                  numpy.array([ex[detector + 'backgroundrateerr'] for ex in lowbackground]))
#    lowbg.SetMarkerStyle(20)
#    lowbg.Fit(ROOT.TF1('pol0','pol0'), 'Q', '', fitmin, fitmax)
#    lowbg.Draw('PSAME')

    canvas.Print(detector + '_background.pdf')

  bgexps = [ex for ex in experiments if detector + 'irradiationrate' in ex]
  if len(bgexps) > 0:
    irrbg = ROOT.TGraphErrors(len(numpy.concatenate([ex['start'] for ex in bgexps])),
                              numpy.concatenate([[float(min(ex['runs'])) for _ in ex['start']] for ex in bgexps]),
                              numpy.concatenate([ex[detector + 'irradiationrate'] for ex in bgexps]),
                              numpy.concatenate([[0. for _ in ex['start']] for ex in bgexps]),
                              numpy.concatenate([ex[detector + 'irradiationrateerr'] for ex in bgexps]))
    irrbg.SetMarkerStyle(20)
    irrbg.GetXaxis().SetTitle('Run')
    irrbg.GetYaxis().SetTitle('Added background rate during irradiation (s^{-1} #muA^{-1})')
    irrbg.Fit(ROOT.TF1('pol0','pol0'), 'Q', '', fitmin, fitmax)
    irrbg.Draw('AP')
  canvas.Print(detector + '_irradiationbackground.pdf')


def PrintMonitorCounts(experiments):
  canvas = ROOT.TCanvas('c', 'c')
  mh = ROOT.TH2I('monitorcounts', 'monitorcounts', 160, 1840., 2000., 200, 0., 1500.)
  for ex in experiments:
    for m in ex['monitorcounts']:
      mh.Fill(float(min(ex['runs'])), m)
  mh.Draw('COL')
  canvas.Print('monitorcounts.pdf(')

  mh.Draw('CANDLE')
  canvas.Print('monitorcounts.pdf)')

def PrintTemperatureVsCycle(ex, pdf):
  canvas = ROOT.TCanvas('cm', 'cm')
  graph = ROOT.TGraphErrors(len(ex['minvaporpressure']), numpy.array(ex['cyclenumber']),
                            numpy.array([(maxvp + minvp)/2 for maxvp, minvp in zip(ex['maxvaporpressure'], ex['minvaporpressure'])]),
                            numpy.array([0. for _ in ex['cyclenumber']]),
                            numpy.array([(maxvp - minvp)/2 for maxvp, minvp in zip(ex['maxvaporpressure'], ex['minvaporpressure'])]))
  graph.GetXaxis().SetTitle('Cycle')
  graph.GetXaxis().SetLimits(0., max(ex['cyclenumber']))
  graph.GetYaxis().SetTitle('Vapor pressure (torr)')
  graph.SetMarkerStyle(20)
  graph.Draw('AP')
  fHeTemperature = ROOT.TF1('HeTemperature', lambda x: HeTemperature(x[0]), HeTemperature(graph.GetHistogram().GetMinimum()), HeTemperature(graph.GetHistogram().GetMaximum()))
  Taxis = ROOT.TGaxis(max(ex['cyclenumber']), graph.GetHistogram().GetMinimum(), max(ex['cyclenumber']), graph.GetHistogram().GetMaximum(), 'HeTemperature', 510, '+')
  Taxis.SetTitle('Temperature (K)')
  Taxis.SetLabelFont(42)
  Taxis.SetLabelSize(0.035)
  Taxis.SetLabelOffset(0.045)
  Taxis.SetTitleSize(0.035)
  Taxis.SetTitleFont(42)
  Taxis.SetTitleOffset(1)
  Taxis.Draw()
  canvas.Print(pdf)

def PrintBackgroundVsCycle(ex, pdf, detector):
  if detector + 'background' not in ex:
    return 0., 0.
  canvas = ROOT.TCanvas('bg', 'bg')
  x = numpy.array([float(c) for c, bg in zip(ex['cyclenumber'], ex[detector + 'background'])])
  y = numpy.array([bg/bd if bd > 0. else 0. for bg, bd in zip(ex[detector + 'background'], ex['backgroundduration'])])
  ye = numpy.array([math.sqrt(bg)/bd for bg, bd in zip(ex[detector + 'background'], ex['backgroundduration']) if bg > 0])
  yeh = numpy.array([0.5*ROOT.TMath.ChisquareQuantile(1 - 0.317/2, 2*bg + 2)/bd - bg/bd for bg, bd in zip(ex[detector + 'background'], ex['backgroundduration'])])
  yel = numpy.array([bg/bd - 0.5*ROOT.TMath.ChisquareQuantile(0.317/2, 2*bg)/bd for bg, bd in zip(ex[detector + 'background'], ex['backgroundduration'])])
  graph = ROOT.TGraphAsymmErrors(len(x), x, y, numpy.array([0. for _ in x]), numpy.array([0. for _ in x]), yel, yeh)
  graph.SetTitle('')
  graph.GetXaxis().SetTitle('Cycle')
  graph.GetYaxis().SetTitle(detector + ' background rate (1/s)')
  graph.SetMarkerStyle(20)
  fit = graph.Fit('pol0', 'SQ')
  graph.Draw('AP')
  canvas.Print(pdf)
  if fit.Ndf() > 0 and fit.Chi2()/fit.Ndf() > 3.:
    print('High background variation!')
  return fit.Parameter(0), fit.ParError(0)

def PrintIrradiationBackgroundVsCycle(ex, pdf, detector):
  if detector + 'irradiationrate' not in ex or sum(ex[detector + 'irradiationrate']) == 0:
    return
  canvas = ROOT.TCanvas('bg', 'bg')
  x = numpy.array([float(c) for c in ex['cyclenumber']])
  y = numpy.array(ex[detector + 'irradiationrate'])
  ye = numpy.array(ex[detector + 'irradiationrateerr'])
  graph = ROOT.TGraphErrors(len(x), x, y, numpy.array([0. for _ in x]), ye)
  graph.SetTitle('')
  graph.GetXaxis().SetTitle('Cycle')
  graph.GetYaxis().SetTitle(detector + ' background rate during irradiation (1/(s #muA))')
  graph.SetMarkerStyle(20)
  fit = graph.Fit('pol0', 'SQ')
  graph.Draw('AP')
  canvas.Print(pdf)
  if fit.Chi2()/fit.Ndf() > 3.:
    print('High irradiation-background variation!')

