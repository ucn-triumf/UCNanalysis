import ROOT
import math
import numpy

DetectorBackground = {'li6': (2.16, 0.03), 'he3': (0.0403, 0.0017)}

def SingleExpo():
  SingleExpo = ROOT.TF1('SingleExpo', '[0]*exp(-x/[1])')
  SingleExpo.SetParameters(10, 10)
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


def SubtractBackgroundAndNormalize(counts, countdurations, detector, normalization, normalizationerr):
  bgsub = [c - DetectorBackground[detector][0]*cd if c > 0 else 0. for c, cd in zip(counts, countdurations)]
  bgsuberr = [math.sqrt(c + DetectorBackground[detector][1]**2*cd**2) if c > 0 else 0. for c,cd in zip(counts, countdurations)]
 
  norm = [bgs/m for bgs, m in zip(bgsub, normalization)]
  normerr = [math.sqrt((bgserr/m)**2 + (dm*bgs/m**2)**2) for bgserr, bgs, m, dm in zip(bgsuberr, bgsub, normalization, normalizationerr)]

  return norm, normerr


def SubtractBackgroundAndNormalizeRate(counts, countdurations, detector, normalization, normalizationerr):
  norm, normerr = SubtractBackgroundAndNormalize(counts, countdurations, detector, normalization, normalizationerr)
  return [n/d for n, d in zip(norm, countdurations)], [ne/d for ne, d in zip(normerr, countdurations)]


def BackgroundRate(counts, durations):
#  if len(counts) == 0 or len(durations) == 0:
#    return 0., 0.
  return float(sum(counts))/sum(durations), math.sqrt(sum(counts))/sum(durations)


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
    bg.SetMarkerColor(ROOT.kRed)
    bg.SetMarkerStyle(20)
    bg.Draw('AP')

  lowbackground = [ex for ex in bgexps if ex[detector + 'backgroundrate'] < 2.7]
  if len(lowbackground) > 0:
    lowbg = ROOT.TGraphErrors(len(lowbackground), numpy.array([float(min(ex['runs'])) for ex in lowbackground]), 
                                                  numpy.array([ex[detector + 'backgroundrate'] for ex in lowbackground]),
                                                  numpy.array([0. for _ in lowbackground]),
                                                  numpy.array([ex[detector + 'backgroundrateerr'] for ex in lowbackground]))
    lowbg.SetMarkerStyle(20)
    lowbg.Fit(ROOT.TF1('pol0','pol0'), 'Q', '', fitmin, fitmax)
    lowbg.Draw('PSAME')

    canvas.Print(detector + '_background.pdf')

  bgexps = [ex for ex in experiments if detector + 'irradiationrate' in ex and ex[detector + 'irradiationrateerr'] > 0]
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
  mh = ROOT.TH2I('monitorcounts', 'monitorcounts', 270, 930., 1200., 200, 0., 5000.)
  for ex in experiments:
    for m in ex['monitorcounts']:
      mh.Fill(float(min(ex['runs'])), m)
  mh.Draw('COL')
  canvas.Print('monitorcounts.pdf(')

  mh.Draw('CANDLE')
  canvas.Print('monitorcounts.pdf)')


