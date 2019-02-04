import ROOT
import math
import numpy

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


def SubtractBackgroundAndNormalizeToBeam(counts, countdurations, backgroundcounts, backgrounddurations, beamcurrent):
  brate = sum(backgroundcounts)/sum(backgrounddurations)
  brateerr = math.sqrt(sum(backgroundcounts))/sum(backgrounddurations)

  bgsub = [c - brate*cd for c, cd in zip(counts, countdurations)]
  bgsuberr = [math.sqrt(c + brateerr**2*cd**2) for c,cd in zip(counts, countdurations)]
 
  meanbeam = [numpy.mean(beam) for beam in beamcurrent]
  beamerr = [numpy.std(beam) for beam in beamcurrent]
  norm = [bgs/cur for bgs, cur in zip(bgsub, meanbeam)]
  normerr = [math.sqrt((bgserr/cur)**2 + (curerr*bgs/cur**2)**2) for bgserr, bgs, cur, curerr in zip(bgsuberr, bgsub, meanbeam, beamerr)]

  return brate, brateerr, norm, normerr

  
def SubtractBackgroundAndNormalizeToMonitor(counts, countdurations, backgroundrate, backgroundrateerr, monitorcounts):
  bgsub = [c - backgroundrate*cd for c, cd in zip(counts, countdurations)]
  bgsuberr = [math.sqrt(c + backgroundrateerr**2*cd**2) for c,cd in zip(counts, countdurations)]
 
  norm = [bgs/m for bgs, m in zip(bgsub, monitorcounts)]
  normerr = [math.sqrt((bgserr/m)**2 + (math.sqrt(m)*bgs/m**2)**2) for bgserr, bgs, m in zip(bgsuberr, bgsub, monitorcounts)]

  return norm, normerr


def PrintBackground(experiments, detector = 'li6', fitmin = 0, fitmax = 0):
  canvas = ROOT.TCanvas('c', 'c')
  bg = ROOT.TGraphErrors(len(experiments), numpy.array([float(min(ex['runs'])) for ex in experiments]), 
                                           numpy.array([ex[detector + 'backgroundrate'] for ex in experiments]),
                                           numpy.array([0. for _ in experiments]),
                                           numpy.array([ex[detector + 'backgroundrateerr'] for ex in experiments]))
  bg.SetTitle(detector + ' background')
  bg.GetXaxis().SetTitle('Run')
  bg.GetYaxis().SetTitle('Background rate (s^{-1})')
  bg.SetMarkerColor(ROOT.kRed)
  bg.SetMarkerStyle(20)
  bg.Draw('AP')

  lowbackground = [ex for ex in experiments if ex[detector + 'backgroundrate'] < 2.7]
  lowbg = ROOT.TGraphErrors(len(lowbackground), numpy.array([float(min(ex['runs'])) for ex in lowbackground]), 
                                                numpy.array([ex[detector + 'backgroundrate'] for ex in lowbackground]),
                                                numpy.array([0. for _ in lowbackground]),
                                                numpy.array([ex[detector + 'backgroundrateerr'] for ex in lowbackground]))
  lowbg.SetMarkerStyle(20)
  lowbg.Fit(ROOT.TF1('pol0','pol0'), 'Q', '', fitmin, fitmax)
  lowbg.Draw('PSAME')

  canvas.Print(detector + '_background.pdf')
