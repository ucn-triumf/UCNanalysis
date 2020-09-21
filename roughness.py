import sys
import openpyxl
import numpy
import math
import ROOT
import scipy.signal
import allantools

# Read profilometer scans from Excel reports and calculate different measures of roughness

ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)

c = ROOT.TCanvas('c', 'c')
for filename in sys.argv[1:]:
  suffix = ''
  if len(sys.argv) > 2:
    if filename == sys.argv[1]:
      suffix = '('
    elif filename == sys.argv[-1]:
      suffix = ')'
      
  print(filename)

  wb = openpyxl.load_workbook(filename)
  sheet = wb['DATA']
  pos = numpy.array([row.value*25.4e-3 for row in sheet['E'] if row.value])
  amp = numpy.array([row.value*25.4e-3*1e-6 for row in sheet['F'] if row.value])
  Ra = float(sheet['A1'].value.split()[1])*0.0254*1e-6 # Ra from report
  Ra2 = sum([abs(a) for a in amp])/len(amp) # Calculate Ra from raw scan data

  c.SetLogx(0)
  c.SetLogy(0)
  gr = ROOT.TGraph(len(amp), pos, amp)
  averoughness = numpy.convolve(amp, numpy.ones(20)/20, mode = 'valid') # Calculate running average over 20 samples
  grave = ROOT.TGraph(len(averoughness), pos[19:-19], averoughness)
  grave.SetLineColor(ROOT.kBlue)
  mg = ROOT.TMultiGraph()
  mg.Add(gr)
  mg.Add(grave) # Plot raw data and running average in one graph
  mg.SetTitle(filename + ' R_{{a}} = {0:.0f} / {1:.0f} nm;Position (m);Amplitude(m)'.format(Ra*1e9, Ra2*1e9))
  mg.GetYaxis().SetRangeUser(-5e-6,5e-6)
  mg.Draw('AL')
  c.Print('roughness.pdf' + suffix)

  c.SetLogx(1)
  c.SetLogy(1)

  spec = abs(numpy.fft.rfft(amp))[1:] # Calculate fourier transform of raw scan data
  wavelength = numpy.array([25.4e-3/2/k for k in range(1, len(spec) + 1)])
  mg = ROOT.TMultiGraph()
  mg.SetTitle(filename + '{0:.3g};Wavelength (m);Amplitude (m)'.format(numpy.mean(spec[1155:1411]))) # print average between 9um and 11um in title
  gr = ROOT.TGraph(len(spec), wavelength, spec)
  mg.Add(gr)

#  avespec = numpy.convolve(spec, numpy.ones(100)/100, mode = 'valid') # Calculate running average of fourier transform over 100 samples
#  grave = ROOT.TGraph(len(avespec), wavelength[99:-99], avespec)
#  grave.SetLineColor(ROOT.kBlue)
#  f = ROOT.TF1('f', '[0]*x^[1]')
#  f.SetParameters(0.1, 2)
#  grave.Fit(f, 'Q', '', 3e-6,1e-5)
#  mg.Add(grave)
#  mg.Draw('AL')
#  mg.GetXaxis().SetRangeUser(1e-6,1e-4)
  mg.GetYaxis().SetRangeUser(1e-8, 1e-3)
  mg.Draw('AL')
  c.Print('spectrum.pdf' + suffix)

  # calculate difference of each sample in a to local average of n neighboring samples
  def LocalDeviation(a, n):
    side = numpy.ones(int((n - 1)/2))/n
    window = numpy.concatenate((-side, [1. - 1./n], -side))
    conv = numpy.convolve(a, window, mode = 'valid')
    return math.sqrt(sum(numpy.power(conv, 2))/(n - 1))
  
  locdev = numpy.array([LocalDeviation(amp, n) for n in range(3,2001,2)]) # calculate difference to local averages with all averaging distances
  avgdist = numpy.array([n*(pos[1] - pos[0]) for n in range(3,2001,2)])
  grlocdev = ROOT.TGraph(len(locdev), avgdist, locdev)
  grlocdev.SetTitle(filename + ' ({0:.2g} #mum);Averaging distance (m);Standard deviation from local average (m)'.format(locdev[7]*1e6))
  grlocdev.GetYaxis().SetRangeUser(1e-7,1e-5)
  grlocdev.Draw('AP')
  c.Print('locdev.pdf' + suffix)

#  c.SetLogy(0)
#  c.SetLogz(1)
#  f, x, S = scipy.signal.spectrogram(amp, fs=1./1.5e-6, mode='magnitude')
#  h2 = ROOT.TH2D('h2', ';Frequency (1/m);Position (m)', len(f), f[0], f[-1], len(x), x[0], x[-1])
#  for i in range(len(f)):
#    for j in range(len(x)):
#      h2.Fill(f[i], x[j], S[i][j])
#  h2.Draw('COLZ')
#  c.Print('spectrogram.pdf' + suffix)

  c.SetLogx(0)
  dev = [math.degrees(math.atan((y2 - y1)/(x2 - x1))) for x1, x2, y1, y2 in zip(pos[:-1], pos[1:], amp[:-1], amp[1:])] # calculate angle to flat surface of all segments between neighboring samples
  hist = ROOT.TH1D('hist', filename+';Deviation from flat surface (deg);', 200, -10, 10)
  for v in dev:
    hist.Fill(v)
  hist.Fit('gaus', 'Q')
  hist.Draw()
  c.Print('deviation.pdf' + suffix)

#  spec2 = abs(numpy.fft.rfft(dev))[1:] # calculate fourier transform of angles calculated above
#  gr = ROOT.TGraph(len(spec2), wavelength, spec2)
#  c.SetLogx(1)
#  c.SetLogy(1)
#  gr.SetTitle(filename + ';Wavelength (m);Amplitude (deg)')
#  gr.GetYaxis().SetRangeUser(0.1, 1000.)
#  gr.Draw('AL')
#  c.Print('spectrum2.pdf' + suffix)
#
#  c.SetLogx(0)
#  c.SetLogy(0)
#  conv = numpy.correlate(amp[2000:-2000], amp, 'valid') # calculate auto-correlation function
#  conv = conv[len(conv)//2:]
#  correllen = next(x for x, y in zip(pos, conv) if y < conv[0]*0.1) # determine correlation length (where auto-correlation function drops below 10% of its maximum)
#  gr = ROOT.TGraph(len(conv), numpy.linspace(0., pos[2000], len(conv)), conv)
#  gr.SetTitle(filename + 'Correlation length (10%): {0:.2} m;Correlation length (m);Correlation'.format(correllen))
#  gr.GetXaxis().SetRangeUser(0., pos[2000])
#  gr.Draw('AL')
#  c.Print('correl.pdf' + suffix)
#
#  c.SetLogx(1)
#  taus2, ad, ade, ns = allantools.oadev(amp, 1./(pos[1] - pos[0])) # calculate allan deviation
#  gr = ROOT.TGraphErrors(len(taus2), taus2, ad, numpy.zeros(len(taus2)), ade)
#  gr.SetTitle(filename + ';Averaging length (m);Allan deviation')
#  gr.Draw('AP')
#  c.Print('allan.pdf' + suffix)

