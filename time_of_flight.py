import ROOT
import sys
import numpy
import math
import UCN


def ReadCycles(infile, experiments):
  countperiod = 1
  monitorperiod = 0
  tofrange = (62, 150)

  for ex in experiments:
    ex['start'] = []
    ex['cyclenumber'] = []
    ex['tofspectrum'] = ROOT.TH1D('TCN{0}'.format(ex['TCN']), 'TCN{0}'.format(ex['TCN']), tofrange[1] - tofrange[0], tofrange[0], tofrange[1])
    ex['tofspectrum'].GetXaxis().SetTitle('Time in cycle (s)')
    ex['tofspectrum'].SetBit(ROOT.TH1.kIsAverage)
    ex['tofspectrum'].SetDirectory(0)
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
    
    for ex in experiments:
      if run not in ex['runs']:
        continue

      ex['start'].append(cycle.start)
      ex['cyclenumber'].append(float(cycle.cyclenumber))
      ex['SCMcurrent'].append([v/250e-6 for v in cycle.SCMVoltages3])
  
      duration = cycle.beamonduration + cycle.beamoffduration
      li6hits, bins = numpy.histogram([hit for hit in getattr(cycle, 'Li6/hits')], int(math.floor(duration)), (0., math.floor(duration)))
      norm, normerr = UCN.SubtractBackgroundAndNormalizeToMonitor(li6hits, [1.0 for _ in li6hits], 2.16, 0.01, [He3[monitorperiod] for _ in li6hits])
      tof = ROOT.TH1D('tofspectrum', 'tofspectrum', tofrange[1] - tofrange[0], tofrange[0], tofrange[1])
      for n, ne, binlo, binhi in zip(norm, normerr, bins[:-1], bins[1:]):
        b = tof.FindBin((binlo + binhi)/2.)
        tof.SetBinContent(b, n)
        tof.SetBinError(b, ne)
      tof.SetBit(ROOT.TH1.kIsAverage)
      ex['tofspectrum'].Add(tof)


# normalize one time-of-flight spectrum to another and print to pdf
def NormalizeTOF(experiments, toftcn, reftcn):
  tof = next((ex for ex in experiments if ex['TCN'] == toftcn), None) # find tof experiment with given TCN number
  ref = next((ex for ex in experiments if ex['TCN'] == reftcn), None) # find reference experiment with given TCN number
  if not tof or not ref:
    return
  tofspec = tof['tofspectrum'].Clone() # make copy of tof spectrum
  tofspec.Divide(ref['tofspectrum']) # normalize to reference spectrum
  c = ROOT.TCanvas('c', 'c')
  tofspec.Draw()
  c.Print('TCN{0}_TCN{1}.pdf'.format(toftcn, reftcn)) # print to pdf
	  

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

ReadCycles(ROOT.TFile(sys.argv[1]), experiments) # read data from cycles

c = ROOT.TCanvas('c', 'c')
for ex in experiments: # draw all tof spectra
  ex['tofspectrum'].Draw()
  c.Print('TCN{0}.pdf'.format(ex['TCN']))

for cur in ['25A', '50A', '75A', '100A', '125A', '150A', '175A', '200A']: # normalize tof spectra of SCM to tof with 0 current in SCM
  NormalizeTOF(experiments, '18-065_{0}'.format(cur), '18-065_0A')
  NormalizeTOF(experiments, '18-265_{0}'.format(cur), '18-265_0A')
  
NormalizeTOF(experiments, '18-215', '18-380') # normalize tof spectra of foil transmission experiments
NormalizeTOF(experiments, '18-115', '18-480')
NormalizeTOF(experiments, '18-240', '18-302')
