import ROOT
import sys
from array import array
import numpy
import math
import itertools
import UCN

#Background Li6: 2.16 +/- 0.03
#           He3: 0.0403 +/- 0.0017

# analyze storage time from list of runs
def ReadCycles(infile, experiments):
  countperiodP1 = 1
  countperiodP2 = 2
  countperiodP3 = 3
  countperiodP4 = 4
  countperiodP5 = 4
  monitorperiod = 0
  backgroundperiod = 0
  storageperiod = 0

  for ex in experiments:
    ex['start'] = []
    ex['beamcurrent'] = []
    ex['li6counts'] = []
    ex['li6countsfull'] = []
    ex['countduration'] = []
    ex['monitorcounts'] = []
    ex['monitorduration'] = []
    ex['monitorcounts2'] = []
    ex['monitorduration2'] = []
    ex['li6background'] = []
    ex['backgroundduration'] = []
    ex['storageduration'] = []
    ex['beamcurrent'] = []
    ex['mintemperature'] = []
    ex['maxtemperature'] = []
    ex['minvaporpressure'] = []
    ex['maxvaporpressure'] = []
    ex['SCMcurrent'] = []
    ex['he3rate'] = []
    ex['li6rate'] = []
    ex['cycle'] =[]
    ex['N00'] = [] 
    ex['N01'] = []
    ex['N10'] = []
    ex['N11'] = []
    ex['monitorN00ir'] = [] 
    ex['monitorN01ir'] = []
    ex['monitorN10ir'] = []
    ex['monitorN11ir'] = []
    ex['monitorN00ct'] = [] 
    ex['monitorN01ct'] = []
    ex['monitorN10ct'] = []
    ex['monitorN11ct'] = []

  for cycle in infile.cycledata:
    run = cycle.runnumber
    if not any(run in ex['runs'] for ex in experiments): # if there is no experiment using this cycle
      continue
    Li6 = cycle.countsLi6
    He3 = cycle.countsHe3
    d = cycle.durations
    # filter useless cycles
    if min(cycle.B1V_KSM_PREDCUR) < 0.1:
      print('SKIPPING cycle {0} in run {1} because beam dropped below 0.1uA ({2}uA)'.format(cycle.cyclenumber, cycle.runnumber, min(cycle.B1V_KSM_PREDCUR)))
      continue
    if numpy.std(cycle.B1V_KSM_PREDCUR) > 0.2:
      print('SKIPPING cycle {0} in run {1} because beam fluctuated by {2}uA'.format(cycle.cyclenumber, cycle.runnumber, numpy.std(cycle.B1V_KSM_PREDCUR)))
      continue
    if Li6[10] == 0:
      print('SKIPPING cycle {0} in run {1} because Li6 does not contain data in all periods'.format(cycle.cyclenumber, cycle.runnumber))
      continue
    if He3[monitorperiod] < 1000:
      print('SKIPPING cycle {0} in run {1} because He3 saw less than 1000 monitor counts ({2})'.format(cycle.cyclenumber, cycle.runnumber, He3[monitorperiod]))
      continue
    if d[backgroundperiod] > 0 and Li6[backgroundperiod]/d[backgroundperiod] > 10:
      print('SKIPPING cycle {0} in run {1} because of high Li6 background ({2}/s)'.format(cycle.cyclenumber, cycle.runnumber, Li6[backgroundperiod]/d[backgroundperiod]))
      continue
    for ex in experiments:
      if run not in ex['runs']:
        continue
      ex['cycle'].append(cycle.cyclenumber)
      ex['start'].append(cycle.start)
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
      ex['li6counts'].append(Li6[countperiodP1]+Li6[countperiodP2]+Li6[countperiodP3]+Li6[countperiodP4])
      #print('li6counts are {0} for cycle {1}'.format(Li6[countperiodP1]+Li6[countperiodP2]+Li6[countperiodP3]+Li6[countperiodP4], cycle.cyclenumber))
      ex['countduration'].append(d[countperiodP1]+d[countperiodP2]+Li6[countperiodP3]+Li6[countperiodP4])
      ex['monitorcounts'].append(He3[monitorperiod])
      ex['monitorduration'].append(d[monitorperiod])
      ex['monitorcounts2'].append(He3[countperiodP1]+He3[countperiodP2]+He3[countperiodP3]+He3[countperiodP4])
      ex['monitorduration2'].append(d[countperiodP1]+d[countperiodP2]+Li6[countperiodP3]+Li6[countperiodP4])
      ex['li6background'].append(Li6[backgroundperiod])
      ex['backgroundduration'].append(d[backgroundperiod])
      ex['storageduration'].append(d[storageperiod])
      #print('Monitor Conuts is {0}').format(He3[monitorperiod],)
      li6rate = ROOT.TH1D('Le6_{0}_{1}'.format(cycle.runnumber, cycle.cyclenumber), 'Li6 detector rate', int(math.floor(sum(d))-d[monitorperiod]), d[monitorperiod], math.floor(sum(d)))
      for c in getattr(cycle, 'Li6/hits'):
        if cycle.cyclenumber%4==0 and c>60.:
          ex['N00'].append(c)
        if cycle.cyclenumber%4==1 and c>60.:
          ex['N10'].append(c)    
        if cycle.cyclenumber%4==2 and c>60.:
          ex['N01'].append(c)
        if cycle.cyclenumber%4==3 and c>60.:
          ex['N11'].append(c)
        if c>60.:
          li6rate.Fill(c)
      for c in getattr(cycle, 'He3/hits'):
        if c<60.:
          if cycle.cyclenumber%4==0:
            ex['monitorN00ir'].append(c)
          if cycle.cyclenumber%4==1:
            ex['monitorN10ir'].append(c)    
          if cycle.cyclenumber%4==2:
            ex['monitorN01ir'].append(c)
          if cycle.cyclenumber%4==3:
            ex['monitorN11ir'].append(c)
        else:
          if cycle.cyclenumber%4==0:
            ex['monitorN00ct'].append(c)
          if cycle.cyclenumber%4==1:
            ex['monitorN10ct'].append(c)    
          if cycle.cyclenumber%4==2:
            ex['monitorN01ct'].append(c)
          if cycle.cyclenumber%4==3:
            ex['monitorN11ct'].append(c)
      norm = 1./float(He3[monitorperiod])#li6rate.GetEntries()
      li6rate.Scale(norm)
      li6rate.GetXaxis().SetTitle('Time (s)')
      li6rate.GetYaxis().SetTitle('Li6 rate (s^{-1})')
      li6rate.SetDirectory(0)
      ex['li6rate'].append(li6rate)
	  
	  
# analyze storage time from list of runs
def PolarizationOvertime(ex):

  #Hardcoded background error
  bg1 =2.16 #n/s
  bge1 = 0.03

  print('\nAnalyzing TCN' + ex['TCN'])
  ex['li6backgroundrate'] = sum(ex['li6background'])/sum(ex['backgroundduration'])
  ex['li6backgroundrateerr'] = math.sqrt(sum(ex['li6background']))/sum(ex['backgroundduration'])
  print('Li6 detector background rate: {0} +/- {1} 1/s'.format(ex['li6backgroundrate'], ex['li6backgroundrateerr']))

  # report average monitor counts, range of beam current, range of He-II temperature
  monitoravg = numpy.average(ex['monitorcounts'], None, [1./m for m in ex['monitorcounts']], True)
  print 'Monitor counts: {0} +/- {1}'.format(monitoravg[0], 1./math.sqrt(monitoravg[1]))
  print('Beam current from {0} to {1} uA'.format(min(min(c) for c in ex['beamcurrent']), max(max(c) for c in ex['beamcurrent'])))
  #print 'Temperatures from {0} to {1} K'.format(min(ex['mintemperature']), max(ex['maxtemperature']))

  #seting up float variable arrays to plot
  x = [float(c) for c in ex['cycle']]
  xerr = [0. for cycle in ex['cycle']]
  y = [float(c) for c in ex['li6counts']]
  yerr = [math.sqrt(float(c)) for c in ex['li6counts']]
  
  # plot uncorrected ucn counts in Li6 det vs cycle
  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr)) 
  graph.SetTitle('TCN{0} Li6 counts per cycle'.format(ex['TCN']))
  graph.GetXaxis().SetTitle('Cycle Number')
  graph.GetYaxis().SetTitle('UCN-count unnormalized')
   
  #create the canvas and print
  canvas = ROOT.TCanvas('c', 'c')
  graph.Draw('AP')
  pdf = 'TCN{0}.pdf'.format(ex['TCN'])
  canvas.Print(pdf + '(')
  # Background uncorrected
  y = []
  yerr = []
  # plot background-uncorrected UCN counts correlated with monitor detector
  for cl,ch in zip(ex['li6counts'],ex['monitorcounts']):
    y.append(float(cl)/float(ch)/2)
    yerr.append(math.sqrt(float(cl)/float(ch)/float(ch)/4+float(cl)*float(cl)/float(ch)/float(ch)/float(ch)/4))

  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph.SetTitle('TCN{0} (UCN Ratio for monitor counts norm to Irradiation He3)'.format(ex['TCN']))
  graph.GetXaxis().SetTitle('Cycle Number')
  graph.GetYaxis().SetTitle('Li6 counts/He3 Counts')
  # Ratio background uncorrected
  #graph.Draw('AP')
  #canvas.Print(pdf)

  #Ratio of counts background corrected
  for cl,ch in zip(ex['li6counts'],ex['monitorcounts2']):
    y.append(float(cl)/float(ch)/2)
    Dcl = math.sqrt(cl/120/120+bge1*bge1)
    yerr.append(math.sqrt(float(Dcl)*float(Dcl)/float(ch)/float(ch)*60*60+float(cl)/120*float(cl)/120/float(ch)/float(ch)/float(ch)*60*60))
  graph = ROOT.TGraphErrors(len(x), numpy.array(x), numpy.array(y), numpy.array(xerr), numpy.array(yerr))
  graph.SetTitle('TCN{0} (UCN Ratio Li6 count/s / He3 count/s)'.format(ex['TCN']))
  graph.GetXaxis().SetTitle('Cycle Number')
  graph.GetYaxis().SetTitle('Li6 counts/He3 counts')
  graph.SetMarkerStyle(21)
  graph.SetMarkerSize(0.5)
  graph.Draw('AP')
  canvas.Print(pdf)

  #bins = [60., 62., 64., 67., 70., 74., 80., 100., 180.] 
  bins = [60., 62., 62.5, 63., 63.5, 64., 67., 70., 74., 80., 100., 180.] 
  #bins size is choosen to roughly get equal statistics in each bin
  times1 = [2. , 0.5, 0.5,  0.5, 0.5,  3.,  3.,  4.,  6., 20.,   80.]
  numbin = 11

  li6n00graph = ROOT.TH1D("Li6_N00","Li N00 Detector Rate",numbin,array('d',bins))
  li6n00graph.GetXaxis().SetTitle('Time (s)')
  li6n00graph.GetYaxis().SetTitle('Li6 rate (s^{-1})')
  li6n00graph.SetDirectory(0)
  for c in ex['N00']:
    li6n00graph.Fill(c)
  li6n00graph.Sumw2()
  li6n01graph = ROOT.TH1D("Li6_N01","Li N01 Detector Rate",numbin,array('d',bins))
  li6n01graph.GetXaxis().SetTitle('Time (s)')
  li6n01graph.GetYaxis().SetTitle('Li6 rate (s^{-1})')
  li6n01graph.SetDirectory(0)
  for c in ex['N01']:
    li6n01graph.Fill(c)
  li6n01graph.Sumw2()
  li6n10graph = ROOT.TH1D("Li6_N10","Li N00 Detector Rate",numbin,array('d',bins))
  li6n10graph.GetXaxis().SetTitle('Time (s)')
  li6n10graph.GetYaxis().SetTitle('Li6 rate (s^{-1})')
  li6n10graph.SetDirectory(0)
  for c in ex['N10']:
    li6n10graph.Fill(c)
  li6n10graph.Sumw2()
  li6n11graph = ROOT.TH1D('Li6_N11','Li N11 Detector Rate',numbin,array('d',bins))
  li6n11graph.GetXaxis().SetTitle('Time (s)')
  li6n11graph.GetYaxis().SetTitle('Li6 rate (s^{-1})')
  li6n11graph.SetDirectory(0)
  for c in ex['N11']:
     li6n11graph.Fill(c)
  li6n11graph.Sumw2()

  ratio00graph = ROOT.TH1D('Ratio_N00','Ratio N00/He3 Detector Ratio',numbin,array('d',bins))
  ratio00graph.GetXaxis().SetTitle('Time (s)')
  ratio00graph.GetYaxis().SetTitle('Li6counts/He3counts')

  ratio01graph = ROOT.TH1D('Ratio_N01','Ratio N01/He3 Detector Ratio',numbin,array('d',bins))
  ratio01graph.GetXaxis().SetTitle('Time (s)')
  ratio01graph.GetYaxis().SetTitle('Li6counts/He3counts')

  ratio10graph = ROOT.TH1D('Ratio_N10','Ratio N10/He3 Detector Ratio',numbin,array('d',bins))
  ratio10graph.GetXaxis().SetTitle('Time (s)')
  ratio10graph.GetYaxis().SetTitle('Li6counts/He3counts')

  ratio11graph = ROOT.TH1D('Ratio_N11','Ratio N11/He3 Detector Ratio',numbin,array('d',bins))
  ratio11graph.GetXaxis().SetTitle('Time (s)')
  ratio11graph.GetYaxis().SetTitle('Li6counts/He3counts')


  P0a = [0.100098677763,
         -0.000965492981703,
         0.0191408900653,
         0.016765387599,
         -0.0172769269617,
         0.00649448022435,
         -0.010230698995,
         0.00120924841799,
         -0.0028494230592,
         0.00876659415867,
         0.010301943088 ]
  P50a = [-0.0627870210314,
          0.0556699893438,
          -0.0011035711578,
          0.0310203835222,
          -0.0248401290515,
          -0.0221099547341,
          -0.0109209848292,
          -0.0230597455152,
          0.00683191708133,
          0.0137940601043,
          -0.00680863333668]
  P100a=[ 0.164848067484,
          0.193924651509,
          0.14962266601,
          0.142157216139,
          0.127758904671,
          0.129688938709,
          0.120240715505,
          0.110587446937,
          0.118071635723,
          0.138187128096,
          0.206303739708]
  P150a=[0.652875112075,
         0.660322663705,
         0.583391140926,
         0.515128786831,
         0.45191641782,
         0.47775630114,
         0.468276617196,
         0.444742604989,
         0.445052277317,
         0.442398381339,
         0.410539695132]
  P175a = [0.697503979791,
           0.710924736323,
           0.661081296954,
           0.573558593248,
           0.555969218126,
           0.54283430417,
           0.520198315635,
           0.452673223767,
           0.459892053948,
           0.447700201553,
           0.417064240845]
  P200a = [0.74380416843,
           0.744316702284,
           0.670309498774,
           0.611842838351,
           0.533765943973,
           0.536430067423,
           0.550625648042,
           0.472490859081,
           0.476603872211,
           0.451244650208,
           0.424388041012]
  #These are the counts per bin for the unquie bin structure of the array bins[] for the SCM to be included to generate a comparison graph of each value.  


  Pol200 = ROOT.TH1D('Pol_SCM','Polarization SCM 200A',numbin,array('d',bins))
  Pol200.GetXaxis().SetTitle('Time (s)')
  Pol200.GetYaxis().SetTitle('SCM Polarization')
  Pol200.SetDirectory(0)
  Pol200.SetLineColor(1)
  Pol200.SetLineWidth(3)
  Pol50 = ROOT.TH1D('Pol_SCM','Polarization SCM 50A',numbin,array('d',bins))
  Pol50.GetXaxis().SetTitle('Time (s)')
  Pol50.SetLineColor(2)
  Pol50.SetLineWidth(3)
  Pol50.GetYaxis().SetTitle('SCM Polarization')
  Pol100 = ROOT.TH1D('Pol_SCM','Polarization SCM 50A',numbin,array('d',bins))
  Pol100.GetXaxis().SetTitle('Time (s)')
  Pol100.GetYaxis().SetTitle('SCM Polarization')
  Pol100.SetLineColor(3)
  Pol100.SetLineWidth(3)
  Pol150 = ROOT.TH1D('Pol_SCM','Polarization SCM 150A',numbin,array('d',bins))
  Pol150.GetXaxis().SetTitle('Time (s)')
  Pol150.GetYaxis().SetTitle('SCM Polarization')
  Pol150.SetLineColor(4)
  Pol150.SetLineWidth(3)
  Pol175 = ROOT.TH1D('Pol_SCM','Polarization SCM 175A',numbin,array('d',bins))
  Pol175.GetXaxis().SetTitle('Time (s)')
  Pol175.GetYaxis().SetTitle('SCM Polarization')
  Pol175.SetLineColor(5)
  Pol175.SetLineWidth(3)
  Pol0 = ROOT.TH1D('Pol_SCM','Polarization SCM 0A',numbin,array('d',bins))
  Pol0.GetXaxis().SetTitle('Time (s)')
  Pol0.GetYaxis().SetTitle('SCM Polarization')
  Pol0.SetLineColor(6)
  Pol0.SetLineWidth(3)

  i=0.0
  polbin=1
  for scm in P200a:
    #print('entering loop {0}'.format(scm))
    if i==0:
      Pol200.SetBinContent(polbin,0.0)
      i=1.0
     # print("skip the first bin")
    else:  
      Pol200.SetBinContent(polbin,scm)
    polbin=polbin+1
 
  Pol200.Draw()
  i=0.0
  polbin=1
  for scm in P175a:
    #print('entering loop {0}'.format(scm))
    if i==0:
      Pol175.SetBinContent(polbin,0.0)
      i=1.0
     # print("skip the first bin")
    else:  
      Pol175.SetBinContent(polbin,scm)
    polbin=polbin+1
  Pol175.Draw("SAME")
  i=0.0
  polbin=1
  for scm in P150a:
    #print('entering loop {0}'.format(scm))
    if i==0:
      Pol150.SetBinContent(polbin,0.0)
      i=1.0
     # print("skip the first bin")
    else:  
      Pol150.SetBinContent(polbin,scm)
    polbin=polbin+1
  Pol150.Draw("SAME")
  i=0.0
  polbin=1
  for scm in P100a:
    #print('entering loop {0}'.format(scm))
    if i==0:
      Pol100.SetBinContent(polbin,0.0)
      i=1.0
     # print("skip the first bin")
    else:  
      Pol100.SetBinContent(polbin,scm)
    polbin=polbin+1
  Pol100.Draw("SAME")
  i=0.0
  polbin=1
  for scm in P50a:
    #print('entering loop {0}'.format(scm))
    if i==0:
      Pol50.SetBinContent(polbin,0.0)
      i=1.0
     # print("skip the first bin")
    else:  
      Pol50.SetBinContent(polbin,scm)
    polbin=polbin+1
  Pol50.Draw("SAME")
  i=0.0
  polbin=1
  for scm in P0a:
    #print('entering loop {0}'.format(scm))
    if i==0:
      Pol0.SetBinContent(polbin,0.0)
      i=1.0
     # print("skip the first bin")
    else:  
      Pol0.SetBinContent(polbin,scm)
    polbin=polbin+1
  Pol0.Draw("SAME")
  canvas.BuildLegend(0.9,.3,.9,.3,"","")
  canvas.Print(pdf)

  
  #Helium counts rate and the nubmer of cycles it was included
  He3N00 = 0. #15 cycles
  He3N10 = 0. #16 cycles
  He3N01 = 0. #15 cycles
  He3N11 = 0. #15 cycles

  for cycle,monitor in zip(ex['cycle'],ex['monitorcounts']):
    if cycle%4==0:
      He3N00 = He3N00 + float(monitor)
    if cycle%4==1:
      He3N10 = He3N10 + float(monitor)
    if cycle%4==2:
      He3N01 = He3N01 + float(monitor)
    if cycle%4==3:
      He3N11 = He3N11 + float(monitor)

 

  print("n00, n10, n01, n11 He counts  {0}, {1}, {2}, {3}".format(He3N00,He3N10,He3N01,He3N11) )

  #i = 1
  # Calculating a background corrected shifted version that accounts for the number of cycles for a experimment
  #for t in times1:
    #LI = (li6n00graph.GetBinContent(i))/15
    #print ("{0}, {1}, {2} ".format(math.sqrt(LI)/LI, math.sqrt(He3N00/15)*15/He3N00, bge1/bg1))
    #RelError = math.sqrt((LI/(15*15)+bge1**2*t**2)/abs(LI/(15)-bg1*t)**2+math.sqrt(15)/(He3N00/120))
    #LI = li6n00graph.GetBinContent(i)/15.-bg1*t
    #print ("counts {0}, scaled counts {1}, bg = {2}, Dt = {3}, BGcorrCounts = {4}".format(li6n00graph.GetBinContent(i),li6n00graph.GetBinContent(i)/15., bg1, t, li6n00graph.GetBinContent(i)/15.-bg1*t))
    #Ratio = LI/((He3N00/15)/120*t)
    #Error = Ratio*RelError
    #li6n00graph.SetBinContent(i, Ratio)
    #li6n00graph.SetBinError(i, Error)
    #error1  = Error
    
    #print('N00 ratio is {0}+/-{1} with rel {2}'.format(Ratio,error1,error1/Ratio))
    #numRun = 16 # 

    #LI = li6n10graph.GetBinContent(i)/numRun
    #RelError = math.sqrt((LI/(numRun**2)+bge1**2*t**2)/abs(LI/(numRun)-bg1*t)**2+math.sqrt(numRun)/(He3N10))
    #LI = li6n10graph.GetBinContent(i)/numRun-bg1*t
    #Ratio = LI/((He3N10/numRun)/120*t)
    #Error = Ratio*RelError
    #li6n10graph.SetBinContent(i,Ratio)
    #li6n10graph.SetBinError(i,Error)
    #error2  = Error

    #numRun = 15 # 
    #LI = li6n01graph.GetBinContent(i)/numRun
    #RelError = math.sqrt((LI/(numRun**2)+bge1**2*t**2)/abs(LI/(numRun)-bg1*t)**2+math.sqrt(numRun)/(He3N01))
    #LI = li6n01graph.GetBinContent(i)/numRun-bg1*t
    #Ratio = LI/((He3N01/numRun)/120*t)
    #Error = Ratio*RelError
    #li6n01graph.SetBinContent(i, Ratio)
    #li6n01graph.SetBinError(i, Error)
    #error3  = Error

    #LI = li6n11graph.GetBinContent(i)/15
    #RelError = math.sqrt((LI/(15*15)+bge1**2*t**2)/abs(LI/(15)-bg1*t)**2+math.sqrt(15)/(He3N11))
    #LI = li6n11graph.GetBinContent(i)/15-bg1*t
    #Ratio = LI/((He3N11/15)/120*t)
    #Error = Ratio*RelError
    #li6n11graph.SetBinContent(i, Ratio)
    #li6n11graph.SetBinError(i,Error)
    #error4 = Error
    
    #i=i+1


  polCheck = li6n00graph.Clone("polCheck")
  polCheck.SetTitle('polCheck Polarization Over Time (normalized)')
  polCheck.GetYaxis().SetTitle('Polarization Power, p')
  polCheck.GetXaxis().SetTitle('Time (s)')

  f = open("TCN{0}-polarization.txt".format(ex['TCN']),"w+")
  i = 1
  # Calculating a background corrected version that uses rates 
  # Note: there are 15 runs in some cases and 16 in others. The background needs to be subtracted that many times. 
  for t in times1:
    LI = li6n00graph.GetBinContent(i)
    Nc = LI-bg1*t*15
    ebg = bge1/15 
    eNc = math.sqrt(math.sqrt(LI**2)+(ebg*t)**2)
    HEc = (He3N00/120)*t #scaling the counts to bin size
    eHEc = (math.sqrt(He3N00)/120)*t #scaling the error to bin size 
    Rc = Nc/HEc
    eRc = math.sqrt((eNc/HEc)*(eNc/HEc)+(eHEc*Nc/(HEc*HEc))*(eHEc*Nc/(HEc*HEc)))
    ratio00graph.SetBinContent(i, Rc)
    ratio00graph.SetBinError(i, eRc)
  
    LI = li6n01graph.GetBinContent(i)
    Nc = LI-bg1*t*16
    ebg = bge1/16
    eNc = math.sqrt(math.sqrt(LI**2)+(ebg*t)**2)
    HEc = (He3N01/120)*t #scaling the counts to bin size
    eHEc = (math.sqrt(He3N01)/120)*t #scaling the error to bin size 
    Rc = Nc/HEc
    eRc = math.sqrt((eNc/HEc)*(eNc/HEc)+(eHEc*Nc/(HEc*HEc))*(eHEc*Nc/(HEc*HEc)))
    ratio01graph.SetBinContent(i, Rc)
    ratio01graph.SetBinError(i, eRc)

    LI = li6n10graph.GetBinContent(i)
    Nc = LI-bg1*t*15
    ebg = bge1/15 
    eNc = math.sqrt(math.sqrt(LI**2)+(ebg*t)**2)
    HEc = (He3N10/120)*t #scaling the counts to bin size
    eHEc = (math.sqrt(He3N10)/120)*t #scaling the error to bin size 
    Rc = Nc/HEc
    eRc = math.sqrt((eNc/HEc)*(eNc/HEc)+(eHEc*Nc/(HEc*HEc))*(eHEc*Nc/(HEc*HEc)))
    ratio10graph.SetBinContent(i, Rc)
    ratio10graph.SetBinError(i, eRc)

    LI = li6n11graph.GetBinContent(i)
    Nc = LI-bg1*t*15
    ebg = bge1/15 
    eNc = math.sqrt(math.sqrt(LI**2)+(ebg*t)**2)
    HEc = (He3N11/120)*t #scaling the counts to bin size
    eHEc = (math.sqrt(He3N11)/120)*t #scaling the error to bin size 
    Rc = Nc/HEc
    eRc = math.sqrt((eNc/HEc)*(eNc/HEc)+(eHEc*Nc/(HEc*HEc))*(eHEc*Nc/(HEc*HEc)))
    ratio11graph.SetBinContent(i, Rc)
    ratio11graph.SetBinError(i, eRc)
  
    gN00 = ratio00graph.GetBinContent(i)
    egN00 = ratio00graph.GetBinError(i)
    gN01 = ratio01graph.GetBinContent(i)
    egN01 = ratio01graph.GetBinError(i)
    gN10 = ratio10graph.GetBinContent(i)
    egN10 = ratio10graph.GetBinError(i)
    gN11 = ratio11graph.GetBinContent(i)
    egN11 = ratio11graph.GetBinError(i)

    
    eA = math.sqrt(2*egN01)
    eB = math.sqrt(egN00*egN00+egN10*egN10)
    eC = abs(gN11*gN00)*math.sqrt(egN11/gN11*egN11/gN11+egN00/gN00*egN00/gN00)
    #p = B^2/(C-A^2)  : B = n00-n01, A = n01, C = n11*n00
    eCA = math.sqrt(eA*eA+eC*eC)
    B = (gN00-gN10)
    BB = B*B
    CA = (gN11*gN00-gN01*gN01)
    esCA = math.sqrt(1/2)*eCA
    sCA = math.sqrt(CA)
    P = BB/CA
    P = math.sqrt(P)
    eP = abs (P)*math.sqrt(eB/B*eB/B+esCA/sCA*esCA/sCA)

    #eNPC = math.sqrt(egN00*egN00+egN10*egN10)
    #eDPC = math.sqrt(egN00*egN00+egN10*egN10)
    #NPC = gN00-gN01
    #DPC = gN00+gN01
    #PC = NPC/DPC
    #ePC = abs(PC)*math.sqrt(eNPC/NPC*eNPC/NPC+eDPC/DPC*eDPC/DPC)
    #print("{0} is the pol").format(PC)

    polCheck.SetBinContent(i,P)
    polCheck.SetBinError(i,eP)

    print("polarization efficiency is {0} +/- {1} for bin {2}").format( polCheck.GetBinContent(i),polCheck.GetBinError(i), i)
    
    f.write("{0} {1}\n".format(polCheck.GetBinContent(i),polCheck.GetBinError(i)))

    i=i+1
  f.close()
  polCheck.Draw('e')
  canvas.Print(pdf)

  numerator2 = ratio00graph.Clone("numerator2")
  numerator2.Add(ratio10graph, -1.)
  numerator2.Multiply(numerator2)
  denom2 = ratio01graph.Clone("denom2")
  denom2.Multiply(denom2)
  denom3 = ratio00graph.Clone("denom3")
  denom3.Multiply(ratio11graph)
  denom3.Add(denom2, -1)
  numerator2.Divide(denom3)
  numerator2.SetTitle('Polarization Over Time (normalized)')
  numerator2.GetYaxis().SetTitle('Polarization Power, p')
  numerator2.GetXaxis().SetTitle('Time (s)')

  for i in range(1,12):
    yvalue=numerator2.GetBinContent(i)
    yerror=numerator2.GetBinError(i)
    if i == 1:
      numerator2.SetBinContent(i,0.)
      numerator2.SetBinError(i,0.)
    else:
      numerator2.SetBinContent(i, math.sqrt(abs(yvalue)))
      numerator2.SetBinError(i,1./2.*yerror/yvalue)

    
      
  numerator2.Draw()
  canvas.Print(pdf)
  #Ratio of Li6/He3 normailized via seciton 3.1 in Fall run 2018 - analysis report https://github.com/ucn-triumf/UCNanalysis2018/blob/master/report/report.pdf

 
  kRed = 632
  kBlue = 600
  kYellow = 400
  kGreen = 416

  ROOT.gStyle.SetOptTitle(0)


  #againt this is not weight by the correct number of counts and cycles for correct comparison
  li6n00g = ROOT.TH1D("Li6_N00_rate","Li N00 Detector Rate",120,60.,180.)
  li6n00g.GetXaxis().SetTitle('Time (s)')
  li6n00g.GetYaxis().SetTitle('UCN Counts Li6 (Neutrons)')
  li6n00graph.SetDirectory(0)
  for c in ex['N00']:
    li6n00g.Fill(c)
  li6n00g.Sumw2()
  li6n00g.SetLineWidth(2)
  li6n00g.SetLineStyle(1)
  li6n00g.SetLineColor(kRed+2)
  li6n00g.SetMarkerStyle(21)
  li6n00g.SetMarkerSize(0.5)
  li6n00g.SetMarkerColor(kRed)
  li6n00g.Draw("lep")
  #canvas.Print(pdf)
  li6n10g = ROOT.TH1D("Li6_N10_rate","Li N10 Detector Rate",120,60.,180.)
  li6n10g.GetXaxis().SetTitle('Time (s)')
  li6n10g.GetYaxis().SetTitle('Li6 rate (s^{-1})')
  li6n10graph.SetDirectory(0)
  for c in ex['N10']:
    li6n10g.Fill(c)
  li6n10g.Sumw2()
  li6n10g.SetLineWidth(2)
  li6n10g.SetLineStyle(1)
  li6n10g.SetLineColor(kBlue+2)
  li6n10g.SetMarkerStyle(22)
  li6n10g.SetMarkerSize(0.5)
  li6n10g.SetMarkerColor(kBlue)
  li6n10g.Draw("lepSAME")
  #canvas.Print(pdf)
  li6n01g = ROOT.TH1D("Li6_N01_rate","Li N01 Detector Rate",120,60.,180.)
  li6n01g.GetXaxis().SetTitle('Time (s)')
  li6n01g.GetYaxis().SetTitle('Li6 rate (s^{-1})')
  li6n01g.SetDirectory(0)
  for c in ex['N01']:
    li6n01g.Fill(c)
  li6n01g.Sumw2()
  li6n01g.SetName("li6n01g")
  li6n01g.SetLineWidth(2)
  li6n01g.SetLineStyle(1)
  li6n01g.SetLineColor(kGreen+2)
  li6n01g.SetMarkerStyle(23)
  li6n01g.SetMarkerSize(0.5)
  li6n01g.SetMarkerColor(kGreen)
  li6n01g.Draw("lepSAME")
  #canvas.Print(pdf)
  li6n11g = ROOT.TH1D("Li6_N11_rate","Li N11 Detector Rate",120,60.,180.)
  li6n11g.GetXaxis().SetTitle('Time (s)')
  li6n11g.GetYaxis().SetTitle('Li6 rate (s^{-1})')
  li6n11g.SetDirectory(0)
  for c in ex['N11']:
    li6n11g.Fill(c)
  li6n11g.Sumw2()
  li6n11g.SetLineWidth(2)
  li6n11g.SetLineStyle(1)
  li6n11g.SetLineColor(kYellow+2)
  li6n11g.SetMarkerStyle(20)
  li6n11g.SetMarkerSize(0.5)
  li6n11g.SetMarkerColor(kYellow)
  li6n11g.Draw("SAME")
  
  canvas.BuildLegend(0.9,.3,.9,.3,"","")
  canvas.Print(pdf)


  ## Here we use the previous graphs and manipulate only the bins in interval of 1 second

  li6n00ge = ROOT.TH1D("Li6_N00_rate_e","Li N00 Detector Rate ",120,60.,180.)
  li6n00ge.GetXaxis().SetTitle('Time (s)')
  li6n00ge.GetYaxis().SetTitle('UCN Count in Li6 (counts)')
  li6n00ge.SetDirectory(0)
  for b in range(120):
    li6n00ge.SetBinContent(b, li6n00g.GetBinContent(b)-bg1)
    li6n00ge.SetBinError(b, math.sqrt(li6n00g.GetBinContent(b)+bge1*bge1))
  li6n00ge.SetLineWidth(2)
  li6n00ge.SetLineStyle(1)
  li6n00ge.SetLineColor(kRed+2)
  li6n00ge.SetMarkerStyle(21)
  li6n00ge.SetMarkerSize(0.5)
  li6n00ge.SetMarkerColor(kRed)
  li6n00ge.Draw("lep")
    #canvas.Print(pdf)
  li6n10ge = ROOT.TH1D("Li6_N10_rate_e","Li N10 Detector Rate ",120,60.,180.)
  li6n10ge.GetXaxis().SetTitle('Time (s)')
  li6n10ge.GetYaxis().SetTitle('Li6 rate (s^{-1})')
  #li6n10ge.SetDirectory(0)
  for b in range(120):
    li6n10ge.SetBinContent(b, li6n10g.GetBinContent(b)-bg1)
    li6n10ge.SetBinError(b, math.sqrt(li6n10g.GetBinContent(b)+bge1*bge1))
  li6n10ge.SetLineWidth(2)
  li6n10ge.SetLineStyle(1)
  li6n10ge.SetLineColor(kBlue+2)
  li6n10ge.SetMarkerStyle(21)
  li6n10ge.SetMarkerSize(0.5)
  li6n10ge.SetMarkerColor(kBlue)
  li6n10ge.Draw("lepSAME")
  li6n01ge = ROOT.TH1D("Li6_N01_rate_e","Li N01 Detector Rate ",120,60.,180.)
  li6n01ge.GetXaxis().SetTitle('Time (s)')
  li6n01ge.GetYaxis().SetTitle('Li6 rate (s^{-1})')
  #li6n01ge.SetDirectory(0)
  for b in range(120):
    li6n01ge.SetBinContent(b, li6n01g.GetBinContent(b)-bg1)
    li6n01ge.SetBinError(b, math.sqrt(li6n01g.GetBinContent(b)+bge1*bge1))
  li6n01ge.SetLineWidth(2)
  li6n01ge.SetLineStyle(1)
  li6n01ge.SetLineColor(kGreen+2)
  li6n01ge.SetMarkerStyle(21)
  li6n01ge.SetMarkerSize(0.5)
  li6n01ge.SetMarkerColor(kGreen)
  li6n01ge.Draw("lepSAME")
  li6n11ge = ROOT.TH1D("Li6_N11_rate_e","Li N11 Detector Rate ",120,60.,180.)
  li6n11ge.GetXaxis().SetTitle('Time (s)')
  li6n11ge.GetYaxis().SetTitle('Li6 rate (s^{-1})')
  #li6n11ge.SetDirectory(0)
  for b in range(120):
    li6n11ge.SetBinContent(b, li6n11g.GetBinContent(b)-bg1)
    li6n11ge.SetBinError(b, math.sqrt(li6n11g.GetBinContent(b)+bge1*bge1))
  li6n11ge.SetLineWidth(2)
  li6n11ge.SetLineStyle(1)
  li6n11ge.SetLineColor(kYellow+2)
  li6n11ge.SetMarkerStyle(21)
  li6n11ge.SetMarkerSize(0.5)
  li6n11ge.SetMarkerColor(kYellow)
  li6n11ge.Draw("lepSAME")

  canvas.BuildLegend(0.9,.3,.9,.3,"","")
  canvas.Print(pdf)

  He3I00g = ROOT.TH1D("He3_N00_irradiation_rate","He N00 Detector Rate irr",60,0.,60.)
  He3I00g.GetXaxis().SetTitle('Time (s)')
  He3I00g.GetYaxis().SetTitle('He3Det neutron rate (s^{-1})')
  He3I00g.SetDirectory(0)
  for c in ex['monitorN00ir']:
    He3I00g.Fill(c)
  He3I00g.Sumw2()
  He3I00g.Draw()
  #canvas.Print(pdf)
  He3I10g = ROOT.TH1D("He3_N10_irradiation_rate","He N10 Detector Rate irr",60,0.,60.)
  He3I10g.GetXaxis().SetTitle('Time (s)')
  He3I10g.GetYaxis().SetTitle('He3Det neutron rate (s^{-1})')
  #He3I10g.SetDirectory(0)
  for c in ex['monitorN10ir']:
    He3I10g.Fill(c)
  He3I10g.Sumw2()
  He3I10g.Draw("same")
  #canvas.Print(pdf)
  He3I01g = ROOT.TH1D("He3_N01_irradiation_rate","He N01 Detector Rate irr",60,0.,60.)
  He3I01g.GetXaxis().SetTitle('Time (s)')
  He3I01g.GetYaxis().SetTitle('He3Det neutron rate (s^{-1})')
  #He3I01g.SetDirectory(0)
  for c in ex['monitorN01ir']:
    He3I01g.Fill(c)
  He3I01g.Sumw2()
  He3I01g.Draw("same")
  #canvas.Print(pdf)
  He3I11g = ROOT.TH1D("He3_N11_irradiation_rate","He N11 Detector Rate irr",60,0.,60.)
  He3I11g.GetXaxis().SetTitle('Time (s)')
  He3I11g.GetYaxis().SetTitle('He3Det neutron rate (s^{-1})')
  #He3I11g.SetDirectory(0)
  for c in ex['monitorN11ir']:
    He3I11g.Fill(c)
  He3I11g.Sumw2()
  He3I11g.Draw("same")
  canvas.Print(pdf)

  He3C00g = ROOT.TH1D("He3_N00_count_period_rate","He N00 Detector Rate ct",120,60.,180.)
  He3C00g.GetXaxis().SetTitle('Time (s)')
  He3C00g.GetYaxis().SetTitle('He3Det neutron rate (s^{-1})')
  He3C00g.SetDirectory(0)
  for c in ex['monitorN00ct']:
    He3C00g.Fill(c)
  He3C00g.Sumw2()
  He3C00g.Draw()
  He3C10g = ROOT.TH1D("He3_N10_count_period_rate","He N10 Detector Rate ct",120,60.,180.)
  He3C10g.GetXaxis().SetTitle('Time (s)')
  He3C10g.GetYaxis().SetTitle('He3Det neutron rate (s^{-1})')
  #He3C10g.SetDirectory(0)
  for c in ex['monitorN10ct']:
    He3C10g.Fill(c)
  He3C10g.Sumw2()
  He3C10g.Draw("same")
  #canvas.Print(pdf)
  He3C01g = ROOT.TH1D("He3_N01_count_period_rate","He N01 Detector Rate ct",120,60.,180.)
  He3C01g.GetXaxis().SetTitle('Time (s)')
  He3C01g.GetYaxis().SetTitle('He3Det neutron rate (s^{-1})')
  #He3C01g.SetDirectory(0)
  for c in ex['monitorN01ct']:
    He3C01g.Fill(c)
  He3C01g.Sumw2()
  He3C01g.Draw("same")
  #canvas.Print(pdf)
  He3C11g = ROOT.TH1D("He3_N11_count_period_rate","He N11 Detector Rate ct",120,60.,180.)
  He3C11g.GetXaxis().SetTitle('Time (s)')
  He3C11g.GetYaxis().SetTitle('He3Det neutron rate (s^{-1})')
  #He3C11g.SetDirectory(0)
  for c in ex['monitorN11ct']:
    He3C11g.Fill(c)
  He3C11g.Sumw2()
  He3C11g.Draw("same")
  canvas.Print(pdf)

  RatioN00ge = ROOT.TH1D("Ratio_n00_rate_e","Detector Ratio e",120,60.,180.)
  RatioN00ge.GetXaxis().SetTitle('Time (s)')
  RatioN00ge.GetYaxis().SetTitle('Ration of UCN Counts Li6/He3')
  RatioN00ge.SetDirectory(0)
  for nbin in range(120):
    RatioN00ge.SetBinContent(nbin,li6n00ge.GetBinContent(nbin)/(He3I00g.Integral(0,120)/120))
    RatioN00ge.SetBinError(nbin, math.sqrt( (li6n00ge.GetBinError(nbin)/li6n00ge.GetBinContent(nbin))**2 + ((li6n00ge.GetBinContent(nbin))**2)/((He3I00g.Integral(0,120)/120)**3) ))
  RatioN00ge.Draw()
  RatioN10ge = ROOT.TH1D("Ratio_n10_rate_e","Detector Ratio n10 e",120,60.,180.)
  for nbin in range(120):
    RatioN10ge.SetBinContent(nbin,li6n10ge.GetBinContent(nbin)/(He3I10g.Integral(0,120)/120))
    RatioN10ge.SetBinError(nbin, math.sqrt( (li6n10ge.GetBinError(nbin)/li6n10ge.GetBinContent(nbin))**2 + ((li6n10ge.GetBinContent(nbin))**2)/((He3I10g.Integral(0,120)/120)**3) ))
  RatioN10ge.Draw('same')
  RatioN01ge = ROOT.TH1D("Ratio_n01_rate_e","Detector Ratio n01 e",120,60.,180.)
  for nbin in range(120):
    RatioN01ge.SetBinContent(nbin,li6n01ge.GetBinContent(nbin)/(He3I01g.Integral(0,120)/120))
    RatioN01ge.SetBinError(nbin, math.sqrt( (li6n01ge.GetBinError(nbin)/li6n01ge.GetBinContent(nbin))**2 + ((li6n01ge.GetBinContent(nbin))**2)/((He3I01g.Integral(0,120)/120)**3) ))
  RatioN01ge.Draw('same')
  RatioN11ge = ROOT.TH1D("Ratio_n11_rate_e","Detector Ratio n11 e",120,60.,180.)
  for nbin in range(120):
    RatioN11ge.SetBinContent(nbin,li6n11ge.GetBinContent(nbin)/(He3I11g.Integral(0,120)/120))
    RatioN11ge.SetBinError(nbin, math.sqrt( (li6n11ge.GetBinError(nbin)/li6n11ge.GetBinContent(nbin))**2 + ((li6n11ge.GetBinContent(nbin))**2)/((He3I11g.Integral(0,120)/120)**3) ))
  RatioN11ge.Draw('same')

  
  gPolNorm = ROOT.TH1D("Polarization n00 and n10","",120,60.,180.)
  gPolNorm.GetXaxis().SetTitle('time (s)')
  gPolNorm.GetYaxis().SetTitle('polarization')
  gPolNorm2 = ROOT.TH1D("Polarization n00 and n01","",120,60.,180.)
  gPolNorm2.GetXaxis().SetTitle('time (s)')
  gPolNorm2.GetYaxis().SetTitle('polarization')
  gPolNorm2.SetLineColor(2)
  gPolNorm3 = ROOT.TH1D("Polarization n11 and n10","",120,60.,180.)
  gPolNorm3.GetXaxis().SetTitle('time (s)')
  gPolNorm3.GetYaxis().SetTitle('polarization')
  gPolNorm3.SetLineColor(3)
  gPolNorm4 = ROOT.TH1D("Polarization n11 and n01","",120,60.,180.)
  gPolNorm4.GetXaxis().SetTitle('time (s)')
  gPolNorm4.GetYaxis().SetTitle('polarization')
  gPolNorm4.SetLineColor(4)
  gPolNorm5 = ROOT.TH1D("Polarization method f1=f2","",120,60.,180.)
  gPolNorm5.GetXaxis().SetTitle('time (s)')
  gPolNorm6 = ROOT.TH1D("Polarization method f1=f2 v2","",120,60.,180.)
  gPolNorm6.GetXaxis().SetTitle('time (s)')

  for nbin in range(120):
    numN00=RatioN00ge.GetBinContent(nbin)
    errN00=RatioN00ge.GetBinError(nbin)
    numN10=RatioN10ge.GetBinContent(nbin)
    errN10=RatioN10ge.GetBinError(nbin)
    numN01=RatioN01ge.GetBinContent(nbin)
    errN01=RatioN01ge.GetBinError(nbin)
    numN11=RatioN11ge.GetBinContent(nbin)
    errN11=RatioN11ge.GetBinError(nbin)

    if ((numN00-numN10)/(numN00+numN10))<0.:
      pol1 = -1*math.sqrt(abs((numN00-numN10)/(numN00+numN10)) )
    else: pol1 = math.sqrt(abs((numN00-numN10)/(numN00+numN10)) )
    gPolNorm.SetBinContent(nbin, pol1)
   # gPolNorm.SetBinError(nbin, math.sqrt( ((numN10/((math.sqrt(abs(pol1))*(numN10+numN00))**2))**2)*(errN00)**2+(((numN00/((math.sqrt(abs(pol1))*(numN10+numN00))**2))**2)*(errN10)**2 ) ))
    if ((numN00-numN01)/(numN00+numN01))<0.:
      pol2 = -1*math.sqrt(abs((numN00-numN01)/(numN00+numN01)) )
    else: pol2 = math.sqrt(abs((numN00-numN01)/(numN00+numN01)) )
    gPolNorm2.SetBinContent(nbin, pol2)
   # gPolNorm2.SetBinError(nbin, math.sqrt( ((numN01/((math.sqrt(abs(pol2))*(numN01+numN00))**2))**2)*(errN00)**2+(((numN00/((math.sqrt(abs(pol2))*(numN01+numN00))**2))**2)*(errN01)**2 ) ))
    if ((numN11-numN10)/(numN11+numN10))<0.:
      pol3 = -1*math.sqrt(abs((numN11-numN10)/(numN11+numN10)) )
    else: pol3 = math.sqrt(abs((numN11-numN10)/(numN11+numN10)) )
    gPolNorm3.SetBinContent(nbin, pol3)
   # gPolNorm3.SetBinError(nbin, math.sqrt( ((numN10/((math.sqrt(abs(pol3))*(numN10+numN11))**2))**2)*(errN11)**2+(((numN11/((math.sqrt(abs(pol3))*(numN10+numN11)**2)))**2)*(errN10)**2 ) ) )
    if ((numN11-numN01)/(numN00+numN01))<0.:
      pol4 = -1*math.sqrt(abs((numN11-numN01)/(numN11+numN01)) )
    else: pol4 = math.sqrt(abs((numN11-numN01)/(numN11+numN01)) )
    gPolNorm4.SetBinContent(nbin, pol4)
    gPolNorm4.SetBinError(nbin, math.sqrt( ((numN01/((math.sqrt(abs(pol4))*(numN01+numN11))**2))**2)*(errN11)**2+(((numN11/((math.sqrt(abs(pol4))*(numN01+numN11)**2)))**2)*(errN01)**2 ) ) )
    

    if ((numN00*numN11-numN10**2))<0.:
      pol5 = -1*abs(numN00-numN10)/math.sqrt(abs(numN00*numN11-numN10**2))
    else: pol5 = abs(numN00-numN10)/math.sqrt(abs(numN00*numN11-numN10**2))
    gPolNorm5.SetBinContent(nbin, pol5)
    gPolNorm5.SetBinError(nbin, 0 )
    if ((numN00*numN11-numN01**2))<0.:
      pol6 = -1*abs(numN00-numN01)/math.sqrt(abs(numN00*numN11-numN01**2))
    else: pol6 = abs(numN00-numN01)/math.sqrt(abs(numN00*numN11-numN01**2))
    gPolNorm6.SetBinContent(nbin, pol6)
    gPolNorm6.SetBinError(nbin, 0 )

  
  gPolNorm.SetDirectory(0)
  gPolNorm.Draw()
  gPolNorm2.Draw('same')
  gPolNorm3.Draw('same')
  gPolNorm4.Draw('same')
  canvas.Print(pdf)
  gPolNorm5.SetDirectory(0)
  gPolNorm5.Draw()
  gPolNorm6.Draw('same')
  #canvas.Print(pdf)


  

  canvas.Print(pdf+ ')')


#ROOT.gStyle.SetOptStat(1001111)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1

# list runs for each experiment
#experiments = [{'TCN': '18-180-v3', 'runs': [1029]}]
#experiments = [{'TCN': '18-070', 'runs': [1029]}]
# 1147, 1148,  1149, 1150
experiments = [{'TCN': '18-280-v3', 'runs': [1150]}]

ReadCycles(ROOT.TFile(sys.argv[1]), experiments)
			  
# loop over experiments
for ex in experiments:
 PolarizationOvertime(ex)



