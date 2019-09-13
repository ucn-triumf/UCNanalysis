import ROOT
import sys
from array import array
import numpy
import math
import itertools
import UCN


def Polarization(infile1,P,eP,graphsID,graph):
   
    SCM_A = []
    SCM_eA = []  

    f = open(infile1,"r")
    for line in f:
        columns = line.split()
        SCM_A.append(float(columns[0]))
        SCM_eA.append(float(columns[1]))
    f.close()

    j=0
    while j<numbin:
        A = SCM_A[j]
        eA= SCM_eA[j]
        p_f = P[j]
        ep_f = eP[j]
        
        #A = p_scm*p_f
        P_scm = A/p_f
                  
        eP_scm  = abs(P_scm)*math.sqrt(((eA*eA)/(A*A))+((ep_f*ep_f)/(p_f*p_f)))
        if graph == "TCN18-180":
            if graphsID ==1:
                gSCM180_0.SetBinContent(j+1, P_scm)
                gSCM180_0.SetBinError(j+1, eP_scm)
            if graphsID ==2:
                gSCM180_50.SetBinContent(j+1, P_scm)
                gSCM180_50.SetBinError(j+1, eP_scm)
            if graphsID ==3:
                gSCM180_100.SetBinContent(j+1, P_scm)
                gSCM180_100.SetBinError(j+1, eP_scm)  
            if graphsID ==4:
                gSCM180_150.SetBinContent(j+1, P_scm)
                gSCM180_150.SetBinError(j+1, eP_scm)
            if graphsID ==5:
                gSCM180_175.SetBinContent(j+1, P_scm)
                gSCM180_175.SetBinError(j+1, eP_scm)
            if graphsID ==6:
                gSCM180_200.SetBinContent(j+1, P_scm)
                gSCM180_200.SetBinError(j+1, eP_scm)
        if graph == "TCN18-280":
            if graphsID ==1:
                gSCM280_0.SetBinContent(j+1, P_scm)
                gSCM280_0.SetBinError(j+1, eP_scm)
            if graphsID ==2:
                gSCM280_50.SetBinContent(j+1, P_scm)
                gSCM280_50.SetBinError(j+1, eP_scm)
            if graphsID ==3:
                gSCM280_100.SetBinContent(j+1, P_scm)
                gSCM280_100.SetBinError(j+1, eP_scm)  
            if graphsID ==4:
                gSCM280_150.SetBinContent(j+1, P_scm)
                gSCM280_150.SetBinError(j+1, eP_scm)
            if graphsID ==5:
                gSCM280_175.SetBinContent(j+1, P_scm)
                gSCM280_175.SetBinError(j+1, eP_scm)
            if graphsID ==6:
                gSCM280_200.SetBinContent(j+1, P_scm)
                gSCM280_200.SetBinError(j+1, eP_scm)
        if graph == "p_avg":
            if graphsID ==1:
                gSCMavg_0.SetBinContent(j+1, P_scm)
                gSCMavg_0.SetBinError(j+1, eP_scm)
            if graphsID ==2:
                gSCMavg_50.SetBinContent(j+1, P_scm)
                gSCMavg_50.SetBinError(j+1, eP_scm)
            if graphsID ==3:
                gSCMavg_100.SetBinContent(j+1, P_scm)
                gSCMavg_100.SetBinError(j+1, eP_scm)  
            if graphsID ==4:
                gSCMavg_150.SetBinContent(j+1, P_scm)
                gSCMavg_150.SetBinError(j+1, eP_scm)
            if graphsID ==5:
                gSCMavg_175.SetBinContent(j+1, P_scm)
                gSCMavg_175.SetBinError(j+1, eP_scm)
            if graphsID ==6:
                gSCMavg_200.SetBinContent(j+1, P_scm)
                gSCMavg_200.SetBinError(j+1, eP_scm)

        #print("{0} +/- {1} avg for the {2}-th bin").format(avgP[i],avgeP[i],i)
        j=j+1
    
    

bins = [60., 62., 62.5, 63., 63.5, 64., 67., 70., 74., 80., 100., 180.] 
numbin = 11

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gROOT.SetBatch(1)
ROOT.gErrorIgnoreLevel = ROOT.kInfo + 1

polA180 = []
poleA180 = []
polA280 = []
poleA280 = []
exSCM= ['TCN18-070-0A-polarization.txt',
        'TCN18-070-50A-polarization.txt',
        'TCN18-070-100A-polarization.txt',
        'TCN18-070-150A-polarization.txt',
        'TCN18-070-175A-polarization.txt',
        'TCN18-070-200A-polarization.txt']

avgP = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
avgeP = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]

canvas = ROOT.TCanvas('c', 'c')

kRed = 632
kYellow = 400
kGreen = 416
kCyan = 432
kBlue = 600
kGreen = 416
kPurple = 616

file180 = 'TCN18-180-v3-polarization.txt'
file280 = 'TCN18-280-v3-polarization.txt'

P180graph = ROOT.TH1D("P180","Polarization over Time, TCN18-180",numbin,array('d',bins))
P180graph.GetXaxis().SetTitle('Time (s)')
P180graph.GetYaxis().SetTitle('Polarizing Power, P_{foil}')

P280graph = ROOT.TH1D("P280","Polarization over Time, TCN18-280",numbin,array('d',bins))
P280graph.GetXaxis().SetTitle('Time (s)')
P280graph.GetYaxis().SetTitle('Polarizing Power, P_{foil}')

avgPgraph = ROOT.TH1D("avgP","Polarization over Time, mean of TCN18-180 and -280",numbin,array('d',bins))
avgPgraph.GetXaxis().SetTitle('Time (s)')
avgPgraph.GetYaxis().SetTitle('Polarizing Power, P_{foil}')

f = open(file180,'r')
for line in f:
    columns = line.split()
    polA180.append(float(columns[0]))
    poleA180.append(float(columns[1]))
f.close()

f = open(file280,'r')
for line in f:
    columns = line.split()
    polA280.append(float(columns[0]))
    poleA280.append(float(columns[1]))
f.close()

i=0

pdf = 'Polarization.pdf'

while i<numbin:
    p1 = polA180[i]
    ep1= poleA180[i]
    p2 = polA280[i]
    ep2 = poleA280[i]
    
    avgP[i]=(p1+p2)/2.
    avgeP[i]= 1./2.*math.sqrt(ep1*ep1+ep2*ep2)
    #print("{0} +/- {1} avg for the {2}-th bin").format(p1,ep1,i)
    #print("{0} +/- {1} avg for the {2}-th bin").format(p2,ep2,i)
    #print("{0} +/- {1} avg for the {2}-th bin").format(avgP[i],avgeP[i],i)
    
    P180graph.SetBinContent(i+1, p1)
    P180graph.SetBinError(i+1, ep1)
    P280graph.SetBinContent(i+1, p2)
    P280graph.SetBinError(i+1, ep2)
    avgPgraph.SetBinContent(i+1, avgP[i])
    avgPgraph.SetBinError(i+1, avgeP[i])

    i=i+1

#graph of polarization
P180graph.SetDirectory(0)
P180graph.SetLineWidth(2)
P180graph.SetLineStyle(1)
P180graph.SetLineColor(kBlue+2)
P180graph.SetMarkerStyle(21)
P180graph.SetMarkerSize(0.5)
P180graph.SetMarkerColor(kBlue)
P180graph.Draw()
P180graph.SetMinimum(0.)
P180graph.SetMaximum(1.)


P280graph.SetLineWidth(2)
P280graph.SetLineStyle(1)
P280graph.SetLineColor(kRed+2)
P280graph.SetMarkerStyle(21)
P280graph.SetMarkerSize(0.5)
P280graph.SetMarkerColor(kRed)
P280graph.Draw('same')


avgPgraph.SetLineWidth(2)
avgPgraph.SetLineStyle(1)
avgPgraph.SetLineColor(kPurple+2)
avgPgraph.SetMarkerStyle(21)
avgPgraph.SetMarkerSize(0.5)
avgPgraph.SetMarkerColor(kPurple)
avgPgraph.Draw('same')

canvas.BuildLegend(0.9,.3,.9,.3,"","")
canvas.Print(pdf +'(')



graph = ["TCN18-180","TCN18-280","p_avg"]

gSCM180_0 = ROOT.TH1D("SCM_"+graph[0]+"_0A","Polarization over Time,  0 A SCM foil using  p_{foil} "+graph[0],numbin,array('d',bins))
gSCM180_50 = ROOT.TH1D("SCM_"+graph[0]+"_50A","Polarization over Time, 50 A SCM foil using  p_{foil} "+graph[0],numbin,array('d',bins))
gSCM180_100 = ROOT.TH1D("SCM_"+graph[0]+"_100A","Polarization over Time, 100 A SCM foil using  p_{foil} "+graph[0],numbin,array('d',bins))
gSCM180_150 = ROOT.TH1D("SCM_"+graph[0]+"_150A","Polarization over Time, 150 A SCM foil using  p_{foil} "+graph[0],numbin,array('d',bins))
gSCM180_175 = ROOT.TH1D("SCM_"+graph[0]+"_175A","Polarization over Time, 175 A SCM foil using  p_{foil} "+graph[0],numbin,array('d',bins))
gSCM180_200 = ROOT.TH1D("SCM_"+graph[0]+"_200A","Polarization over Time, 200 A SCM foil using  p_{foil} "+graph[0],numbin,array('d',bins))

gSCM180_0.GetXaxis().SetTitle('Time (s)')
gSCM180_0.GetYaxis().SetTitle('Polarizing Power, P_{SCM}')
gSCM180_0.SetDirectory(0)
gSCM180_0.SetMinimum(0.)
gSCM180_0.SetMaximum(1.)

gSCM280_0 = ROOT.TH1D("SCM_"+graph[1]+"_0A","Polarization over Time, 0 A SCM foil using  p_{foil} "+graph[1],numbin,array('d',bins))
gSCM280_50 = ROOT.TH1D("SCM_"+graph[1]+"_50A","Polarization over Time, 50 A SCM foil using  p_{foil} "+graph[1],numbin,array('d',bins))
gSCM280_100 = ROOT.TH1D("SCM_"+graph[1]+"_100A","Polarization over Time, 100 A SCM foil using  p_{foil} "+graph[1],numbin,array('d',bins))
gSCM280_150 = ROOT.TH1D("SCM_"+graph[1]+"_150A","Polarization over Time, 150 A SCM foil using  p_{foil} "+graph[1],numbin,array('d',bins))
gSCM280_175 = ROOT.TH1D("SCM_"+graph[1]+"_175A","Polarization over Time, 175 A SCM foil using  p_{foil} "+graph[1],numbin,array('d',bins))
gSCM280_200 = ROOT.TH1D("SCM_"+graph[1]+"_200A","Polarization over Time, 200 A SCM foil using  p_{foil} "+graph[1],numbin,array('d',bins))

gSCM280_0.GetXaxis().SetTitle('Time (s)')
gSCM280_0.GetYaxis().SetTitle('Polarizing Power, P_{SCM}')
gSCM280_0.SetMinimum(0.)
gSCM280_0.SetMaximum(1.)

gSCMavg_0 = ROOT.TH1D("SCM_"+graph[2]+"_0A","Polarization over Time, 0 A SCM foil using  p_{foil} "+graph[2],numbin,array('d',bins))
gSCMavg_50 = ROOT.TH1D("SCM_"+graph[2]+"_50A","Polarization over Time, 50 A SCM foil using  p_{foil} "+graph[2],numbin,array('d',bins))
gSCMavg_100 = ROOT.TH1D("SCM_"+graph[2]+"_100A","Polarization over Time, 100 A SCM foil using  p_{foil} "+graph[2],numbin,array('d',bins))
gSCMavg_150 = ROOT.TH1D("SCM_"+graph[2]+"_150A","Polarization over Time, 150 A SCM foil using  p_{foil} "+graph[2],numbin,array('d',bins))
gSCMavg_175 = ROOT.TH1D("SCM_"+graph[2]+"175A","Polarization over Time, 175 A SCM foil using  p_{foil} "+graph[2],numbin,array('d',bins))
gSCMavg_200 = ROOT.TH1D("SCM_"+graph[2]+"_200A","Polarization over Time, 200 A SCM foil using  p_{foil} "+graph[2],numbin,array('d',bins))

gSCMavg_0.GetXaxis().SetTitle('Time (s)')
gSCMavg_0.GetYaxis().SetTitle('Polarizing Power, P_{SCM}')

gSCMavg_0.SetMinimum(0.)
gSCMavg_0.SetMaximum(1.)

kRed = 632







#Style and Colour for 0 A SCM current
gSCM180_0.SetLineWidth(2)
gSCM180_0.SetLineStyle(1)
gSCM180_0.SetLineColor(kPurple+2)
gSCM180_0.SetMarkerStyle(21)
gSCM180_0.SetMarkerSize(0.5)
gSCM180_0.SetMarkerColor(kPurple)

gSCM280_0.SetLineWidth(2)
gSCM280_0.SetLineStyle(1)
gSCM280_0.SetLineColor(kPurple+2)
gSCM280_0.SetMarkerStyle(21)
gSCM280_0.SetMarkerSize(0.5)
gSCM280_0.SetMarkerColor(kPurple)

gSCMavg_0.SetLineWidth(2)
gSCMavg_0.SetLineStyle(1)
gSCMavg_0.SetLineColor(kPurple+2)
gSCMavg_0.SetMarkerStyle(21)
gSCMavg_0.SetMarkerSize(0.5)
gSCMavg_0.SetMarkerColor(kPurple)

#Style and Colour for 50 A SCM current
gSCM180_50.SetLineWidth(2)
gSCM180_50.SetLineStyle(1)
gSCM180_50.SetLineColor(kBlue+2)
gSCM180_50.SetMarkerStyle(21)
gSCM180_50.SetMarkerSize(0.5)
gSCM180_50.SetMarkerColor(kBlue)

gSCM280_50.SetLineWidth(2)
gSCM280_50.SetLineStyle(1)
gSCM280_50.SetLineColor(kBlue+2)
gSCM280_50.SetMarkerStyle(21)
gSCM280_50.SetMarkerSize(0.5)
gSCM280_50.SetMarkerColor(kBlue)

gSCMavg_50.SetLineWidth(2)
gSCMavg_50.SetLineStyle(1)
gSCMavg_50.SetLineColor(kBlue+2)
gSCMavg_50.SetMarkerStyle(21)
gSCMavg_50.SetMarkerSize(0.5)
gSCMavg_50.SetMarkerColor(kBlue)

#Style and Colour for 100 A SCM current
gSCM180_100.SetLineWidth(2)
gSCM180_100.SetLineStyle(1)
gSCM180_100.SetLineColor(kCyan+2)
gSCM180_100.SetMarkerStyle(21)
gSCM180_100.SetMarkerSize(0.5)
gSCM180_100.SetMarkerColor(kCyan)

gSCM280_100.SetLineWidth(2)
gSCM280_100.SetLineStyle(1)
gSCM280_100.SetLineColor(kCyan+2)
gSCM280_100.SetMarkerStyle(21)
gSCM280_100.SetMarkerSize(0.5)
gSCM280_100.SetMarkerColor(kCyan)

gSCMavg_100.SetLineWidth(2)
gSCMavg_100.SetLineStyle(1)
gSCMavg_100.SetLineColor(kCyan+2)
gSCMavg_100.SetMarkerStyle(21)
gSCMavg_100.SetMarkerSize(0.5)
gSCMavg_100.SetMarkerColor(kCyan)

#Style and Colour for 150 A SCM current
gSCM180_150.SetLineWidth(2)
gSCM180_150.SetLineStyle(1)
gSCM180_150.SetLineColor(kGreen+2)
gSCM180_150.SetMarkerStyle(21)
gSCM180_150.SetMarkerSize(0.5)
gSCM180_150.SetMarkerColor(kGreen)

gSCM280_150.SetLineWidth(2)
gSCM280_150.SetLineStyle(1)
gSCM280_150.SetLineColor(kGreen+2)
gSCM280_150.SetMarkerStyle(21)
gSCM280_150.SetMarkerSize(0.5)
gSCM280_150.SetMarkerColor(kGreen)

gSCMavg_150.SetLineWidth(2)
gSCMavg_150.SetLineStyle(1)
gSCMavg_150.SetLineColor(kGreen+2)
gSCMavg_150.SetMarkerStyle(21)
gSCMavg_150.SetMarkerSize(0.5)
gSCMavg_150.SetMarkerColor(kGreen)

#Style and Colour for 175 A SCM current
gSCM180_175.SetLineWidth(2)
gSCM180_175.SetLineStyle(1)
gSCM180_175.SetLineColor(kYellow+2)
gSCM180_175.SetMarkerStyle(21)
gSCM180_175.SetMarkerSize(0.5)
gSCM180_175.SetMarkerColor(kYellow)

gSCM280_175.SetLineWidth(2)
gSCM280_175.SetLineStyle(1)
gSCM280_175.SetLineColor(kYellow+2)
gSCM280_175.SetMarkerStyle(21)
gSCM280_175.SetMarkerSize(0.5)
gSCM280_175.SetMarkerColor(kYellow)

gSCMavg_175.SetLineWidth(2)
gSCMavg_175.SetLineStyle(1)
gSCMavg_175.SetLineColor(kYellow+2)
gSCMavg_175.SetMarkerStyle(21)
gSCMavg_175.SetMarkerSize(0.5)
gSCMavg_175.SetMarkerColor(kYellow)

#Style and Colour for 200 A SCM current
gSCM180_200.SetLineWidth(2)
gSCM180_200.SetLineStyle(1)
gSCM180_200.SetLineColor(kRed+2)
gSCM180_200.SetMarkerStyle(21)
gSCM180_200.SetMarkerSize(0.5)
gSCM180_200.SetMarkerColor(kRed)

gSCM280_200.SetLineWidth(2)
gSCM280_200.SetLineStyle(1)
gSCM280_200.SetLineColor(kRed+2)
gSCM280_200.SetMarkerStyle(21)
gSCM280_200.SetMarkerSize(0.5)
gSCM280_200.SetMarkerColor(kRed)

gSCMavg_200.SetLineWidth(2)
gSCMavg_200.SetLineStyle(1)
gSCMavg_200.SetLineColor(kRed+2)
gSCMavg_200.SetMarkerStyle(21)
gSCMavg_200.SetMarkerSize(0.5)
gSCMavg_200.SetMarkerColor(kRed)


for i in range(6):
    A = exSCM[i]
    #graph of SCM versus p_foil from 180
    Polarization(A,polA180,poleA180,i+1,graph[0])
gSCM180_0.Draw()
gSCM180_50.Draw('same')
gSCM180_100.Draw('same')
gSCM180_150.Draw('same')
gSCM180_175.Draw('same')
gSCM180_200.Draw('same')
canvas.BuildLegend(0.9,.3,.9,.3,"","")
canvas.Print(pdf)

gSCM280_0.SetDirectory(0)
for i in range(6):
    A = exSCM[i]
    #graph of SCM versus p_foil from 180
    Polarization(A,polA280,poleA280,i+1,graph[1])
gSCM280_0.Draw()
gSCM280_50.Draw('same')
gSCM280_100.Draw('same')
gSCM280_150.Draw('same')
gSCM280_175.Draw('same')
gSCM280_200.Draw('same')
canvas.BuildLegend(0.9,.3,.9,.3,"","")
canvas.Print(pdf)

gSCMavg_0.SetDirectory(0)
for i in range(6):
    A = exSCM[i]
    #graph of SCM versus p_foil from avg
    Polarization(A,avgP,avgeP,i+1,graph[2])
gSCMavg_0.Draw()
gSCMavg_50.Draw('same')
gSCMavg_100.Draw('same')
gSCMavg_150.Draw('same')
gSCMavg_175.Draw('same')
gSCMavg_200.Draw('same')

canvas.BuildLegend(0.9,.3,.9,.3,"","")
canvas.Print(pdf+')')



