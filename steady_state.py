import ROOT
import sys
import math
import numpy
import scipy.optimize
import inspect
import UCN

f = ROOT.TFile(sys.argv[1])

## Pass the run numbers to identify which cycles are to be used

listOfRuns = [1162, 1163]

beamOffBuffer = 100		# How long to wait after a beam-off period before data is considered valid agian (seconds)

beamOffThreshold = 0.9		# The threshold beam current which triggers an interval of beam-off (uA)

beamOnThreshold = 0.9		# Once this is reached during a beam-off period, after $beamOffBuffer seconds data is valid again (uA)

Li6Background = 4.84		# The background rate calculated by Wolfgang in the Li6 detector ( s^-1 uA^-1 )

totalBackground = 2.16		# The constant overall background rate

## Also create an overall plot, to which all of the different count vs temperature plots will be added

canvas = ROOT.TCanvas('name1', 'name2')

## Loop through all cycles, which creates a histogram of counts
## and rates for relevant cycles

def makePlots( detector, combined ):

	## If a combined plot is desired, then the TMultiGraph is create now. Will
	## end the function call if neither 'combined' nor 'separate' is returned

	if combined == 'combined' or combined == 'Combined':
		counts = ROOT.TMultiGraph()
	elif combined != 'separate' and combined != 'Separate':
		print('The plots must be either "combined", or "separate" ', combined, ' is not an option.')
		return None

	for cycle in f.cycledata:
	
		if cycle.runnumber in listOfRuns:
	
		 	## Retrieve relevant data
	
			reftime = cycle.start
			li6hits = numpy.array([h for h in getattr(cycle, 'Li6/hits')])
			he3hits = numpy.array([h for h in getattr(cycle, 'He3/hits')])
			Ttime = numpy.array([t for t in getattr(cycle, 'Source/timestamp')]) - reftime
			beam = numpy.array([c*k for c, k in zip(getattr(cycle, 'Beamline/B1V_KSM_PREDCUR'), getattr(cycle, 'Beamline/B1V_KSM_BONPRD'))])
			kicker = numpy.array([k for k in getattr(cycle, 'Beamline/B1V_KSM_BONPRD')])
			Btime = numpy.array([t for t in getattr(cycle, 'Beamline/timestamp')]) - reftime
			iv1 = numpy.array([v for v in getattr(cycle, 'Source/UCN_UGD_IV1_STATON')])
			ts11 = numpy.array([t for t in getattr(cycle, 'Source/UCN_ISO_TS11_RDTEMP')])
			ts12 = numpy.array([t for t in getattr(cycle, 'Source/UCN_ISO_TS12_RDTEMP')])
			ts14 = numpy.array([t for t in getattr(cycle, 'Source/UCN_ISO_TS14_RDTEMP')])
			ts16 = numpy.array([t for t in getattr(cycle, 'Source/UCN_ISO_TS16_RDTEMP')])
			pg9l = numpy.array([p for p in getattr(cycle, 'Source/UCN_ISO_PG9L_RDPRESS')])
			pg9h = numpy.array([p for p in getattr(cycle, 'Source/UCN_ISO_PG9H_RDPRESS')])

			## Depending on the detector to be used as specified in the function call
			## set hits to refer to either Li6 or He3 hits

			if detector == 'Li6' or detector == 'li6':
				hits = li6hits
			elif detector == 'He3' or detector == 'he3':
				hits = he3hits
			else:
				print('There is no ', detector, ' detector!')
				return None

			## Firstly, let's plot the UCN count rate over time to determine where the beam current dies
	
			## Bin counts to determine count rates

			## Adding additional background if using the Li6 detector:

			if detector == 'li6' or detector == 'Li6':
				background = totalBackground
			else:
				background = 0

			hist = numpy.histogram(hits, Ttime)
			rate = hist[0]/numpy.diff(Ttime) - background
			Ttimestamps = ((Ttime[:-1] + Ttime[1:]) / 2)
	
			## Make nice plots of count rates vs time and save
	
			ratevstime = ROOT.TGraph(len(beam), Btime, beam)
			ratevstime.GetXaxis().SetTitle('Time ( s )')
			ratevstime.GetYaxis().SetTitle('Beam Current ( #muA )')
			ratevstime.SetTitle('Beam Current vs Time in Run {0}'.format(cycle.runnumber))
			ratevstime.Draw('AP')
			canvas.Print('steady_state/beamVsTimeRun{0}.pdf'.format(cycle.runnumber))
	
			## Now let's try to identify points where the beam current is too low to provide
			## reliable information, and filter out all count rates based on that
	
			## To identify intervals in which the beam is off, every time the it drops below
			## 0.1uA, all measurements between that timestamp and 60 seconds after the first
			## time in which it recovers to at least 0.9uA are thrown out
	
			beamOffIntervals = []
	
			## Iterate through the beam current measurements, create a new interval within
			## beamOffIntervals which bookends a period in which the beam is off
	
			i = 0
	
			beamOff = False
	
			while i < len(beam):
				if beam[i] < beamOffThreshold and beamOff == False:
					beamOff = True
					interval = [Btime[i]]
				elif beam[i] > beamOnThreshold and beamOff == True:
					beamOff = False
					interval.append(Btime[i]+ beamOffBuffer)
					beamOffIntervals.append(interval)
				elif beam[i] < 1 and beamOff == True and i == len(beam) - 1:	# Special case if IV1 is still closed at end of cycle
					interval.append(Btime[i])
					beamOffIntervals.append(interval)
				i += 1
	
			## Looping through the array of timestamps for source measurements, discard those which
			## fall within any of the intervals of beam-off in beamOffIntervals
	
			beamOffIndices = []
	
			for i in range(len(Ttimestamps)):
				for interval in beamOffIntervals:
					if Ttimestamps[i] >= interval[0] and Ttimestamps[i] < interval[1]:
						beamOffIndices.append(i)
	
			beamOffIndices2 = []
	
			for i in range(len(Btime)):
				for interval in beamOffIntervals:
					if Btime[i] >= interval[0] and Btime[i] < interval[1]:
						beamOffIndices2.append(i)
	
			## Now, remove all of the measurements at those indices so that only
			## data points in which the beam is turned on are included
	
			listOfParams = [Ttimestamps, rate, iv1, ts11, ts12, ts14, ts16, pg9l, pg9h]
	
			Ttimestamps = numpy.delete(Ttimestamps, beamOffIndices)
			rate = numpy.delete(rate, beamOffIndices)
			iv1 = numpy.delete(iv1, beamOffIndices)
			ts11 = numpy.delete(ts11, beamOffIndices)
			ts12 = numpy.delete(ts12, beamOffIndices)
			ts14 = numpy.delete(ts14, beamOffIndices)
			ts16 = numpy.delete(ts16, beamOffIndices)
			pg9l = numpy.delete(pg9l, beamOffIndices)
			pg9h = numpy.delete(pg9h, beamOffIndices)
	
			Btime = numpy.delete(Btime, beamOffIndices2)
			beam = numpy.delete(beam, beamOffIndices2)
	
			## Repeat this same process for the instances in which IV1 is closed
	
			iv1ClosedIntervals = []
	
			i = 0
	
			iv1Closed = False
	
			while i < len(iv1):
				if iv1[i] < 1 and iv1Closed == False:
					iv1Closed = True
					interval = [Ttimestamps[i]]
				elif iv1[i] == 1 and iv1Closed == True:
					iv1Closed = False
					interval.append(Ttimestamps[i] + beamOffBuffer)
					iv1ClosedIntervals.append(interval)
				elif iv1[i] < 1 and iv1Closed == True and i == len(iv1) - 1:	# Special case if IV1 is still closed at end of cycle
					interval.append(Ttimestamps[i])
					iv1ClosedIntervals.append(interval)
				i += 1
	
			iv1ClosedIndices = []
	
			for i in range(len(Ttimestamps)):
				for interval in iv1ClosedIntervals:
					if Ttimestamps[i] >= interval[0] and Ttimestamps[i] < interval[1]:
						iv1ClosedIndices.append(i)
	
			iv1ClosedIndices2 = []
	
			for i in range(len(Btime)):
				for interval in iv1ClosedIntervals:
					if Btime[i] >= interval[0] and Btime[i] < interval[1]:
						iv1ClosedIndices2.append(i)
	
			Ttimestamps = numpy.delete(Ttimestamps, iv1ClosedIndices)
			rate = numpy.delete(rate, iv1ClosedIndices)
			ts11 = numpy.delete(ts11, iv1ClosedIndices)
			ts12 = numpy.delete(ts12, iv1ClosedIndices)
			ts14 = numpy.delete(ts14, iv1ClosedIndices)
			ts16 = numpy.delete(ts16, iv1ClosedIndices)
			pg9l = numpy.delete(pg9l, iv1ClosedIndices)
			pg9h = numpy.delete(pg9h, iv1ClosedIndices)
	
			Btime = numpy.delete(Btime, iv1ClosedIndices2)
			beam = numpy.delete(beam, iv1ClosedIndices2)
	
			## This part creates plots of the beam current and count rate as functions
			## of time, so that the success of the beam-off period removal can be evaluated
	
			## ratevstime = ROOT.TGraph(len(rate), Ttimestamps, rate)
			## ratevstime.GetXaxis().SetTitle('Time (s)')
			## ratevstime.GetYaxis().SetTitle('UCN count rate ( s^{-1} )')
			## ratevstime.SetTitle('UCN count rate in Run {0}'.format(cycle.runnumber))
			## ratevstime.Draw('AP')
			## canvas.Print('steady_state/VsTimeRun{0}.pdf'.format(cycle.runnumber))
	
			## beamvstime = ROOT.TGraph(len(beam), Btime, beam)
			## beamvstime.GetXaxis().SetTitle('Time (s)')
			## beamvstime.GetYaxis().SetTitle('Beam current ( #muA )')
			## beamvstime.SetTitle('Beam current in Run {0}'.format(cycle.runnumber))
			## beamvstime.Draw('AP')
			## canvas.Print('steady_state/beamOnVsTimeRun{0}.pdf'.format(cycle.runnumber))
			
			## Now, combine the measurements from PG9L and PG9H by selecting from PG9L
			## where PG9L < 2 Torr, and PG9H >= 2 Torr
	
			lowIndices = numpy.where(pg9l < 2)[0]
			highIndices = numpy.where(pg9l >= 2)[0]
	
			pg9 = numpy.empty(len(pg9l))
			pg9[lowIndices] = pg9l[lowIndices]
			pg9[highIndices] = pg9h[highIndices]
	
			## Using a vapour pressure-temperature correlation and then inverting
			## Correlation from https://doi.org/10.1103/PhysRev.100.743
	
			def HeVaporPressure(T):
		
				I = 4.6202
				A = 6.399
				B = 2.541
				C = 0.00612
				D = 0.5197
				a = 7.
				b = 14.14
	
				lnP = I - A/T + B*math.log(T) + C/2*T**2 - D*(a*b/(b**2 + 1) - 1./T)*math.atan(a*T - b)
				- a*D/2/(b**2 + 1)*math.log(T**2/(1 + (a*T - b)**2))
	
				return math.exp(lnP)
	
			def HeTemperature(P):
				
				return scipy.optimize.brentq(lambda T: HeVaporPressure(T) - P, 0.1, 4)
	
			## Create an array to fill with temperatures calculated from PG9
			## and then fill by looping through combined PG9 readings
			
			pg9Temps = numpy.empty(len(pg9))
	
			for i in range(len(pg9)):
				pg9Temps[i] = HeTemperature(pg9[i])
	
			## Finally, normalize the UCN count rate to the beam current
	
			beamCur = numpy.interp(Ttimestamps, Btime, beam)
	
			rate = rate/beamCur

			## Now that beam current is normalized, can subtract Li6 background

			if detector == 'li6' or detector == 'Li6':
				rate = rate - Li6Background
	
			## Now, make plots of the UCN count rate vs temperature using TS11, 
			## TS12, TS14, TS16, and the temperature calculated from PG9L & PG9H

			## As the count rate in the He3 detector is lower, the plot dimensions will be different			

			if detector == 'li6' or detector == 'Li6':
				canvas.DrawFrame(0.8, 0, 2.0, 4000)
			else:
				canvas.DrawFrame(0.8, 0, 2.0, 70)

			ratevsts11 = ROOT.TGraph(len(rate), ts11, rate)
			ratevsts11.GetXaxis().SetTitle('TS11 (K)')
			ratevsts11.GetYaxis().SetTitle('UCN count rate ( s^{-1} #muA^{-1})')
			ratevsts11.SetTitle('UCN count rate vs TS11 in Run {0}'.format(cycle.runnumber))
			ratevsts11.SetMarkerColor(ROOT.kPink)
			ratevsts11.SetLineColor(ROOT.kPink)
	
			ratevsts12 = ROOT.TGraph(len(rate), ts12, rate)
			ratevsts12.GetXaxis().SetTitle('TS12 (K)')
			ratevsts12.GetYaxis().SetTitle('UCN count rate ( s^{-1} #muA^{-1})')
			ratevsts12.SetTitle('UCN count rate vs TS12 in Run {0}'.format(cycle.runnumber))
			ratevsts12.SetMarkerColor(ROOT.kRed)
			ratevsts12.SetLineColor(ROOT.kRed)
	
			ratevsts14 = ROOT.TGraph(len(rate), ts14, rate)
			ratevsts14.GetXaxis().SetTitle('TS14 (K)')
			ratevsts14.GetYaxis().SetTitle('UCN count rate ( s^{-1} #muA^{-1})')
			ratevsts14.SetTitle('UCN count rate vs TS14 in Run {0}'.format(cycle.runnumber))
			ratevsts14.SetMarkerColor(ROOT.kGreen)
			ratevsts14.SetLineColor(ROOT.kGreen)
	
			ratevsts16 = ROOT.TGraph(len(rate), ts16, rate)
			ratevsts16.GetXaxis().SetTitle('TS16 (K)')
			ratevsts16.GetYaxis().SetTitle('UCN count rate ( s^{-1} #muA^{-1})')
			ratevsts16.SetTitle('UCN count rate vs TS16 in Run {0}'.format(cycle.runnumber))
			ratevsts16.SetMarkerColor(ROOT.kBlue)
			ratevsts16.SetLineColor(ROOT.kBlue)
	
			ratevspg9T = ROOT.TGraph(len(rate), pg9Temps, rate)
			ratevspg9T.GetXaxis().SetTitle('PG9 VP-Temperature (K)')
			ratevspg9T.GetYaxis().SetTitle('UCN count rate ( s^{-1} #muA^{-1})')
			ratevspg9T.SetTitle('UCN count rate vs PG9 VP-Temperature in Run {0}'.format(cycle.runnumber))
			ratevspg9T.SetMarkerColor(ROOT.kBlack)
			ratevspg9T.SetLineColor(ROOT.kBlack)

			## If the plots are desired to be separate, a new plot is created for every cycle

			if combined == 'separate' or combined == 'Separate':
				ratevstemp = ROOT.TMultiGraph()
				
				ratevstemp.Add(ratevsts11)
				ratevstemp.Add(ratevsts12)
				ratevstemp.Add(ratevsts14)
				ratevstemp.Add(ratevsts16)
				ratevstemp.Add(ratevspg9T)
	
				legend = ROOT.TLegend(0.1, 0.7, 0.48, 0.9)
				legend.AddEntry(ratevsts11, 'TS11', 'L')
				legend.AddEntry(ratevsts12, 'TS12', 'L')
				legend.AddEntry(ratevsts14, 'TS14', 'L')
				legend.AddEntry(ratevsts16, 'TS16', 'L')
				legend.AddEntry(ratevspg9T, 'PG9', 'L')
				legend.Draw('P')
	
				ratevstemp.GetXaxis().SetTitle('Temperature ( K )')
				ratevstemp.GetYaxis().SetTitle('UCN count rate ( s^{-1} #muA^{-1} )')
				ratevstemp.SetTitle('Li6 count rates by different temperature measuremnets')
				ratevstemp.Draw('P')
	
				canvas.Print('steady_state/{0}RateVsTempRun{1}.pdf'.format(detector, cycle.runnumber))
			else:
				counts.Add(ratevsts11)
				counts.Add(ratevsts12)
				counts.Add(ratevsts14)
				counts.Add(ratevsts16)
				counts.Add(ratevspg9T)

	if combined == 'combined' or combined == 'Combined':
		
		legend = ROOT.TLegend(0.1, 0.7, 0.48, 0.9)
		legend.AddEntry(ratevsts11, 'TS11', 'L')
		legend.AddEntry(ratevsts12, 'TS12', 'L')
		legend.AddEntry(ratevsts14, 'TS14', 'L')
		legend.AddEntry(ratevsts16, 'TS16', 'L')
		legend.AddEntry(ratevspg9T, 'PG9', 'L')
		legend.Draw('P')

		counts.GetXaxis().SetTitle('Temperature ( K )')
		counts.GetYaxis().SetTitle('UCN count rate ( s^{-1} #muA^{-1} )')
		counts.SetTitle('Li6 count rates by different temperature measuremnets')
		counts.Draw('P')
	
		canvas.Print('steady_state/{0}RateVsTempCombinedRuns.pdf'.format(detector))

## Now call this four times to get separate and combined plots for the Li6 and He3 detectors

makePlots('li6', 'separate')
makePlots('li6', 'combined')
makePlots('he3', 'separate')
makePlots('he3', 'combined')
