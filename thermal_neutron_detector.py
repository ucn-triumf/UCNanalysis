import ROOT
import sys
import math
import numpy
import inspect
import UCN
import csv

f = ROOT.TFile(sys.argv[1])

canvas = ROOT.TCanvas('name1', 'name2')

firstRun = [cycle.runnumber for cycle in f.cycledata][0]

startupTime = 30				# The time from the start of the cycle, after irradiation has begun, until the LND reading stabilizes (s)

beamDropThreshold = 0.1				# The threshold value which will make any cycle in which the beam drops to this value invalid (uA)

beamDeviationThreshold = 0.02			# The threshold value of std. dev. of the beam which, if exceeded, will discard the cycle (uA)

lndPlateaus = []				# This will be populated with the average LND reading at all plateaus

beamPlateaus = []				# This will be populated with the average beam current at all plateaus

skippedCycles = []				# This will be populated with the run number of every skipped cycle

lastRun = -999					# Placeholder value

avgReading = {}					# Will be populated by the average LND readings for each run and then averaged for each TCN number

lowReading = -0.3e-6				# The value for which, when the average LND reading of a TCN experiment is below this, is considered low

## Taking the first to columns from spreadsheet which have a TCN number for each run and exporting as CSV - reading in

with open('RunToTCN.csv', mode='r') as infile:
	reader = csv.reader(infile)
	runToTCN = {rows[0]:rows[1] for rows in reader}

for TCN in runToTCN.values():
	avgReading[TCN] = []

for cycle in f.cycledata:

	## Retrieve relevant data

	lndReading = numpy.array([r for r in getattr(cycle, 'LND/LND_Reading')])
	Ttime = numpy.array([t for t in getattr(cycle, 'LND/timestamp')]) - cycle.start
	beam = numpy.array([b for b in getattr(cycle, 'Beamline/B1V_KSM_PREDCUR')])
	Btime = numpy.array([t for t in getattr(cycle, 'Beamline/timestamp')]) - cycle.start

	irradiationTime = cycle.beamonduration 

	if len(lndReading) == 0:
		continue
	elif Ttime[-1] < startupTime:
		continue

	## Discard cycles in which the beam drops below some lower threshold, $beamDropThreshold, or
	## the standard deviation of the beam is greater than some value, $beamDeviationThreshold

	if min(beam) < beamDropThreshold:
		skippedCycles.append(cycle.runnumber)
		continue
	if numpy.std(beam) > beamDeviationThreshold:
		skippedCycles.append(cycle.runnumber)
		continue

	# Find the beam at every instance of time at the LND refresh rate by interpolation

	beamReading = numpy.interp(Ttime, Btime, beam)

	## Make plots of LND reading vs time for all runs

	runTCN = runToTCN[str(cycle.runnumber)]

	if lastRun != -999:
		lastRunTCN = runToTCN[str(lastRun)]
	else:
		lastRunTCN = -999

	if lastRunTCN != runTCN and lastRunTCN != -999:
		canvas.Print('thermal_neutron_detector/lndReadingVsTime{0}.pdf)'.format(lastRunTCN))

	lndVsTime = ROOT.TGraph(len(lndReading), Ttime, lndReading)
	lndVsTime.GetXaxis().SetTitle('Time ( s )')
	lndVsTime.GetYaxis().SetTitle('LND Reading')
	lndVsTime.SetTitle('Normalized Thermal Neutron Detector Reading')
	lndVsTime.Draw('AP')
	
	if lastRunTCN != runTCN:
		canvas.Print('thermal_neutron_detector/lndReadingVsTime{0}.pdf('.format(runTCN))
	else:
		canvas.Print('thermal_neutron_detector/lndReadingVsTime{0}.pdf'.format(runTCN))

	lastRun = cycle.runnumber

	## Now want to extract the plateaus for each run. Do this by looping through all values of the LND
	## between $startupTime and either the end of the cyle, or $irradiationTime and then averaging them

	plateauVals = []
	plateauBeam = []

	irradiationOverIndex = numpy.where(Ttime > irradiationTime)[0]

	irradiationStartIndex = numpy.where(Ttime > startupTime)[0]

	## Ensure that only plateaus during the time in which beam is on are considered

	if len(irradiationStartIndex) == 0:
		continue
	elif len(irradiationOverIndex) == 0:			# The last index where the LND reading is during irradiation	
		duration = len(lndReading[irradiationStartIndex[0]:])
		init = irradiationStartIndex[0]
	else:
		duration = len(lndReading[irradiationStartIndex[0]:irradiationOverIndex[0]])
		init = irradiationStartIndex[0]

	if duration == 0:
		continue

	atPlateau = False

	for i in range(init, duration + init ):

		plateauVals.append(lndReading[i])
		plateauBeam.append(beamReading[i])
	
	lndPlateaus.append( sum(plateauVals)/len(plateauVals) )
	beamPlateaus.append( sum(plateauBeam)/len(plateauBeam) )

	normalizedReading = (sum(plateauVals)/len(plateauVals))/(sum(plateauBeam)/len(plateauBeam))

	avgReading[runTCN].append(normalizedReading)
	

## Finish .pdf of last run

lastRunTCN = runToTCN[str(lastRun)]

canvas.Print('thermal_neutron_detector/lndReadingVsTime{0}.pdf)'.format(lastRun))

lowReadings = open("thermal_neutron_detector/lowReadings.txt", "w+")

## Write the TCN numbers of those which have an average LND reading lower than $lowReading to a text file

for TCN in runToTCN.values():
	if len(avgReading[TCN]) != 0:
		averageVal = sum(avgReading[TCN])/len(avgReading[TCN])
		if averageVal > lowReading:
			lowReadings.write(TCN + "\n")

lowReadings.close()

## Create a new file that has no duplicate TCN numbers, and remove the old file

linesSeen = set()

badTCNs = open("thermal_neutron_detector/lowReadingTCNs.txt", "w+")

for line in open("thermal_neutron_detector/lowReadings.txt", "r"):
    if line not in linesSeen:
        badTCNs.write(line)
        linesSeen.add(line)

badTCNs.close()

## Now make plots of the LND plateaus vs the beam current

lndVsBeam = ROOT.TGraph(len(beamPlateaus), numpy.array(beamPlateaus), numpy.array(lndPlateaus))
lndVsBeam.GetXaxis().SetTitle('Beam Current ( #muA )')
lndVsBeam.GetYaxis().SetTitle('LND plateau reading')
lndVsBeam.SetTitle('')
lndVsBeam.Draw('AP')

canvas.Print('thermal_neutron_detector/lndReadingVsBeam.pdf')
