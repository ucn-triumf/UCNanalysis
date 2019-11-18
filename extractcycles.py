import ROOT
import multiprocessing
import re
import sys
import math
import os
from array import array
import numpy
import itertools


# This script reads Li6-detector data, He3-detector data, SourceEpics data, and BeamlineEpics data for each cycle in a set of Run-data files
# and maps it into a python dictionary
# - key1: [array with single entry]
# - key2: [array2 with single entry]
# - Source: { - channel: [array3 with N entries]
#             - channel2: [array4 with N entries] }
# - Beamline: { - channel: [array5 with N entries]
#               - channel2: [array6 with N entries] }
# - key3: [array7 with single entry]
# ... and so on
# 
# All data is stored in arrays, to make it compatible with ROOT trees.
# The script can be easily extended by adding new keys, dictionaries, and channels to function ReadUCNTree.
# This data structure is dynamically mapped into a ROOT tree and written to file ucn_output.root
# Make sure that all entries in the dictionary are properly initialized, otherwise the output tree might contain garbage


# pair up cycles from both detector frontends with smallest differences in start time
def MatchTransitions(runnumber, he3cyclestart, li6cyclestart):
  pairs = sorted(itertools.product(he3cyclestart, li6cyclestart), key = lambda t: abs(t[0] - t[1])) # get all possible pairs and sort by time difference
  matchedhe3 = []
  matchedli6 = []
  for pair in pairs: # go through all possible pairs
    if pair[0] not in matchedhe3 and pair[1] not in matchedli6: # if none of the two start times are already in the matched list, add the pair to the matched list
      offset = pair[0] - pair[1]
      if abs(offset) > 20: # discard if time difference is too large
        print(' Removing cycle at {0} from run {1} because offset is {2}'.format(pair[0], runnumber, offset))
      else:
        matchedhe3.append(pair[0])
        matchedli6.append(pair[1])
  matchedhe3.sort()
  matchedli6.sort()
  # report on any cycles that could not be matched
  if any(t not in matchedhe3 for t in he3cyclestart):
    print(' Found no match to He3 cycles at {0} in run {1}'.format([t for t in he3cyclestart if t not in matchedhe3], runnumber))
  if any(t not in matchedli6 for t in li6cyclestart):
    print(' Found no match to Li6 cycles at {0} in run {1}'.format([t for t in li6cyclestart if t not in matchedli6], runnumber))
  return matchedhe3, matchedli6


# filter spikes in UCN rate within one second after period transitions
def FilterPileup(hits, runnumber, transition):
  if runnumber not in range(2023, 2036): # filter only runs known to be affected
    return hits
  pileupbins = []
  for periodstart in transition: # go through period transitions
    start = periodstart
    end = periodstart + 1.
    rate, edges = numpy.histogram([h[0] for h in hits], (end - start)/0.001, (start, end)) # histogram detector hits within one second of period transition with 1ms resolution
    pileupbins = pileupbins + [(edges[i], edges[i + 1]) for i, r in enumerate(rate) if r > 3] # find bins where rate is above threshold
  filtered = [h for h in hits if not any([b[0] <= h[0] < b[1] for b in pileupbins])] # keep only detector hits in bins where rate is smaller than threshold
  if len(filtered) != len(hits):
    print(' Removed {0} suspected pileup events ({1}%) from run {2}'.format(len(hits) - len(filtered), (len(hits) - len(filtered))*100./len(hits), runnumber))
  return filtered


# plot detector rate by histogramming hits within cycle with given bin width
def RatePlot(trans, hits, resolution):
  bins = math.ceil(trans[-1] - trans[0])
  rateplot = ROOT.TH1I('rate', 'rate', int(bins/resolution), 0., float(bins))
  for hit in hits:
    if hit[0] >= trans[0] and hit[0] < trans[-1]:
      rateplot.Fill(hit[0] - trans[0])
  rateplot.SetDirectory(0)

  return rateplot


### Read data from Epics channels between timestamps start and end
### Parameters: tree to read from, list of channel names, start and end timestamps
### Returns: dictionary containing arrays with data for each channel
def ReadEpicsData(EpicsTree, channels, start, end):
  match = [re.match('(\w+)(?:\[(\d+)\])?', c) for c in channels]
  types = [EpicsTree.GetBranch(m.group(1)).GetTitle()[-1].swapcase() for m in match]
  data = dict(zip((m.group(1) + (m.group(2) if m.group(2) else '') for m in match), (array(t) for t in types)))
  # loop over Epics tree and copy data if timestamp between start and end
  for entry in EpicsTree:
    if start <= entry.timestamp <= end:
      for m, t in zip(match, types):
        if m.group(2):
          data[m.group(1) + m.group(2)].append(getattr(entry, m.group(1))[int(m.group(2))])
        else:
          data[m.group(0)].append(getattr(entry, m.group(1)))
  return data


### Read a raw run file potentially containing several cycles
### Parameters: filename
### Returns: run number, list of dictionaries for each cycle as explained in file header, lists of Li6- and He3-detector-rate histograms for each cycle
def ReadUCNTree(fn):
  if not os.path.exists(fn):
    print(' Skipping non-existing file {0}'.format(fn))
    return 0, {}, [], []
  # disable graphical output
  ROOT.gROOT.SetBatch(1)

  # open file and get all relevant trees
  runfile = ROOT.TFile(fn)
  ttranshe3 = runfile.Get('RunTransitions_He3')
  ttransli6 = runfile.Get('RunTransitions_Li-6')
  tbeam = runfile.Get('BeamlineEpicsTree')
  tli6 = runfile.Get('UCNHits_Li-6')
  the3 = runfile.Get('UCNHits_He3')
  tsource = runfile.Get('SourceEpicsTree')
  tscm = runfile.Get('SCMTree')
  tsequencer = runfile.Get('SequencerTree')
  tlnd = runfile.Get('LNDDetectorTree')
  if not ttranshe3 or not ttransli6 or not tbeam or not tli6 or not the3 or not tsource or not tscm or not tsequencer or not tlnd:
    print(' Skipping file {0} because it does not contain all necessary data'.format(fn))
    return 0, {}, [], []
  
  # extract run number from filename
  match = re.match('ucn_tree_(\d+).root',os.path.split(fn)[1])
  runnumber = int(match.group(1))

  # check if there is beam data available
  if tbeam.GetEntries() == 0:
    print(' Skipping run {0} because it contains no beam data to estimate beam current!'.format(runnumber))
    return 0, {}, [], []

  UCNhits = {'He3': [(hit.tUnixTime, hit.tChannel) for hit in the3 if hit.tIsUCN == 1],
             'Li6': [(hit.tUnixTime, hit.tChannel) for hit in tli6 if hit.tIsUCN == 1]}
  if len(UCNhits['He3']) == 0 and len(UCNhits['Li6']) == 0:
    print(' Skipping run {0} because no UCN were detected'.format(runnumber))
    return 0, {}, [], []

  sequencerenabled = max([0] + [s.sequencerEnabled for s in tsequencer])
  if ttranshe3.GetEntries() == 0 and ttransli6.GetEntries() == 0 and sequencerenabled == 0:
    # if the sequencer was not enabled we add a single cycle with start and end based on the available slowcontrol-data timestamps
    timestamps = [t.timestamp for t in tbeam] + [t.timestamp for t in tsource] + [t.timestamp for t in tscm] + [t.timestamp for t in tlnd] + [t.timestamp for t in tsequencer]
    print(' Sequencer was disabled in run {0}, trying to add artifical cycle from {1} to {2}'.format(runnumber,min(timestamps),max(timestamps)))
    startsHe3 = startsLi6 = [min(timestamps)]
  elif ttranshe3.GetEntries() == 0 or ttransli6.GetEntries() == 0:
    # skip run if the sequencer was enabled but there are no run transitions
    print(' Skipping run {0} because He3 ({1}) or Li6 ({2}) detector does not contain any cycles'.format(runnumber, ttranshe3.GetEntries(), ttransli6.GetEntries()))
    return 0, {}, [], []
  else:
    # if there are run transitions for both detectors try to match them based on their start time
    print('Reading {0} and {1} cycles from run {2} in {3}'.format(ttranshe3.GetEntries(), ttransli6.GetEntries(), runnumber, fn))
    startsHe3, startsLi6 = MatchTransitions(runnumber, [t.cycleStartTime for t in ttranshe3], [t.cycleStartTime for t in ttransli6])

  # initialize lists containing data dictionaries and detector-rate histograms for each cycle
  rundata = []
  rateplots = {'He3': [], 'Li6': []}

  # loop over the cycle start times
  for cyclenumber, (startHe3, startLi6) in enumerate(zip(startsHe3, startsLi6)):
    # initialize data dictionary for single cycle and read basic run parameters
    # all data is stored in arrays, to make it compatible with ROOT trees
    data = {}
    data['runnumber'] = array('i', [runnumber])
    data['cyclenumber'] = array('i', [cyclenumber])
    data['supercyclenumber'] = array('i', [0])
    data['start'] = array('d', [startHe3])
    data['timingoffset'] = array('d', [startHe3 - startLi6])
    data['beamonduration'] = array('d', [0.])
    data['beamoffduration'] = array('d', [0.])
    data['periods'] = {}
    data['periods']['durations'] = array('d', [0]*11)
    data['periods']['countsLi6'] = array('i', [0]*11)
    data['periods']['countsHe3'] = array('i', [0]*11)
    data['periods'] = dict(('valve{0}state'.format(valve), array('i', [0]*11)) for valve in range(8))

    timestamp = 0
    # get beam-on- and -off-durations from beamline Epics data closest to cycle start time
    beamstart = min([(abs(b.timestamp - data['start'][0]), b.B1V_KSM_RDBEAMON_VAL1*0.00088801, b.B1V_KSM_RDBEAMOFF_VAL1*0.00088801) for b in tbeam])
    data['beamonduration'] = array('d', [beamstart[1]])
    data['beamoffduration'] = array('d', [beamstart[2]])
    
    # get run transitions belonging to cycle start times
    transition = {'He3': next((t for t in ttranshe3 if t.cycleStartTime == startHe3), None),
                  'Li6': next((t for t in ttransli6 if t.cycleStartTime == startLi6), None)}

    periods = {}
    if not transition['He3'] and not transition['Li6']:
      # if there are no run transitions, set up cycle with single period with start and end based on slowcontrol-data timestamps
      lasthit = max(timestamps)
      periods = {'He3': [startHe3] + [lasthit]*11,
                 'Li6': [startLi6] + [lasthit]*11}
    else:
      # if there are run transitions get period timing from sequencer settings
      beaminterval = data['beamonduration'][0] + data['beamoffduration'][0]
      periods = {'He3': [startHe3] + [getattr(transition['He3'], 'cyclePeriod{0}EndTime'.format(p)) for p in range(10)] + [startHe3 + beaminterval],
                 'Li6': [startLi6] + [getattr(transition['Li6'], 'cyclePeriod{0}EndTime'.format(p)) for p in range(10)] + [startLi6 + beaminterval]}
      seqduration = periods['He3'][-2] - periods['He3'][0]
      if seqduration > beaminterval:
        print(' Skipping cycle {0} in run {1} because cycle sequence ({2}s) is longer than irradiation interval ({3}s+{4}s)'.format(cyclenumber, runnumber, seqduration, data['beamonduration'][0],
data['beamoffduration'][0]))
        continue

    # skip cycles that have total duration of 0s
    if periods['He3'][0] == periods['He3'][-1]:
      print(' Skipping cycle {0} in run {1} because it has a duration of 0s'.format(cyclenumber, runnumber))
      continue

    # extract period durations and valve states
    data['periods']['durations'] = array('d', numpy.diff(periods['He3']))
    if any(data['periods']['durations'] != numpy.diff(periods['Li6'])):
      # check that period durations match between detector frontends
      print(' Skipping cycle {0} in run {1} because there is a mismatch between period durations in He3 {2} and Li6 {3} detectors'.format(cyclenumber, runnumber, data['periods']['durations'], numpy.diff(periods['Li6'])))
      continue
    if transition['He3']:
      data['supercyclenumber'][0] = transition['He3'].superCycleIndex
      for valve in range(8):
        data['periods']['valve{0}state'.format(valve)] = array('i', [getattr(transition['He3'], 'valveStatePeriod{0}'.format(period))[valve] for period in range(10)] + [0])
    
    UCNhits['Li6'] = FilterPileup(UCNhits['Li6'], runnumber, periods['Li6'])

    # collect hits during cycle and plot rates for both detectors
    for det, resolution in zip(['He3', 'Li6'], [1., 1.]):
      data['periods']['counts' + det] = array('i', numpy.histogram([h[0] for h in UCNhits[det]], periods[det])[0])
      data[det] = {'hits': array('d', [h[0] - periods[det][0] for h in UCNhits[det] if periods[det][0] <= h[0] < periods[det][-1]]), 
                   'channel': array('i', [h[1] for h in UCNhits[det] if periods[det][0] <= h[0] < periods[det][-1]])}
      rateplots[det].append(RatePlot(periods[det], UCNhits[det], resolution))

    ### add entries to data dictionary here if you want it to be written to the output file ###

    # get relevant SourceEpics data during cycle and store it as sub-dictionary (add more channels if needed)
    data['Source'] = ReadEpicsData(tsource, ['timestamp', 'UCN_ISO_TS11_RDTEMP', 'UCN_ISO_TS12_RDTEMP', 'UCN_ISO_TS14_RDTEMP', 'UCN_ISO_TS16_RDTEMP',\
                                             'UCN_ISO_PG9L_RDPRESS', 'UCN_ISO_PG9H_RDPRESS', 'UCN_UGD_IV1_STATON', 'UCN_EXP_IG5_RDVAC', 'UCN_EXP_IG6_RDVAC'], \
                                   periods['He3'][0], periods['He3'][-1])

    # get relevant BeamlineEpics data during beam-on time and store it as sub-dictionary (add more channels if needed)
    beamstop = periods['He3'][-1] if not transition['He3'] and not transition['Li6'] else data['start'][0] + data['beamonduration'][0]
    data['Beamline'] = ReadEpicsData(tbeam, ['timestamp', 'B1V_KSM_PREDCUR', 'B1V_KSM_BONPRD'], data['start'][0], beamstop)

    # get relevant SCM data during cycle and store it as sub-dictionary
    data['SCM'] = ReadEpicsData(tscm, ['timestamp', 'SCMVoltages[3]'], periods['He3'][0], periods['He3'][-1])

    # get relevant LND data during cycle and store it as sub-dictionary
    data['LND'] = ReadEpicsData(tlnd, ['timestamp', 'LND_Reading'], periods['He3'][0], periods['He3'][-1])

    # when done, add cycle data to list
    rundata.append(data)

  # return run number, list of data dictionaries for each cycle, and lists of detector rate histograms for each cycle
  return runnumber,rundata,rateplots['Li6'],rateplots['He3']


### Write dictionary as explained in file header to ROOT tree
### Parameters: ROOT tree to fill, dictionary containing data for single cycle
def WriteCycleData(otree, cycledata):
  # initialize dictionary storing lengths of channels that have more than one entry
  l = {}
  # loop over data of cycle
  for name in cycledata:
    data = cycledata[name]
    # if the data is another dictionary: create sub-branches for each entry in dictionary
    if isinstance(data, dict):
      # sub-branches contain arrays, each with the same length, store length
      l[name] = array( 'i', [ len(next(iter(data.values()))) ] )
      br = otree.GetBranch(name)
      if not br:
        # if it does not exist create new branch with dictionary entries as leafs, including leaf containing the length of the arrays
        leaflist = name + 'Len/I:' + ':'.join([channel + '[' + name + 'Len]/' + data[channel].typecode.swapcase() for channel in data])
        br = otree.Branch(name, 0, leaflist)

      # set read address for each leaf
      br.GetLeaf(name + 'Len').SetAddress(l[name])
      for channel in data:
        if l[name][0] > 0:
          br.GetLeaf(channel).SetAddress(data[channel])

    # if the data is an array: create top-level branch
    elif isinstance(data, array):
      # top level branches may only contain a single value, no arrays
      assert(len(data) == 1)
      # if it does not exist create new branch
      if not otree.GetBranch(name):
        otree.Branch(name, data, name + '/' + data.typecode.swapcase())
      # if it already does exist just set read address
      else:
        otree.SetBranchAddress(name, data)
    else:
      # only dictionaries and arrays supported at the moment
      assert(True)

  # fill tree from set branch/leaf addresses
  otree.Fill()



### MAIN PROGRAM STARTS HERE ###


# open output file
ofile = ROOT.TFile('ucn_output.root', 'RECREATE')
# create output tree
otree = ROOT.TTree('cycledata', 'cycledata')

# set up pool of threads to handle several files in parallel
pool = multiprocessing.Pool()
# loop over files given as command line parameters and call ReadUCNTree for each
for runnumber,rundata,li6,he3 in pool.imap_unordered(ReadUCNTree, sys.argv[1:]):
#for runnumber,rundata,li6,he3 in itertools.imap(ReadUCNTree, sys.argv[1:]): # replace line above with this one to turn off multithreading for debugging

  if not rundata:
    continue

  # set current file to output file
  ofile.cd()

  # loop over cycles in the data read from file
  print('Writing {0} cycles from run {1}'.format(len(rundata), runnumber))
  for cycledata in rundata:
    WriteCycleData(otree, cycledata)

  # write rate histograms for each detector and cycle
  if not ofile.Get('DetectorRates'):
    ofile.mkdir('DetectorRates')
  if not ofile.Get('run{0}'.format(runnumber)):
    ofile.mkdir('DetectorRates/run{0}'.format(runnumber))
  ROOT.gDirectory.cd('DetectorRates/run{0}'.format(runnumber))
  for cyclenumber,rate in enumerate(li6):
    if rate:
      rate.Write('li6_{0}_{1}'.format(runnumber, cyclenumber))
  for cyclenumber,rate in enumerate(he3):
    if rate:
      rate.Write('he3_{0}_{1}'.format(runnumber, cyclenumber))

ofile.Write()
