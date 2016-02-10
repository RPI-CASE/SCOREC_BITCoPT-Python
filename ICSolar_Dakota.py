import src.icsolar as icsolar
import src.casegeom as casegeom
import src.shading as shading
import src.weather as weather
import src.solar as solar
import src.icsolar_support as support

import test
import numpy as np
import time as cputime
import matplotlib.pyplot as plt
import os
import sys
from joblib import Parallel, delayed
import time as cputime
import src.dakota as dakota
"""
Each day will get its own dakota input file
This file acts as a driver, calling dakota on itself.

When called with no inputs, this file starts the driving process,
If inputs are provided, it is assumed dakota has called it, and a seperate 
routine is called. This file uses a reduced version of the ICSolar.py 
whole-building model, with less input and output and two versions of the solve,

its a bit of an ad-hoc file, but conceptually the idea is,
for each day, DNI = DNI(1+stddev)
where stddev is a gaussian process variable, effectively a %
that should be between [-1,1], and is computed by dakota
Each day has a fixed stddev, based on TMY day in obtained from readTMYUQ

solveUQ below deals with all of this, creating mean and std deviations to pass into 
dakota by creating dakota input files and calling dakota on each, with each day
treated independently, 

Finally, the run function puts everything together once dakota is done,
and files are left in Results/UQ

the quadOrder variable controls the number of coefficients in the PCE expansion

"""


"""
This is the driver of the system, calling dakota over and over
for each day. This is only called once, and organizes everything
"""
def solveUQ(problemInputs,solverInputs):
  # get DNI and DNI standard deviation
  (DNI,DNIUQ,extAirT,DHI,timezone,lat,lon,city) = \
    weather.readTMYUQ(problemInputs['TMY'])
  start = problemInputs['range'][0]
  end = problemInputs['range'][1]

  # start collecting mean data and std deviations, and coefficients
  meanArray = np.zeros(end-start)
  stddevArray = np.zeros(end-start)
  coeffArray = []
  # for each day we are interested in
  for day in range(start,end):
    # the uncertainty on the day is the max seen in the TMY data
    UQ = max(DNIUQ[(day*24):((day+1)*24)]) 
    # write the dakota input file for that day
    dakotafile = 'ICS'+str(day)+'.in'
    dakota.writeICSDakotaInputFile('Results/UQ/'+problemInputs['directory']+'/'+dakotafile,
      UQ,day,problemInputs['directory'],problemInputs['quadOrder'])
    # call dakota
    os.system('dakota --output=Results/UQ/'+problemInputs['directory']+'/'+'dakota'+str(day)+'.out '+
      '--error=Results/UQ/'+problemInputs['directory']+'/dakota'+str(day)+'.err Results/UQ/'
      +problemInputs['directory']+'/'+dakotafile)
    # read the output file and collect output
    (mean,stddev,coeff) = dakota.readOutputFile('Results/UQ/'+problemInputs['directory']+\
      '/dakota'+str(day)+'.out',problemInputs['quadOrder'])
    meanArray[day-start] = mean
    stddevArray[day-start] = stddev
    coeffArray.append(coeff)
  return (meanArray,stddevArray,coeffArray)

"""
This takes in two sets of inputs (see __main__ at the bottom)
and is used per process, in parallel

it also calls icsolar.solve and handles
solving at each timestep, and keeping track of data

This function is called by dakota, over and over for each day

"""
def solve(problemInputs,solverInputs):
  ############################################### take data from problemInputs
  # set up geometry, unlike ICSolar.py geometry and TMY data is re-created on each call
  # to this, since dakota is calling this almost directly
  geometry = casegeom.readFile(problemInputs['geometry'],"ICSolar",
    init['useSunlitFraction'])
  geometry.computeBlockCounts(icsolar.moduleHeight,icsolar.moduleWidth)

  geometry.initializeData(solverInputs['moduleDataNames'])
  casegeom.tiltGeometry(geometry,init['tilt'],init['useSunlitFraction'])
  useSunlitFraction = problemInputs['useSunlitFraction']
  directory = problemInputs['directory']
  startDay = problemInputs['range'][0]
  endDay = problemInputs['range'][1]

  (DNI,exteriorAirTemp,DHI,timezone,lat,lon,city) = weather.readTMY(problemInputs['TMY'])
  solar.setTimezone(timezone)
  solar.setLocation(lat,lon)

  stepsPerDay = problemInputs['stepsPerDay']

  days = endDay - startDay
  timesteps = days*stepsPerDay
  startStep = startDay*stepsPerDay
  endStep = (days+startDay)*stepsPerDay

  ############################################### set up results 

  # this bins things by direction
  bins = set([g.dir for g in geometry]);
  bins = list(bins)
  bins.append('total')

  directionData = {name:{b:np.zeros(timesteps) for b in bins} \
    for name in solverInputs['directionDataNames']}

  ############################################### pre-computed information
  # this is a lookup table to connect each module to a shading table
  
  shadingIndices = []
  for i in range(len(geometry)):
    shadingIndices.append([])
  for g in geometry:
    for i in range(g.nY):
      shadingIndices[geometry.index(g)].append \
        (shading.getShadingIndex(0,i,1,g.nY))

  ############################################### daytime setup
  # idea is that we don't solve when its not daytime

  daytime = False
  previousDayTime = False

  ############################################### solver
  clockStart = cputime.time()
  for ts in range(startStep,endStep):
    # its daytime if DNI > 0 for three straight hours
    # this logic is suspect at best

    if(1 < ts and ts < endStep-1):
      daytime = (DNI[ts-startStep-1]+DNI[ts-startStep]+DNI[ts-startStep+1] > 0)

    # Do the DNI calculation for that hour, and zero it if the standard deviation is < -1
    DNI[ts-startStep] = DNI[ts-startStep]*(1+solverInputs['DNISTD'])
    if(DNI[ts-startStep] < 1e-10):
      DNI[ts-startStep] = 0.

    # once the first daytime hits, we are done initializing
    time = solverInputs['dt']*ts
    solverInputs['exteriorAirTemp'] = exteriorAirTemp[ts-startStep]

    sunPosition = solar.getSunPosition(time)

    for g in geometry:
      index = geometry.index(g)
      matches = geometry.getMatches(g)
      # if there are no matches, or the first match index is greater than
      # this index, it is the first one that needs to be solved
      # otherwise, don't waste time solving this
      if not matches or geometry.index(matches[0]) > index:

        # averaged -> shading vector applied uniformly, not to top or bottom

        if useSunlitFraction is True:
          (sunlit,averaged) = shading.getMeanSunlitFraction(geometry,g,time,solverInputs['dt'],5)
        else:
          # if not using sunlit fraction, then just pick 0 or 1
          averaged = True
          if(DNI[ts-startStep] == 0 and DHI[ts-startStep] == 0):
            sunlit = 0.
          else:
            sunlit = 1.0

        # iterate over modules, from bottom to top,
        # computing local shading fractions
        shade = np.ones(g.nY)
        for m in range(g.nY):
          shade[m] = shading.getUnshadedFraction(sunPosition,g,
            shadingIndices[index][m])

        glaze = max(solar.getGlazingTransmitted(sunPosition,g,1),0)

        # set up air temperature
        g.data['externalAirT'] = np.ones(g.nY)*exteriorAirTemp[ts-startStep]

        # don't solve this facade if the sun can't see it, but zero 
        # all data
        if(not daytime):
          for name in g.data:
            g.data[name] = np.zeros(g.nY)
          continue

        shadedVector = shading.applySunlitFraction(sunlit,g,averaged)

        g.data['DNI'] = DNI[ts-startStep]*shadedVector
        g.data['DHI'] = DHI[ts-startStep]*shadedVector
        g.data['Glazing'] = glaze*shadedVector

        g.data['DNIatModule'] = g.data['DNI']*glaze*shade
        g.data['DHIatModule'] = g.data['DHI']*glaze*shade*0.5*(1.+np.sin(g.tilt))

        # solverInputs['Q_w'] = 0.024801*(0.66)*1e-3*g.data['DNIatModule']
        solverInputs['Q_d'] = 0.866*625.5*0.0001*g.data['DNIatModule']*(1.-solverInputs['eta'])
        solverInputs['Q_a'] = np.zeros(g.nY)
        # (Q_c,Q_I) = support.getRadiativeGain(sunPosition,g,g.data['DNIatModule'],g.data['DHIatModule'])
        # solverInputs['Q_c'] = g.data['DHIatModule']*icsolar.moduleHeight*icsolar.moduleWidth
        solverInputs['numModules'] = g.nY
        
        # set up previous temperature
        if (not previousDayTime):
          solverInputs['previousWaterModuleT'] = solverInputs['inletWaterTemp']*np.ones(g.nY)
          solverInputs['previousWaterTubeT'] = solverInputs['inletWaterTemp']*np.ones(g.nY)
        else:
          solverInputs['previousWaterModuleT'] = g.data['waterModuleT']
          solverInputs['previousWaterTubeT'] = g.data['waterTubeT']

        # solve the problem
        results = icsolar.solve(solverInputs)

        # process results for storage 
        g.data['waterModuleT'] = results['waterModule']
        g.data['waterTubeT'] = results['waterTube']
        g.data['inletAirTemp'] = results['airModule']
        g.data['receiverT'] = results['receiver']

        # electrical calculations
        electData = np.sum(solverInputs['eta'] \
          *1e-3*0.866*625.5*0.0001*g.data['DNIatModule'])*g.nX

        directionData['elect'][g.dir][ts-startStep] += electData*(1+len(matches))

        # thermal calculations
        thermalData = np.sum(solverInputs['waterFlowRate']*4.218 \
          *(g.data['waterModuleT']-g.data['waterTubeT']))*g.nX

        g.data['thermal'] = solverInputs['waterFlowRate']*4.218 \
          *(g.data['waterModuleT']-g.data['waterTubeT'])

        # only add the contribution if its positive
        if(thermalData > 0):
          directionData['thermal'][g.dir][ts-startStep] += thermalData*(1+len(matches))
        # fill in matching data for thermal, electrical, and epc


    # done with the step        
    previousDayTime = daytime

    # set up directional data by summing over things
    directionData['thermal']['total'][ts-startStep] = \
      sum([directionData['thermal'][b][ts-startStep] for b in bins if b is not 'total'])
    directionData['elect']['total'][ts-startStep] = \
      sum([directionData['elect'][b][ts-startStep] for b in bins if b is not 'total'])
    # post processing and cleanup
    
  return directionData
"""
This is a load balancing step, because dividing up the days by the number of processors
isn't ideal, we can do a lot better with a little bit of work. The assumption is that
the workload/runtime is proportional to the number of hours of sunlight,
so we divide everything into days, getting a rough balance of work. As 
there are more hours in the summer, this generally leads to an uneven split,
for example, NYC with 4 processors splits as [(0, 102), (102, 181), (181, 261), (261, 365)]
with the first 102 days on proc, the next 79 on one, the next 80 on another, and the last 104
on the final processor

"""
def loadBalance(days,startDay,numProc,DNI):

  hoursPerProc = np.count_nonzero(DNI[startDay*24:((startDay+days)*24)])/numProc
  pairs = []

  procStartDay = startDay
  endDay = startDay
  for i in range(numProc-1):
    hours = 0
    while hours < hoursPerProc and endDay < startDay+days-1:
      hours += np.count_nonzero(DNI[endDay*24:(endDay+1)*24])
      endDay = endDay+1
    pairs.append((procStartDay,endDay))
    procStartDay = endDay
  pairs.append((procStartDay,startDay+days))
  return pairs

"""
This function gets called once, which starts solveUQ in parallel, 
each day can still be done in parallel, with dakota running in serial on each day
"""
def run(init,solverInputs):
  stepsPerDay = 24
  timesteps = init['days']*stepsPerDay
  startStep = init['startDay']*stepsPerDay
  endStep = (init['days']+init['startDay'])*stepsPerDay

  procSplit = loadBalance(init['days'],init['startDay'],init['numProcs'],DNI)
  clockStart = cputime.time()

  ############################################### set up results
  # create the directories
  if not os.path.exists('Results'):
    os.makedirs('Results')
  if not os.path.exists('Results/UQ'):
    os.makedirs('Results/UQ')
  if not os.path.exists('Results/UQ/'+init['directory']):
    os.makedirs('Results/UQ/'+init['directory'])

  ############################################### run in parallel

  """
  This is the tricky pair, the processor splits are in procSplit, and each portion of the weather data
  thats needed is passed into the solver itself, so theres a lot of offsetting to be done,

  This sets up the inputs to pass into the solver itself
  """
  problemInputs = []
  for (start,end) in procSplit:
    stepRange = slice(start*stepsPerDay,end*stepsPerDay)
    inputDict = {
    'range':(start,end),
    'directory':init['directory'],
    'geometry':init['geometry'],
    'useSunlitFraction':init['useSunlitFraction'],
    'TMY':init['TMY'],
    'DNI':DNI[stepRange],
    'exteriorAirTemp':exteriorAirTemp[stepRange],
    'DHI':DHI[stepRange],
    'stepsPerDay':init['stepsPerDay'],
    'quadOrder':init['quadOrder'],
    }
    problemInputs.append(inputDict)

  # sets up the data to plot, calls solveUQ
  mean = np.zeros(init['days'])
  stddev = np.zeros(init['days'])
  coeffArray = []
  for j in range(init['quadOrder']):
    coeffArray.append(np.zeros(init['days']))

  results = Parallel(n_jobs = init['numProcs'])(delayed(solveUQ)(problemInputs[i],solverInputs) \
    for i in range(init['numProcs']))

  for i in range(init['numProcs']):
    (start,end) = procSplit[i]
    mean[start:end] = results[i][0] # mean
    stddev[start:end] = results[i][1] # stddev
    for d in range(start,end):
      for j in range(init['quadOrder']):
      	coeffArray[j][d] = results[i][2][d-start][j]
  meansum = np.cumsum(mean)
  stddevsum = np.sqrt(np.cumsum(stddev*stddev))

  plt.subplot(2,1,1)
  plt.errorbar(range(init['startDay'],init['startDay']+init['days']),mean,yerr=2.0*stddev,fmt='o')
  plt.ylabel('Daily Output (kWh per m^2)'), plt.xlabel('Time (Days)')
  plt.xlim(init['startDay'],init['startDay']+init['days'])
  plt.subplot(2,1,2)
  plt.plot(range(init['startDay'],init['startDay']+init['days']),meansum,linewidth=2.0)
  plt.fill_between(range(init['startDay'],init['startDay']+init['days']), meansum-2.0*stddevsum, 
    meansum+2.0*stddevsum, facecolor='yellow', alpha=0.5)
  plt.ylabel('Total Output (kWh)'), plt.xlabel('Time (Days)')
  plt.xlim(init['startDay'],init['startDay']+init['days']-1)

  plt.gcf().set_size_inches(15,5)
  plt.tight_layout()
  plt.savefig('Results/UQ/'+init['directory']+'/'+str(init['quadOrder'])+'_outputs.png')
  plt.close()

  for i in range(init['quadOrder']):
    plt.semilogy(range(init['startDay'],init['startDay']+init['days']),abs(coeffArray[i])/abs(coeffArray[0]),'-',label=str(i),linewidth=2)
  plt.ylabel('|PCE Coefficient|'), plt.xlabel('Time (Days)')
  plt.xlim(init['startDay'],init['startDay']+init['days']-1)
  plt.legend(loc=0)
  plt.gcf().set_size_inches(5,5)
  plt.tight_layout()
  plt.savefig('Results/UQ/'+init['directory']+'/'+str(init['quadOrder'])+'_expansion.png')
  plt.close()

  # lets try and cleanup
  fileList = os.listdir('.')
  for f in fileList:
    if f.startswith('LHS') or '.rst' in f or 'S4' in f:
      os.rename(f,'Results/UQ/'+f)

if __name__ == "__main__":
  tilt = 0

  # these are general variables, days should be greater than 1 for meaningful plots
  init = { 
  'numProcs':8,
  'tilt':tilt,
  'startDay':0,
  'days':365,
  'directory':'NYC'+str(tilt),
  'TMY':'data/TMY/NYC.csv',
  'geometry':'data/geometry/whole-building.txt',
  'useSunlitFraction':True,
  'stepsPerDay':24,
  'quadOrder':2,
  }

  init['range'] = (init['startDay'],init['startDay']+init['days'])
  # these are run specific variables
  solverInputs = {
  'dt':init['stepsPerDay']/24.*3600.,
  'eta':0.33, # electrical efficiency, constant for now.
  'length':icsolar.moduleHeight,
  'inletWaterTemp':20.0,
  'inletAirTemp':20.0,
  'interiorAirTemp':20.0,
  'waterFlowRate':1.6*1000.*1e-6, # 1.6 mL/s in kg/s
  # 2.0 m/s * 0.16, cross sectional area, times density in kg/s
  'airFlowRate':2.0*0.16*1.200,
  # data for each facade (one number for each facade)
  'facadeDataNames':[],
  # data on modules (a number for each module)
  'moduleDataNames':['DNI','DNIatModule','DHI','DHIatModule', 'Glazing', \
                     'waterModuleT','waterTubeT','inletAirTemp','externalAirT',
                    'receiverT','thermal','electrical'],
  # data per direction
  'directionDataNames':['thermal','elect'],
  }

  # if there are arguments, dakota is calling this function,
  # bypass solveUQ here
  if (len(sys.argv)) > 1:
    paramfile = sys.argv[1]
    # read the parameter file, and start the run again
    paramsdict = dakota.readInputFile(paramfile)

    init['startDay'] = int(paramsdict['day'])
    solverInputs['DNISTD'] = float(paramsdict['DNISTD'])
    directionData = solve(init,solverInputs)
    outfile = open(sys.argv[2],'w')
    outfile.write(str(sum(directionData['thermal']['total'] + \
      directionData['elect']['total'])) +' f0\n')
    outfile.close()
  else:
    run(init,solverInputs)
