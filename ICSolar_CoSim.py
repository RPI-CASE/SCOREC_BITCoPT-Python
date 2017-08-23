import src.icsolar as icsolar
import src.casegeom as casegeom
import src.shading as shading
import src.weather as weather
import src.solar as solar
import src.icsolar_support as support
import src.base.geom as sgeom
import test
import numpy as np
import time as cputime
import matplotlib.pyplot as plt
import os
import sys
from joblib import Parallel, delayed
from datetime import timedelta
from pyfmi import load_fmu
"""
This file contains an example run for ICSolar with geometry
It relies on the icsolar problem setup in src.icsolar,
as well as support functions in other files.

All results will end up in Results/<directory name>

There are three ways to get data out of the simulation.
The names of this data is in the solverInputs  <name>DataNames

There are Paraview VTK files for each module on the facade,
the data of which is dictated in the geometry, g.data

There are facadeData with facadeDataNames which are per facade

There is directionData with directionDataNames which are per direction

Each of these is defaulted to zero and needs to be
filled inside the solve function

"""



"""
This is used to disable the warnings which occur when 
solving an ill-conditioned problem
such as DNI ~ 0. 

The warning:

minpack.py:236: RuntimeWarning: The iteration is not making 
  good progress, as measured by the 
  improvement from the last ten iterations.

is fairly common, and as such is filtered out
"""

import warnings
warnings.filterwarnings("ignore", \
  message = 'The iteration is not making good progress')

"""
This takes in two sets of inputs (see __main__ at the bottom)
and is used per process, in parallel

it also calls icsolar.solve and handles
solving at each timestep, and keeping track of data

This is where data gets created and stored

"""
def solve(init,problemInputs,solverInputs):

  ############################################### take data from problemInputs
  geometry = problemInputs['geometry']
  useSunlitFraction = problemInputs['useSunlitFraction']
  directory = problemInputs['directory']
  startDay = problemInputs['range'][0]
  endDay = problemInputs['range'][1]
  DNI = problemInputs['DNI']
  exteriorAirTemp = problemInputs['exteriorAirTemp']
  DHI = problemInputs['DHI']

  dt = solverInputs['dt']

  conversionFactor = solverInputs['lensEfficiency']*\
    solverInputs['concentrationRatio']*solverInputs['cellArea']

  stepsPerDay = problemInputs['stepsPerDay']

  interpTime = range(startDay*24,endDay*24) # hourly value from 0 to 8760, needs to be int
  # print "interpTime length:", len(interpTime)
  # print "DNI length:", len(DNI)
  # print "exteriorAirTemp length:", len(exteriorAirTemp)
  days = endDay - startDay
  timesteps = days*stepsPerDay
  startStep = startDay*stepsPerDay
  endStep = endDay*stepsPerDay
  # print "startStep", startStep, "endStep", endStep

  ############################################### load EnergyPlus FMu for co-simulation
  # fill in load code here 

  # load FMU
  fmuModel = load_fmu(init['fmuModelName'])

  # load options
  opts = fmuModel.simulate_options()
  # set number of communication points equal to timesteps
  opts['ncp'] = timesteps
  # manual initialize
  opts['initialize'] = False
  fmuModel.initialize(startDay*86400,endDay*86400)

  ############################################### set up results 

  # this bins things by facade direction
  bins = set([g.dir for g in geometry]);
  bins = list(bins)
  bins.append('total')
  print 'bins =', bins

  directionData = {name:{b:np.zeros(timesteps) for b in bins} \
    for name in solverInputs['directionDataNames']}

  # this collects things by facade number
  facadeBins = range(len(geometry))
  facadeData = {name:[] for name in solverInputs['facadeDataNames']}
  for name in solverInputs['facadeDataNames']:
    for i in range(len(geometry)):
      facadeData[name].append(np.zeros(timesteps))

  ############################################### pre-computed information
  # this is a lookup table to connect each module to a shading table
  
  shadingIndices = []
  for i in range(len(geometry)):
    shadingIndices.append([])
  for g in geometry:
    if g.facadeType == 'window':
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
    # print "ts:", ts
    # Change the timestep (ts) to time (in seconds) that the model can work with
    tsToSec = ts*dt
    tsToHour = ts*dt/3600.0 

    # interpolate weather data values onces as the beginning of the timestep
    # weather data is hourly, simulation runs sub-hourly
    dniForTS = np.interp(tsToHour,interpTime,DNI)
    dhiForTS = np.interp(tsToHour,interpTime,DHI) # Consider adding some heat to cavity air from DHI
    exteriorAirTempForTS = np.interp(tsToHour,interpTime,exteriorAirTemp)
    # print "exterior air temperature =", exteriorAirTempForTS

    # define inlet cavity air temperature to use the outdoor drybulb (or indoor drybulb, later)
    solverInputs['inletAirTemp'] = exteriorAirTempForTS

    # its daytime if DNI > 0 for three straight hours
    # this logic is suspect at best
    if(1 < ts and ts < endStep-1):
      # daytime = (DNI[ts-startStep-1]+DNI[ts-startStep]+DNI[ts-startStep+1] > 0)
      daytime = (np.interp(tsToHour-dt,interpTime,DNI)+dniForTS+np.interp(tsToHour+dt,interpTime,DNI) > 0)

    # once the first daytime hits, we are done initializing
    # time = solverInputs['dt']*ts
    solverInputs['exteriorAirTemp'] = exteriorAirTempForTS

    sunPosition = solar.getSunPosition(tsToSec)

    for g in geometry:
      if g.facadeType != 'window':
        continue
      index = geometry.index(g)
      matches = geometry.getMatches(g)

      # if there are no matches, or the first match index is greater than
      # this index, it is the first one that needs to be solved
      # otherwise, don't waste time solving this
      if not matches or geometry.index(matches[0]) > index:

        # averaged -> shading vector applied uniformly, not to top or bottom

        if useSunlitFraction is True:
          (sunlit,averaged) = shading.getMeanSunlitFraction(geometry,g,tsToSec,solverInputs['dt'],5)
        else:
          # if not using sunlit fraction, then just pick 0 or 1
          averaged = True
          if(dniForTS == 0 and dhiForTS == 0):
            sunlit = 0.
          else:
            sunlit = 1.0

        # iterate over modules, from bottom to top,
        # computing local shading fractions
        shade = np.ones(g.nY)
        for m in range(g.nY):
          shade[m] = shading.getUnshadedFraction(sunPosition,g,
            shadingIndices[index][m])

        # calculate a few variables 
        glaze = max(solar.getGlazingTransmitted(sunPosition,g,1),0)
        AOI = solar.getAOI(sunPosition,g)
        (pitch,yaw) = solar.getArrayAngle(sunPosition,g)
        # lets collect some data
        facadeData['glaze'][index][ts-startStep] = glaze
        facadeData['aoi'][index][ts-startStep] = AOI
        facadeData['yaw'][index][ts-startStep] = yaw
        facadeData['pitch'][index][ts-startStep] = pitch
        facadeData['shade'][index][ts-startStep] = shade[int(g.nY/2)] # shading in the middle
        facadeData['dni'][index][ts-startStep] = dniForTS

        # handle the matching facades, the ones we wont solve for, but make sure their data is
        # properly stored
        for name in ['glaze','aoi','yaw','pitch','shade','dni']:
          for match in matches:
            facadeData[name][geometry.index(match)] = facadeData[name][index]

        # set up air temperature
        g.data['externalAirT'] = np.ones(g.nY)*exteriorAirTempForTS

        # don't solve this facade if the sun can't see it, but zero 
        # all data
        if(not daytime):
          for name in g.data:
            g.data[name] = np.zeros(g.nY)
          # continue # commented out so that the model still provides values for all timestep

        shadedVector = shading.applySunlitFraction(sunlit,g,averaged)
        facadeData['shadedVector'] = shadedVector

        # solar energy after shading from other ICSF arrays 
        g.data['DNI'] = dniForTS*shadedVector
        g.data['DHI'] = dhiForTS*shadedVector
        g.data['Glazing'] = glaze*shadedVector

        # solar energy after the exterior glazing
        g.data['DNIafterExtGlass'] = g.data['DNI']*glaze
        g.data['DHIafterExtGlass'] = g.data['DHI']*glaze*0.5*(1.+np.sin(g.tilt))

        # solar energy at the module
        g.data['DNIatModule'] = g.data['DNIafterExtGlass']*shade
        g.data['DHIatModule'] = g.data['DHIafterExtGlass']*shade

        # energy per cell in the middle module
        facadeData['epc'][index][ts-startStep] = conversionFactor*g.data['DNIatModule'][int(g.nY/2)]
        
        # solverInputs['Q_w'] = 0.024801*(0.66)*1e-3*g.data['DNIatModule']
        solverInputs['Q_d'] = conversionFactor*g.data['DNIatModule']*(1.-solverInputs['eta'])
        
        radGain = support.getRadiativeGain(sunPosition,g,solverInputs['lensEfficiency'])
        solverInputs['Q_c'] = radGain['cavRadHeatGain']*icsolar.moduleHeight*icsolar.moduleWidth
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

        if g.dir == 'south':
          avgCavityAirTemp = np.average(results['airTube'])
          avgModAirTemp = np.average(results['airModule'])
          intRadHeatGain = np.average(radGain['intRadHeatGain'])
          fmuModel.set('SouthYaw',yaw)
          fmuModel.set('SouthPitch',pitch)
          fmuModel.set('SouthShadingVector',np.average(shadedVector))
          fmuModel.set('BEECavAirTemp',avgCavityAirTemp) # Will need to get smarter with this when working with multiple windows
          fmuModel.set('BEEWindowSolarInside',intRadHeatGain)
          if init['SHW'] == True:
            fmuModel.set('ICSFSouthThermalFluidFlow_kgps',g.nX*solverInputs['waterFlowRate'])
            fmuModel.set('ICSFSouthThermalFluidTemperature',results['waterModule'][-1])

        if g.dir == 'roof': # this needs to be updated to be more intelligent with how it sets output values to EnergyPlus inputs
          if init['SHW'] == True:
            fmuModel.set('ICSFRoofThermalFluidFlow_kgps',g.nX*solverInputs['waterFlowRate'])
            fmuModel.set('ICSFRoofThermalFluidTemperature',results['waterModule'][-1])
          fmuModel.set('RoofYaw',yaw)
          fmuModel.set('RoofPitch',pitch)
          fmuModel.set('RoofShadingVector',np.average(shadedVector))


        # co-simulation variables
        # init['interiorAirTemp']
        # init['exteriorAirTemp'] # can be pulled from the weather file
        # results['airModule'] and results['airTube']

        # electrical calculations
        eGenPerMod = solverInputs['eta']*conversionFactor*g.data['DNIatModule'] # Watts (W)
        electData = np.sum(solverInputs['eta']*1e-3\
          *conversionFactor*g.data['DNIatModule'])*g.nX # kilowatts (kW)

        directionData['elect'][g.dir][ts-startStep] += electData*(1+len(matches))
        facadeData['elect'][index][ts-startStep] = electData


        # thermal calculations
        thermalData = np.sum(solverInputs['waterFlowRate']*4.218 \
          *(g.data['waterModuleT']-g.data['waterTubeT']))*g.nX

        g.data['thermal'] = solverInputs['waterFlowRate']*4.218 \
          *(g.data['waterModuleT']-g.data['waterTubeT'])

        # only add the contribution if its positive
        if(thermalData > 0):
          directionData['thermal'][g.dir][ts-startStep] += thermalData*(1+len(matches))
          facadeData['thermal'][index][ts-startStep] = thermalData
        # fill in matching data for thermal, electrical, and epc
        for name in ['thermal','elect','epc']:
          for match in matches:
            facadeData[name][geometry.index(match)] = facadeData[name][index]


        if init['printToCMD'] == True:
          if tsToHour%1.0==0.0 and np.average(g.data['DNI']) > 0:
            disDate = timedelta(seconds=int(tsToSec))
            print ('{};  Q_c={:.1f}; ExtT [C]={:.1f}; ModT={:.1f}; CavT={:.1f}'.format(disDate,np.average(solverInputs['Q_c']),exteriorAirTempForTS,np.average(results['airModule']),np.average(results['airTube'])))
            print ('                     Alt={:.0f};  Azi={:.0f};  DNI={:.0f};  DHI={:.0f};  EPC={:.1f}'.format(sunPosition['altitude']*180/np.pi,sunPosition['azimuth']*180/np.pi,dniForTS,dhiForTS,facadeData['epc'][index][ts-startStep]))
            print ('                     ShadeVector={:.0f};  Sunlit={:.0f}').format(np.average(shadedVector), np.average(sunlit))
            print ('                     dniAfterExtGlass={:.0f};  dhiAfterExtGlass={:.0f}'.format(np.average(g.data['DNIafterExtGlass']),np.average(g.data['DHIafterExtGlass'])))
            print ('                     dniAtMod={:.0f};  dniToInterior={:.0f};  dhiToInterior={:.0f}'.format(np.average(g.data['DNIatModule']),np.average(radGain['dniToInterior']),np.average(radGain['dhiToInterior'])))
            print ('                     radEgen={:.0f};  radQgen={:.0f}'.format(np.average(radGain['Egen']),np.average(radGain['Qgen'])))
            print ('                     electData={:.0f};  thermalData={:.0f}'.format(electData*1000,thermalData*1000))

    # done with the step        
    previousDayTime = daytime

    # set up directional data by summing over things
    directionData['thermal']['total'][ts-startStep] = \
      sum([directionData['thermal'][b][ts-startStep] for b in bins if b is not 'total'])
    directionData['elect']['total'][ts-startStep] = \
      sum([directionData['elect'][b][ts-startStep] for b in bins if b is not 'total'])
    # post processing and cleanup

    # co-simulation
    fmuModel.set('ICSFElectricGeneration',directionData['elect']['total'][ts-startStep]*1000)
    fmuModel.set('ICSFThermalGeneration',directionData['thermal']['total'][ts-startStep]*1000)
    cosimRes = fmuModel.do_step(current_t=tsToSec,step_size=dt, new_step=True)
    
    if(daytime and problemInputs['writeVTKFiles']):
      casegeom.writeVTKfile(geometry,'Results/'+directory+'/VTK/geom'+'0'*(4-len(str(ts)))+str(ts)+'.vtk','')
  
  # terminate co-simulation
  fmuModel.terminate()
  
  print "runtime for day",startDay,"to end of day",startDay+days-1,":",'%.2f' % (cputime.time()-clockStart)
  return (directionData,facadeData)

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
This is the base program, after main is called, that organizes all the parallel data
collects information, makes a few plots, and sets up parameters for each run,
and saves all the files

This also does the parallelization, load balancing, time zones, etc

"""
def run(init,solverInputs):
  ############################################### data loading
  if init['TMY'].endswith('.csv'):
    data = weather.readTMY(init['TMY'])
  elif init['TMY'].endswith('.epw'):
    data = weather.readEPW(init['TMY'])
  else:
    sys.exit('The weather file was not loaded correctly. Please provide either an .epw file from \
      https://energyplus.net/weather or a .csv file from http://rredc.nrel.gov/solar/old_data/nsrdb/1991-2005/tmy3/by_state_and_city.html') 

  saveVars = zip(data['airTemp'])
  # head = ['extAirTemp']

  with open('WeatherDataCheck.csv', 'wb') as myfile:
    np.savetxt(myfile, saveVars, delimiter=',')
    myfile.close()

  DNI = data['DNI']
  DHI = data['DHI']
  exteriorAirTemp = data['airTemp']

  solar.setTimezone(data['timezone'])
  solar.setLocation(data['lat'],data['lon'])

  # set up geometry
  if init['useIDFGeom'] == True:
    geometry = casegeom.readNewFile(init['geomfile'],"ICSolar",
      init['useSunlitFraction'])
  else:
    geometry = casegeom.readFile(init['geomfile'],"ICSolar",
      init['useSunlitFraction'])
  geometry.computeBlockCounts(icsolar.moduleHeight,icsolar.moduleWidth)

  geometry.initializeData(solverInputs['moduleDataNames'])
  casegeom.tiltGeometry(geometry,init['tilt'],init['useSunlitFraction'])

  stepsPerDay = init['stepsPerDay']
  timesteps = init['days']*stepsPerDay
  startStep = init['startDay']*stepsPerDay
  endStep = (init['days']+init['startDay'])*stepsPerDay
  print "Run startStep:", startStep
  print "Run endStep:", endStep
  print "Run timesteps:", timesteps

  procSplit = loadBalance(init['days'],init['startDay'],init['numProcs'],DNI)
  clockStart = cputime.time()

  ############################################### set up results
  # create the directories
  if not os.path.exists('Results/'+init['directory']):
    os.makedirs('Results/'+init['directory'])

  if ( init['writeVTKFiles'] ):
    if not os.path.exists('Results/'+init['directory']+'/VTK'):
      os.makedirs('Results/'+init['directory']+'/VTK')    

  # this bins things by facade direction, 
  # it is directionData['variablename']['direction'][data]
  # and for facadeData['variablename'][index][data]
  bins = list(set([g.dir for g in geometry]));
  bins.append('total')
  directionData = {name:{b:np.zeros(timesteps) for b in bins} \
    for name in solverInputs['directionDataNames']}

  # this collects things by facade number
  facadeData = {name:[] for name in solverInputs['facadeDataNames']}
  for name in solverInputs['facadeDataNames']:
    for i in range(len(geometry)):
      facadeData[name].append(np.zeros(timesteps))

  area = geometry.getAreas()

  ############################################### run in parallel

  """
  This is the tricky pair, the processor splits are in procSplit, and each portion of the weather data
  thats needed is passed into the solver itself, so theres a lot of offsetting to be done,

  This sets up the inputs to pass into the solver itself
  """
  problemInputs = []
  for (start,end) in procSplit:
    stepRange = slice(start*24,end*24)
    # print "stepRange ", stepRange
    inputDict = {
    'range':(start,end),
    'directory':init['directory'],
    'geometry':geometry,
    'useSunlitFraction':init['useSunlitFraction'],
    'writeVTKFiles':init['writeVTKFiles'],
    'DNI':DNI[stepRange],
    'exteriorAirTemp':exteriorAirTemp[stepRange],
    'DHI':DHI[stepRange],
    'stepsPerDay':init['stepsPerDay'],
    }
    problemInputs.append(inputDict)

  # print ('The number of modules in a stack is {} and the number of stacks is {}').format()
  # results = solve(problemInputs[0],solverInputs)
  results = Parallel(n_jobs = init['numProcs'])(delayed(solve)(init,problemInputs[i],solverInputs) for i in range(init['numProcs']))
  print data['city'],'final runtime is','%.2f' % (cputime.time()-clockStart)

  ############################################### collapse results

  for i in range(init['numProcs']):
    (start,end) = procSplit[i]
    resultRange = slice((start-init['startDay'])*init['stepsPerDay'],(end-init['startDay'])*init['stepsPerDay'])
    # print "resultRange:",resultRange
    for b in bins:
      for name in solverInputs['directionDataNames']:
        directionData[name][b][resultRange] = results[i][0][name][b]
    for j in range(len(geometry)):
      for name in solverInputs['facadeDataNames']:
        facadeData[name][j][resultRange] = results[i][1][name][j]

  ############################################### write Data Files if requested

  if init['writeDataFiles']:
    outputfiles = {}
    outputFacadeFiles = {}

    # put headers into files
    for b in bins:
      outputfiles[b] = (open('Results/'+init['directory']+'/output_'+b+'.txt','w'))
      outputfiles[b].write('# time,'+','.join(solverInputs['directionDataNames'])+'\n')

    for j in range(len(geometry)):
      outputFacadeFiles[j] = (open('Results/'+init['directory']+'/output_'+str(j)+'.txt','w'))
      outputFacadeFiles[j].write('# time,'+','.join(solverInputs['facadeDataNames'])+'\n')

    # write the data
    for ts in range(startStep,endStep):
      for b in bins:
        outputfiles[b].write(str(ts*solverInputs['dt'])+','+','.join(['%.3e' % directionData[name][b][ts-startStep] \
          for name in solverInputs['directionDataNames']])+'\n')
      for j in range(len(geometry)):
        outputFacadeFiles[j].write(str(ts*solverInputs['dt'])+','+','.join(['%.3e' % facadeData[name][j][ts-startStep] \
          for name in solverInputs['facadeDataNames']])+'\n')

    # close the files
    for b in bins:
      outputfiles[b].close()
    for j in range(len(geometry)):
      outputFacadeFiles[j].close()

  # generate some plots in the Results/ directory
  # check the data is available
  # if ('thermal' in directionData and 'elect' in directionData):

  #   linestyles = {'south':'--sk', 'roof':'--k', 'east':'--xk', 'west':'--ok', 'total':'-k'}

  #   plt.subplot(1,3,1)
  #   for b in bins:
  #     plt.plot(np.arange(startStep,endStep)/stepsPerDay,np.cumsum(directionData['thermal'][b])/area[b],linestyles[b],
  #       linewidth=2.0,label=b,markevery=30*stepsPerDay,markersize=7,fillstyle='none',markeredgewidth=2)

  #   axes = plt.gca()
  #   axes.tick_params(axis='both', which='major', labelsize=18)

  #   plt.legend(loc='upper left',fontsize=18)
  #   plt.ylabel('Thermal Output ($kWh/m^2$)',fontsize=18), plt.xlabel('Time (Days)',fontsize=18)
  #   plt.xlim(init['startDay'],init['startDay']+init['days'])

  #   plt.subplot(1,3,2)
  #   for b in bins:
  #     plt.plot(np.arange(startStep,endStep)/stepsPerDay,np.cumsum(directionData['elect'][b])/area[b],linestyles[b],
  #       linewidth=2.0,label=b,markevery=30*stepsPerDay,markersize=7,fillstyle='none',markeredgewidth=2)

  #   plt.legend(loc='upper left',fontsize=18)
  #   plt.ylabel('Electrical Output ($kWh/m^2$)',fontsize=18), plt.xlabel('Time (Days)',fontsize=18)
  #   plt.xlim(init['startDay'],init['startDay']+init['days'])
  #   axes = plt.gca()
  #   axes.tick_params(axis='both', which='major', labelsize=18)

  #   plt.subplot(1,3,3)
  #   plt.plot(np.arange(startStep,endStep)/stepsPerDay,
  #     np.cumsum(directionData['thermal']['total']+directionData['elect']['total'])/area['total'],
  #     linewidth=2.0,label='total',linestyle='-',color='0.01')
  #   plt.plot(np.arange(startStep,endStep)/stepsPerDay,np.cumsum(directionData['thermal']['total'])/area['total'],
  #     linewidth=2.0,label='thermal',linestyle='--',color='0.20')
  #   plt.plot(np.arange(startStep,endStep)/stepsPerDay,np.cumsum(directionData['elect']['total'])/area['total'],
  #     linewidth=2.0,label='electrical',linestyle='-.',color='0.4')
  #   plt.ylabel('Total Output ($kWh/m^2$)',fontsize=18), plt.xlabel('Time (Days)',fontsize=18)
  #   plt.xlim(init['startDay'],init['startDay']+init['days'])
  #   axes = plt.gca()
  #   plt.legend(loc='upper left',fontsize=18)
  #   axes.tick_params(axis='both', which='major', labelsize=18)

  #   plt.gcf().set_size_inches(15,5)
  #   plt.tight_layout()
  #   plt.savefig('Results/'+init['directory']+'_outputs.png')
  #   plt.close()

if __name__ == "__main__":
  test.runTests()

  if not os.path.exists('Results'):
    os.makedirs('Results')

  tilt = 0
  stepsPerHour = 6

  # these are general variables
  init = { 
  'numProcs':1, # Do not change this value for now
  'tilt':tilt,
  'startDay':0,
  'days':365,
  'directory':'BEEUnit'+str(tilt),
  'TMY':'data/TMY/USA_NY_CentralPark.epw',
  'geomfile':'data/geometry/BEEWindowIDF.txt',
  'fmuModelName':'./data/fmu/in_cosim_NYC.fmu',
  'useIDFGeom':True,
  'SHW':True,
  'useSunlitFraction':True,
  'writeDataFiles':True,
  'writeVTKFiles':False,
  'printToCMD':False,
  'stepsPerHour':stepsPerHour,
  'stepsPerDay':24*stepsPerHour,
  }

  # these are run specific variables
  solverInputs = {
  'dt':24./init['stepsPerDay']*3600.,
  'lensEfficiency':0.886, # transmittance efficiency of the concentrating lens
  'concentrationRatio':625.5, # concentration of DNI by the fresnel lens
  'cellArea':0.0001, # area of the photovoltaic cell in m^2
  'eta':0.33, # electrical efficiency, constant for now.
  'length':icsolar.moduleHeight,
  'inletWaterTemp':20.0,
  'inletAirTemp':15.0,
  'interiorAirTemp':20.0,
  'waterFlowRate':1.6*1000.*1e-6, # 1.6 mL/s in kg/s
  # 2.0 m/s * 0.16, cross sectional area, times density in kg/s
  'airFlowRate':1.0*0.16*1.200,
  # data for each facade (one number for each facade)
  'facadeDataNames':['epc','glaze','aoi','yaw','pitch',
               'shade','shadedVector','dni','thermal','elect'],
  # data on modules (a number for each module)
  'moduleDataNames':['DNI','DNIafterExtGlass','DNIatModule','DHI','DHIafterExtGlass','DHIatModule', 'Glazing', \
                     'waterModuleT','waterTubeT','inletAirTemp','externalAirT',
                    'receiverT','thermal','electrical'],
  # data per direction
  'directionDataNames':['thermal','elect'],
  }

  print 'dt', solverInputs['dt']

  run(init,solverInputs)