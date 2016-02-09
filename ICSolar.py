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
from joblib import Parallel, delayed

"""
This file contains an example run for ICSolar with geometry
It relies on the icsolar problem setup in src.icsolar,
as well as support functions in other files.

All results will end up in Results/<directory name>

"""

def solve(inputs):

  ############################################### input parameters
  geometry = inputs['geometry']
  useSunlitFraction = inputs['useSunlitFraction']
  directory = inputs['directory']
  startDay = inputs['range'][0]
  endDay = inputs['range'][1]
  DNI = inputs['DNI']
  extAirT = inputs['extAirT']
  DHI = inputs['DHI']
  
  days = endDay - startDay
  stepsPerDay = 24
  timesteps = days*stepsPerDay
  startStep = startDay*stepsPerDay
  endStep = (days+startDay)*stepsPerDay
  eta = 0.33 # electrical efficiency, constant for now.

  inletWaterTemp = 20.0 # inlet water temperature, 
  inletAirTemp = 20.0 # inlet air temperature
  interiorAirTemp = 20.0 # interior air temperature
  # 1.6 mL/s 
  waterFlowRate = 1.6*1000.*1e-6 # kg/s
  # 2.0 m/s * 0.16, cross sectional area, times density
  airFlowRate = 2.0*0.16*1.200 # kg/s

  init = {'dt':3600.,'length':icsolar.moduleH,'waterFlowRate':waterFlowRate,
  'airFlowRate':airFlowRate,'waterT':inletWaterTemp,'airT':inletAirTemp,
  'airInterior':interiorAirTemp}

  ############################################### set up results 
  if not os.path.exists('Results/'+directory):
    os.makedirs('Results/'+directory)

  # this bins things by facade direction
  bins = set([g.dir for g in geometry]);
  bins = list(bins)
  bins.append('total')
  thermal = dict([(bb,np.zeros(timesteps)) for bb in bins])
  elect = dict([(bb,np.zeros(timesteps)) for bb in bins])

  # this collects things by facade number
  facadeBins = range(len(geometry))
  epcFacade = []
  glazeFacade = []
  AOIFacade = []
  yawFacade = []
  pitchFacade = []
  shadeFacade = []
  thermalFacade = []
  electFacade = []
  dniFacade = []

  for i in facadeBins:
    epcFacade.append(np.zeros(timesteps))
    glazeFacade.append(np.zeros(timesteps))
    AOIFacade.append(np.zeros(timesteps))
    yawFacade.append(np.zeros(timesteps))
    pitchFacade.append(np.zeros(timesteps))
    shadeFacade.append(np.zeros(timesteps))
    thermalFacade.append(np.zeros(timesteps))
    electFacade.append(np.zeros(timesteps))
    dniFacade.append(np.zeros(timesteps))

  area = {bb:np.sum([g.height()*g.width()
   for g in geometry if g.dir == bb]) for bb in bins}
  area['total'] = np.sum([area[bb] for bb in bins])
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

    # once the first daytime hits, we are done initializing
    time = init['dt']*ts
    init['airExterior'] = extAirT[ts-startStep]

    # solarVector
    sunPosition = solar.getSunPosition(time)
    # print time/3600.,"---------------------------------"
    for g in geometry:
      index = geometry.index(g)
      matches = geometry.getMatches(g)
      # if there are no matches, or the first match index is greater than
      # this index, it is the first one that needs to be solved
      # otherwise, don't waste time solving this
      if not matches or geometry.index(matches[0]) > index:
        g.data['airTExternal'] = np.ones(g.nY)*extAirT[ts-startStep]
        if(DNI[ts-startStep] == 0 and DHI[ts-startStep] == 0):
          sunlit = 0.
          averaged = True
        elif useSunlitFraction is True:
          (sunlit,averaged) = shading.getMeanSunlitFraction(geometry,g,time,init['dt'],5)
        else:
          sunlit = 1.0
          averaged = True
        # iterate over modules, from bottom to top,
        # computing local shading fractions
        shade = np.ones(g.nY)
        for m in range(g.nY):
          shade[m] = shading.getUnshadedFraction(sunPosition,g,
            shadingIndices[index][m])


        glaze = max(solar.getGlazingTransmitted(sunPosition,g,1),0)
        AOI = solar.getAOI(sunPosition,g)
        (pitch,yaw) = solar.getArrayAngle(sunPosition,g)

        glazeFacade[geometry.index(g)][ts-startStep] = glaze
        AOIFacade[geometry.index(g)][ts-startStep] = AOI
        yawFacade[geometry.index(g)][ts-startStep] = yaw
        pitchFacade[geometry.index(g)][ts-startStep] = pitch
        shadeFacade[geometry.index(g)][ts-startStep] = shade[int(g.nY/2)]
        dniFacade[geometry.index(g)][ts-startStep] = DNI[ts-startStep]

        for match in matches:
          glazeFacade[geometry.index(match)][ts-startStep] = glazeFacade[geometry.index(g)][ts-startStep]
          AOIFacade[geometry.index(match)][ts-startStep] = AOIFacade[geometry.index(g)][ts-startStep]
          yawFacade[geometry.index(match)][ts-startStep] = yawFacade[geometry.index(g)][ts-startStep]
          pitchFacade[geometry.index(match)][ts-startStep] = pitchFacade[geometry.index(g)][ts-startStep]
          shadeFacade[geometry.index(match)][ts-startStep] = shadeFacade[geometry.index(g)][ts-startStep]
          dniFacade[geometry.index(match)][ts-startStep] = dniFacade[geometry.index(g)][ts-startStep]

        # don't solve this facade if the sun can't see it, but zero 
        # all data

        if(not daytime):
          g.data['DNI'] = np.zeros(g.nY)
          g.data['DHI'] = np.zeros(g.nY)
          g.data['DNIatModule'] = np.zeros(g.nY)
          g.data['DHIatModule'] = np.zeros(g.nY)
          g.data['receiverT'] = np.zeros(g.nY)
          g.data['airT'] = np.ones(g.nY)*inletAirTemp
          g.data['waterModuleT'] = np.ones(g.nY)*inletWaterTemp
          g.data['waterTubeT'] = np.ones(g.nY)*inletWaterTemp
          g.data['thermal'] = np.zeros(g.nY)
          g.data['electrical'] = np.zeros(g.nY)
          continue

        shadedVector = shading.applySunlitFraction(sunlit,g,averaged)

        g.data['DNI'] = DNI[ts-startStep]*shadedVector
        g.data['DHI'] = DHI[ts-startStep]*shadedVector
        g.data['Glazing'] = glaze*shadedVector

        g.data['DNIatModule'] = g.data['DNI']*glaze*shade
        # I think this simplifies your math, cos(pi/2-tilt) = sin(tilt)
        g.data['DHIatModule'] = g.data['DHI']*glaze*shade*0.5*(1.+np.sin(g.tilt))

        # init['Q_w'] = 0.024801*(0.66)*1e-3*g.data['DNIatModule']
        init['Q_d'] = 0.866*625.5*0.0001*g.data['DNIatModule']*(1.-eta)
        init['Q_a'] = np.zeros(g.nY)
        # (Q_c,Q_I) = support.getRadiativeGain(sunPosition,g,g.data['DNIatModule'],g.data['DHIatModule'])
        # init['Q_c'] = g.data['DHIatModule']*icsolar.moduleH*icsolar.moduleW
        init['numModules'] = g.nY
        
        # set up previous temperature
        if (not previousDayTime):
          init['previousWaterModuleT'] = init['waterT']*np.ones(g.nY)
          init['previousWaterTubeT'] = init['waterT']*np.ones(g.nY)
        else:
          init['previousWaterModuleT'] = g.data['waterModuleT']
          init['previousWaterTubeT'] = g.data['waterTubeT']

        # solve the problem
        results = icsolar.solve(init)

        # process results for storage 
        g.data['waterModuleT'] = results['waterModule']
        g.data['waterTubeT'] = results['waterTube']
        g.data['airT'] = results['airModule']
        g.data['receiverT'] = results['receiver']

        # compute the electrical and thermal energy
        electData = np.sum(eta \
          *1e-3*0.866*625.5*0.0001*g.data['DNIatModule'])*g.nX

        elect[g.dir][ts-startStep] += electData*(1+len(matches))
        electFacade[geometry.index(g)][ts-startStep] = electData
        # g.data['electrical'] = eta*1e-3*0.866*625.5*0.0001*g.data['DNIatModule']

        thermalData = np.sum(init['waterFlowRate']*4.218 \
          *(g.data['waterModuleT']-g.data['waterTubeT']))*g.nX
        g.data['thermal'] = init['waterFlowRate']*4.218 \
          *(g.data['waterModuleT']-g.data['waterTubeT'])

        epcFacade[geometry.index(g)][ts-startStep] = 0.866*625.5*0.0001*g.data['DNIatModule'][int(g.nY/2)]
        for match in matches:
          epcFacade[geometry.index(match)][ts-startStep] = epcFacade[geometry.index(g)][ts-startStep]
          electFacade[geometry.index(match)][ts-startStep] = electData
  
        # only add the contribution if its positive
        if(thermalData > 0):
          thermal[g.dir][ts-startStep] += thermalData*(1+len(matches))
          thermalFacade[index][ts-startStep] = thermalData
          for match in matches:
            thermalFacade[geometry.index(match)][ts-startStep] = thermalData
    previousDayTime = daytime
    thermal['total'][ts-startStep] = sum([thermal[b][ts-startStep] for b in bins])
    elect['total'][ts-startStep] = sum([elect[b][ts-startStep] for b in bins])
    # post processing and cleanup
    
    if(daytime and inputs['writeVTKFiles']):
      casegeom.writeVTKfile(geometry,'Results/'+directory+'/geom'+'0'*(4-len(str(ts)))+str(ts)+'.vtk','')
 
  print "runtime for days",startDay,"to",startDay+days-1,":",'%.2f' % (cputime.time()-clockStart)
  return (elect,thermal,electFacade,thermalFacade,epcFacade,glazeFacade,AOIFacade,shadeFacade,yawFacade,pitchFacade,dniFacade)
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

def run(init):
  ############################################### data loading 
  (DNI,extAirT,DHI,timezone,lat,lon,city) = weather.readTMY(init['TMY'])
  solar.setTimezone(timezone)
  solar.setLocation(lat,lon)
  geometry = casegeom.readFile(init['geomfile'],"ICSolar",
    init['useSunlitFraction'])
  geometry.computeBlockCounts(icsolar.moduleH,icsolar.moduleW)

  dataNames = ['DNI','DNIatModule','DHI','DHIatModule', 'Glazing', \
               'waterModuleT','waterTubeT','airT','airTExternal',
               'receiverT','thermal','electrical']

  geometry.initializeData(dataNames)
  casegeom.tiltGeometry(geometry,init['tilt'],init['useSunlitFraction'])

  stepsPerDay = 24
  timesteps = init['days']*stepsPerDay
  startStep = init['startDay']*stepsPerDay
  endStep = (init['days']+init['startDay'])*stepsPerDay

  procSplit = loadBalance(init['days'],init['startDay'],init['numProcs'],DNI)
  clockStart = cputime.time()

  ############################################### set up results
  if not os.path.exists(init['directory']):
    os.makedirs(init['directory'])

  # this bins things by facade direction
  bins = list(set([g.dir for g in geometry]));
  bins.append('total')
  thermal = dict([(bb,np.zeros(timesteps)) for bb in bins])
  elect = dict([(bb,np.zeros(timesteps)) for bb in bins])

  # this collects things by facade number
  facadeBins = range(len(geometry))
  epcFacade = []
  thermalFacade = []
  electFacade = []
  glazeFacade = []
  AOIFacade = []
  yawFacade = []
  pitchFacade = []
  shadeFacade = []
  dniFacade = []

  for i in facadeBins:
    epcFacade.append(np.zeros(timesteps))
    glazeFacade.append(np.zeros(timesteps))
    AOIFacade.append(np.zeros(timesteps))
    yawFacade.append(np.zeros(timesteps))
    pitchFacade.append(np.zeros(timesteps))
    shadeFacade.append(np.zeros(timesteps))
    thermalFacade.append(np.zeros(timesteps))
    electFacade.append(np.zeros(timesteps))
    dniFacade.append(np.zeros(timesteps))

  area = {bb:np.sum([g.height()*g.width()
   for g in geometry if g.dir == bb]) for bb in bins}
  area['total'] = np.sum([area[bb] for bb in bins])

  ############################################### run in parallel

  """
  This is the tricky pair, the processor splits are in procSplit, and each portion of the weather data
  thats needed is passed into the solver itself, so theres a lot of offsetting to be done,
  """
  inputs = []
  for (start,end) in procSplit:
    stepRange = slice(start*stepsPerDay,end*stepsPerDay)
    inputDict = {
    'range':(start,end),
    'directory':init['directory'],
    'geometry':geometry,
    'useSunlitFraction':init['useSunlitFraction'],
    'writeVTKFiles':init['writeVTKFiles'],
    'DNI':DNI[stepRange],
    'extAirT':extAirT[stepRange],
    'DHI':DHI[stepRange],
    }
    inputs.append(inputDict)

  results = Parallel(n_jobs = init['numProcs'])(delayed(solve)(inputs[i]) for i in range(init['numProcs']))
  print city,'final runtime is','%.2f' % (cputime.time()-clockStart)

  ############################################### collapse results

  for i in range(init['numProcs']):
    (start,end) = procSplit[i]
    resultRange = slice((start-init['startDay'])*stepsPerDay,(end-init['startDay'])*stepsPerDay)
    for b in bins:
      elect[b][resultRange] = results[i][0][b]
      thermal[b][resultRange] = results[i][1][b]
    for j in facadeBins:
      electFacade[j][resultRange] = results[i][2][j]
      thermalFacade[j][resultRange] = results[i][3][j]
      epcFacade[j][resultRange] = results[i][4][j]
      glazeFacade[j][resultRange] = results[i][5][j]
      AOIFacade[j][resultRange] = results[i][6][j]
      shadeFacade[j][resultRange] = results[i][7][j]
      yawFacade[j][resultRange] = results[i][8][j]
      pitchFacade[j][resultRange] = results[i][9][j]
      dniFacade[j][resultRange] = results[i][10][j]

  if init['writeDataFiles']:
    outputfiles = {}
    for b in bins:
      outputfiles[b] = (open('Results/'+init['directory']+'/output_'+b[0]+'.txt','w'))
      outputfiles[b].write('#time,thermal,electrical\n')

    outputFacadeFiles = {}
    for b in facadeBins:
      outputFacadeFiles[b] = (open('Results/'+init['directory']+'/output_'+str(b)+'.txt','w'))
      outputFacadeFiles[b].write('#time,thermal,electrical,epc,glaze,AOI,shade,yaw,pitch\n')

    for ts in range(startStep,endStep):

      for b in bins:
        outputfiles[b].write(','.join([str(ts),'%.12e' % thermal[b][ts-startStep], \
          '%.12e' % elect[b][ts-startStep]])+'\n')
      for b in facadeBins:
        outputFacadeFiles[b].write(','.join([str(ts),'%.3e' % thermalFacade[b][ts-startStep],\
          '%.3e' % electFacade[b][ts-startStep],'%.3e' % epcFacade[b][ts-startStep],\
          '%.3e' % glazeFacade[b][ts-startStep], '%.3e' % AOIFacade[b][ts-startStep],
          '%.3e' % shadeFacade[b][ts-startStep], '%.3e' % yawFacade[b][ts-startStep],
          '%.3e' % pitchFacade[b][ts-startStep],'%.3e' % dniFacade[b][ts-startStep]])+'\n')

    for b in bins:
      outputfiles[b].close()
    for b in facadeBins:
      outputFacadeFiles[b].close()
    ############################################### plot outputs

  linestyles = {'south':'--sk', 'roof':'--k', 'east':'--xk', 'west':'--ok', 'total':'-k'}

  ################# Cumulative Sum per Facade Area

  plt.subplot(1,3,1)
  for b in bins:
    plt.plot(np.arange(startStep,endStep)/stepsPerDay,np.cumsum(thermal[b])/area[b],linestyles[b],
      linewidth=2.0,label=b,markevery=30*stepsPerDay,markersize=7,fillstyle='none',markeredgewidth=2)

  axes = plt.gca()
  axes.tick_params(axis='both', which='major', labelsize=18)

  plt.legend(loc='upper left',fontsize=18)
  plt.ylabel('Thermal Output ($kWh/m^2$)',fontsize=18), plt.xlabel('Time (Days)',fontsize=18)
  plt.xlim(init['startDay'],init['startDay']+init['days'])

  plt.subplot(1,3,2)
  for b in bins:
    plt.plot(np.arange(startStep,endStep)/stepsPerDay,np.cumsum(elect[b])/area[b],linestyles[b],
      linewidth=2.0,label=b,markevery=30*stepsPerDay,markersize=7,fillstyle='none',markeredgewidth=2)

  plt.legend(loc='upper left',fontsize=18)
  plt.ylabel('Electrical Output ($kWh/m^2$)',fontsize=18), plt.xlabel('Time (Days)',fontsize=18)
  plt.xlim(init['startDay'],init['startDay']+init['days'])
  axes = plt.gca()
  axes.tick_params(axis='both', which='major', labelsize=18)

  plt.subplot(1,3,3)
  plt.plot(np.arange(startStep,endStep)/stepsPerDay,
    np.cumsum(elect['total']+thermal['total'])/area['total'],
    linewidth=2.0,label='total',linestyle='-',color='0.01')
  plt.plot(np.arange(startStep,endStep)/stepsPerDay,np.cumsum(thermal['total'])/area['total'],
    linewidth=2.0,label='thermal',linestyle='--',color='0.20')
  plt.plot(np.arange(startStep,endStep)/stepsPerDay,np.cumsum(elect['total'])/area['total'],
    linewidth=2.0,label='electrical',linestyle='-.',color='0.4')
  plt.ylabel('Total Output ($kWh/m^2$)',fontsize=18), plt.xlabel('Time (Days)',fontsize=18)
  plt.xlim(init['startDay'],init['startDay']+init['days'])
  axes = plt.gca()
  plt.legend(loc='upper left',fontsize=18)
  axes.tick_params(axis='both', which='major', labelsize=18)

  plt.gcf().set_size_inches(15,5)
  plt.tight_layout()
  plt.savefig('Results/'+init['directory']+'_outputsSum.png')
  plt.close()

if __name__ == "__main__":
  test.runTests()

  if not os.path.exists('Results'):
    os.makedirs('Results')

  tilt = 0

  init = { 
  'numProcs':8,
  'tilt':tilt,
  'startDay':0,
  'days':365,
  'directory':'NYC'+str(tilt),
  'TMY':'data/TMY/NYC.csv',
  'geomfile':'data/geometry/whole-building.txt',
  'useSunlitFraction':True,
  'writeDataFiles':True,
  'writeVTKFiles':True,
  }
  run(init)
