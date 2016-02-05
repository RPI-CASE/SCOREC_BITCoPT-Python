"""
Document here, for DRIPS model
"""

""" Required Modules """
import src.drips as drips
""" Optional Modules """
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import test
from src.casegeom import readValidationFile
from scipy.interpolate import interp1d

var = ['Ta1','Ta2','mvdot','mdot','wdot','Tw']
quantities = ['humidity','water','air','vapor','velocity','airT1','airT2','waterT','pressure']
qlabels = ['Relative Humidity','Water Mass','Air Flow Rate','Vapor Flow Rate','Velocity (m/s)',
  'Air Temp (Region)','Air Temp (Des.)','Water Temp (Des.)','pressure']
  
area = drips.W*drips.H
volume = area*drips.L


if __name__ == "__main__":
  # test.runTests()
  data = readValidationFile('data/drips/drips.csv')
  
  # startHumidity = 0.75
  # finalHumidity = 0.75
  # lagTime = 100.0
  for i in range(len(data['time'])):
    data['time'][i] *= 60.

  for i in range(len(data['inletT'])):
    data['inletT'][i] = (data['inletT'][i] -32.)/1.8
    data['outletT'][i] = (data['outletT'][i] -32.)/1.8

  humidityFunc = interp1d(data['time'],data['inletRH'])
  temperatureFunc = interp1d(data['time'],data['inletT'])
  humidityFuncOut = interp1d(data['time'],data['outletRH'])
  temperatureFuncOut = interp1d(data['time'],data['outletT'])
  temperature = data['inletT'][0] # degreees C
  pressure = 101325 # Pa
  velocity = 1.39 # m/s
  # humidity = startHumidity # fraction
  humidity = data['inletRH'][0]*0.01
  numRegions = 1
  finalTime = data['time'][-1]
  drips.setAdsorptionHalflife(60*60)
  nSteps = int(finalTime)/60


  weightFraction = drips.weightFraction(humidity,pressure,temperature)
  density = drips.density(weightFraction,pressure,temperature)
  inletMassFlowRate = area*density*velocity
  init = {}
  # initialize everything
  init['n'] = numRegions
  init['mw'] = np.ones(init['n']+1)*0.
  init['airFlowRate'] = np.ones(init['n']+1)*inletMassFlowRate
  init['airT1'] = np.ones(init['n']+1)*temperature
  init['airT2'] = np.ones(init['n']+1)*temperature
  init['vaporFlowRate'] = np.ones(init['n']+1)*weightFraction*inletMassFlowRate
  init['waterT'] = np.ones(init['n']+1)*temperature
  init['waterFlowRate'] = np.ones(init['n']+1)*0.
  init['pressure'] = np.ones(init['n']+1)*pressure
  init['dt'] = finalTime/nSteps
  init['airInterior'] = 24
  # initialize solution to plot
  solution = {}
  for v in quantities:
    solution[v] = np.zeros([nSteps+1,init['n']+1])
  solution['velocity'][0,:] = np.ones(init['n']+1)*velocity
  solution['humidity'][0,:] = np.ones(init['n']+1)*humidity
  solution['air'][0,:] = init['airFlowRate']
  solution['vapor'][0,:] = init['vaporFlowRate']
  solution['water'][0,:] = init['mw']
  solution['airT1'][0,:] = init['airT1']
  solution['airT2'][0,:] = init['airT2']
  solution['waterT'][0,:] = init['waterT']
  solution['pressure'][0,:] = init['pressure']

  solution['velocity'][:,0] = np.ones(nSteps+1)*velocity
  solution['humidity'][:,0] = np.ones(nSteps+1)*humidity
  solution['air'][:,0] = np.ones(nSteps+1)*init['airFlowRate'][0]
  solution['vapor'][:,0] = np.ones(nSteps+1)*init['vaporFlowRate'][0]
  solution['water'][:,0] = np.ones(nSteps+1)*init['mw'][0]
  solution['airT1'][:,0] = np.ones(nSteps+1)*init['airT1'][0]
  solution['airT2'][:,0] = np.ones(nSteps+1)*init['airT2'][0]
  solution['waterT'][:,0] = np.ones(nSteps+1)*init['waterT'][0]
  solution['pressure'][:,0] = np.ones(nSteps+1)*init['pressure'][0]

  # start stepping
  for i in range(1,nSteps+1):

    results = drips.solve(init)
    time = init['dt']*i

    for j in range(1,init['n']+1):
      T = results['Ta2'][j-1]
      W = results['mvdot'][j-1]/results['mdot'][j-1]
      init['mw'][j] = init['mw'][j] + init['dt']*(results['mwdot'][j-1])

      solution['humidity'][i,j] = drips.relHumidity(pressure,W,T)
      density = drips.density(humidity,pressure,T)
      solution['pressure'][i,j] = pressure
      solution['velocity'][i,j] = results['mdot'][j-1]/(density*area)
      solution['air'][i,j] = results['mdot'][j-1]
      solution['vapor'][i,j] = results['mvdot'][j-1]
      solution['water'][i,j] = init['mw'][j]
      solution['airT1'][i,j] = results['Ta1'][j-1]
      solution['airT2'][i,j] = results['Ta2'][j-1]
      solution['waterT'][i,j] = results['Tw'][j-1]

      init['airFlowRate'][j] = results['mdot'][j-1]
      init['vaporFlowRate'][j] = results['mvdot'][j-1]
      init['waterFlowRate'][j] = results['mwdot'][j-1]
      init['airT1'][j] = results['Ta1'][j-1]
      init['airT2'][j] = results['Ta2'][j-1]
      init['waterT'][j] = results['Tw'][j-1]


    temperature = temperatureFunc(time)
    humidity = humidityFunc(time)*0.01
    if( i % 100 == 0 ):
      print time,humidity
    weightFraction = drips.weightFraction(humidity,pressure,temperature)
    density = drips.density(weightFraction,pressure,temperature)
    inletMassFlowRate = area*density*velocity

    init['airFlowRate'][0] = inletMassFlowRate
    init['vaporFlowRate'][0] = weightFraction*inletMassFlowRate
    init['airT1'][0] = temperature
    solution['humidity'][i,0] = humidity
    solution['vapor'][i,0] = init['vaporFlowRate'][0]
    solution['air'][i,0] = inletMassFlowRate
    solution['airT1'][i,0] = temperature

  fig = plt.figure()
  plt.gcf().set_size_inches(15,8)

  timesteps = np.arange(nSteps+1)*init['dt']
  for i in range(8):
    plt.subplot(2,4,i+1)
    for j in range(init['n']+1):
      plt.plot(timesteps,solution[quantities[i]][:,j],'-',linewidth=2.0,
        markersize=10,label=str(j))
    plt.title(qlabels[i])

    plt.xlabel('time (s)')
    plt.xlim([0,finalTime])
    plt.ticklabel_format(useOffset=False)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    if(quantities[i] is 'humidity'):
      plt.ylim([0,max(data['inletRH'])])
      plt.legend(loc=0)
  plt.tight_layout()
  plt.show()

  fig = plt.figure()
  plt.gcf().set_size_inches(8,4)

  timesteps = np.arange(nSteps+1)*init['dt']
  plt.subplot(1,2,1)
  plt.plot(np.array(data['time'])/60.,data['inletRH'],'o',linewidth=2.0,markersize=5,label='inlet')
  plt.plot(np.array(data['time'])/60.,data['outletRH'],'o',linewidth=2.0,markersize=5,label='expt out')
  plt.plot(timesteps/60,solution['humidity'][:,1]*100.,'-',linewidth=2.0,markersize=10,label='model out')

  plt.title('Relative Humidity (%)')
  plt.xlabel('time (mins)')
  plt.xlim([0,finalTime/60])
  plt.ticklabel_format(useOffset=False)
  plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
  plt.legend(loc=0)

  plt.subplot(1,2,2)
  plt.plot(np.array(data['time'])/60.,data['inletT'],'o',linewidth=2.0,markersize=5,label='inlet')
  plt.plot(np.array(data['time'])/60.,data['outletT'],'o',linewidth=2.0,markersize=5,label='expt out')
  plt.plot(timesteps/60,solution['airT1'][:,1],'-',linewidth=2.0,markersize=2,label='model out')

  plt.title('airTemperature (C)')
  plt.xlabel('time (mins)')
  plt.xlim([0,finalTime/60])
  plt.ticklabel_format(useOffset=False)
  plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
  plt.tight_layout()
  plt.show()