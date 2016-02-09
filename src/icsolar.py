"""
icsolar.py

This defines the ICSolar model
where we our problem is similar to this

          ||  |\   |w|   |
          ||  | \  |w|   |
          ||--| 4|-|w|---|
          ||  | /  |w|   |
          ||  |/   |w|   |
          ||       |w|   |
exterior  ||    3  |w|   |  interior 
          ||       |w|   |
          ||  |\   |w|   |
          ||  | \  |w|   |
          ||--| 2|-|w|---|
          ||  | /  |w|   |
          ||  |/   |w|   |
          ||    1  |w|   |
          ||-------|w|---|
"""

import base.blocks as b
import base.flux as f
import base.problem as p
import base.source as s
import base.geom as g

import time as cputime
import numpy as np
from sys import exit
from icsolar_support import *

""" module constants here, for lack of better place """

moduleHeight = 0.3429
moduleWidth = 0.2794

_solverTime = 0.
"""
solve:          Builds the problem

input(s):       init, dictionary of physical parameters
                  numModules,   number of modules
                  Q_a,          heat added to air (optional)
                  Q_w or Q_d, heat added to water or device
                    one or the other needs to be present
                  exteriorAirTemp,  outside air temperature
                  interiorAirTemp,  inside air temperature
                  inletWaterTemp,       inlet water temperature
                  inletAirTemp,         inlet air temperature
                  waterFlowRate,  mass flow rate
                  airFlowRate,    mass flow rate
                  previousWaterModuleT,  water temperature in last step
                  previousWaterTubeT,    water temperature in last step
                  previousWaterFlowRate, water flow rate in last step
                  length,       length of tube between modules
                  dt,           timestep size
                  Q_c, heat added to cavity in air
output(s):      dictionary of results for four zones
"""

def solve(init):

  clockTime = cputime.time()

  dt = init['dt']

  if('previousWaterFlowRate' not in init.keys()):
    previousWaterFlowRate = init['waterFlowRate']
  else:
    previousWaterFlowRate = init['previousWaterFlowRate']

  # Inlet blocks, the material is a dictionary of parameter functions,
  # and mass flow rate is just a material
  waterInlet = b.Block('waterInlet',{'T':init['inletWaterTemp']},constWater,\
    {'mdot':init['waterFlowRate']})
  airInlet = b.Block('airInlet',{'T':init['inletAirTemp']},constAir,\
    {'mdot':init['airFlowRate']})

  # exterior and interior regions
  airExt = b.Block('Exterior',{'T':init['exteriorAirTemp']})
  airInt = b.Block('Interior',{'T':init['interiorAirTemp']})

  waterModule = []
  airModule = []
  waterTube = []
  airTube = []

  # initialize all regions
  for i in range(init['numModules']):
    waterModule.append(b.Block('waterModule_'+str(i),
      {'T':init['inletWaterTemp']},constWater,{'mdot':init['waterFlowRate']}))
    airModule.append(b.Block('airModule_'+str(i),
      {'T':init['inletAirTemp']},constAir,{'mdot':init['airFlowRate']}))   
    waterTube.append(b.Block('waterTube_'+str(i),
      {'T':init['inletWaterTemp']},constWater,{'mdot':init['waterFlowRate']}))
    airTube.append(b.Block('airTube_'+str(i),
      {'T':init['inletAirTemp']},constAir,{'mdot':init['airFlowRate']}))
  
  # now set up tube regions
  for i in range(init['numModules']):
    # water-air 
    waterTube[i].addFlux(f.Flux(airTube[i],heatCondWaterAir,{'L':init['length']}))
    airTube[i].addFlux(f.Flux(waterTube[i],heatCondWaterAir,{'L':init['length']}))

    # interior and exterior heat transfer
    airTube[i].addFlux(f.Flux(airInt,heatCondInterior,{'L':init['length']}))
    airTube[i].addFlux(f.Flux(airExt,heatCondExterior,{'L':init['length']}))

    if('Q_c' in init.keys()):
      airModule[i].addSource(s.Source('const',T = -init['Q_c'][i]))

  # now set the module regions
  for i in range(init['numModules']):

    # Energy source
    if('Q_w' in init.keys() and 'Q_d' not in init.keys()):
     waterModule[i].addSource(s.Source(s.constant,{'T':-init['Q_w'][i]}))
    elif('Q_d' in init.keys() and 'Q_w' not in init.keys()):
     waterModule[i].addSource(s.Source(s.constant,{'T':-init['Q_d'][i]
      *receiver['h'](waterModule[i])/receiver['mCp'](waterModule[i])}))
    else:
      exit("need only one energy source into the receiver or water\n")

    #  previous temp change due to receiver
    waterModule[i].addSource(s.Source(s.constant,{'T':-receiver['mCp'](waterModule[i])/\
      receiver['h'](waterModule[i])/dt*init['waterFlowRate']*waterCp(waterModule[i])*\
      (init['previousWaterModuleT'][i]-init['previousWaterTubeT'][i])}))

    tempChangeInTimeCoeff = receiver['mCp'](waterModule[i])/dt
    waterFlowRateCoeff = -(init['waterFlowRate']-previousWaterFlowRate)/dt\
    *waterCp(waterModule[i])*receiver['mCp'](waterModule[i])/receiver['h'](waterModule[i])

    # Sources due to temp change in time and derivative of water flow rate

    waterModule[i].addSource(s.Source(s.linear,{'T':tempChangeInTimeCoeff+waterFlowRateCoeff}))
    waterModule[i].addSource(s.Source(s.constant,{'T':-(tempChangeInTimeCoeff+waterFlowRateCoeff)*\
      init['previousWaterModuleT'][i]}))

    # part of the flux
    waterModule[i].addFlux(f.Flux(waterTube[i],heatConvDevice,{'dt':dt}))

    # Q_a from DNI
    if('Q_w' in init.keys()):
      airModule[i].addSource(s.Source(s.constant,{'T':-init['Q_a'][i]}))


  # connect up the regions
  for i in range(init['numModules']):
    airModule[i].addFlux(f.Flux(airTube[i],heatConv))
    waterModule[i].addFlux(f.Flux(waterTube[i],heatConv))

  airTube[0].addFlux(f.Flux(airInlet,heatConv))
  waterTube[0].addFlux(f.Flux(waterInlet,heatConv))
  for i in range(1,init['numModules']):
    airTube[i].addFlux(f.Flux(airModule[i-1],heatConv))
    waterTube[i].addFlux(f.Flux(waterModule[i-1],heatConv))

  # arrange them in a connected order
  blocksToSolve = []
  for i in range(init['numModules']):
    blocksToSolve.append(airTube[i])
    blocksToSolve.append(waterTube[i])
    blocksToSolve.append(airModule[i])
    blocksToSolve.append(waterModule[i])
  ICSolar = p.Problem(blocksToSolve)
  ICSolar.setBand((2,2))
  ICSolar.solve()

  results = {}
  results['waterModule'] = np.array([waterModule[i]['T'] 
    for i in range(init['numModules'])])
  results['waterTube'] = np.array([waterTube[i]['T'] 
    for i in range(init['numModules'])])
  results['airModule'] = np.array([airModule[i]['T'] 
    for i in range(init['numModules'])])
  results['airTube'] = np.array([airTube[i]['T'] 
    for i in range(init['numModules'])])
  results['receiver'] = np.array([init['waterFlowRate']*waterCp(waterModule[i]) \
    /receiver['h'](waterInlet)*(results['waterModule'][i]-results['waterTube'][i]) \
    + results['waterModule'][i] for i in range(init['numModules'])])
  global _solverTime
  _solverTime += cputime.time() - clockTime 

  return results

def getSolverTime():
  return _solverTime
