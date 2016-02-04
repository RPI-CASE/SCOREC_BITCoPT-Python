"""

drips.py contains a model for the drips system based on Mae-Ling's work

the functions are pretty straight forward, as the model is simple,
with support functions in dripssupport.py







"""



import base.blocks as b
import base.flux as f
import base.problem as p
import base.source as s

import time as cputime
import numpy as np
from drips_support import *


_solverTime = 0.

L = 5*0.3048
W = 3*0.3048
H = 3.1666*0.3048

def setAdsorptionHalflife(thalf):
  supportSetAdsorptionHalflife(thalf)

def relHumidity(pA,W,T):
  return pA/pvsat(T)*(W/(0.62198+W))

def weightFraction(relativeHumidity,pA,T):
  return 0.62198*relativeHumidity*pvsat(T)/(pA-relativeHumidity*pvsat(T))

def density(w,p,T):
  rho = p/(286.9*(T+273.15))
  return rho*(1+w)/(1+1.609*w)

def solve(init):

  clockTime = cputime.time()

  initial = {'mdot':init['airFlowRate'][0], 
            'mvdot':init['vaporFlowRate'][0],
            'Tw':init['waterT'][0],
            'Ta1':init['airT1'][0], 
            'Ta2':init['airT2'][0],
            'mwdot':0.}
  inlet = b.Block('inlet',desiccant,initial,{'mw':0.})

  airInt = b.Block('Interior','air',{Ta1:init['airInterior']})

  dripsBlocks = []

  # initialize air regions
  for i in range(1,init['n']+1):
    dripsBlocks.append(b.Block('dripsBlocks'+str(i),desiccant,
      initial,{'mw':init['mw'][i]}))

  for i in range(1,init['n']):
    dripsBlocks[i].addFlux(f.Flux(dripsBlocks[i-1],airMassFlux))
    dripsBlocks[i].addFlux(f.Flux(dripsBlocks[i-1],heatConvDesiccant))
    dripsBlocks[i].addSource(s.Source(heatConvAirVapor)) 
    dripsBlocks[i].addSource(s.Source(heatConvWater))
    dripsBlocks[i].addFlux(f.Flux(airInt,heatCondInterior,{'L':L}))

  dripsBlocks[0].addFlux(f.Flux(inlet,airMassFlux))
  dripsBlocks[0].addFlux(f.Flux(inlet,heatConvdripsBlocksiccant))
  dripsBlocks[0].addSource(s.Source(heatConvAirVapor))
  dripsBlocks[0].addSource(s.Source(heatConvWater))
  dripsBlocks[0].addFlux(f.Flux(airInt,heatCondInterior,{'L':L}))

  drips = p.Problem(dripsBlocks)

  drips.solve()
  results = {}

  var = ['Ta1','Ta2','mvdot','mdot','mwdot','Tw']
  for v in var:
    results[v] = np.array([dripsBlocks[i][v] for i in range(init['n'])])

  global _solverTime
  _solverTime += cputime.time() - clockTime 

  return results

def getSolverTime():
  return _solverTime
