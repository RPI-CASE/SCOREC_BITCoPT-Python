"""
icsolarsupport.py contains the physics needed for ICsolar problems

This is independent of the geometry

functions:

"""
import base.blocks as b
import base.flux as f
import base.problem as p
import base.source as s

from math import pi
import numpy as np

""" Materials """
""" everything should in base units, J, g, m, etc """
def waterRho(B):
  return 1000.0 # kg/m^3
def waterCp(B):
  return 4218.0 # J/kg/K
def airRho(B):
  return 1.200 # kg/m^3 
def airCp(B):
  return 1005.0 # J/kg/K
def vaporCp(B):
  return 1864.0 # J/kg/K
def hVap(B):
  return 2.260e6 # J/kg
desiccant = \
{
  'CpAir': airCp,
  'CpVapor': vaporCp,
  'hVap': hVap,
  'airRho': airRho,
  'waterRho': waterRho,
}

# Pa
def pvsat(T):
  return 10**(8.07131-1730.63/(233.426+T))*133.322368
""" 
Fluxes

Where B is a block, N is the block its connected to
and P is a dictionary of parameters
""" 
def enth(B,T):
  mf = B['mvdot']/B['mdot']
  return (1.0-mf)*B.P['CpAir'](B)*B[T]+ \
  mf*(B.P['hVap'](B)+B.P['CpVapor'](B)*B[T])

# initial value
_thalf = 1.0
def supportSetAdsorptionHalflife(thalf):
  global _thalf
  _thalf = thalf

def qw(B):
  global _thalf
  W = B['mvdot']/B['mdot']
  C = 101325/pvsat(B['Ta2'])*(W/(0.62198+W))
  if (C < 0.44):
    return 0.
  mwsat = 0.37*(C-0.44)**0.82*75.0
  # if(mwsat-B.mw) < 0.:
  #   return 0.
  return (mwsat-B.P['mw'])*0.69314718056/(_thalf)

def heatConvDesiccant(B,N,P):

  q = qw(B)
  h = enth(B,'Ta2')
  hN = enth(N,'Ta1')
  # print "Desiccant",B.name,N.name,B['mdot']*h-N['mdot']*hN
  return {'Ta2':B['mdot']*h-N['mdot']*hN+(q*B.P['CpVapor'](B)*B['Ta2']-q*B.P['hVap'](B))}

def heatConvAirVapor(B,P):

  h = enth(B,'Ta1')
  hN = enth(B,'Ta2')
  return {'Ta1':B['mdot']*(h-hN)}

def heatConvWater(B,P):
  q = qw(B)
  return {'Tw':B['mwdot']*B['Tw']-q*B['Ta2']}
def heatCondInterior(B,N,P):
  return {'Ta1':(B['Ta1']-N['Ta1'])*250.0*P['L']} # W/K
def heatCondExterior(B,N,P):
  return {'Ta1':(B['Ta1']-N['Ta1'])*0.5233*P['L']} # W/K

def airMassFlux(B,N,P):
  q = qw(B)
  return {'mdot':B['mdot']-(N['mdot']-q),
    'mvdot':B['mvdot']-(N['mvdot']-q), 
    'mwdot':B['mwdot']-q}

