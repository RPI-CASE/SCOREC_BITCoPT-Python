"""
icsolarsupport.py contains the physics needed for ICsolar problems

This is independent of the geometry

functions:

"""
import base.blocks as b
import base.flux as f
import base.problem as p
import base.source as s

import solar as solar
from math import pi

""" Materials 

these are just dictionaries of parameter functions used in blocks
since these functions are called often in flux and such calculations

"""
""" everything should in base units, J, g, m, etc """
def waterRho(B):
  return 1000.0 # kg/m^3
def waterCp(B):
  return 4218.0 # J/kg/K
def airRho(B):
  return 1.200 # kg/m^3 
def airCp(B):
  return 1005.0 # J/kg/K
def receiverh(B): # W/K
  return 1./(0.2 + 0.2487*(B.p['mdot']/(waterRho(B)*1e-6))**(-0.773))
def receivermCp(B):
  return 4.14 # J/K

constWater = \
{
'rho': waterRho,
'Cp': waterCp
}
constAir = \
{
'Cp': airCp,
'rho': airRho
}
receiver = \
{
'h': receiverh,
'mCp': receivermCp
}

""" 
Fluxes

Where B is a block, N is the block its connected to
and P is a dictionary of parameters
""" 
def heatCondWaterAir(B,N,P):
  return {'T':(B['T']-N['T'])*0.2333*P['L']} # W/K
def heatCondInterior(B,N,P):
  return {'T':(B['T']-N['T'])*1.6124*P['L']} # W/K
def heatCondExterior(B,N,P):
  return {'T':(B['T']-N['T'])*0.5233*P['L']} # W/K

def heatConv(B,N,P):
  return {'T':B.p['mdot']*B.P['Cp'](B)*(B['T']-N['T'])}
def heatConvDevice(B,N,P):
  return {'T':receiver['mCp'](B)/receiver['h'](B)/P['dt'] \
    *B.p['mdot']*B.P['Cp'](B)*(B['T']-N['T'])}

"""
getGapTransmittance:  Fraction of system area that is open for radiation to 
                      transmit through without being intercepted by any 
                      ICSF modules

input(s):             yaw, pitch
output(s):            fraction
"""
def getGapTransmittance(yaw, pitch):
  yaw *= 180.0/pi 
  pitch *= 180.0/pi 
  yawHalfRange = 34 
  pitchRange = 25 
  coverageYawN = 10.0/12.0 # yaw neutral coverage
  coveragePitchN = 10/11.5 # pitch neutral coverage

  #Calculate the fraction of glazing area which has a gap for radiation to transmit
  return 1.-min(coverageYawN+(1-coverageYawN)*abs(yaw/yawHalfRange),1)\
    *min(coveragePitchN+(1-coveragePitchN)*abs(pitch/pitchRange),1)

"""
radiativeGain:  Cavity Radiation Gain
                The Egen and Qgen constants are expected efficiency 
                based on W/m2 of DNI after the glazing, shading 
                and optical losses.

                It would be better to use the simulated electrical 
                and thermal generation as inputs and solve
                but it wouldn't be too innacurate to just assume 
                constant conversion percentages
                20% conversion of DNI at the receiver to Egen (W/m2)
                28% conversion of DNI at the receiver to Qgen (W/m2)
                Could use the below equations for better accuracy 
                but they account for array of modules
                output of this should be W/m2 for now but can change 
                (if geometry is fed in) to W for specific volumes
                Egen=((55)*-2.28e-4+0.1838)*DNIafterExt/1000,
                Qgen=0.024801*(0.66)*1e-3*DNIafterExt

input(s):       position of sun, geometry (and its data), and lens efficiency
output(s):      Solar heat gain to cavity, Solar heat gain to interior
"""
def getRadiativeGainAtTime(time, geometry, lensEfficiency):
  return getRadiativeGain(getSunPosition(time),geometry,lensEfficiency)

def getRadiativeGain(sunPositions, geometry, lensEfficiency):
  Egen = geometry.data['DNIatModule']*lensEfficiency*0.2
  Qgen = geometry.data['DNIatModule']*lensEfficiency*0.28
  (pitch, yaw) = solar.getArrayAngle(sunPositions,geometry)

  # transmittance of double pane interior window
  # will need to update this function with window properties from EP
  intGlazing = solar.getGlazingTransmitted(sunPositions,geometry,2)

  # Calculates a fraction of area not occupied by a module face
  # In other words, 0.1 means 10% of the area is open for light to 
  # freely pass by the module. Geometric gap calculation
  gapTranFactor = getGapTransmittance(yaw, pitch)
  
  # Calculate the DNI and DHI to the building interior using 
  # geometric gap calculations and double pane glazing transmittance
  dniToInterior = geometry.data['DNIafterExtGlass'] * gapTranFactor * intGlazing
  # Calculate the DHI to interior 
  dhiToInterior = geometry.data['DHIafterExtGlass'] * gapTranFactor * intGlazing

  # DNI that is not converted to electrical energy and thermal energy
  # or transmitted to interior is captured as heat inside the cavity
  dniHeat = geometry.data['DNIafterExtGlass'] - Egen - Qgen - dniToInterior
  # DHI that is not transmitted to the interior is captured at heat inside cavity
  dhiHeat = geometry.data['DHIafterExtGlass'] - dhiToInterior

  # Radiative heat gain to add to the cavity air temperature
  cavRadHeatGain = dniHeat + dhiHeat
  # Radiative heat gain to add to the interior heat transfer
  intRadHeatGain = dniToInterior + dhiToInterior

  # Units are W/m2, will need to multiply by 
  # surface area to determine the heat added to a specific volume.
  return {'cavRadHeatGain':cavRadHeatGain, 'intRadHeatGain':intRadHeatGain,\
    'dniToInterior':dniToInterior,'dhiToInterior':dhiToInterior,\
    'dniHeat':dniHeat,'dhiHeat':dhiHeat,'Egen':Egen,'Qgen':Qgen}

