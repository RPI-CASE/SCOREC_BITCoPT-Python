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
getGapTransmittance:

input(s):             yaw, pitch
output(s):
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
                based on W/m2 of DNI after the glazing and optical losses.

                It would be better to use the simulated electrical 
                and thermal generation as inputs and solve
                but it wouldn't be too innacurate to just assume 
                constant conversion percentages
                Egen 10.5% conversion of DNI (W/m2) to Egen (W/m2)
                Qgen 13.0% conversion of DNI (W/m2) to Qgen (W/m2)
                Could use the below equations for better accuracy 
                but they account for array of modules
                output of this should be W/m2 for now but can change 
                (if geometry is fed in) to W for specific volumes
                Egen=((55)*-2.28e-4+0.1838)*DNIafterExt/1000,
                Qgen=0.024801*(0.66)*1e-3*DNIafterExt

input(s):   
output(s):  
"""
def getRadiativeGainAtTime(time, geometry, DNI, DHI):
  return getRadiativeGain(getSunPosition(time),geometry,DNI,DHI)

def getRadiativeGain(sunPositions, geometry, DNI, DHI):
  Egen = 0.105*DNI
  Qgen = 0.13*DNI
  (pitch, yaw) = solar.getArrayAngle(sunPositions,geometry)

  glazing = solar.getGlazingTransmitted(sunPositions,geometry,2)
  # Calculate the DNI to interior using gap and double pane 
  # glazing transmittance
  DNItoInterior = DNI * glazing * getGapTransmittance(yaw, pitch)
  # DNI that is not converted to electricity energy and thermal energy
  # or transmitted to interior is captured as heat inside the cavity
  DNIheat = DNI - Egen - Qgen - DNItoInterior

  # Calculate the DHI to interior 
  # Assumption: Diffuse radiation is not intercepted by ICSolar like direct
  # because it comes at different angles, allowing passage 
  # through transparent materials
  DHItoInterior = DHI * glazing
  # DHI that is not transmitted to the interior is captured at heat inside cavity
  DHIheat = DHI - DHItoInterior

  # Radiative heat gain to add to the cavity air temperature
  CavRadHeatGain = DNIheat + DHIheat
  # Radiative heat gain to add to the interior heat transfer
  IntRadHeatGain = DNItoInterior + DHItoInterior

  # Units are W/m2, will need to multiply by 
  # surface area to determine the heat added to a specific volume.
  return (CavRadHeatGain, IntRadHeatGain)

