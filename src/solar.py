"""
solar.py contains solar calculations required for facade calculations

latitude, longitude, and timezone can be set once for all calculations
using setTimezone and setLocation

Two versions for most functions, to avoid computing sun positions if 
we already have them at a specified time


functions:
secondsToDatetime,  convert seconds since beginning of 2015 to datetime
getSunPosition,     time, lat, long to sun positions
getAOI,             AOI based on sun positions and surface orientation and tilt
getSolarVector,     convert altitude and azimuth to vector
getArrayAngle,      system pitch and yaw response to sunPositions and surface
getGlazingTransmitted,  transmitted solar radiation based on AOI

"""

import base.geom as g
import numpy as np
from sys import exit
from datetime import datetime, timedelta

try:
  import Pysolar.solar as solar
except ImportError:
  exit("pysolar library required for calculations\n \
    https://github.com/pingswept/pysolar\n")
try:
  from pytz import timezone, utc
except ImportError:
  exit("pytz library required for timezone calculations\n")

# store these constants here
# Projective window properites:
xLite = 6 #Thickness of lites (layers) in mm
nAir = 1.0 #Optical index, air
nGlass = 1.53 #Optical index, glass
Rsfc = 0.00001 #Surface soiling coefficient, lower number = less dirty
cDisp = 0.0075 #Coefficient of dispersion, Ultrawhite=0.0075.

# these are global variables, set to NYC studio by default
_latitude = 40.707583
_longitude = -74.010828
_utctimezone = utc
_timezonehours = -5

"""
  allows every function to use the same timezone and
  latitude and longitude
"""
def setTimezone(tz):
  global _timezonehours
  # Manually set timezone hours to bypass daylight savings time
  # DST messes up solar tracking for the time being
  _timezonehours = tz

def setLocation(lat,lon):
  global _latitude
  global _longitude
  _latitude = lat
  _longitude = lon

"""
secondsToDatetime:  Converts seconds into datetime starting at 2015.
                    Datetime needed for solar calculations.

input(s):   time (seconds), since beginning of 2015
output(s):  datetime
"""

def secondsToDatetime(secondsOfYear):
  ################################## Daylight Savings Off
  # Daylight savings does not have an effect on this section. 
  # This is manual at the time because 
  # I'm still debugging why DST doesnt help the results. 
  localDatetime = datetime(2015,1,1,0,30,0)+timedelta(seconds=secondsOfYear)
  globalDatetime = localDatetime - timedelta(hours=_timezonehours)
  ################################## Daylight Savings On
  # # This time is daylight savings aware and will automatically adjust
  # dt = datetime(2015,1,1,0,30,0) + timedelta(seconds=secondsOfYear)
  # localDatetime = _timezone.localize(dt)
  # globalDatetime = localDatetime.astimezone(_utctimezone)
  
  return (localDatetime, globalDatetime)

"""
getSunPosition:   Take seconds of the year, latitude, and longitude,
                  convert to solar positions for given location.

input(s):   time (seconds),   since beginning of 2015
            latitude (degrees)
            longitude (degrees)

output(s):  altitude (rads),  height of sun from horizon
            azimuth (rads),   position of sun (east to west) with West positive, South=0
            declination (rads), angular difference north or south of the equator
            hourAngle (rads), angular distance between celestial sphere and a point
"""
def getSunPosition(time):
  # Only issue at the moment with this function is daylight savings 
  # but this will be addressed later.
  (localDatetime, globalDatetime) = secondsToDatetime(time)
  # Use pysolar to calculate altitude and azimuth
  altitude = solar.GetAltitude(_latitude, _longitude, globalDatetime) * np.pi/180
  # Pysolar considers east positive and west negative
  # reversed to follow the opposite convention.
  azimuth = -solar.GetAzimuth(_latitude, _longitude, globalDatetime) * np.pi/180
  if(azimuth < -np.pi):
    azimuth += 2*np.pi
  elif (azimuth > np.pi):
    azimuth -= 2*np.pi
  # Calculate the solar declination from the day of the year
  declination = solar.GetDeclination(solar.GetDayOfYear(globalDatetime)) * np.pi/180
  hourAngle = -solar.GetHourAngle(globalDatetime, _longitude) * np.pi/180

  return {'altitude':altitude, 'azimuth':azimuth, \
    'declination':declination, 'hourAngle':hourAngle}

"""
getAOI: Convert time, geometry of surfaces, and latitude to
        angle of incidence used for determining glazing losses.

input(s):   time (seconds)
            geometry (Geometry Object)
            latitude (rads)

output(s):  AOI (rads), Angle of Incidence, approach of a ray to a surface
"""
def getAOIAtTime(time,geometry):
  return getAOI(getSunPosition(time),geometry)

def getAOI(sunPositions,geometry):
  lat = _latitude*np.pi/180.

  del_c = np.cos(sunPositions['declination'])
  del_s = np.sin(sunPositions['declination'])
  phi_c = np.cos(lat)
  phi_s = np.sin(lat)
  omg_c = np.cos(sunPositions['hourAngle'])
  omg_s = np.sin(sunPositions['hourAngle'])
  bet_c = np.cos(np.pi / 2. - geometry.tilt)
  bet_s = np.sin(np.pi / 2. - geometry.tilt)
  gam_c = np.cos(-geometry.orient)
  gam_s = np.sin(-geometry.orient)

  return np.arccos(del_s * phi_s * bet_c - \
      del_s * phi_c * bet_s * gam_c + \
      del_c * phi_c * bet_c * omg_c + \
      del_c * phi_s * bet_s * gam_c * omg_c + \
      del_c * bet_s * gam_s * omg_s)

def getAOICheck(sunPositions,orient,tilt):
  lat = _latitude*np.pi/180.

  del_c = np.cos(sunPositions['declination'])
  del_s = np.sin(sunPositions['declination'])
  phi_c = np.cos(lat)
  phi_s = np.sin(lat)
  omg_c = np.cos(sunPositions['hourAngle'])
  omg_s = np.sin(sunPositions['hourAngle'])
  bet_c = np.cos(np.pi / 2. - tilt)
  bet_s = np.sin(np.pi / 2. - tilt)
  gam_c = np.cos(orient)
  gam_s = np.sin(orient)

  return np.arccos(del_s * phi_s * bet_c - \
      del_s * phi_c * bet_s * gam_c + \
      del_c * phi_c * bet_c * omg_c + \
      del_c * phi_s * bet_s * gam_c * omg_c + \
      del_c * bet_s * gam_s * omg_s)

"""
getSolarVector: Convert sun positions into vectors

input(s):   time (seconds)
            geometry (Geometry Object)

output(s):  Solar vector
"""
def getSolarVectorAtTime(time):
  return getSolarVector(getSunPosition(time))

def getSolarVector(sunPositions):
  # Transformation Matrix
  phi = 3.*np.pi/2. - sunPositions['azimuth']
  # Solar zenith angle
  theta = np.pi/2. - sunPositions['altitude'] 

  # Result of Cartesian to Spherical
  return np.array([np.sin(theta)*np.cos(phi),
                   np.sin(theta)*np.sin(phi),
                   np.cos(theta)])
"""
arrayAngle: Convert time, surface geometry, and sun positions to
            system pitch and yaw for solar tracking and self shading.

input(s):   time (seconds)
            geometry (Geometry Object)
            sunPositions (dictionary)

output(s):  pitch (rads), array's north and south response to altitude, orientation, and tilt
            yaw (rads), array's east and west response to azimuth, orientation, and tilt
"""
def getArrayAngleAtTime(time,geometry):
  return getArrayAngle(getSunPosition(time),geometry)

def getArrayAngle(sunPositions,geometry):
  s = getSolarVector(sunPositions)
  sr = s*geometry.Rsolar
  sr = np.array([sr[0,0],sr[0,1],sr[0,2]])
  sr = sr/np.linalg.norm(sr)
  yaw = -(np.arctan2(sr[1],sr[0])+np.pi/2.)
  if(yaw < -np.pi):
    yaw += 2*np.pi
  elif (yaw > np.pi):
    yaw -= 2*np.pi
  pitch = np.pi/2.-np.arctan2(np.sqrt(sr[1]*sr[1]+sr[0]*sr[0]),sr[2])
  # print geometry.dir,geometry._index, pitch*180/np.pi, yaw*180/np.pi
  
  return (pitch, yaw)

"""
getMeanArrayAngle: gets an averaged pitch/yaw based on equally spaced sampling

input(s):   time (seconds)
            geometry (Geometry Object)
            timestep (seconds)
            n (number of samples (actually 2*n+1 samples))

output(s):  (meanPitch,meanYaw)
"""
def getMeanArrayAngle(time,geometry,timestep,n):
  pitchArray = np.zeros(2*n+1)
  yawArray = np.zeros(2*n+1)
  dt = timestep/(2*n)
  for i in range(-n,n+1):
    (pitch,yaw) = getArrayAngleAtTime(time+i*dt,geometry)
    pitchArray[i+n] = pitch
    yawArray[i+n] = yaw

  return (np.mean(pitchArray),np.mean(yawArray))
"""
getGlazingTransmitted:  Convert sunPositions, geometry, and layers of glazing
                        through the Schlick glazing losses calculation
                        to determine the fraction of radiation that passes 
                        through the glazing. 

input(s):   time (seconds)
            sunPositions (dictionary)
            geometry (Geometry Object)
            layers (integer), 1 layer for exterior, 2 layers for interior

output(s):  transGlazingLosses (fraction), what amount of the sun passes through glazing
"""
def getGlazingTransmittedAtTime(time,geometry,layers):
  return getGlazingTransmitted(getSunPosition(time),geometry,layers)

def getGlazingTransmitted(sunPositions,geometry,layers):
  AOI = getAOI(sunPositions,geometry)
  # Normal reflection (Snell's)
  Ro = ((nAir-nGlass) / (nAir+nGlass))**2 
  # Fresnel reflectance from one surface
  RFres = Ro + (1-Ro) * (1-np.cos(AOI))**5

  # Schlick glazing losses calculation
  transGlazingLosses = (1 - Rsfc / np.cos(AOI)) * \
    ((1 - cDisp * xLite / np.cos(AOI)) * \
    (1 - 2 * RFres / (1 + RFres)))**layers
  
  # this can return negative numbers, but in every case
  # that occurs, its multiplied by zero later on, so who cares
  return transGlazingLosses