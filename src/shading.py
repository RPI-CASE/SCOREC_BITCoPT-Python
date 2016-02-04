"""
shading.py contains the shading matrices for modules shading modules
as well as other shading support

functions:
getShadingIndex
getUnshadedFraction
getStudioShadingIndex
getStudioUnshadedFraction
applySunlitFraction
getMeanSunlitFraction
"""

import base.geom as g
import numpy as np
import time as cputime
from solar import *
from os.path import isfile
from sys import exit

indices = [11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,\
           44,45,51,52,53,54,55]

for index in indices:
  if (not isfile('data/shading/5x5/'+str(index)+'.txt')):
    exit('warning: shading tables not found in data/shading/5x5\n')

shading = [np.loadtxt('data/shading/5x5/'+str(index)+'.txt') \
  for index in indices]

"""
Two different sets of studio shading tables
"""
# studioShading = [np.loadtxt('data/shading/validation/'+str(index)+'.txt') \
#   for index in range(12)] 
for index in range(0,6):
  if (not isfile('data/shading/studio/module'+str(index+1)+'.txt')):
    exit('error shading tables not found in data/shading/studio\n')

studioShading = [np.loadtxt('data/shading/studio/module'+str(index+1)+'.txt') 
  for index in range(0,6)] 

_shadingTime = 0.

"""
getShadingIndex: 
            given a module on a facade with nX by nY modules
            figure out what shading matrix to use

55 54 53 52 51
45 44 43 42 41
35 34 33 32 31
25 24 23 22 21
15 14 13 12 11


input(s):   i,j, location of the module (from bottom left)
            nX,nY, geometry specific
output(s):  integer
"""
def _getShadingIndex(i,j,nX,nY):
  # interior
  m = nX - 1
  n = nY - 1

  if (nX == 1):
    ind = 3
  elif(i == 0 or i == 1):
    ind = i+1
  elif(i == m-1 or i == m):
    ind = 5-(m-i)
  else:
    ind = 3
  
  if(j == 0):
    return 10+ind
  elif(j == 1):
    return 20+ind
  elif(j == n-1):
    return 40+ind
  elif(j == n): 
    return 50+ind
  else:
    return 30+ind
  return 33

def getShadingIndex(i,j,nX,nY):
  return indices.index(_getShadingIndex(i,j,nX,nY))
"""
getUnshadedFraction: 
            This uses the generic 5x5 shading matrices
            and bilinearly interpolates the data to get a shading
            factor at each point

input(s):   index, integer from 0 to 24
            geometry, geometry object 
output(s):  float, between 0 and 1
"""
def getUnshadedFractionAtTime(time,geometry,index):
  return getUnshadedFraction(getSunPosition(time),geometry,index)

def getUnshadedFraction(sunPositions,geometry,index):
  (pitch, yaw) = getArrayAngle(sunPositions,geometry)
  yaw *= 180.0/np.pi
  pitch *= 180.0/np.pi

  if (abs(yaw) > 71.99 or abs(pitch) > 71.99):
    return 0.0

  iPitch = 24-pitch/3.
  iYaw = 24-yaw/3.

  i = np.floor(iPitch)
  j = np.floor(iYaw)
  x = iPitch - i
  y = iYaw - j
  # x in [-72,72], y in [-72,72] 
  # top right is [72,72]
  return shading[index][i][j]*(1-x)*(1-y) + \
  shading[index][i+1][j]*x*(1-y) + \
  shading[index][i][j+1]*(1-x)*y + \
  shading[index][i+1][j+1]*x*y

"""
These are the same as above, for validation, so duplicate code
is acceptable in this case, and easier

11 5
10 4
9 3
8 2
7 1
6 0
input(s):         i = 0,1, j = 0,5
output(s):        float, between 0 and 1
"""
def getStudioShadingIndex(i,j):
  return (1-i)*6+j

def getStudioUnshadedFractionAtTime(time,geometry,index):
  return getStudioUnshadedFraction(getSunPosition(time),geometry,index)

def getStudioUnshadedFraction(sunPositions,geometry,index):
  (pitch, yaw) = getArrayAngle(sunPositions,geometry)

  yaw *= 180.0/np.pi
  pitch *= 180.0/np.pi

  if (abs(yaw) > 71.99 or abs(pitch) > 71.99):
    return 0.0
  iPitch = 24-pitch/3.
  iYaw = 24-yaw/3.

  i = np.floor(iPitch)
  j = np.floor(iYaw)
  x = iPitch - i
  y = iYaw - j

  # x in [-72,72], y in [-72,72] 
  # top right is [72,72]
  return studioShading[index][i][j]*(1-x)*(1-y) + \
  studioShading[index][i+1][j]*x*(1-y) + \
  studioShading[index][i][j+1]*(1-x)*y + \
  studioShading[index][i+1][j+1]*x*y

"""
applySunlitFraction:  applies the sunlit fraction to an array 
                      in the geometry

input(s):   sunlit fraction
            geometry object
            averaged, a boolean describing how to apply it 
            (whether to apply it from one end or to average across)

output(s):  array of length g.nY
"""
def applySunlitFraction(sunlit,geometry,averaged):
  nY = geometry.nY
  if sunlit > 1.-1e-10:
    return np.ones(nY)
  if sunlit < 1e-10:
    return np.zeros(nY)
  if averaged:
    return np.ones(nY)*sunlit
  numShaded = int((1.-sunlit)*nY)
  fracRemaining = (1.-sunlit)*nY-numShaded
  shading = np.ones(nY)

  if(geometry.tilt > 2.*np.pi/7.): #bottom shaded
    for j in range(numShaded):     
      shading[j] = 0.
    shading[numShaded] *= (1.-fracRemaining)
  else: #top shaded
    for j in range(nY-numShaded,nY):
      shading[j] = 0.
    shading[nY-numShaded-1] *= (1.-fracRemaining)

  return shading

"""
calculateAverageSunlitFraction: 
            calculates the average sunlit fraction over a time,
            based on 2n+1 points

input(s):   geometrySet, complete set of geometries
            geometry, individual geometry
            time, exact time
            timestep, length of period to average over
            n, number of samples to average (done over timestep/2 on either 
              side of the time)
output(s):  float between 0 and 1
            boolean telling whether its less than 1 because of shading, 
              or because there was a portion of time where it wasnt seen
"""  
def getMeanSunlitFraction(geometrySet,
  geometry,time,timestep,n):
  clockTime = cputime.time()
  sunlitArray = np.zeros(2*n+1)
  dt = timestep/(2*n)
  for i in range(-n,n+1):
    nS = -getSolarVectorAtTime(time+i*dt)
    sunlitArray[i+n] = geometrySet.getUnshadedFraction(geometry,nS)
  
  count = 0;
  if(not geometrySet.s[geometrySet.index(geometry)]):
    averaged = False    
  else:
    averaged = True
  global _shadingTime
  _shadingTime += cputime.time() - clockTime

  return (np.mean(sunlitArray),averaged)

"""
  returns time in seconds used for shading
"""
def getShadingTime():
  return _shadingTime
