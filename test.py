"""
This will serve as a place to run unit tests, 
as this is becoming more necessary,

Basic components of the code can be tested here
The idea is that this should be very quick to run, 
and can be included the main run files, and run 
prior to every simulation
"""

import poisson2D
import diffusion2D
import src.base.geom as g
import src.solar as solar
import src.casegeom as cg
import src.shading as shading
import numpy as np
import time as clocktime
import src.icsolar as icsolar

def runGeometryTest():

  latitude = 40.707583
  longitude = -74.010828
  timezone = -5
  solar.setTimezone(timezone)
  solar.setLocation(latitude,longitude)
  geometry = cg.readFile('data/geometry/whole-building.txt',"whole-building")
  geom = geometry[13]
  nS = -shading.getSolarVectorAtTime(11.5*3600)
  assert abs(geometry.getUnshadedFraction(geom,nS) - 0.5097854) < 1e-6, \
    "failed sunlit fraction test"
  geometrySet2 = g.GeometrySet('tilt')
  geometrySet2.append(geom)
  oldg = g.Geometry([np.copy(c) for c in geom])
  geom.rotate(np.array([1,1,1]),geom[0],0)
  for i in range(4):
    assert np.abs(np.linalg.norm(oldg[i]-geom[i])) < 1e-6, \
      "failed rotation test"
  geom.rotate(geom[2]-geom[1],geom[1],-90./180.*np.pi)
  geometrySet2.append(oldg)

  assert np.abs(np.linalg.norm(oldg[1]-geom[1])) < 1e-6, \
    "failed rotation test"
  assert np.abs(np.linalg.norm(oldg[2]-geom[2])) < 1e-6, \
    "failed rotation test"  
  assert np.abs(np.dot(oldg[1]-oldg[0],geom[1]-geom[0])) < 1e-6, \
    "failed rotation test"      

  # geometry2 = cg.readFile('data/geometry/tilt20.txt',"test20")
  # cg.writeVTKfile(geometry2,'tilt20.vtk','tilt20')

def runSolverTests():
  rate = poisson2D.test()
  assert (rate[0] > 2.3 and rate[1] > 2.3), \
    "testing failed poisson2D: solver error"
  rate = diffusion2D.test()
  assert (rate[0] > 2.3 and rate[1] > 2.3), \
    "testing failed diffusion2D: solver error"

def runSolarTest():
  latitude = 40.707583
  longitude = -74.010828
  timezone = -5
  solar.setTimezone(timezone)
  solar.setLocation(latitude,longitude)
  #(localTime, globalTime) = solar.secondsToDatetime(12*3600)
  earlyAngles = solar.getSunPosition(10*3600)
  noonAngles = solar.getSunPosition(12*3600) 
  lateAngles = solar.getSunPosition(14*3600) 
  summerAngle = solar.getSunPosition(171*24*3600+12*3600)
  winterAngle = solar.getSunPosition(1*24*3600+12*3600)
  assert (earlyAngles['azimuth'] < noonAngles['azimuth'] < lateAngles['azimuth']), \
    "Unexpected: testing azimuth angles - found azimuth angles to be reverse/incorrect"
  assert (earlyAngles['altitude'] < noonAngles['altitude'] > lateAngles['altitude']), \
    "Unexpected: testing altitude - the sun is not highest at noon"
  assert (winterAngle['declination'] < summerAngle['declination']), \
    "Unexpected: testing delication - sun is not higher in summer than winter"

  geometry = cg.readFile('data/geometry/whole-building.txt',"whole-building")
  geom = geometry[13]
  # JS: I do not like this test because as we change geometry it will break but it's here. 
  # assert abs(solar.getAOI(noonAngles,geom) - solar.getAOICheck(noonAngles,0.0,(60*np.pi/180.))) < 1e-6, \
  #   "Unexpected - Angle of incidence on roof has changed"
  (pitchSummer, yawSummer) = solar.getArrayAngle(summerAngle,geom)
  (pitchWinter, yawWinter) = solar.getArrayAngle(winterAngle,geom)
  assert (pitchWinter < pitchSummer), \
    "Unexpected: testing pitch - the winter pitch is higher than summer"

  assert(solar.getGlazingTransmitted(earlyAngles,geom,1) < \
    solar.getGlazingTransmitted(noonAngles,geom,1) > \
      solar.getGlazingTransmitted(lateAngles,geom,1)), \
    "Unexpected: testing glazing transmittance - was not greatest at noon for the roof"

  

def runTests():
  start = clocktime.time()
  print "running tests ...",
  runGeometryTest()
  runSolverTests()
  runSolarTest()
  print "all passed in", '%.2f' % (clocktime.time()-start),"seconds"
if __name__ == '__main__':
  runTests()