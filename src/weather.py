import numpy as np

"""
weather.py contains the reader for TMY files and other
weather related functions

functions:
getData
"""

"""
getData:  gets DNI, air temperature, and DHI

input(s):   filename
output(s):  DNI, air temperature, DHI
"""
def readTMY(filename):
  f = open(filename,'r')
  lines = f.readlines()
  f.close()

  data = {}
  data['city'] = lines[0].split(',')[1]
  data['timezone'] = float(lines[0].split(',')[3])
  data['lat'] = float(lines[0].split(',')[4])
  data['lon'] = float(lines[0].split(',')[5])

  # Get index number in second row of the weather file 
  # that corresponds to the header string
  AIRindex = lines[1].split(',').index('Dry-bulb (C)')
  DNIindex = lines[1].split(',').index('DNI (W/m^2)')
  DHIindex = lines[1].split(',').index('DHI (W/m^2)')
  # Direct Normal Irradiance
  data['DNI'] = [float(lines[i].split(',')[DNIindex]) 
    for i in range(2,len(lines))]
  # Dry-Bulb Air Temperature
  data['airTemp'] = [float(lines[i].split(',')[AIRindex]) 
    for i in range(2,len(lines))]
  # Diffuse Horizontal Irradiance
  data['DHI'] = [float(lines[i].split(',')[DHIindex]) 
    for i in range(2,len(lines))]

  return data

def readEPW(filename):
  # Same as read TMY but for an EPW format
  # Read header for location info
  f = open(filename,'r')
  firstLine = f.readline()
  f.close()

  data = {}
  data['city'] = firstLine.split(',')[1]
  data['lat'] = float(firstLine.split(',')[6])
  data['lon'] = float(firstLine.split(',')[7])
  data['timezone'] = float(firstLine.split(',')[8])

  # Read columns for weather data and assign dictionary names
  wd = np.genfromtxt(filename, delimiter=',', skip_header=8, names=('year','month','day','hour','x1','x2','drybulb','dewpoint','relhum','atmospress','exthorzrad','extdirrad','horzIRsky','globhorzrad','dirnormrad','difhorzrad','globhorzillum','dirnormillum','difhorzillum','zenlum','windir','windspd','totskycvr','opaqskycvr','visibility','ceilinghgt','presweatobs','preswethcodes','precipwtr','aerosol','snowdepth','daylastsnow','albedo','rain','rainquant'))

  # Direct Normal Irradiance
  data['DNI'] = wd['dirnormrad']
  # Dry-Bulb Air Temperature
  data['airTemp'] = wd['drybulb']
  # Diffuse Horizontal Irradiance
  data['DHI'] = wd['difhorzrad']

  return data

def readTMYUQ(filename):
  f = open(filename,'r')
  lines = f.readlines()
  f.close()
  data = {}
  data['city'] = lines[0].split(',')[1]
  data['timezone'] = float(lines[0].split(',')[3])
  data['lat'] = float(lines[0].split(',')[4])
  data['lon'] = float(lines[0].split(',')[5])

  AIRindex = lines[1].split(',').index('Dry-bulb (C)')
  DNIindex = lines[1].split(',').index('DNI (W/m^2)')
  DNIUQindex = lines[1].split(',').index('DNI uncert (%)')

  DHIindex = lines[1].split(',').index('DHI (W/m^2)')
  # Direct Normal Irradiance
  data['DNI'] = [float(lines[i].split(',')[DNIindex]) 
    for i in range(2,len(lines))]
  data['DNIUQ'] = [float(lines[i].split(',')[DNIUQindex]) 
    for i in range(2,len(lines))]
  # Dry-Bulb Air Temperature
  data['airTemp'] = [float(lines[i].split(',')[AIRindex]) 
    for i in range(2,len(lines))]
  # Diffuse Horizontal Irradiance
  data['DHI'] = [float(lines[i].split(',')[DHIindex]) 
    for i in range(2,len(lines))]

  return data
