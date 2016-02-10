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
