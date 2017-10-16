def readLightingSchedule(cosimDirection, filename):
  f = open(filename,'r')
  lines = f.readlines()
  f.close()

  data = {}
  # data['city'] = lines[0].split(',')[1]
  # data['timezone'] = float(lines[0].split(',')[3])
  data['lat'] = float(lines[1].split(',')[15])

  # Load lighting fraction for each co-simulation direction
  for direction in cosimDirection:
    if str(direction).lower() == 'core':
      direction = 'roof'
    header = str(direction).lower()+'.InteriorLightingFraction'

    # Get index number in first row of the lighting fraction file 
    # that corresponds to the header string
    fractionIndex = lines[0].split(',').index(header)

    # Lighting Fractions
    data[direction] = [float(lines[i].split(',')[fractionIndex]) 
    for i in range(1,len(lines))]

  # southIndex = lines[0].split(',').index('south.LightingFraction')
  # eastIndex = lines[0].split(',').index('east.LightingFraction')
  # westIndex = lines[0].split(',').index('west.LightingFraction')
  # roofIndex = lines[0].split(',').index('roof.LightingFraction')

  # # South Lighting Fractions
  # data['southLightingFraction'] = [float(lines[i].split(',')[southIndex]) 
  #   for i in range(1,len(lines))]
  # # East Lighting Fractions
  # data['eastLightingFraction'] = [float(lines[i].split(',')[eastIndex]) 
  #   for i in range(1,len(lines))]
  # # West Lighting Fractions
  # data['westLightingFraction'] = [float(lines[i].split(',')[westIndex]) 
  #   for i in range(1,len(lines))]
  # # Roof Lighting Fractions
  # data['roofLightingFraction'] = [float(lines[i].split(',')[roofIndex]) 
  #   for i in range(1,len(lines))]

  return data