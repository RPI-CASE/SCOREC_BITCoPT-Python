def readLightingSchedule(cosimDirection, filename):
  f = open(filename,'r')
  lines = f.readlines()
  f.close()

  data = {}
  # data['city'] = lines[0].split(',')[1]
  # data['timezone'] = float(lines[0].split(',')[3])
  data['lat'] = float(lines[1].split(',')[8])

  # Load lighting fraction for each co-simulation direction
  for direction in cosimDirection:
    # String name for header
    header = str(direction).lower()+'. LightingFraction'

    # Get index number in first row of the lighting fraction file 
    # that corresponds to the header string
    fractionIndex = lines[0].split(',').index(header)

    # Building energy model uses the name core but the surface is called roof
    if str(direction).lower() == 'roof':
      direction = 'Core'
    # Lighting Fractions
    data[direction] = [float(lines[i].split(',')[fractionIndex]) 
    for i in range(1,len(lines))]

  return data