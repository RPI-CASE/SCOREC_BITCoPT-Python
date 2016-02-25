"""
casegeom.py contains case specific functions to build off
of the base infrastructure, such as reading file types or
writing problem specific output, and any other case specific
geometry handling

functions:

readFile
writeVTKFile
readValidationFile
createSimpleGeometry
"""
import sys
import csv
import numpy as np
import base.geom as g
from collections import defaultdict

"""
readFile:       reads from the original file format,
                with coordinates in {x0,y0,z0} {x1,y1,z1}  ...
                on a single line

                each line is  a single facade, see data/geometry/whole-building.txt

input(s):       filename, geometry name
output(s):      GeometrySet object
"""
def readFile(filename,name,useSunlit = True):
  geometrySet = g.GeometrySet(name)
  f = open(filename,'r')
  surfaces = [line.rstrip('\n') for line in f.readlines() if line[0] is '{']
  f.close()
  for s in surfaces:
    c = [d.lstrip('{').lstrip(' {').split(',') \
      for d in s.split('}') if d is not '']
    if(len(c) != 4):
      sys.exit("attempted to import surface with "+str(len(c))+" coordinates\n")
    coords = [np.array(coord,float) for coord in c]
    for coord in coords:
      if(len(coord) != 3):
        sys.exit("coordinate is invalid\n")
    geom = g.Geometry(coords)
    geometrySet.append(geom)

  geometrySet.computeSubsets(useSunlit)
  geometrySet.computeMatches(useSunlit)

  return geometrySet

"""
readNewFile:    reads from the new file format,
                with much more information, including walls

                each line is  a single facade, see data/geometry/whole-building.txt

input(s):       filename, geometry name
output(s):      GeometrySet object
"""
def readNewFile(filename,name,useSunlit = True):
  geometrySet = g.GeometrySet('name')

  f = open(filename,'r')
  lines = f.readlines()
  f.close()
  index = 0
  data = {}
  # loop through all the data
  while index < len(lines):
    header = lines[index][0:lines[index].find(',')]
    # if the header is already found (main line) 
    # then add another entry, otherwise initialize it
    if header in data:
      data[header].append([])
    else:
      data[header] = [[]]
    # loop until you hit an empty line, denoting a new entry
    while index < len(lines) and lines[index].strip():
      # remove the endline and initial tabs that are present
      data[header][-1].append(lines[index].rstrip('\n').lstrip('\t').lstrip(' '))
      index = index+1
    # skip the blank line
    index = index+1

  # we now have all the data in a nice header
  # first lets process the geometry info
  # this is kind of crude, but seems to work
  geometry_list = data['GlobalGeometryRules'][0]
  for item in geometry_list:
    if ( 'Starting Vertex Position' in item):
      starting_vertex_position = item[0:item.find(',')]
    if ( 'Vertex Entry Direction' in item):
      vertex_entry_direction = item[0:item.find(',')]

  if starting_vertex_position == 'UpperLeftCorner':
    coord_offset = 2 # my starting corner is bottom right

  if vertex_entry_direction == 'Counterclockwise':
    pass # this is fine, its the same orientation I use

  building_surface_list = data['BuildingSurface:Detailed']
  walls = {}
  for surface in building_surface_list:

    name = surface[1].split(',')[0]
    surface_type = surface[2].split(',')[0]
    # ignore all surfaces not named "Wall"
    if surface_type != 'Wall':
      continue
    coords = [[a[:-2] for a in c.split(' ')[0:3]] for c in surface[-4:]]
    coords = [coords[(i + coord_offset) % 4] for i in range(len(coords))]
    coords = [np.array(coord,float) for coord in coords]
    geom = g.Geometry(coords,'wall')
    walls[name] = geom
    geometrySet.append(geom)

  window_list = data['FenestrationSurface:Detailed']
  for surface in window_list:
    name = surface[1].split(',')[0]
    surface_type = surface[2].split(',')[0]
    wall_name = surface[4].split(',')[0]
    coords = [[a[:-2] for a in c.split(' ')[0:3]] for c in surface[-4:]]
    coords = [coords[(i + coord_offset) % 4] for i in range(len(coords))]
    coords = [np.array(coord,float) for coord in coords]
    geom = g.Geometry(coords)

    if wall_name in walls:
      geom.setWall(walls[wall_name])
    else:
      geom = g.Geometry(coords)
    geometrySet.append(geom)

  geometrySet.computeSubsets(useSunlit)
  geometrySet.computeMatches(useSunlit)

  # for i in range(len(geometrySet)):
  #   print i, [geometrySet.index(s) for s in geometrySet.s[i]], geometrySet[i].dir

  # g.plot(geometrySet,numbers = True, normals = True)
  return geometrySet
"""
writeVTKFile:   writes a paraview VTK file of the data stored in the geometry

input(s):       geometrySet:
                filename:
                note:         a short string to describe the file
                nX:           the number of cells in the x-direction
output(s):      none
"""

def writeVTKfile(geometrySet, filename, note, nX=1):
  numCells = sum([g.nY*nX for g in geometrySet])
  numPoints = 4*numCells
  f = open(filename,'w')
  f.write("# vtk DataFile Version 3.0\n")
  f.write("CASE " + note+"\n")
  f.write("ASCII\n")
  f.write("DATASET UNSTRUCTURED_GRID\n")
  f.write("POINTS "+str(numPoints)+str(" float\n"))

  c = [[],[],[],[]]
  for g in geometrySet:
    dX = (g[3]-g[0])/nX
    dY = (g[1]-g[0])/g.nY
    for x in range(nX):
      for y in range(g.nY):
        c[0] = dX*x    +dY*y
        c[1] = dX*(x+1)+dY*y
        c[2] = dX*(x+1)+dY*(y+1)
        c[3] = dX*x    +dY*(y+1)
        for p in c:
          f.write(str(p[0]+g[0][0])+" "\
            +str(p[1]+g[0][1])+" "\
            +str(p[2]+g[0][2])+"\n")

  f.write("CELLS "+str(numCells)+" "+str(5*numCells)+"\n");
  for i in range(numCells):
    f.write("4 "+" ".join([str(j) for j in range(4*i,4*(i+1))])+"\n")

  f.write("CELL_TYPES "+str(numCells)+"\n")
  for i in range(numCells):
    f.write("9\n")

  # this is all the data
  geometrySet.collectDataNames()
  if(len(geometrySet.dataNames)) > 0:
    f.write("CELL_DATA "+str(numCells)+"\n")
  for dataName in geometrySet.dataNames:
    f.write("SCALARS "+dataName+" float 1\n")
    f.write("LOOKUP_TABLE default\n")
    for A in geometrySet:
      m = geometrySet.getMatches(A)
      if(len(m) > 0 and geometrySet.index(m[0]) < geometrySet.index(A)):
        for j in range(nX):
          for i in range(A.nY):
            f.write(str(geometrySet[geometrySet.index(m[0])].data[dataName][i])+"\n")
      else:
        for j in range(nX):     
          for i in range(A.nY):
            f.write(str(A.data[dataName][i])+"\n")
  f.close()
"""
readValidationFile: reads a csv file with validation data

input(s):       filename
output(s):      data, an ordered dictionary of data
"""
def readValidationFile(filename):
  csvfile = open(filename,'rU')
  cr = csv.DictReader(csvfile)
  # read in all the data
  data = defaultdict(list) 
  for row in cr: # read a row as {column1: value1, column2: value2,...}
    for (k,v) in row.items(): # go over each column name and value 
      data[k].append(float(v))

  return data

def createSimpleGeometry(orientation,tilt):
  return g.SimpleGeometry(orientation,tilt)

def tiltGeometry(geometry,tilt,useSunlit):
  if type(tilt) is int or type(tilt) is float:
    for geom in geometry:
      if geom.dir is not 'roof':
        geom.rotate(geom[2]-geom[1],geom[1],tilt*np.pi/180.)
  else:
    for i in range(len(tilt)):
      geom = geometry[i]
      geom.rotate(geom[2]-geom[1],geom[1],tilt[i]*np.pi/180.)
  geometry.computeSubsets(useSunlit)
  geometry.computeMatches(useSunlit)

