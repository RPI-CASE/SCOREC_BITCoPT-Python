""" 
geom.py contains the Geometry and GeemetrySet classes

Each geometry represents a planar, rectangular facade,
given by a set of coordinates in counter-clockwise order,
such that normals point outwards and the first point is on 
the bottom right.

The class itself acts as a wrapper around the coordinates,
and coordinates can be accessed simply by iterating over the
object itself, and standard get/set/len type operations.

Quantities that are used often are stored, and quantities that
are not required often, are calculations for now.

Requires the SHAPELY library to handle facades shading other facades

Each Geometry Object has:
  (.n) normal vector, as an np.array
  (.orient) orientation angle, in radians
  (.tilt) tilt angle, in radians
  (.coords) coordinates
  (._index) a number representing its location in the GeometrySet
            let the GeometrySet deal with this, dont access directly
  (.nX) the number of blocks in the X direction
  (.nY) the number of blocks in the Y direction
  (.dir) (NSEW/Roof)
  (.data) a dictionary of data
  (.R) rotation matrix, transposed for ease of use,
       transpose(R(tilt)*R(orient))
  (.poly) Polygon Representation of the surface in x-z plane

Each GeometrySet Object is a collection of Geometry Objects,
in an easily accessible way, allowing for queries
it is more or less, an overloaded list

geometry objects do not know what 'index' they have
or what other geometry objects shade them, this is all
stored in the GeometrySet

Each GeometrySet Object has:
  (.g) a list of geometry objects
  (.s) a list of subsets for each geometry object
  (.m) a list of matches, similar geometries, 
       for each geometry object
  (.name) a name
  (.dataNames) a list of data stored in the geometry

There is also a plot function, to plot geometrySets, see below

There is also a SimpleGeometry object,
which only has orientation and tilt, and really acts
as a way to use functions that only need these two things, 
without needing coordinates. 

Rather than have a base class
for geometry that inherits, this, its just easier to have this trivial
class created, and pass it in the place of geometry objects when
coordinates are not needed
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sys import exit
try:
  from shapely.geometry import Polygon
  from shapely.geos import TopologicalError
  SHAPELY = True
except ImportError:
  SHAPELY = False

class SimpleGeometry(object):
  """ 
  SimpleGeometry Class

  __init__:   Object Constructor

  input(s):   orientation (radians)
              tilt (radians)
  output(s):  None
  """
  def __init__(self,orient,tilt):
    self.orient = orient
    self.tilt = tilt
    self.Rx = np.matrix([[1,0,0],[0,np.cos(self.tilt), -np.sin(self.tilt)],
      [0,np.sin(self.tilt), np.cos(self.tilt)]],float)
    self.Rz = np.matrix([[np.cos(self.orient), -np.sin(self.orient), 0], 
      [np.sin(self.orient), np.cos(self.orient), 0], [0,0,1]],float)
    self.R = np.transpose(self.Rx*self.Rz)
    self.Rsolar = np.transpose(self.Rx)*self.Rz

    if self.tilt < -1e10:
      print "warning, tilt is less than 0 for simple geometry object\n"
  def __repr__(self):
    return 'SimpleGeometry Object with orient '+str(self.orient) \
      +' and tilt '+str(self.tilt)

class Geometry(object):
  """ 
  Geometry Class

  __init__:   Object Constructor

  input(s):   coordinates
  output(s):  None
  """

  def __init__(self,coords,wall_coords = None):
    if(len(coords) != 4):
      exit("surface with more than 4 coordinates.\n \
        only quadrilateral facades supported for now\n")
    if wall_coords:
    if(len(wall_coords) != 4):
      exit("wall with more than 4 coordinates.\n \
        only quadrilateral facades supported for now\n")      
    self.coords = coords
    self._setGeometricData()

    # default information
    self._index = 0
    self.data = {}
    self.nX = 1
    self.nY = 1
  """ 
  _setGeometricData:  Internal function used if coordinates have changed,
                      or when initializing

  input(s):   None
  output(s):  None
  """   
  def _setGeometricData(self):
    # normal vector
    self.n = np.cross(self.coords[2]-self.coords[0],
      self.coords[1]-self.coords[0])
    self.n = self.n/np.linalg.norm(self.n)
 
    self.orient = np.pi/2.+np.arctan2(self.n[1],self.n[0])
    self.tilt = np.pi/2.- np.arctan2(np.sqrt(self.n[0]*self.n[0]
      +self.n[1]*self.n[1]),self.n[2])
    if self.orient > 0.999999999*np.pi:
      self.orient -= 2.*np.pi
    """
    This is part where directions are defined

    correct directions are critical for shading calculations,
    they are used as an easy way to check if facades will interact
    """
    # its a roof if the tilt is greater than 51 degrees
    if(abs(self.orient) <= np.pi/4. and self.tilt < 2.*np.pi/7.):
      self.dir = 'south'
    elif(abs(self.orient-np.pi/2.) <= np.pi/4. and self.tilt < 2.*np.pi/7.):
      self.dir = 'east'
    elif(abs(self.orient+np.pi) <= np.pi/4. and self.tilt < 2.*np.pi/7.):
      self.dir = 'north'
    elif(abs(self.orient+np.pi/2.) <= np.pi/4. and self.tilt < 2.*np.pi/7.):
      self.dir = 'west'
    elif(self.tilt > 2.*np.pi/7.):
      self.dir = 'roof'
    else:
      self.dir = 'unknown'
      print 'warning: facade with tilt of '+str(self.tilt)+' and orient of '+str(self.orient) \
        + ' has unknown direction'

    self.Rx = np.matrix([[1,0,0],[0,np.cos(self.tilt), -np.sin(self.tilt)],
      [0,np.sin(self.tilt), np.cos(self.tilt)]],float)
    self.Rz = np.matrix([[np.cos(self.orient), -np.sin(self.orient), 0], 
      [np.sin(self.orient), np.cos(self.orient), 0], [0,0,1]],float)
    self.R = np.transpose(self.Rx*self.Rz)
    self.Rsolar = np.transpose(self.Rx)*self.Rz

    # precompute the planar polygot used for intersection in the x-z plane
    if SHAPELY:
      self.poly = Polygon([(c[:,0],c[:,2]) for c in 
        [np.asmatrix(self[i])*self.R
        for i in range(len(self))]])
    else:
      print "\n!!! Warning: cannot compute shading \
      without shapely package!!! \n \
      https://pypi.python.org/pypi/Shapely\n"

    if self.tilt < -1e10:
      print "warning, tilt is less than 0 for geometry with "
      print self
      print "\n"
  """
              overloading of container functions

  input(s):   None
  output(s):  
  """

  def __len__(self):
    return 4
  def __iter__(self):
    return iter(self.coords)
  def __getitem__(self,key):
    return self.coords[key]
  def __setitem__(self,key,val):
    if(len(val) != 3):
      exit("coordinates must arrays of length 3")
    self.coords[key] = val
  def __repr__(self):
    return "Facade with height "+str(self.height())+" width " \
      + str(self.width()) + " in direction "+self.dir

  """
  height:     computes height and width of block
  width:
  
  input(s):   None
  output(s):  height, width
  """
  def height(self):
    return np.linalg.norm(self.coords[1]-self.coords[0])

  def width(self):
    return np.linalg.norm(self.coords[2]-self.coords[1])

  """
  rotate:     rotates all coordinates about an axis
  
  input(s):   axis (np array of size 3)
              point (np array of size 3)
              angle (radians)
  output(s):  height, width
  """
  def rotate(self, axis, point, angle):
    point = np.copy(point)
    R = rotationMatrix(axis,angle)

    for i in range(4):
      self.coords[i] -= point

    rotatedCoords = [np.squeeze(np.asarray(R*np.transpose(np.asmatrix(c)))) 
      for c in self.coords]
    
    for i in range(4):
      self.coords[i] = rotatedCoords[i] + point

    self._setGeometricData()

class GeometrySet(object):
  """ 
  GeometrySet Class

  __init__:   Object Constructor

  input(s):   (name) string corresponding to name
  output(s):  None
  """
  def __init__(self,name = ''):
    self.g = []
    self.s = None
    self.m = None
    self.name = name
    self.dataNames = None
    if not SHAPELY:
      print "\n!!! Warning: cannot compute shading \
      without shapely package!!! \n \
      https://pypi.python.org/pypi/Shapely\n"

  """
  len,iter,       these allow the GeometrySet to
  get/set item    be treated as a container, and
  append          iterated over

  input(s):       key, integer
                  val, Geometry Object
  output(s):      as shown
  """

  def __len__(self):
    return len(self.g)
  def __iter__(self):
    return iter(self.g)
  def __getitem__(self,key):
    return self.g[key]
  def __setitem__(self,key,val):
    self.g[key] = val

  """
  index:          This is a useful function, 
                  returning the index of a Geometry Object
                  in the container

  input(s):       Geometry Object
  output(s):      integer indicating location
  """  
  def index(self,geom):
    return geom._index

  def append(self,geom):
    geom._index = len(self)
    self.g.append(geom)

  """
  __repr__    overloading for print command

  input(s):   None
  output(s):  name
  """
  def __repr__(self):
    return "GeometrySet " + self.name +" with "+str(len(self))+" facades"

  """
  getMatch:   returns the match or subset list for particular
  getSubset   surface

  input(s):   Geometry Object
  output(s):  List of Geometry Objects
  """
  def getMatches(self,A):
    return self.m[self.index(A)]

  def getSubset(self,A):
    return self.s[self.index(A)]

  """
  computeBlockCounts:   sets nX, nY for all geometries

  input(s):   height and width
  output(s):  none
  """
  def computeBlockCounts(self,blockHeight,blockWidth):
    for g in self.g:
      g.nX = int(g.width()/blockWidth)
      g.nY = int(g.height()/blockHeight)

  """
  getAreas:   gets Areas, computes them if they haven't already been

  input(s):   none
  output(s):  dictionary of areas
  """
  def getAreas(self):
    bins = list(set([g.dir for g in self.g]));
    area = {b:0. for b in bins}
    for g in self.g:
      area[g.dir] += area[g.dir] + g.height()*g.width()

    area['total'] = np.sum([area[b] for b in bins])
    return area

  """
  collectDataNames:   sets the list of strings corresponding
                      to all the data stored in the geometry
                      really should only be used for output

  input(s):           none
  output(s):          none
  """
  def collectDataNames(self):
    if self.dataNames is None:
        # collect DataNames
      dataNames = set()
      for g in self.g:
        for k in g.data.keys():
          dataNames.add(k)
      self.dataNames = list(dataNames)

  """
  initializeData:   initialize data arrays as arrays of length nY

  input(s):         (names) list of data to initialize
  output(s):        none
  """      
  def initializeData(self,names):
    for g in self.g:
      for name in names:
        g.data[name] = np.zeros(g.nY)

  """
  computeSubsets:   computes the subset of surfaces that
                    could shade another surface, based
                    on solar rays. See doc/shaded.tex

                    self.s[index] is a list of indices of facades
                    in the geometry set that could shade it

  input(s):   useSunlit, boolean determining whether this even matters
  output(s):  none
  """
  def computeSubsets(self,useSunlit = False):
    self.s = []
    for i in range(len(self.g)):
      self.s.append([])  
    if not useSunlit:
      return
    if not SHAPELY:
      return  
    for A in self.g:
      # minimum z coordinate
      min_z = min([c[2] for c in A])

      # two passes are done to build a set before actual shading is checked
      firstpass = [B for B in self.g 
      if self.g.index(A) != self.g.index(B) 
      and A.dir == B.dir 
      and np.dot(A.n,B.n) > 1e-5
      and max([c[2] for c in B]) >= min_z]

      secondpass = []
      for B in firstpass:
        for i in range(4):
          p = (B[i]-A[i])/np.linalg.norm(B[i]-A[i])
          if np.dot(p,A.n) > 1e-5:
            secondpass.append(B)
            break
      # check if C blocks B from shading A
      for B in secondpass:
        polygonB = B.poly
        for C in secondpass:
          if(self.g.index(B) == self.g.index(C)): continue

          pC = [] # projected coordinates
          for i in range(4):
            n = B[i]-A[i]
            # project point from A through C onto B in direction of B-A
            s = np.dot(B.n,(B[0]-C[i]))/np.dot(B.n,n)
            if s < 1e-10:
              continue
            pC.append(C[i]+n*s) # projected points

          polygonB = self.getIntersection(polygonB,pC,B.R)

        # if they intersect a bit, add this to the set
        if(polygonB.area/B.poly.area > 1e-5):
          self.s[self.index(A)].append(B)

  """
  similar:    determines whether two objects are geometrically similar
              same normal, height, width, and shading
  input(s):   two geometry objects
  output(s):  True/False
  """
  def similar(self,A,B):
    tol = 1e-4
    if(abs(A.height() - B.height()) > tol):
      return False
    if(abs(A.width() - B.width()) > tol):
      return False
    if(np.dot(A.n,B.n) < 1.-tol):
      return False
    return True

  """
  computeMatches:   computes the subset of surfaces that are
                    geometrically similar, see above

                    self.m[index] contains a list of matches, these
                    are duplicated such that if self.m[1] = 2,
                    self.m[2] = 1, both contain the matching information
                    and only the lower indexed facade should be solved

  input(s):         useSunlit, boolean determining whether 
                    geometric shading occurs
  output(s):        none
  """
  
  def computeMatches(self,useSunlit = False):
    self.m = []
    for i in range(len(self.g)):
      self.m.append([])
    if useSunlit is False:
      for A in self.g:
        for B in self.g:
          if (B not in self.m[self.index(A)] and 
            self.similar(A,B) and self.index(A) != self.index(B)):
              self.m[self.index(A)].append(B)
              self.m[self.index(B)].append(A)
    else:
      if not SHAPELY:
        return      
      for A in self.g:
        # B is in the match of A if:
        # it hasnt been added, its similar, its not exactly the same one
        # and for now, its subset is either empty or size 1
        similar = [B for B in self.g 
          if (B not in self.m[self.index(A)] and 
            self.similar(A,B) and self.index(A) != self.index(B)
            and len(self.s[self.index(A)]) == len(self.s[self.index(B)])
            and len(self.s[self.index(A)]) <= 1)]

        # now if they both have sizes of 1 or less
        for B in similar:
          if len(self.s[self.index(A)]) == 1:
            AA = self.s[self.index(A)][0]
            BB = self.s[self.index(B)][0]
            if self.similar(AA,BB) and \
              sum([np.dot(AA[i]-A[i],BB[i]-B[i]) for i in range(4)]):
              self.m[self.index(A)].append(B)
              self.m[self.index(B)].append(A)             
          else:
            self.m[self.index(A)].append(B)
            self.m[self.index(B)].append(A)

  """
  getIntersection:  returns the remaining part of polygonA
                    not intersected by the coordinates in B
                    A - intersection(A,B)

  input(s):         polygonA, Polygon from geometry object A
                    B, list of coordinates
                    R, rotation matrix
                    plot (optional argument to plot these)
  output(s):        Polygon object

  This input format is used since we know what A is
  on the x-z plane, and then make a bunch of comparisons with other 
  geometry objects, B, and this reduces the number of times we create
  polygons, while simplifying the user API.
  """
  def getIntersection(self,polygonA,B,R,plot=False):
    if not SHAPELY:
      exit("unable to compute intersections without shapely package")
    if(len(B) < 3):
      return polygonA
    RB = (B*R)
    polygonB = Polygon([(RB[i,0],RB[i,2]) for i in range(len(B))])
    try:
      intersection = polygonA-polygonB.intersection(polygonA)
    except TopologicalError:
      intersection = polygonA

    if(plot):
      plotIntersections(polygonA,polygonB,intersection)

    return intersection

  """
  getUnshadedFraction:  gets the sunlit fraction of facade A

  inputs(s):            A, Geometry Object
                        n, vector
  outputs(s):           float, between 0 and 1
  """
  def getUnshadedFraction(self,A,n):
    if not SHAPELY:
      return 1.0
    # is the vector pointing into the surface
    dotAn = np.dot(A.n,n);
    if(dotAn) > 1e-10:
      return 0.

    polygonA = A.poly
    for B in self.getSubset(A):
      pPB = []
      for c in B:
        s = np.dot(A.n,(A[0]-c))/dotAn
        if s < 1e-10:
          continue          
        pPB.append(c+n*s)
      polygonA = self.getIntersection(polygonA,pPB,A.R)
    return polygonA.area/A.poly.area

  def printMatches(self):
    for A in self.g:
      print self.index(A),[self.index(B) for B in self.m[self.index(A)]]

  def printSubsets(self):
    for A in self.g:
      print self.index(A),[self.index(B) for B in self.s[self.index(A)]]    

"""
plotIntersections:    support for getIntersections
                      plots the intersection in the x-z plane
                      purely a debugging function, should never be
                      directly called, but using plot argument of 
                      getIntersections

input(s):       pA,pB,pI (Polygon Objects)
output(s):      none
"""      
def plotIntersections(pA,pB,pI,iA = 0,iB = 0):
  plt.plot([x for (x,y) in pI.exterior.coords],
    [y for (x,y) in pI.exterior.coords],
    label='intersection',linewidth=5.0)
  plt.plot([x for (x,y) in pA.exterior.coords],
    [y for (x,y) in pA.exterior.coords],
    label='A',linewidth=2.0,color='red')
  plt.plot([x for (x,y) in pB.exterior.coords],
    [y for (x,y) in pB.exterior.coords],
    label='B',linewidth=2.0,color='black')
  plt.legend(loc=0)
  plt.axis('equal')
  plt.show()

"""
plot:           support for GeometrySet, does some plotting 

input(s):       GeometrySet Object
                options, dictionary with True/False
                  modules, plot individual modules (very slow)
                  numbers, plot index numbers
                  normals, plot facade normals

output(s):      none
""" 
def plot(GeometrySet,**options):
  ax = plt.gcf().gca(projection='3d')
  ax.view_init(60,210)
  ax.set_xlabel('South',fontsize=18)
  ax.set_ylabel('West',fontsize=18)
  # ax.set_zlabel('Height',fontsize=18)

  for g in GeometrySet:
    xyz = []
    xyz_ = []
    for i in range(3):
      xyz.append([c[i] for c in g])
      xyz_.append(np.mean(xyz[i]))
      xyz[i].append(xyz[i][0])

    ax.plot(xyz[0],xyz[1],zs=xyz[2],color='black')
    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_zticks([])
    # Get rid of the panes
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # Get rid of the spines
    ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    if(options.get('normals',False)):
      ax.quiver(xyz_[0],xyz_[1],xyz_[2],g.n[0],g.n[1],g.n[2],
        length=5.0,linewidth=2.0,color='red')
    if(options.get('numbers',False)):
      ax.text(xyz_[0],xyz_[1],xyz_[2],str(GeometrySet.index(g)))
    if(options.get('modules',True) and g.nX*g.nY > 1):

      dX = (g[3]-g[0])/g.nX
      dY = (g[1]-g[0])/g.nY
      for x in range(g.nX):
        for y in range(g.nY):
          c = [[],[],[],[]]
          c[0] = dX*x    +dY*y
          c[1] = dX*(x+1)+dY*y
          c[2] = dX*(x+1)+dY*(y+1)
          c[3] = dX*x    +dY*(y+1)
          c.append(c[0])
          xx = np.zeros(5)
          yy = np.zeros(5)
          zz = np.zeros(5)
          for i in range(5):
            xx[i] = c[i][0] + g[0][0]
            yy[i] = c[i][1] + g[0][1]
            zz[i] = c[i][2] + g[0][2]
          ax.plot(xx,yy,zs=zz,color='black')
  Xmax = max([x[0] for x in g for g in GeometrySet])
  Ymax = max([x[1] for x in g for g in GeometrySet])
  Zmax = max([x[2] for x in g for g in GeometrySet])
  Xmin = min([x[0] for x in g for g in GeometrySet])
  Ymin = min([x[1] for x in g for g in GeometrySet])
  Zmin = min([x[2] for x in g for g in GeometrySet])

  max_range = np.array([Xmax-Xmin, Ymax-Ymin, Zmax-Zmin]).max()
  Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(Xmax+Xmin)
  Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Ymax+Ymin)
  Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Zmax+Zmin)

  for xb, yb, zb in zip(Xb, Yb, Zb):
     ax.plot([xb], [yb], [zb], 'w')
  plt.show()

"""
rotationMatrix: returns a rotation matrix

input(s):       axis (arbitrary), angle (radians)
output(s):      numpy matrix
""" 
def rotationMatrix(axis, angle):
  a = axis/np.linalg.norm(axis)
  return np.transpose(np.cos(angle)*np.identity(3)+\
    np.sin(angle)*np.matrix([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])+\
    (1.-np.cos(angle))*np.outer(a,a))
