""" 
blocks.py contains the Block Class

Each Block has:
  (.name) string identifying the block's name (or ID)
  (.P) parameters, a dictionary of parameters or functions
     which depend on the states inside the block
     These are common and to be used in fluxes, sources, etc
     (optional)
  (.state) physical state
  (.F) list of fluxes, which return a dictionary
  (.S) list of sources, which return a dictionary
    Sources and Fluxes do not need to be ordered 
    since they are never explicitly globally unwrapped
  (.t) time
  (.T) Time function, returns dict of coefficient on time terms

The equation for the block is
R(state) = Sum(Fluxes(state)) + Sum(Sources(state)) = 0
In the unsteady case
d/d(state) = R(state)

All blocks are connected through fluxes, defined in flux.py


"""
from collections import OrderedDict

class Block(object):
  """ 
  Block Class

  __init__:   Object Constructor

  input(s):   (name) string corresponding to block name
              (initial) dictionary of initial conditions
              (parameterFunctions) dictionary of parameter functions,
                functions that are shared among fluxes/sources,
                such as material properties
              (parameters) dictionary of parameters
              (t) time

  output(s):  None
  """

  def __init__(self,name,initial,parameterFunctions=None,\
    parameters=None,t=0):
    self.name = name
    self.state = OrderedDict(initial)
    self.P = parameterFunctions
    self.p = parameters
    self.F = []
    self.S = [] 
    self.t = t
    self.T = lambda B : dict([(s,1) for s in B.state])

  """
  These overload the [] operator, such that the 
  states can be set much more easily

  """
  def __getitem__(self,key):
    return self.state[key]

  def __setitem__(self,key,val):
    self.state[key] = val
  """
  addFlux:
  addSource: Block Setup Functions

  input(s):  (F,S) Flux objects, Source objects
  output(s): None

  these functions don't add anything new, but make
  building blocks easier using addFlux/addSource instead
  of directly modifying the lists
  """
  def addFlux(self,F):
    self.F.append(F)
    F.B = self

  def removeFlux(self,name):
    self.F[:] = [F for F in self.F if F.name != name]

  def addSource(self,S):
    self.S.append(S)
    S.B = self

  def removeSource(self,name):
    self.S[:] = [S for S in self.S if S.name != name]

  """
  R: Residual function

  input(s):  None
  output(s): dict of states and residuals corresponding to state

  sums over sources and fluxes to calculate the residual
  """
  def R(self):
    R = OrderedDict([(s,0) for s in self.state])
    for d in [F.flux() for F in self.F]+[S.source() for S in self.S]:
      for s in d:
        R[s] += d[s]
    return R

  """
  __repr__

  input(s):  None
  output(s): str representation of the object. 

  allows for using print direction on an object
  """ 
  def __repr__(self):
    return "block named " +self.name+" with "+str(len(self.state))+"\n"  \
      + ",".join([s + '=' + str(self[s]) for s in self.state])
