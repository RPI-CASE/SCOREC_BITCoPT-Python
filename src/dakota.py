"""
dakota.py contains support to interact with dakota,
specifically the file readers and writers. This is
copied largely off of their example file

This really isnt needed now, but deleting it seems foolish, 
so lets let it exist

"""
import os
import re
import numpy as np

"""
writes a dakota input file where the standard deviation of DNI is the variable
of interest

https://dakota.sandia.gov/

"""
def writeICSDakotaInputFile(filename,DNISTD,day,directory,quadOrder):
  f = open(filename,'w')
  f.write("# Dakota Input File: "+filename+"\n")
  f.write("environment\n")
  f.write("\n")
  f.write("method\n")
  f.write("  polynomial_chaos\n")
  f.write("    quadrature_order   = "+str(quadOrder)+"\n")
  f.write("    dimension_preference = 1\n")
  f.write("    samples = 1000\n")
  f.write("    seed = 12347 rng rnum2\n")
  f.write("    response_levels =\n")
  f.write("    .1 1. 50.\n")
  f.write("    variance_based_decomp\n")
  f.write("\n")
  f.write("variables\n")
  f.write("  normal_uncertain = 1\n")
  f.write("    means             = 0\n")
  f.write("    std_deviations    = "+str(0.01*DNISTD)+"\n")
  f.write("    descriptors       = 'DNISTD'\n")
  f.write("  discrete_state_set\n")
  f.write("  string 1\n")
  f.write("    descriptors 'day'\n")
  f.write("    elements_per_variable 1\n")
  f.write("    elements '"+str(day)+"'\n")
  f.write("\n")
  f.write("interface\n")
  f.write("  system\n")
  f.write("  analysis_drivers = 'python ICSolar_Dakota.py'\n")
  f.write("  parameters_file = 'Results/UQ/"+directory+'/'+"params"+str(day)+".in'\n")
  f.write("  results_file = 'Results/UQ/"+directory+'/'+"results"+str(day)+".out'\n")
  f.write("  file_tag\n  file_save\n")
  f.write("\n")
  f.write("responses\n")
  f.write("  response_functions = 1\n")
  f.write("  no_gradients\n")
  f.write("  no_hessian\n")
  f.close()

"""
based this function off an example in dakota
"""
def readInputFile(inputFile):
  # ----------------------------
  # Parse DAKOTA parameters file
  # ----------------------------

  # setup regular expressions for parameter/label matching
  e = '-?(?:\\d+\\.?\\d*|\\.\\d+)[eEdD](?:\\+|-)?\\d+' # exponential notation
  f = '-?\\d+\\.\\d*|-?\\.\\d+'                        # floating point
  i = '-?\\d+'                                         # integer
  value = e+'|'+f+'|'+i                                # numeric field
  tag = '\\w+(?::\\w+)*'                               # text tag field

  # regular expression for aprepro parameters format
  aprepro_regex = re.compile('^\s*\{\s*(' + tag + ')\s*=\s*(' + value +')\s*\}$')
  # regular expression for standard parameters format
  standard_regex = re.compile('^\s*(' + value +')\s+(' + tag + ')$')

  # open DAKOTA parameters file for reading
  paramsfile = open(inputFile, 'r')

  # extract the parameters from the file and store in a dictionary
  paramsdict = {}
  for line in paramsfile:
    m = aprepro_regex.match(line)
    if m:
      paramsdict[m.group(1)] = m.group(2)
    else:
      m = standard_regex.match(line)
    if m:
      paramsdict[m.group(2)] = m.group(1)

  paramsfile.close()

  # crude error checking; handle both standard and aprepro cases
  num_vars = 0
  if ('variables' in paramsdict):
    num_vars = int(paramsdict['variables'])
  elif ('DAKOTA_VARS' in paramsdict):
    num_vars = int(paramsdict['DAKOTA_VARS'])

  num_fns = 0
  if ('functions' in paramsdict):
    num_fns = int(paramsdict['functions'])
  elif ('DAKOTA_FNS' in paramsdict):
    num_fns = int(paramsdict['DAKOTA_FNS'])

  # -------------------------------
  # Convert and send to application
  # -------------------------------
  return paramsdict

"""
based this function off an example in dakota
"""
def readOutputFile(outputFile,quadOrder):
  f = open(outputFile,'r')
  coeff = np.zeros(quadOrder)
  for line in f:
    index = line.find('He')
    if (index >= 0):
      offset = 3 - (line[2] == '-')
      value = float(line.split(' ')[offset])
      coeff[int(re.search(r'\d+', line[index:-1]).group())] = value
    if line.startswith('  expansion:'):
      mean = float(line.split()[1])
      stddev = float(line.split()[2])
      return (mean,stddev,coeff)
