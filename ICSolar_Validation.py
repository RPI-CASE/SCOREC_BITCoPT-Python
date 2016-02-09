import src.icsolar as icsolar
import src.casegeom as casegeom
import src.shading as shading
import numpy as np
import time as cputime
import matplotlib.pyplot as plt
import os
import test

from matplotlib import gridspec
from joblib import Parallel, delayed

"""
This file will perform the validation runs for ICSolar using the data
in data/ICSolar/


"""
def run(name,nicename):
  ############################################### timer

  clockStart = cputime.time()

  ############################################### data loading 

  data = casegeom.readValidationFile('data/ICSolar/'+name+'.csv')
  n = len(data['Timestamp'])

  print 'loading data has taken','%.2f' % (cputime.time()-clockStart)

  # create a geometry for this case to get solar information off of
  # with orientation -40, and tilt of 0.
  geometry = casegeom.createSimpleGeometry(-40*np.pi/180.,0)
  ############################################### set up results 

  # initialize results, list of numpy arrays
  moduleTemp = [np.zeros(n),np.zeros(n),np.zeros(n),
                np.zeros(n),np.zeros(n),np.zeros(n)]
  tubeTemp = [np.zeros(n),np.zeros(n),np.zeros(n),
              np.zeros(n),np.zeros(n),np.zeros(n)]
  receiverTemp = [np.zeros(n),np.zeros(n),np.zeros(n),
                  np.zeros(n),np.zeros(n),np.zeros(n)]

  ############################################### pre-computed information
  # this is a lookup table to connect each module to a shading table
  shadingIndices = [0,1,2,3,4,5]

  ############################################### input parameters 
  interiorAirTemp = 20.0 # Celcius

  # 2.0 m/s * 0.16, cross sectional area, times density
  airFlowRate = 2.0*0.16*1.200

  # timestep is ten seconds
  init = {'dt':10.,'length':icsolar.moduleHeight,'interiorAirTemp':interiorAirTemp, 
    'inletAirTemp':interiorAirTemp,'airFlowRate':airFlowRate,
    'numModules':6}

  ############################################### set up results 
  if not os.path.exists('ValidationResults'):
    os.makedirs('ValidationResults')
  ############################################### solver

  for ts in range(n):
    # times the length of the step
    clockStepStart = cputime.time()

    init['exteriorAirTemp'] = data['Tamb'][ts] 
    init['inletWaterTemp'] = data['exp_inlet'][ts]
    init['waterFlowRate'] = data['exp_flowrate'][ts]*1000.*1e-6 # kg/s
    time = float(data['Timestamp'][ts])
    shade = np.ones(6)

    for i in range(6):
      shade[i] = shading.getStudioUnshadedFractionAtTime(time,geometry,shadingIndices[i])

    # Uncomment this to use exact energy in as the input, for debugging
    # init['Q_w'] = [data['heatgen_m'+str(7-i)][ts] for i in range(1,7)]
    # assume 30% electrical efficiency
    init['Q_d'] = 0.57*625.5*0.0001*data['DNI'][ts]*shade*(1.-0.30)
    # set up previous temperature, if we can
    if (ts == 0):
      init['previousWaterModuleT'] = init['inletWaterTemp']*np.ones(6)
      init['previousWaterTubeT'] = init['inletWaterTemp']*np.ones(6)
      init['previousWaterFlowRate'] = init['waterFlowRate']
    else:
      init['dt'] = time - data['Timestamp'][ts-1]
      init['previousWaterModuleT'] = results['waterModule']
      init['previousWaterTubeT'] = results['waterTube']
      init['previousWaterFlowRate'] = data['exp_flowrate'][ts-1]*1000.*1e-6

    # solve the problem
    results = icsolar.solve(init)

    for i in range(6):
      moduleTemp[i][ts] = results['waterModule'][i]
      tubeTemp[i][ts] = results['waterTube'][i]
      receiverTemp[i][ts] = results['receiver'][i]

  for i in range(6):
    np.savetxt('ValidationResults/'+name+'_'+str(i)+'.txt',moduleTemp[i],
      fmt='%2.2f',delimiter=',')

  ############################################### cleanup
  print nicename,'final runtime is','%.2f' % (cputime.time()-clockStart),'for',ts,'steps'
  print nicename,'solving ODES:','%.2f' % icsolar.getSolverTime()
  ############################################### lets put a measure of accuracy
  T_RMSE = np.zeros(6)
  for i in range(6):
    T_RMSE[i] = np.sqrt(np.mean((moduleTemp[i]-data['m'+str(7-(i+1))+'_out'])**2))
    print nicename,'RMS of Temperature for module ',i,'is','%.2f' % T_RMSE[i]

  ############################################### plot outputs
  t = np.arange(n)/6. # minutes

  for i in range(6):

    gs = gridspec.GridSpec(2,1, height_ratios=[3, 1]) 

    ax1 = plt.subplot(gs[0])
    plt.title('Module '+str(i+1),fontsize=22)

    ax1.plot(t,data['m'+str(7-(i+1))+'_out'],'--k',linewidth=2.)
    ax1.plot(t,moduleTemp[i],'-k',linewidth=2.)

    ax1.set_ylabel(r'T ($^\circ$C)',fontsize=20)
    ax1.set_ylim([min(data['m1_out'])*0.8,max(data['m1_out'])*1.2])
    ax1.get_xaxis().set_ticks(np.arange(0,int(max(t)/10)*10+10,int(max(t)/6/10)*10))

    ax1.get_xaxis().set_ticklabels([])
    ax1.tick_params(axis='both', which='major', labelsize=14)

    plt.gca().set_xlim([t[0],int(max(t)/10)*10])

    ax2 = plt.subplot(gs[1])
    ax2.plot(t,data['DNI'],'-k',linewidth=2.)
    ax2.set_ylabel('DNI',fontsize=20)
    ax2.set_xlabel(r'Time (minutes)',fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax2.set_ylim([0,max(data['DNI'])*1.2])
    ax2.get_yaxis().set_ticks([0,int(max(data['DNI']))])
    ax2.get_xaxis().set_ticks(np.arange(0,int(max(t)/10)*10+10,int(max(t)/6/10)*10))
    ax2.yaxis.grid(True,linewidth=2.)
    plt.gcf().set_size_inches(4,4)
    plt.gca().set_xlim([t[0],int(max(t)/10)*10])
    plt.tight_layout()

    plt.savefig('ValidationResults/' + name+'_module_'+str(i+1)+'.png')

    plt.close()

if __name__ == "__main__":
  test.runTests()
  # actual filename
  basenames = ['Jan28','Jan31','Feb6','Feb11','Feb27','Feb28','Mar06','Mar09','Mar19']
  # nice string for printing
  nicenames = ['Jan 28','Jan 31','Feb 6','Feb 11','Feb 27','Feb 28','Mar 6','Mar 9','Mar 19']
  n = len(basenames)
  Parallel(n_jobs=n)(delayed(run)(basenames[i],nicenames[i]) \
    for i in range(n))
