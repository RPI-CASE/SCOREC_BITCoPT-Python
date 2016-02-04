import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
indices = [11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,\
           44,45,51,52,53,54,55]
indices = [53,33,13]
titles = ['Bottom','Middle','Top']
for i in range(3):

  x = np.arange(-72,73, 3)
  y = np.arange(-72,73, 3)

  Y, P = np.meshgrid(x, y)
  z = np.loadtxt(str(indices[i])+'.txt')
  plt.subplot(1,3,i+1)
  pc = plt.pcolor(Y,P,z,cmap = cm.Greys_r)

  plt.xlabel('Yaw',fontsize=14)
  if i == 0:
    plt.ylabel('Pitch',fontsize=14)
  else:
    plt.gca().get_yaxis().set_ticklabels([])

  plt.xlim([-72,72])
  plt.ylim([-72,72])
  plt.title(titles[i],fontsize=14)
  plt.gca().tick_params(axis='both', which='major', labelsize=12,direction='out')
plt.gcf().subplots_adjust(right=0.88,bottom=0.15)
#plt.tight_layout()
cbar_ax = plt.gcf().add_axes([0.90, 0.15, 0.02, 0.75])
plt.gcf().colorbar(pc, cax=cbar_ax)
plt.gcf().set_size_inches(10,3.5)

plt.savefig('composite.png')
plt.close()
