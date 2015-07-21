

import numpy as np
import matplotlib.pyplot as pl
import sys, glob
#from matplotlib import pyplot as plt
#import matplotlib.mlab as mlab
#import pylab as pl
	
if len(sys.argv) <2:
      print "No data file name given. Please enter"
      datafile = raw_input("-> ")
      if len(glob.glob(datafile))==0:
            print "Data file %s not found. Exiting" % datafile
            sys.exit() 

for x in range(1,len(sys.argv)):
    datafile=sys.argv[x]
    if len(glob.glob(datafile))==0:
        print "Data file %s not found. Exiting" % datafile
        sys.exit() 
    data=np.loadtxt(datafile)
#    pl.plot(data[:,0], data[:,1],marker='o',markersize=3,label= datafile)



plotname = 'plot.png'
      
v_hist = np.ravel(data[:,1])   # 'flatten' van klein naar groot zetten
fig = pl.figure()
ax1 = fig.add_subplot(111)

n, bins, patches = ax1.hist(v_hist, bins = 5, normed=1, facecolor='green')
#bincenters = 0.5*(bins[1:]+bins[:-1])
#y = mlab.normpdf( bincenters)
#l = ax.plot(bincenters, y, 'r--', linewidth=1)
ax1.set_ylabel(r'$P \ (\Delta x) $', fontsize=15)
ax1.set_xlabel(r'$\frac{\Delta x}{\sigma} $', fontsize=15)
ax1.set_xlim(0.95, 1.05)
ax1.set_ylim(0, 20)
#pl.savefig(plotname, dpi=100) 
pl.show()
