

import numpy as np
import matplotlib.pyplot as pl
import sys, glob





	
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
    pl.plot(data[:,0], data[:,1],marker='o',markersize=3,label= datafile)

xmin =0
xmax = 100
ymin = 0
ymax = 1
pl.axis([xmin, xmax, ymin, ymax])
pl.xlabel('t', fontsize = 16)

pl.ylabel(r'$\mathbb{E}(\rho_N(t))$', fontsize = 20)
#pl.legend(bbox_to_anchor=best, numpoints = 1 )
pl.legend( loc='best', numpoints = 1 )
pl.show()

