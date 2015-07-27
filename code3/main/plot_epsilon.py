

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
    pl.plot(data[:,0], data[:,1],marker='o',markersize=3,label= datafile)



plotname = 'plot.png'
      
xmin =1e-18
xmax = 1
ymin = 1e-18
ymax = 1
pl.axis([xmin, xmax, ymin, ymax])

pl.yscale('log')
pl.xscale('log')

pl.xlabel(r'$\epsilon $', fontsize=15)
pl.ylabel(r'$\phi_T(U +\epsilon) -\phi_T(U) $', fontsize=15)
#ax1.set_xlim(-4, 4)
#ax1.set_ylim(0, 1)
#plt.savefig(plotname, dpi=100) 

#pl.legend([plot, plot2], ('add', 'mul'), 'best', numpoints=1)
#pl.legend([line1, line2, line3],  'best')
#pl.legend( loc='best', numpoints = 1 )
#first_legend = pl.legend([line1], loc=1)
pl.legend(bbox_to_anchor=(1, 0.96), numpoints = 1 )
pl.show()
