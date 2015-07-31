

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

# for x in range(1,len(sys.argv)):
#     datafile=sys.argv[x]
#     if len(glob.glob(datafile))==0:
#         print "Data file %s not found. Exiting" % datafile
#         sys.exit() 
#     data=np.loadtxt(datafile)
#     pl.plot(data[:,0], data[:,1],marker='o',markersize=3,label= datafile)

if len(sys.argv) >1:
      datafile = sys.argv[1]
      if len(glob.glob(datafile))==0:
            print "Data file %s not found. Exiting" % datafile
            sys.exit() 
      data = np.loadtxt(datafile)
      line1 = pl.plot(data[:,0], data[:,1] , 'r',marker='o', markersize=3, label=r'$\mathcal{L}$')
if len(sys.argv) >2:
      datafile2 = sys.argv[2]
      data2 = np.loadtxt(datafile2)
      line2 = pl.plot(data2[:,0], data2[:,1], 'b',marker='o',  markersize=3, label=r'$\mathcal{L_W}$' )
if len(sys.argv) >3:
      datafile3 = sys.argv[3]
      data3 = np.loadtxt(datafile3)
      line3 = pl.plot(data3[:,0], data3[:,1], 'g',marker='o',markersize=3,  label= r'$\bar{\nu}=0.5$'  )

plotname = 'plot.png'
      
xmin =1e-18
xmax = 1
ymin = 1e-18
ymax = 1
pl.axis([xmin, xmax, ymin, ymax])

pl.yscale('log')
pl.xscale('log')

pl.xlabel(r'$\epsilon $', fontsize=15)
pl.ylabel(r'$\epsilon -\phi_T(U +\epsilon) + \phi_T(U) $', fontsize=15)
#ax1.set_xlim(-4, 4)
#ax1.set_ylim(0, 1)
#plt.savefig(plotname, dpi=100) 

#pl.legend([plot, plot2], ('add', 'mul'), 'best', numpoints=1)
#pl.legend([line1, line2, line3],  'best')
#pl.legend( loc='best', numpoints = 1 )
#first_legend = pl.legend([line1], loc=1)
pl.legend(bbox_to_anchor=(1, 0.2), numpoints = 1 )
pl.show()
