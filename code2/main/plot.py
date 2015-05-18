#readFileAndPlot.py
import numpy as np
import pylab as pl
import matplotlib as mpl
#import matplotlib.pyplot as plt
# Use numpy to load the data contained in the file
#from matplotlib.ticker import formatter
import sys, glob
# Use numpy to load the data contained in the file
#from matplotlib.ticker import formatter

	
if len(sys.argv) > 1:
      datafile = sys.argv[1]
      datafile2 = sys.argv[2]
      datafile3 = sys.argv[3]
     # datafile4 = sys.argv[4]
else:
      print "No data file name given. Please enter"
      datafile = raw_input("-> ")
if len(glob.glob(datafile))==0:
      print "Data file %s not found. Exiting" % datafile
      sys.exit() 
data = np.loadtxt(datafile)
data2 = np.loadtxt(datafile2)
data3 = np.loadtxt(datafile3)

plot= pl.plot(data[:,0], data[:,1] , 'r', label=r'$\bar{\nu}=0.5$')
plot2= pl.plot(data2[:,0], data2[:,1], 'b', label=r'$\bar{\nu}=0.05$')
plot3= pl.plot(data2[:,0], data3[:,1], 'g', label= r'$\bar{\nu}=0.5$')

xmin =0
xmax = 30
ymin = 0
ymax = 1
pl.axis([xmin, xmax, ymin, ymax])
#pl.yscale('log')
pl.xlabel('t', fontsize = 16)
#pl.ylabel(r'$T^*$', fontsize = 20)
pl.ylabel(r'$\mathbb{E}(\rho_N(t))$', fontsize = 20)
#pl.legend([plot, plot2], ('add', 'mul'), 'best', numpoints=1)
# pl.legend([plot, plot2, plot3], (r'$\nu}=0.5$',r'$\nu=0.05$',r'$\nu=0.5$'), 'best')
pl.legend( loc='best', numpoints = 1 )
pl.show()
