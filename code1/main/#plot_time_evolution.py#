#readFileAndPlot.py
import numpy as np
#import pylab as pl
import matplotlib.pyplot as pl
# Use numpy to load the data contained in the file
#from matplotlib.ticker import formatter
import sys, glob
# Use numpy to load the data contained in the file
#from matplotlib.ticker import formatter




	
if len(sys.argv) <2:
      print "No data file name given. Please enter"
      datafile = raw_input("-> ")
      if len(glob.glob(datafile))==0:
            print "Data file %s not found. Exiting" % datafile
            sys.exit() 
if len(sys.argv) >1:
      datafile = sys.argv[1]
      if len(glob.glob(datafile))==0:
            print "Data file %s not found. Exiting" % datafile
            sys.exit() 
      data = np.loadtxt(datafile)
      line1 = pl.plot(data[:,0], data[:,1] , 'r',marker='o', markersize=3, label=r'$\bar{\nu}=0.5$')
if len(sys.argv) >2:
      datafile2 = sys.argv[2]
      data2 = np.loadtxt(datafile2)
      line2 = pl.plot(data2[:,0], data2[:,1], 'b',marker='o',  markersize=3, label=r'$\bar{\nu}=0.05$'  )
if len(sys.argv) >3:
      datafile3 = sys.argv[3]
      data3 = np.loadtxt(datafile3)
      line3 = pl.plot(data3[:,0], data3[:,1], 'g',marker='o',markersize=3,  label= r'$\bar{\nu}=0.5$'  )


# if len(sys.argv) >2:
#       datafile2 = sys.argv[2]
#       datafile3 = sys.argv[3]
#      # datafile4 = sys.argv[4]
# else:
#       print "No data file name given. Please enter"
#       datafile = raw_input("-> ")
# if len(glob.glob(datafile))==0:
#       print "Data file %s not found. Exiting" % datafile
#       sys.exit() 
# data = np.loadtxt(datafile)
# data2 = np.loadtxt(datafile2)
# data3 = np.loadtxt(datafile3)

# line1 = pl.plot(data[:,0], data[:,1] , 'r',marker='o', markersize=3, label=r'$\bar{\nu}=0.5$')
# line2 = pl.plot(data2[:,0], data2[:,1], 'b',marker='o',  markersize=3, label=r'$\bar{\nu}=0.05$'  )
# line3 = pl.plot(data3[:,0], data3[:,1], 'g',marker='o',markersize=3,  label= r'$\bar{\nu}=0.5$'  )


xmin =0
xmax = 100
ymin = 0
ymax = 1
pl.axis([xmin, xmax, ymin, ymax])
#pl.yscale('log')
pl.xlabel('t', fontsize = 16)
#pl.ylabel(r'$T^*$', fontsize = 20)
pl.ylabel(r'$\mathbb{E}(\rho_N(t))$', fontsize = 20)
#pl.legend([plot, plot2], ('add', 'mul'), 'best', numpoints=1)
#pl.legend([line1, line2, line3],  'best')
#pl.legend( loc='best', numpoints = 1 )
#first_legend = pl.legend([line1], loc=1)
pl.legend(bbox_to_anchor=(1, 0.96), numpoints = 1 )
#pl.legend([line3], bbox_to_anchor=(0.3, 0.6))
#ax = pl.gca().add_artist(first_legend)




pl.show()
