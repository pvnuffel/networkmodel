

import numpy as np
import matplotlib.pyplot as pl
import sys, glob
from scipy.stats import norm
#from matplotlib import pyplot as plt
#import matplotlib.mlab as mlab
#import pylab as pl
	
if len(sys.argv) <2:
      print "No data file name given. Please enter"
      datafile = raw_input("-> ")
      if len(glob.glob(datafile))==0:
            print "Data file %s not found. Exiting" % datafile
            sys.exit() 


fig = pl.figure()
ax1 = fig.add_subplot(111)

for x in range(1,len(sys.argv)):
    datafile=sys.argv[x]
    if len(glob.glob(datafile))==0:
        print "Data file %s not found. Exiting" % datafile
        sys.exit() 
    data=np.loadtxt(datafile)
    mu, std = norm.fit(data[:,1])
    print mu
    print std
    x = np.linspace(0,2, 100)
    p = norm.pdf(x, mu, std)
    pl.plot(x, p, 'k', linewidth=1)
#    pl.plot(data[:,0], data[:,1],marker='o',markersize=3,label= datafile)

#testgaussdata = norm.rvs(1, 0.05, size=4000)

plotname = 'plot.png'
      
#v_hist = np.ravel(data[:,1])   # 'flatten' van klein naar groot zetten


#n, bins, patches = ax1.hist(v_hist, normed=True, bins = 60, facecolor='green')
#bincenters = 0.5*(bins[1:]+bins[:-1])
#y = mlab.normpdf( bincenters)
#l = ax.plot(bincenters, y, 'r--', linewidth=1)

#pl.hist(data1, bins=100, normed=1, alpha=0.6, color='g')


ax1.annotate('M', xy=(1.5, 0.5),  xycoords='data',
                xytext=(-40, -30), textcoords='offset points',
                arrowprops=dict(arrowstyle="<-",
                                connectionstyle="arc,angleA=0,armA=20,angleB=-90,armB=15,rad=7"),
                )

#pl.text(0.9, 2.7, r'$M=10^3$', fontsize=15)

ax1.set_ylabel(r'$P \ (w) $', fontsize=15)
ax1.set_xlabel(r'$w$', fontsize=15)
ax1.set_xlim(0, 2)
ax1.set_ylim(0, 8)
#pl.savefig(plotname, dpi=100) 
pl.show()


