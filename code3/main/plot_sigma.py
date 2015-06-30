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
# if len(sys.argv) >1:
#       datafile = sys.argv[1]
#       if len(glob.glob(datafile))==0:
#             print "Data file %s not found. Exiting" % datafile
#             sys.exit() 
#       data = np.loadtxt(datafile)
#       line1 = pl.plot(data[:,1], data[:,2], 'r',marker='o', markersize=3, label=r'$\bar{\nu}=0.5$')
# if len(sys.argv) >2:
#       datafile2 = sys.argv[2]
#       data2 = np.loadtxt(datafile2)
#       line2 = pl.plot(data2[:,1], data2[:,2], 'b',marker='o',  markersize=3, label=r'$\bar{\nu}=0.05$'  )
# if len(sys.argv) >3:
#       datafile3 = sys.argv[3]
#       data3 = np.loadtxt(datafile3)
#       line3 = pl.plot(data3[:,1], data3[:,2], 'g',marker='o',markersize=3,  label= r'$\bar{\nu}=0.5$'  )



parabola_file = open("data/endpoints.txt", "w")

 

for x in range(1,len(sys.argv)):
    datafile=sys.argv[x]
    if len(glob.glob(datafile))==0:
        print "Data file %s not found. Exiting" % datafile
        sys.exit() 
    for endline in open(datafile):pass                 #to save point (mean,var) at last timestep 
    parabola_file.write(endline)        
    data=np.loadtxt(datafile)        
    pl.plot(data[:,1], data[:,2], marker='o', color='black', markersize=3,label= datafile)

parabola_file.close();

endpoints= np.loadtxt("data/endpoints.txt");
pl.plot(endpoints[:,1], endpoints[:,2]); 

# parabola = np.polyfit(endpoints[:,1], endpoints[:,2],2)  #degree 2 = parabola       #this method for fitting parabola did not work for some reason
# plotpar = np.poly1d(parabola)
# pl.plot(plotpar,'--')


#q = plt.quiver(X,Y,u,v,angles='xy',scale=1000,color='r')  #check this for plotting arrows


xmin =0
xmax = 1
ymin = 0
ymax = 0.5#dependent on number of realizations 
pl.axis([xmin, xmax, ymin, ymax]);
#pl.yscale('log')
pl.xlabel(r'$\mathbb{E}(\rho_N(t))$', fontsize = 16)
#pl.ylabel(r'$T^*$', fontsize = 20)
pl.ylabel(r'$\sigma(\rho_N(t))$', fontsize = 20)
#pl.legend([plot, plot2], ('add', 'mul'), 'best', numpoints=1)
#pl.legend([line1, line2, line3],  'best')
#pl.legend( loc='best', numpoints = 1 )
#first_legend = pl.legend([line1], loc=1)
#pl.legend(bbox_to_anchor=(1, 0.96), numpoints = 1 )
#pl.legend([line3], bbox_to_anchor=(0.3, 0.6))
#ax = pl.gca().add_artist(first_legend)




pl.show()