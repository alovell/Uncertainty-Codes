import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker

inputfile = sys.argv[1]

#load the 2d data, usecols tells you which columns in your file to import
f = open(inputfile, 'r')
x, y, z = np.loadtxt(f, usecols=(0,1,2), unpack = True)

#define x and y grid from unique points in your x,y sets
xset = set(x)
yset = set(y)
xval = list(xset)
yval = list(yset)
xval.sort()
yval.sort()
z2d = np.transpose(np.reshape(z, (len(xval),len(yval))))


#finding the absolute minimum
zmin = z2d.argmin()
idx = np.unravel_index(zmin, z2d.shape)
xmax = xval[idx[1]]
ymax = yval[idx[0]]
#print z2d.min()
#print xmax, ymax

#normalise the countour so it doesn't look weird, pick cmap
mnl = matplotlib.ticker.MaxNLocator(nbins=30)
levels = mnl.bin_boundaries(z2d.min(),z2d.max())
cmap = plt.get_cmap('Spectral')
#cmap = plt.get_cmap('gray')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

#set the axis labels and title with LaTeX
plt.rc('text',usetex=True)
plt.rc('font',family='serif')
strx = r'r [fm]' #set xlabel
#stry = r'$\delta_v$ [fm]' #set ylabel
stry = r'as [fm]'
strt = r'V-Ws contour' #set title

#plot the graphs and min
plt.contour(xval,yval,z2d, 20, cmap = cmap)
plt.plot(xmax,ymax,'k*', markersize=10)
cbar = plt.colorbar()
cbar.set_label('$\sqrt{\chi ^2}$',fontsize=20)
plt.xlim(min(xval),max(xval))
plt.ylim(min(yval),max(yval))
plt.xlabel(strx,fontsize=20)
plt.ylabel(stry,fontsize=20)
#plt.title(strt)
plt.show()
