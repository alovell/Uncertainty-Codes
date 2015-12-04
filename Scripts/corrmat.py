# I should really make this more general

import matplotlib.pyplot as plt
import numpy as np
import sys

# input file is given in the command line
inputfile = sys.argv[1]

# load the data into the array corr
corr = np.loadtxt(inputfile)
correl = np.zeros((15,200))
corrinel = np.zeros((15,200))

# get rid of the first two columns (data values and errors)
for x in range(0,15):
	for y in range(0,200):
		#print x,y
		correl[x][y] = corr[x,y+2]
		corrinel[x][y] = corr[x+15,y+2]
		
# calculate the correlation coefficients
covel = np.cov(correl)
covinel = np.cov(corrinel)

# write array to file
#np.savetxt("angularcorr_el.txt",correl)
#np.savetxt("angularcorr_inel.txt",corrinel)
np.savetxt("covel.txt",covel)
np.savetxt("covinel.txt",covinel)