import numpy as np
from scipy import linalg as la
import sys

#input file
wijfile = sys.argv[1]
errfile = sys.argv[2]

# load the data
wij = np.loadtxt(wijfile)
err = np.loadtxt(errfile)

errt = err.transpose()
print err
print errt

#la.matmul(wij,err,temp)
#la.matmul(temp,errt,chi)

print chi