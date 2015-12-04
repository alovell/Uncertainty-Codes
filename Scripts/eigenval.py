import numpy as np
from scipy import linalg as la
import sys

# input file given in command line
inputfile = sys.argv[1]
sigfile = sys.argv[2]

# load the data
corr = np.loadtxt(inputfile)
sig = np.loadtxt(sigfile)

# add sigma squared to the diagonal correlation elements
for x in range(0,17):
	#print corr[x][x],sig[x]
	corr[x][x] = corr[x][x] + sig[x]**2
	#print corr[x][x]

# calculate the eigenvalue and eigenvectors
evals,evecs = la.eig(corr)

# print eigenvalues
print evals

# compute the inverse of corr and find eigenvealues
icorr = la.inv(corr)

ievals,ievecs = la.eig(icorr)

print ievals