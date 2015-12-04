import matplotlib.pyplot as plt
import numpy as np
import sys

# input file is given in the command line
inputfile = sys.argv[1]

# load the data into the array corr
corr = np.loadtxt(inputfile)

colors = ['red','green']
for x in range(2,202):
	colors.append( 'blue' )

# for elastic scattering (specific to 12C @ 17 MeV)	
fig, arr = plt.subplots(16,16)

for x in range(0,16):
	for y in range(0,16):
		#strx = r'$\left (\frac{d\sigma}{d\Omega} \right )$ (angle ' + str(x+1) + ')'# set x label
		#stry = r'$\left (\frac{d\sigma}{d\Omega} \right )$ (angle ' + str(x+2) + ')'# set y label
		if x == y:
			# put the histograms along the diagonal elements of the array (for plotting)
			arr[y][x].hist(corr[x,:],bins=20)
			#arr[y][x].axis((0,180,0,200))
			#arr[y][x].axes.set_xlim([0,180])
			#arr[y][x].locator_params(nbins=2)
		else:
			#plt.scatter(corr[x,:],corr[y,:],c=colors,s=10)
			# off diagonal elements are the scatter plots
			arr[y][x].scatter(corr[x,:],corr[y,:],c=colors,s=1)
			#arr[y][x].axis((0,200,0,200))
			#arr[y][x].locator_params(nbins=2)
		# remove the white spaces between the plots
		fig.subplots_adjust(hspace=0,wspace=0)
# remove the axis marks from all of the subplots
plt.setp([a.get_yticklabels() for a in fig.axes[:]],visible=False)
plt.setp([a.get_xticklabels() for a in fig.axes[:]],visible=False)
plt.setp([a.get_xticklabels() for a in arr[15,:]],visible=True)
plt.setp([a.get_yticklabels() for a in arr[:,0]],visible=True)
#plt.set_title('Elastic Correlations in Parameter Space')
plt.show()

# for inelastic scattering (specific to 12C @ 17 MeV)
# same as above
fig, arr2 = plt.subplots(16,16)

for x in range(17,33):
	for y in range(17,33):
		if x == y:
			arr2[y-17][x-17].hist(corr[x,:],bins=20)
		else:
			arr2[y-17][x-17].scatter(corr[x,:],corr[y,:],c=colors,s=1)
	fig.subplots_adjust(hspace=0,wspace=0)
plt.setp([a.get_yticklabels() for a in fig.axes[:]],visible=False)
plt.setp([a.get_xticklabels() for a in fig.axes[:]],visible=False)
plt.setp([a.get_xticklabels() for a in arr2[15,:]],visible=True)
plt.setp([a.get_yticklabels() for a in arr2[:,0]],visible=True)
#plt.set_title('Inelastic Correlations in Parameter Space')
plt.show()