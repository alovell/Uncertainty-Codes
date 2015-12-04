from numpy.random import uniform, seed
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np

#this needs to be fixed to account for the new layout of the data
#want to read in the two numbers and then only save those two columns
#read in data instead
#cd = np.loadtxt('chi1.txt')
#a = cd[:,0]
#b = cd[:,1]
#c = cd[:,2]
#d = cd[:,3]
#e = cd[:,4]
#f = cd[:,5]
#z = cd[:,6]

#load files (defined by parse.pl)
g12 = np.loadtxt('fi12.txt') 
'''
g13 = np.loadtxt('file13.txt')
g14 = np.loadtxt('file14.txt')
g15 = np.loadtxt('file15.txt')
g16 = np.loadtxt('file16.txt')
g23 = np.loadtxt('file23.txt')
g24 = np.loadtxt('file24.txt')
g25 = np.loadtxt('file25.txt')
g26 = np.loadtxt('file26.txt')
g34 = np.loadtxt('file34.txt')
g35 = np.loadtxt('file35.txt')
g36 = np.loadtxt('file36.txt')
g45 = np.loadtxt('file45.txt')
g46 = np.loadtxt('file46.txt')
g56 = np.loadtxt('file56.txt')
'''

#create arrays from files
a12 = g12[:,0]
b12 = g12[:,1]
z12 = g12[:,6]
'''
a13 = g13[:,0]
c13 = g13[:,1]
z13 = g13[:,2]
a14 = g14[:,0]
d14 = g14[:,1]
z14 = g14[:,2]
a15 = g15[:,0]
e15 = g15[:,1]
z15 = g15[:,2]
a16 = g16[:,0]
f16 = g16[:,1]
z16 = g16[:,2]
b23 = g23[:,0]
c23 = g23[:,1]
z23 = g23[:,2]
b24 = g24[:,0]
d24 = g24[:,1]
z24 = g24[:,2]
b25 = g25[:,0]
e25 = g25[:,1]
z25 = g25[:,2]
b26 = g26[:,0]
f26 = g26[:,1]
z26 = g26[:,2]
c34 = g34[:,0]
d34 = g34[:,1]
z34 = g34[:,2]
c35 = g35[:,0]
e35 = g35[:,1]
z35 = g35[:,2]
c36 = g36[:,0]
f36 = g36[:,1]
z36 = g36[:,2]
d45 = g45[:,0]
e45 = g45[:,1]
z45 = g45[:,2]
d46 = g46[:,0]
f46 = g46[:,1]
z46 = g46[:,2]
e56 = g56[:,0]
f56 = g56[:,1]
z56 = g56[:,2]
'''
bf = [162.39593115407664,0.97968354662442858,48.716900173448359,1.1106510638766456,0.68617337297244707,0.30365660123208432]

#define grid
ai = np.linspace(160.197788788350,164.497113660652,11)
bi = np.linspace(0.970854922929466, 0.970854922929466,11)
ci = np.linspace(47.579443640415100,49.854356706481617,11)
di = np.linspace(1.1084824639436155,1.1130365238029787,11)
ei = np.linspace(0.67702749114427185,0.69440466661780498,11)
fi = np.linspace(0.29934955424817056,0.30839435291438888,11)

#ai = np.linspace(160.78368432696311,164.54559359022787,8)
#bi = np.linspace(0.96838871070340948,0.98871941536124386,10)
#ci = np.linspace(47.579443640415100,49.854356706481617,9)
#di = np.linspace(1.1084824639436155,1.1128196638096757,9)
#ei = np.linspace(0.67702749114427185,0.69531925480062229,9)
#fi = np.linspace(0.29934955424817056,0.30796364821599809,9)


#grid data
zab = griddata(a12,b12,z12,ai,bi,interp='linear')
'''
zac = griddata(a13,c13,z13,ai,ci,interp='linear')
zad = griddata(a14,d14,z14,ai,di,interp='linear')
zae = griddata(a15,e15,z15,ai,ei,interp='linear')
zaf = griddata(a16,f16,z16,ai,fi,interp='linear')
zbc = griddata(b23,c23,z23,bi,ci,interp='linear')
zbd = griddata(b24,d24,z24,bi,di,interp='linear')
zbe = griddata(b25,e25,z25,bi,ei,interp='linear')
zbf = griddata(b26,f26,z26,bi,fi,interp='linear')
zcd = griddata(c34,d34,z34,ci,di,interp='linear')
zce = griddata(c35,e35,z35,ci,ei,interp='linear')
zcf = griddata(c36,f36,z36,ci,fi,interp='linear')
zde = griddata(d45,e45,z45,di,ei,interp='linear')
zdf = griddata(d46,f46,z46,di,fi,interp='linear')
zef = griddata(e56,f56,z56,ei,fi,interp='linear')
'''
#my contour plots
abCS = plt.contour(ai,bi,zab,9,linewidths=0.5)
abCB = plt.colorbar(abCS,shrink=0.8,extend='both')
plt.title('V-r contour')
plt.xlabel('V')
plt.ylabel('r')
plt.plot(bf[0],bf[1],'ko')
plt.show()
'''
acCS = plt.contour(ai,ci,zac,9,linewidths=0.5)
acCB = plt.colorbar(acCS,shrink=0.8,extend='both')
plt.title('V-Ws contour')
plt.xlabel('V')
plt.ylabel('Ws')
plt.plot(bf[0],bf[2],'ko')
plt.show()
adCS = plt.contour(ai,di,zad,9,linewidths=0.5)
adCB = plt.colorbar(adCS,shrink=0.8,extend='both')
plt.title('V-rs contour')
plt.xlabel('V')
plt.ylabel('rs')
plt.plot(bf[0],bf[3],'ko')
plt.show()
aeCS = plt.contour(ai,ei,zae,9,linewidths=0.5)
aeCB = plt.colorbar(aeCS,shrink=0.8,extend='both')
plt.title('V-a contour')
plt.xlabel('V')
plt.ylabel('a')
plt.plot(bf[0],bf[4],'ko')
plt.show()
afCS = plt.contour(ai,fi,zaf,9,linewidths=0.5)
afCB = plt.colorbar(afCS,shrink=0.8,extend='both')
plt.title('V-as contour')
plt.xlabel('V')
plt.ylabel('as')
plt.plot(bf[0],bf[5],'ko')
plt.show()
bcCS = plt.contour(bi,ci,zbc,9,linewidths=0.5)
bcCB = plt.colorbar(bcCS,shrink=0.8,extend='both')
plt.title('r-Ws contour')
plt.xlabel('r')
plt.ylabel('Ws')
plt.plot(bf[1],bf[2],'ko')
plt.show()
bdCS = plt.contour(bi,di,zbd,9,linewidths=0.5)
bdCB = plt.colorbar(bdCS,shrink=0.8,extend='both')
plt.title('r-rs contour')
plt.xlabel('r')
plt.ylabel('rs')
plt.plot(bf[1],bf[3],'ko')
plt.show()
beCS = plt.contour(bi,ei,zbe,9,linewidths=0.5)
beCB = plt.colorbar(beCS,shrink=0.8,extend='both')
plt.title('r-a contour')
plt.xlabel('r')
plt.ylabel('a')
plt.plot(bf[1],bf[4],'ko')
plt.show()
bfCS = plt.contour(bi,fi,zbf,9,linewidths=0.5)
bfCB = plt.colorbar(bfCS,shrink=0.8,extend='both')
plt.title('r-as contour')
plt.xlabel('r')
plt.ylabel('as')
plt.plot(bf[1],bf[5],'ko')
plt.show()
cdCS = plt.contour(ci,di,zcd,9,linewidths=0.5)
cdCB = plt.colorbar(cdCS,shrink=0.8,extend='both')
plt.title('Ws-rs contour')
plt.xlabel('Ws')
plt.ylabel('rs')
plt.plot(bf[2],bf[3],'ko')
plt.show()
ceCS = plt.contour(ci,ei,zce,9,linewidths=0.5)
ceCB = plt.colorbar(ceCS,shrink=0.8,extend='both')
plt.title('Ws-a contour')
plt.xlabel('Ws')
plt.ylabel('a')
plt.plot(bf[2],bf[4],'ko')
plt.show()
cfCS = plt.contour(ci,fi,zcf,9,linewidths=0.5)
cfCB = plt.colorbar(cfCS,shrink=0.8,extend='both')
plt.title('Ws-as contour')
plt.xlabel('Ws')
plt.ylabel('as')
plt.plot(bf[2],bf[5],'ko')
plt.show()
deCS = plt.contour(di,ei,zde,9,linewidths=0.5)
deCB = plt.colorbar(deCS,shrink=0.8,extend='both')
plt.title('rs-a contour')
plt.xlabel('rs')
plt.ylabel('a')
plt.plot(bf[3],bf[4],'ko')
plt.show()
dfCS = plt.contour(di,fi,zdf,9,linewidths=0.5)
dfCB = plt.colorbar(dfCS,shrink=0.8,extend='both')
plt.title('rs-as contour')
plt.xlabel('rs')
plt.ylabel('a')
plt.plot(bf[3],bf[5],'ko')
plt.show()
efCS = plt.contour(ei,fi,zef,9,linewidths=0.5)
efCB = plt.colorbar(efCS,shrink=0.8,extend='both')
plt.title('a-as contour')
plt.xlabel('a')
plt.ylabel('as')
plt.plot(bf[4],bf[5],'ko')
plt.show()
'''