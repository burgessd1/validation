#########################################################################################################################################

# This program combines ADCP data and gets rid of the gaps and bad surf data points

############################################################################################################################################

# Combining ADCP files

from pylab import *
import os
from pyseidon_dvt import *
import numpy as np
from utide import  *
from scipy.interpolate import interp1d
import csv
import matplotlib.pyplot as plt
import copy as cp

# Initialize ADCP data

dir = '/EcoII/acadia_uni/workspace/observed/GP/ADCP/GPcabled/'
list = os.listdir(dir)
list.sort()
listADCP = []

for i in list:
        if ".mat" in i:
                 listADCP.append(dir+i)

for index in range(len(listADCP)):
        adcp=ADCP(listADCP[index])

        if index == 0:
                adcps=cp.deepcopy(adcp)
                adcp12=cp.deepcopy(adcp)

        elif index == 1:
                adcp12.Variables.surf=np.hstack((adcp12.Variables.surf,adcp.Variables.surf))
                adcp12.Variables.matlabTime=np.hstack((adcp12.Variables.matlabTime,adcp.Variables.matlabTime))

        if index != 0:
                adcps.Variables.surf=np.hstack((adcps.Variables.surf,adcp.Variables.surf))
                adcps.Variables.matlabTime=np.hstack((adcps.Variables.matlabTime,adcp.Variables.matlabTime))


modifiedtime=np.arange(adcps.Variables.matlabTime[0],adcps.Variables.matlabTime[-1],(adcps.Variables.matlabTime[1]-adcps.Variables.matlabTime[0]))
newadcp=cp.deepcopy(adcp)
newadcp.Variables.matlabTime=modifiedtime

############################################################################################################################################

# Surf
# Harmonic Analysis does not work with anything other than elevation and velocity :((((((((((( 

# Perform Harmonic Analysis of Elevation
harmo1=adcp12.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)

# Reconstruct with original measured data
velos1=adcps.Utils.Harmonic_reconstruction(harmo1)

# Difference
diff_surf = abs(adcps.Variables.surf - velos1['h'])


harmo2=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)
velos2=adcps.Utils.Harmonic_reconstruction(harmo2)

# The following while loop takes mins to run and removes 'bad' data points
while (diff_surf > (np.mean(diff_surf)+2.5*np.std(diff_surf))).any():
        harmo2=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)
        velos2=adcps.Utils.Harmonic_reconstruction(harmo2)
        diff_surf = abs(adcps.Variables.surf - velos2['h'])
        x = np.where(diff_surf>(np.mean(diff_surf)+2.5*np.std(diff_surf)))
        diff_surf=np.delete(diff_surf,x)
        adcps.Variables.surf = np.delete(adcps.Variables.surf,x)
        adcps.Variables.matlabTime = np.delete(adcps.Variables.matlabTime,x)
        print shape(diff_surf)
        print shape(x)

harmo2=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)
velos2=adcps.Utils.Harmonic_reconstruction(harmo2)
diff_surf=abs(adcps.Variables.surf-velos2['h'])

#plt.plot(adcps.Variables.matlabTime,diff_surf,'.')
#show()


#Normalizing residuals

residual=(adcps.Variables.surf-velos2['h'])/(.5*abs(adcps.Variables.surf)+abs(velos2['h']))

mr=np.mean(residual)
sdr=np.std(residual)
plt.plot(residual)
show()
# Hist
plt.hist(residual,bins=100)
show()

# Filling in the gaps
velos=newadcp.Utils.Harmonic_reconstruction(harmo2)
newadcp.Variables.surf=velos['h']
#plt.plot(adcps.Variables.matlabTime,adcps.Variables.surf,'.')
#plt.plot(modifiedtime,velos['h'])
#show()

# find filled in points and add noise


for i in modifiedtime:
        index = np.where((adcps.Variables.matlabTime == modifiedtime))
        if index == []:
                newadcp.Variables.surf[i]=newadcp.Variables.surf[i]*(1+np.random.normal(loc=mr,scale=sdr))

#test final reconstruct
velos=newadcp.Utils.Harmonic_reconstruction(harmo2)
newadcp.Variables.surf=velos['h']
plt.plot(modifiedtime,velos['h'])
show()


