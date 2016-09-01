#######################################################################################################################################

# This program combines ADCP data and gets rid of the gaps and bad velocity data points

# USING VA VELOCITY IN LOOP

######################################################################################################################################

#Combining ADCP files

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
                adcp12.Variables.el=np.hstack((adcp12.Variables.el,adcp.Variables.el))
                adcp12.Variables.ua=np.hstack((adcp12.Variables.ua,adcp.Variables.ua))
                adcp12.Variables.va=np.hstack((adcp12.Variables.va,adcp.Variables.va))
                adcp12.Variables.matlabTime=np.hstack((adcp12.Variables.matlabTime,adcp.Variables.matlabTime))
		
	if index != 0:
                adcps.Variables.el=np.hstack((adcps.Variables.el,adcp.Variables.el))
                adcps.Variables.ua=np.hstack((adcps.Variables.ua,adcp.Variables.ua))
                adcps.Variables.va=np.hstack((adcps.Variables.va,adcp.Variables.va))
                adcps.Variables.matlabTime=np.hstack((adcps.Variables.matlabTime,adcp.Variables.matlabTime))

modifiedtime=np.arange(adcps.Variables.matlabTime[0],adcps.Variables.matlabTime[-1],(adcps.Variables.matlabTime[1]-adcps.Variables.matlabTime[0]))
newadcp=cp.deepcopy(adcp)
newadcp.Variables.matlabTime=modifiedtime


# Velocity

# Do u and v velocities together in same loop

harmo3=adcp12.Utils.Harmonic_analysis(elevation=False,velocity=True)
velos3=adcps.Utils.Harmonic_reconstruction(harmo3)
diff_ua=abs(adcps.Variables.ua - velos3['u'])
diff_va=abs(adcps.Variables.va - velos3['v'])


# The following loop takes 110 mins to run and removes 7522 'bad' data points 
while (diff_va > (np.mean(diff_va)+2.5*np.std(diff_va))).any():
        harmo3=(adcps.Utils.Harmonic_analysis(elevation=False,velocity=True))
        velos3=adcps.Utils.Harmonic_reconstruction(harmo3)
        diff_ua = abs(adcps.Variables.ua - velos3['u'])
        diff_va = abs(adcps.Variables.va - velos3['v'])
        y = np.where(diff_va>(np.mean(diff_va)+2.5*np.std(diff_va)))
        diff_ua=np.delete(diff_ua,y)
        diff_va=np.delete(diff_va,y)
        adcps.Variables.ua = np.delete(adcps.Variables.ua,y)
        adcps.Variables.va = np.delete(adcps.Variables.va,y)
        adcps.Variables.matlabTime = np.delete(adcps.Variables.matlabTime,y)
        print shape(diff_ua)
        print shape(diff_va)
        print shape(y)

harmo3=adcps.Utils.Harmonic_analysis(elevation=False,velocity=True)
velos3=adcps.Utils.Harmonic_reconstruction(harmo3)
diff_ua=abs(adcps.Variables.ua-velos3['u'])
diff_va=abs(adcps.Variables.va-velos3['v'])
plt.plot(adcps.Variables.matlabTime,diff_ua,'.')
show()
plt.plot(adcps.Variables.matlabTime,diff_va,'.')
show()

#Normalizing residuals

residual=(adcps.Variables.ua-velos3['u'])/(.5*abs(adcps.Variables.ua)+abs(velos3['u']))

mr=np.mean(residual)
sdr=np.std(residual)
#plt.plot(residual)
#show()
# Hist
#plt.hist(residual,bins=100)
#show()

# Filling in the gaps

velos=newadcp.Utils.Harmonic_reconstruction(harmo3)
newadcp.Variables.ua=velos['u']
newadcp.Variables.va=velos['v']
#plt.plot(adcps.Variables.matlabTime,adcps.Variables.ua,'.')
#show()
plt.plot(adcps.Variables.matlabTime,adcps.Variables.va,'.')
show()

plt.plot(modifiedtime,velos['u'])
show()
plt.plot(modifiedtime,velos['v'])
show()

# find filled in points and add noise
# might be easier to find where matlabTime != modifiedtime

for i in modifiedtime:
        index = np.where((adcps.Variables.matlabTime == modifiedtime))
        if index == []:
                newadcp.Variables.ua[i]=newadcp.Variables.ua[i]*(1+np.random.normal(loc=mr,scale=sdr))
#plt.plot(modifiedtime,velos['h'])
#show()

############################################################################################################################

