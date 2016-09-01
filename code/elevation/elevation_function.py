
# To execute this function in command line:
# 1. ipython
# 2. from elevation_function.py import elevation
# 3. elevation('path',num)


def elevation(dir, num, Graphics=True, debug = False):
	"""
	This function is for ADCP data that has bad data and gaps in the time series. It combines all ADCP data from a given path, finds and removes the bad elevation data, and then fills in the gaps in the time series.
	
	Inputs:
	  -  dir = directory/path of ADCP file(s)
	  -  num = number of standard deviations from mean requested when removing bad data points
	  -  Graphics=True means that the final figures and residuals will be displayed

	Output:
	  -  Final figure with only good data and no missing gaps

	"""

# Standard deviation messages
	
	if num > 4:
		print "Entered standard deviation not recommended - few data points will be removed"
	if 4 > num > 3:
		print "Program will take 5-50 mins to run"
	if 3 > num > 2:
		print "Program will take ~ 30 mins to run"
	if 2 > num:
		print "Entered standard deviation not recommended - chance of infinite loop" 
	
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

#dir = '/EcoII/acadia_uni/workspace/observed/GP/ADCP/GPcabled/'
	
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
                      	adcp12.Variables.matlabTime=np.hstack((adcp12.Variables.matlabTime,adcp.Variables.matlabTime))

        	if index != 0:
                	adcps.Variables.el=np.hstack((adcps.Variables.el,adcp.Variables.el))
                       	adcps.Variables.matlabTime=np.hstack((adcps.Variables.matlabTime,adcp.Variables.matlabTime))


	
	modifiedtime=np.arange(adcps.Variables.matlabTime[0],adcps.Variables.matlabTime[-1],(adcps.Variables.matlabTime[1]-adcps.Variables.matlabTime[0]))
	newadcp=cp.deepcopy(adcp)
	newadcp.Variables.matlabTime=modifiedtime

############################################################################################################################################

# Elevation


# Perform Harmonic Analysis of original elevation data
	
	harmo1=adcp12.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)

# Reconstruct with original measured data
	
	velos1=adcps.Utils.Harmonic_reconstruction(harmo1)

# Difference
	
	diff_el = abs(adcps.Variables.el - velos1['h'])


	harmo2=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)
	velos2=adcps.Utils.Harmonic_reconstruction(harmo2)
		

# The following while loop removes 'bad' data points
	
	while (diff_el > (np.mean(diff_el)+num*np.std(diff_el))).any():
        	harmo2=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)
        	velos2=adcps.Utils.Harmonic_reconstruction(harmo2)
        	diff_el = abs(adcps.Variables.el - velos2['h'])
        	x = np.where(diff_el>(np.mean(diff_el)+num*np.std(diff_el)))
        	diff_el=np.delete(diff_el,x)
        	adcps.Variables.el = np.delete(adcps.Variables.el,x)
        	adcps.Variables.matlabTime = np.delete(adcps.Variables.matlabTime,x)
        	print shape(diff_el)
        	print shape(x)
	print "Loop finished"
	
	harmo2=adcps.Utils.Harmonic_analysis(elevation=True,velocity=False,trend=False)
	velos2=adcps.Utils.Harmonic_reconstruction(harmo2)
	diff_el=abs(adcps.Variables.el-velos2['h'])

	#plt.plot(adcps.Variables.matlabTime,diff_el,'.')
	#show()


#Normalizing residuals

	residual=(adcps.Variables.el-velos2['h'])/(.5*abs(adcps.Variables.el)+abs(velos2['h']))

	mr=np.mean(residual)
	sdr=np.std(residual)
	#plt.plot(residual)
	#show()
# Hist
	#plt.hist(residual,bins=100)
	#show()

# Filling in the gaps
	
	velos=newadcp.Utils.Harmonic_reconstruction(harmo2)
	newadcp.Variables.el=velos['h']
	#plt.plot(adcps.Variables.matlabTime,adcps.Variables.el,'.')
	#plt.plot(modifiedtime,velos['h'])
	#show()

# find filled in points and add noise


	#for i in modifiedtime:
        	#index = np.where((adcps.Variables.matlabTime == modifiedtime))
        	#if index == []:
                	#newadcp.Variables.el[i]=newadcp.Variables.el[i]*(1+np.random.normal(loc=mr,scale=sdr))

#test final reconstruct
	
	#velos=newadcp.Utils.Harmonic_reconstruction(harmo2)
	#newadcp.Variables.el=velos['h']
	#plt.plot(modifiedtime,velos['h'])
	#show()

	if Graphics:
		plt.plot(residual)
		show()
		plt.hist(residual,bins=100)
		show()
		plt.plot(modifiedtime,velos['h'])
		show()
