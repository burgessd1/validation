##################################################################################################################################

# FUNCTIONS


# Table of Contents

# 1. elevation

# 2. velocity


##################################################################################################################################


# To execute functions in command line:
# 1. ipython
# 2. from *function*.py import function
# 3. function('path',num)

# 1. 










































# 2.

def velocity(dir, num, Graphics=True, debug = False):
        """
        This function is for ADCP data that has bad data and gaps in the time series. It combines all ADCP data from a given path, finds and removes the bad velocity data, and then fills in the gaps in the time series.

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
                print "Program will take less than 85 mins to run"
        if 3 > num > 2:
                print "Program will take ~ 85 mins to run"
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
                        adcp12.Variables.ua=np.hstack((adcp12.Variables.ua,adcp.Variables.ua))
                        adcp12.Variables.va=np.hstack((adcp12.Variables.va,adcp.Variables.va))
                        adcp12.Variables.matlabTime=np.hstack((adcp12.Variables.matlabTime,adcp.Variables.matlabTime))

                if index != 0:
                        adcps.Variables.ua=np.hstack((adcps.Variables.ua,adcp.Variables.ua))
                        adcps.Variables.va=np.hstack((adcps.Variables.va,adcp.Variables.va))
                        adcps.Variables.matlabTime=np.hstack((adcps.Variables.matlabTime,adcp.Variables.matlabTime))
        modifiedtime=np.arange(adcps.Variables.matlabTime[0],adcps.Variables.matlabTime[-1],(adcps.Variables.matlabTime[1]-adcps.Variables.matlabTime[0]))
        newadcp=cp.deepcopy(adcp)
        newadcp.Variables.matlabTime=modifiedtime

############################################################################################################################################

# Velocity


# Perform Harmonic Analysis of original elevation data

        harmo3=adcp12.Utils.Harmonic_analysis(elevation=False,velocity=True)

# Reconstruct with original measured data

        velos3=adcps.Utils.Harmonic_reconstruction(harmo3)

# Difference

        diff_ua = abs(adcps.Variables.ua - velos3['u'])
        diff_va = abs(adcps.Variables.va - velos3['v'])
        diff = sqrt((diff_ua)**2+(diff_va)**2)

        #plt.plot(adcps.Variables.matlabTime,diff,'.')
        #show()

        harmo3=adcps.Utils.Harmonic_analysis(elevation=False,velocity=True)
        velos3=adcps.Utils.Harmonic_reconstruction(harmo3)


# The following loop removes 'bad' data points

        while (diff > (np.mean(diff)+num*np.std(diff))).any():
                harmo3=(adcps.Utils.Harmonic_analysis(elevation=False,velocity=True))
                velos3=adcps.Utils.Harmonic_reconstruction(harmo3)
                diff_ua = abs(adcps.Variables.ua - velos3['u'])
                diff_va = abs(adcps.Variables.va - velos3['v'])
                diff = sqrt(((diff_va)*(diff_va))+((diff_ua)*(diff_ua)))
                y = np.where((diff > (np.mean(diff)+num*np.std(diff))))
                diff_ua=np.delete(diff_ua,y)
                diff_va=np.delete(diff_va,y)
                diff=np.delete(diff,y)
                adcps.Variables.ua = np.delete(adcps.Variables.ua,y)
                adcps.Variables.va = np.delete(adcps.Variables.va,y)
                adcps.Variables.matlabTime = np.delete(adcps.Variables.matlabTime,y)
                print shape(diff)
                print shape(y)
        print "Loop finished"

        harmo3=adcps.Utils.Harmonic_analysis(elevation=False,velocity=True)
        velos3=adcps.Utils.Harmonic_reconstruction(harmo3)
        diff_ua=abs(adcps.Variables.ua-velos3['u'])
        diff_va=abs(adcps.Variables.va-velos3['v'])
        diff=sqrt((diff_ua)**2+(diff_va)**2)

        #plt.plot(adcps.Variables.matlabTime,diff,'.')
        #show()
        #plt.plot(adcps.Variables.matlabTime,diff_ua,'.')
        #show()
        #plt.plot(adcps.Variables.matlabTime,diff_va,'.')
        #show()
#Normalizing residuals

        uresidual=(adcps.Variables.ua-velos3['u'])/(.5*abs(adcps.Variables.ua)+abs(velos3['u']))
        vresidual=(adcps.Variables.va-velos3['v'])/(.5*abs(adcps.Variables.va)+abs(velos3['v']))
        mu=np.mean(uresidual)
        mv=np.mean(vresidual)
        sdu=np.std(uresidual)
        sdv=np.std(vresidual)

# Filling in the gaps

        velos=newadcp.Utils.Harmonic_reconstruction(harmo3)
        newadcp.Variables.ua=velos['u']
        newadcp.Variables.va=velos['v']
        #plt.plot(adcps.Variables.matlabTime,adcps.Variables.ua,'.')
        #show()
        #plt.plot(adcps.Variables.matlabTime,adcps.Variables.va,'.')
        #show()

# find filled in points and add noise


        for i in modifiedtime:
                index = np.where((adcps.Variables.matlabTime[i] == modifiedtime[i]))
                if index == []:
                        newadcp.Variables.va[i]=newadcp.Variables.va[i]*(1+np.random.normal(loc=mv,scale=sdv))
                        newadcp.Variables.ua[i]=newadcp.Variables.ua[i]*(1+np.random.normal(loc=mu,scale=sdu))

        velos=newadcp.Utils.Harmonic_reconstruction(harmo3)
        newadcp.Variables.ua=velos['u']
        newadcp.Variables.va=velos['v']

        if Graphics:
                plt.plot(uresidual)
                show()
                plt.plot(vresidual)
                show()
                plt.hist(uresidual,bins=100)
                show()
                plt.hist(vresidual,bins=100)
                plt.plot(modifiedtime,velos['u'])
                show()
                plt.plot(modifiedtime,velos['v'])
                show()

