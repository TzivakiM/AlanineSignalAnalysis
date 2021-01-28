#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 19:00:33 2019

@author: margarita
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import glob
from bisect import bisect_left, bisect_right
import os

#note: will need to include also height for standard...
#this is the file where everything gets written into
ofile=open("signalheights.csv","w+")
ofile.write('File name, standard pp-half-height b.n. [au], standard peak width b.n.[G], signal pp-height [au], signal peak width [G], normalized signal pp-height [au], normalized signal peak width [G]\n')
ofile.write('name, Mpph, MppwG, Spph, SppwG, NSpph, NSppwG\n')
ofold = "graphs"
if not os.path.isdir(ofold):
    os.mkdir(ofold)

standard = 'yes' #was a standard used in the measurements?
############# STANDARD DATA ####################
#set the values between which the magnetic field corresponding to the standard peak
# is expected
high =  3560.#B-field high value of standard
low = 3538.#B-field low value of standard
gstand = 1.9800 #gfactor of the standard
gstandB = 3540. #Gauss Bfield corresponding to the g-factor in Gauss
#TODO: find the correct value

#set the values between which the magnetic field is expected to be for the signal
s_high = 3510. 
s_low = 3475.

color1 = '#720000'
color2 = '#214f21'

for file in glob.glob("*.asc"):
    #use to look at only one file
    #if file != "AL-01-1.asc":
       #continue
    if standard == 'yes':
    ##################DEAL WITH THE STANDARD FIRST#########################
    #find the maximum between Magnetic field values of high and low
    #those are the spectra files
        data = pd.read_csv(file, skiprows=2, sep='\t')
        #find the index of the bounds
        indlow = data[data['X [G]'].gt(low)].index[0] 
        indhigh = data[data['X [G]'].gt(high)].index[0] 
        #crop dataset to the bounds
        crdata = data[indlow:indhigh]
        
        #find max of the intensity value
        maxI= crdata['Intensity'].max()
        #find the corresponding magnetic field value
        maxB = crdata.loc[crdata['Intensity'].idxmax(),'X [G]']
    
        #find min of the intensity value
        minI= crdata['Intensity'].min()
        #find the corresponding magnetic field value
        minB = crdata.loc[crdata['Intensity'].idxmin(),'X [G]']
        #print(minI)
        #print(minB)
        

        #calculate the pp signal width in Gauss and mT
        ppwG = abs(minB-maxB)#abs needed here because the standard is 180degrees flipped
        ppwmT = 0.1*(ppwG)
        
        #finding the center of the standard
        standC=(maxB+minB)/2
        #shifting everything to the standard g-value
        shift = standC - gstandB
        data_transp = pd.DataFrame()
        data_transp['X [G]'] = data['X [G]'] - shift
        data_transp['Intensity'] = data['Intensity']
        #also shift the Centre value of the standard width!
        standC=standC-shift
        
        #calculate the 1/2 pp signal height
        #get the closest points
        '''    
        Function of bisect:
    - In the case where the value val is in the column, 
    bisect_left will return the precise index of the value in 
    the list and bisect_right will return the index of the next position.
    - In the case where the value is not in the list,
    both bisect_left and bisect_right will return the same index: 
    the one where to insert the value to keep the list sorted.
    '''
        def get_closests(df, col, val):
            lower_idx = bisect_left(df[col].values, val)
            higher_idx = bisect_right(df[col].values, val)
            if higher_idx == lower_idx:      #val is not in the list
                return lower_idx - 1, lower_idx
            else:                            #val is in the list
                return lower_idx, lower_idx
        
        halfwidth_ind= get_closests(data_transp, 'X [G]', standC)
        halfwidth_Xrange=(data_transp.iloc[halfwidth_ind[0],0], data_transp.iloc[halfwidth_ind[1],0])
        halfwidth_Intrange=(data_transp.iloc[halfwidth_ind[0],1], data_transp.iloc[halfwidth_ind[1],1])
        
        #interpolate between values:
        interp=pd.Series([halfwidth_Intrange[0],np.nan, halfwidth_Intrange[1]], index=[halfwidth_Xrange[0], standC, halfwidth_Xrange[1]])
        
        interpolated_halfwidthInt=(interp.interpolate(method='index')).iloc[1]
        print(interp)
        print(interpolated_halfwidthInt)
        
        pphh = maxI-interpolated_halfwidthInt #this is NOT peak-to peak, this is half height to positive peak
        print(pphh)

        
        #define the hight of the standard to be 1 and scale everything based on this
        data_scaled = pd.DataFrame()
        data_scaled['X [G]'] = data_transp['X [G]']
        data_scaled['Intensity'] = data_transp['Intensity']*(1/(2*pphh))*1000000
        #print(data_scaled[data['Intensity']==maxI])
        #print(data_scaled[data['Intensity']==minI])
        
        
        plt.plot(data['X [G]'],data['Intensity'], 'k')
        #plt.vlines([minB,maxB,high,low],-1500000,1000000)
        plt.plot(data_transp['X [G]'],data_transp['Intensity'],color=color1)
        
        plt.plot(data_scaled['X [G]'],data_scaled['Intensity'],color=color2)
        #plt.vlines([s_high,s_low],-1500000,1000000)
        #plt.show()
        
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        #plt.title(filename)
        #plt.vlines([standC],-6e5,6e5)
        plt.savefig(os.path.join('graphs',file[:-4]+'-scaling.svg'))
        plt.show()

        
        #workingdata = data_transp        
        
    
    elif standard == 'no':
        print('We will continue without a standard')
        data = pd.read_csv(file, skiprows=2, sep='\t')
        
        plt.plot(data['X [G]'],data['Intensity'], 'k')
        #plt.vlines([minB,maxB,high,low],-1500000,1000000)

        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        #plt.title(filename)
        plt.savefig(file[:-4]+'-scaling.svg')
        plt.show()
        pphh='nan'
        ppwG='nan'
    
    ####################### SIGNAL ANALYSIS ##############################
    
    #now find the index of the location of the peaks with and without normalization
    def peakHeights(dataframe):
        
        s_indlow = dataframe[dataframe['X [G]'].gt(s_low)].index[0] 
        s_indhigh = dataframe[dataframe['X [G]'].gt(s_high)].index[0] 
        #print(dataframe)
        #crop dataset to the bounds
        crdataframe = dataframe[s_indlow:s_indhigh]
        s_maxI= crdataframe['Intensity'].max()
        #find max of the intensity value
        #find the corresponding magnetic field value
        s_maxB = crdataframe.loc[crdataframe['Intensity'].idxmax(),'X [G]']
        #find min of the intensity value
        s_minI= crdataframe['Intensity'].min()
        #find the corresponding magnetic field value
        s_minB = crdataframe.loc[crdataframe['Intensity'].idxmin(),'X [G]']
        #print(s_minI,s_maxI)
        #print(s_minB,s_maxB)
        #plt.plot(crdataframe['X [G]'],crdataframe['Intensity'])
        #plt.vlines([s_minB,s_maxB],-1500000,1000000)
        
        #calculate the pp signal height
        s_pph = abs(s_maxI-s_minI)
        #calculate the pp signal width in Gauss and mT
        s_ppwG = s_minB-s_maxB
        s_ppwmT = 0.1*(s_ppwG)
        return s_pph,s_ppwG,s_ppwmT
    
    if standard=='no':
        s_pphN='nan'
        s_ppwGN='nan'
        s_ppwmTN='nan'
        s_pph,s_ppwG,s_ppwmT=peakHeights(data)
    
    elif standard == 'yes':
        s_pphN,s_ppwGN,s_ppwmTN=peakHeights(data_scaled)
        s_pph,s_ppwG,s_ppwmT=peakHeights(data_transp)
    
    
    #ofile.write('{:},{:},{:}\n'.format(file,s_pph,s_ppwG))
    ofile.write('{:},{:},{:},{:},{:},{:},{:}\n'.format(file,pphh,ppwG,s_pph,s_ppwG,s_pphN,s_ppwGN))
    
    print(file)
    
ofile.close()





