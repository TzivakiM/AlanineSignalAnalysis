#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 13:02:46 2021

@author: margarita
"""

'''
what will this code do?

1) import both files (p-p height and density)
2) normalize the p-p height by packing density, calculate errors.
3) subtract ppheight 0 from ppheight 3
4) assign a dose to each filename 
5) output: name,dose,ppheightdiff
6) plot dose-ppheigt, dose-pp-heightdiff
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import glob
import os
from scipy.optimize import curve_fit


ofile=open("alanine_alysis.csv","w+")
ofile.write('Sample name, Dose [Gy], p-p heught [a.u], Signal height Difference [a.u.]\n')
#ofile.write('name, Mpph, MppwG, Spph, SppwG, NSpph, NSppwG\n')
ofold = "graphs"
if not os.path.isdir(ofold):
    os.mkdir(ofold)

heightfile = "signalheights_edit.csv"
densityfile = "DensityFile_edit.csv"
#add some sort of Renaming, do manually now
'''
files=[heightfile,densityfile]
for file in files:
    filename=os.path.splitext(file)[0]
    print(filename)
'''

df_density = pd.read_csv(densityfile, sep=',',index_col='Name')
#in mg/mm^3
df_density.sort_index(inplace=True)
densitycols = df_density.columns

df_height = pd.read_csv(heightfile,skiprows=[1], sep=',',index_col='File name')
#in mg/mm^3
df_height.sort_index(inplace=True)

#that's where everything goes
df_master=pd.DataFrame()
#set the error for the intensity measurements, will have to be independently estimated
#deltaI = 0.01
deltaI=0.0

dose= [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,9.0,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]

colcount=0
df_master['Dose'] = dose
pph=df_height[' signal pp-height [au]']
pph_n=df_height[' normalized signal pp-height [au]']
df_master['pph']=pph.values
df_master['I_dn']=df_master['pph'].values/df_density[densitycols[colcount]].values
df_master['deltaI_dn']=np.sqrt((deltaI/df_master['pph'].values)**2 + (df_density[densitycols[colcount+1]].values/df_density[densitycols[colcount]].values)**2)*df_master['pph'].values

df_master['pph_norm']=pph_n.values
df_master['I_dn_norm']=df_master['pph_norm'].values/df_density[densitycols[colcount]].values
df_master['deltaI_dn_norm']=np.sqrt((deltaI/df_master['pph_norm'].values)**2 + (df_density[densitycols[colcount+1]].values/df_density[densitycols[colcount]].values)**2)*df_master['pph_norm'].values
    
#dump everything into an output file
df_master.to_csv('alanineresults.csv', index = False)


#PLOT SOME STUFF
##################################################
dosegroups=df_master.groupby('Dose')
averages=[]
doses=[]
stdevs=[]
df_analysis=pd.DataFrame()

for d in np.unique(df_master['Dose']):
    doses.append(d)
    av = (dosegroups.get_group(d)['I_dn'].mean())
    std = (dosegroups.get_group(d)['I_dn'].std())
    averages.append(av)
    stdevs.append(std)

df_analysis['Dose']=doses
df_analysis['Average']=averages
df_analysis['STDEV']=stdevs  

averages=np.array(averages)
doses=np.array(doses)
stdevs=np.array(stdevs)

def linear(x,m,t):
    return x*m+t

popt, pcov = curve_fit(linear,doses,averages)
plt.errorbar(doses,averages,yerr=stdevs,ls=' ',marker='o',ms=0.1,capsize=4.)
plt.xlabel('Dose [Gy]')
plt.ylabel('Pp-height difference [a.u.]')
plt.plot(doses, linear(doses,popt[0],popt[1]),'-')
plt.title('Peak-to-peak height averages\n all samples')
plt.ticklabel_format(style='sci', axis='y', scilimits=(6,6))
#plt.xlim(0.18,1.2)
#plt.ylim(-1,1000000)
plt.savefig('alanine.png',dpi=600)
plt.show()
