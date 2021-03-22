#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 13:02:46 2021

@author: margarita
"""

'''
what will this code do?
file 0 is the background
file 3 is the measurement
1) import both files
2) use the normalized ppheight
3) subtract ppheight 0 from ppheight 3
4) assign a dose to each filename 1-6:0.2, etc up to 37-42:5Gy
5) output: name,dose,ppheightdiff
6) plot dose-ppheigt, dose-pp-heightdiff
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import glob
import os
from scipy.optimize import curve_fit

mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = "cm"


ofold = "graphs"
if not os.path.isdir(ofold):
    os.mkdir(ofold)

heightfile = "signalheights_edit.csv"
densityfile = "DensityFile_edit.csv"

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

number=[]


colcount=0
df_master['Sample Identifier']=df_density.index
for n in range(len(df_master['Sample Identifier'])):
    num=df_master['Sample Identifier'][n][3:]
    num=num[:2]
    number.append(num)
df_master['Number']=number
df_master['Dose'] = dose
pph=df_height[' signal pp-height [au]']
pph_n=df_height[' normalized signal pp-height [au]']
df_master['pph']=pph.values
df_master['I_dn']=df_master['pph'].values/df_density[densitycols[colcount]].values
df_master['deltaI_dn']=np.sqrt((deltaI/df_master['pph'].values)**2 + (df_density[densitycols[colcount+1]].values/df_density[densitycols[colcount]].values)**2)*df_master['I_dn'].values

df_master['pph_norm']=pph_n.values
df_master['I_dn_norm']=df_master['pph_norm'].values/df_density[densitycols[colcount]].values
df_master['deltaI_dn_norm']=np.sqrt((deltaI/df_master['pph_norm'].values)**2 + (df_density[densitycols[colcount+1]].values/df_density[densitycols[colcount]].values)**2)*df_master['I_dn_norm'].values
    
#dump everything into an output file
df_master.to_csv('alanineresults.csv', index = False)


#PLOT SOME STUFF
##################################################
dosegroups=df_master.groupby('Dose')
averages=[]
doses=[]
stdevs=[]
averages_norm=[]
stdevs_norm=[]
df_analysis=pd.DataFrame()

for d in np.unique(df_master['Dose']):
    doses.append(d)
    av = (dosegroups.get_group(d)['I_dn'].mean())
    av_norm = (dosegroups.get_group(d)['I_dn_norm'].mean())
    std = (dosegroups.get_group(d)['I_dn'].std())
    std_norm = (dosegroups.get_group(d)['I_dn_norm'].std())
    averages.append(av)
    stdevs.append(std)
    averages_norm.append(av_norm)
    stdevs_norm.append(std_norm)

df_analysis['Dose']=doses
df_analysis['Average']=averages
df_analysis['STDEV']=stdevs  
df_analysis['Average_norm']=averages_norm
df_analysis['STDEV_norm']=stdevs_norm

averages=np.array(averages)
doses=np.array(doses)
stdevs=np.array(stdevs)

samplegroups = df_master.groupby('Number')
sampleaverages=[]
sampleaverages_norm=[]
samplenumber=[]
sampledose=[]
samplestdevs=[]
samplestdevs_norm=[]
df_sampleanalysis=pd.DataFrame()

for s in np.unique(df_master['Number']):
    samplenumber.append(s)
    sdose=df_master.loc[df_master['Number']==s,'Dose'].tolist()[0]
    av = (samplegroups.get_group(s)['I_dn'].mean())
    av_norm=(samplegroups.get_group(s)['I_dn_norm'].mean())
    std = (samplegroups.get_group(s)['I_dn'].std())
    std_norm = (samplegroups.get_group(s)['I_dn_norm'].std())
    sampleaverages.append(av)
    samplestdevs.append(std)
    sampleaverages_norm.append(av_norm)
    samplestdevs_norm.append(std_norm)
    sampledose.append(sdose)
df_sampleanalysis['Number']=samplenumber
df_sampleanalysis['Dose']=sampledose
df_sampleanalysis['Average']=sampleaverages
df_sampleanalysis['STDEV']=samplestdevs
df_sampleanalysis['Average Normalized']=sampleaverages_norm
df_sampleanalysis['STDEV Normalized']=samplestdevs_norm

writer=pd.ExcelWriter('alanineAnalysis.xlsx', engine='xlsxwriter')
df_analysis.to_excel(writer, sheet_name='By Dose', index=False)
df_sampleanalysis.to_excel(writer, sheet_name='By Sample', index=False)
writer.save()

def linear(x,m,t):
    return x*m+t

popt, pcov = curve_fit(linear,doses,averages)
plt.figure(figsize=(7.25,5))
plt.errorbar(doses,averages,yerr=stdevs,ls=' ',marker='o',ms=3,capsize=5.,color='black')
plt.xlabel('Dose [Gy]')
plt.ylabel('pp-height difference [a.u.]')
plt.plot(doses, linear(doses,popt[0],popt[1]),linestyle='dotted',color='black')
plt.title('Peak-to-peak height average, all samples')
plt.ticklabel_format(style='sci', axis='y', scilimits=(6,6))
plt.grid(b=True, which='major')
#plt.xlim(0.18,1.2)
#plt.ylim(-1,1000000)
boxprops = {"boxstyle": "square, pad=0.8", "facecolor": "xkcd:light gray", "alpha": 0.9}
plt.text(0.7, 4.0e6, "$y=("+'{:0.01f}'.format(popt[0]/100000) +'x+''{:0.01f}'.format(popt[1]/100000)+') \\times 10^5$',color='black', bbox=boxprops, verticalalignment='center')
plt.savefig('alanine.png',dpi=600)
plt.show()

popt, pcov = curve_fit(linear,doses,averages_norm)
plt.figure(figsize=(7.25,5))
plt.errorbar(doses,averages_norm,yerr=stdevs_norm,ls=' ',marker='o',ms=3,capsize=5.,color='black')
plt.xlabel('Dose [Gy]')
plt.ylabel('pp-height difference [a.u.]')
plt.plot(doses, linear(doses,popt[0],popt[1]),linestyle='dotted',color='black')
plt.title('Peak-to-peak height average, all samples')
plt.ticklabel_format(style='sci', axis='y', scilimits=(6,6))
plt.grid(b=True, which='major')
#plt.xlim(0.18,1.2)
#plt.ylim(-1,1000000)
boxprops = {"boxstyle": "square, pad=0.8", "facecolor": "xkcd:light gray", "alpha": 0.9}
plt.text(0.7, 2.75e6, "$y=("+'{:0.01f}'.format(popt[0]/100000) +'x+''{:0.01f}'.format(popt[1]/100000)+') \\times 10^5$',color='black', bbox=boxprops, verticalalignment='center')
plt.savefig('alanine_norm.png',dpi=600)
plt.show()

#write everything in results file
writer=pd.ExcelWriter('alanineAnalysis.xlsx', engine='xlsxwriter')
df_master.to_excel(writer, sheet_name='All Data', index=False)
df_analysis.to_excel(writer, sheet_name='By Dose', index=False)
df_sampleanalysis.to_excel(writer, sheet_name='By Sample', index=False)
writer.save()
