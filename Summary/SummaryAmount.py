#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This code generates the plot in Figure 13. In order to run the code, create 
a folder name it as "summarydata" and import/copy the csv files from 
"*dataSummary" folders obtained from "summaryPeriodAmt*" where * denotes 
'TGFInh', 'IFN', 'Con', and 'Alt'. 

@authors: Aurelio A. de los Reyes V and Yangjin Kim
"""

from __future__ import division     #floating point division
import matplotlib.pyplot as plt     #plotting library
import os as os
import pandas as pd
import matplotlib.patches as mpatches

#==============================================================================

def read_dataTGF(path,filename):
    workingdir = os.getcwd()

    os.chdir(path)
    dfTGF = pd.read_csv(filename)
    dfTGF.columns = ['C0','I0','TreatmentPeriod0','TreatmentPeriod1','TGFAmountTreatment',\
                  'RelapsePeriod0','RelapsePeriod1','TGFAmountRelapse',\
                  'RelapsePeriod','frequency','totTGFAmount'] 
    os.chdir(workingdir)

    return dfTGF

#==============================================================================

def read_dataIFN(path,filename):
    workingdir = os.getcwd()

    os.chdir(path)
    dfIFN = pd.read_csv(filename)
    dfIFN.columns = ['C0','I0','TreatmentPeriod0','TreatmentPeriod1','IFNAmountTreatment',\
                  'CritPeriod0','CritPeriod1','IFNAmountCrit',\
                  'MainPeriod0','MainPeriod1','IFNAmountMain',\
                  'totIFNAmount','RelapsePeriod','frequency'] 
    os.chdir(workingdir)

    return dfIFN

#==============================================================================

def read_dataCon(path,filename):
    workingdir = os.getcwd()

    os.chdir(path)
    dfCon = pd.read_csv(filename)
    dfCon.columns = ['C0','I0','TreatmentPeriod0','TreatmentPeriod1', 'TGFAmountTreatment',\
                  'IFNAmountTreatment', 'CancerFreePeriod0',\
                  'CancerFreePeriod1','RelapsePeriod0','TGFAmountRelapse',\
                  'IFNAmountRelapse','totTGFAmount','totIFNAmount',\
                  'RelapsePeriod','frequency'] 
    os.chdir(workingdir)

    return dfCon

#==============================================================================

def read_dataAlt(path,filename):
    workingdir = os.getcwd()

    os.chdir(path)
    dfAlt = pd.read_csv(filename)
    dfAlt.columns = ['C0','I0','TreatmentPeriod0','TreatmentPeriod1', 'TGFAmountTreatment',\
                  'IFNAmountTreatment', 'CancerFreePeriod0',\
                  'CancerFreePeriod1','RelapsePeriod0','RelapsePeriod1',\
                  'TGFAmountRelapse','IFNAmountRelapse',\
                  'totTGFAmount','totIFNAmount' ,'RelapsePeriod',\
                  'frequency'] 
    os.chdir(workingdir)

    return dfAlt

#==============================================================================

#accessing the folder with the csv files
path = os.getcwd()+'/'+'summarydata'

filename = ['TGFInhsummary.csv','IFNsummary.csv','Consummary.csv','Altsummary.csv']

dfTGF = read_dataTGF(path,filename[0])
dfIFN = read_dataIFN(path,filename[1])
dfCon = read_dataCon(path,filename[2])
dfAlt = read_dataAlt(path,filename[3])

dfnew = pd.DataFrame({'init_cond': [r'$Q_1$', r'$Q_2$',r'$Q_3$',r'$Q_4$']})
dfnew['totTGFAmount_TGFonly'] = dfTGF['totTGFAmount']/dfTGF['totTGFAmount'][3]
dfnew['totIFNAmount_TGFonly'] = 0
dfnew['totTGFAmount_IFNonly'] = 0
dfnew['totIFNAmount_IFNonly'] = dfIFN['totIFNAmount']/dfIFN['totIFNAmount'][3]
dfnew['totTGFAmount_con'] = dfCon['totTGFAmount']/dfTGF['totTGFAmount'][3]
dfnew['totIFNAmount_con'] = dfCon['totIFNAmount']/dfIFN['totIFNAmount'][3]
dfnew['totTGFAmount_alt'] = dfAlt['totTGFAmount']/dfTGF['totTGFAmount'][3]
dfnew['totIFNAmount_alt'] = dfAlt['totIFNAmount']/dfIFN['totIFNAmount'][3]



plt.close('all')

# Setting the positions and width for the bars
pos = list(range(len(dfnew['init_cond']))) 
width = 0.1 
    
# Plotting the bars
fig, ax = plt.subplots(figsize=(15,10))

TGF_patch = mpatches.Patch(color='lightcoral', label=r'TGF-$\beta$ Inh')
IFN_patch = mpatches.Patch(color='mediumseagreen', label=r'IFN-$\beta$')

# Create a bar with pre_score data,
# in position pos,
plt.bar(pos, dfnew['totTGFAmount_TGFonly'], width, alpha=0.5,
        label=r'TGF-$\beta$ Inh') 
plt.bar([p + width for p in pos], dfnew['totTGFAmount_con'], width, alpha=0.5,
        label='concomitant')
plt.bar([p + width*2 for p in pos], dfnew['totTGFAmount_alt'], width, alpha=0.5,
        label='alternating')

plt.bar([p + width*3.75 for p in pos], dfnew['totIFNAmount_IFNonly'], width, alpha=0.5,
        label=r'IFN-$\beta$')
plt.bar([p + width*4.75 for p in pos], dfnew['totIFNAmount_con'], width, alpha=0.5,
        label='concomitant')
plt.bar([p + width*5.75 for p in pos], dfnew['totIFNAmount_alt'], width, alpha=0.5,
        label='alternating')

plt.legend(bbox_to_anchor=(0, 1, 1., .102), loc=1, ncol=6,borderaxespad=0,\
           prop={'size': 14},frameon=False)

plt.ylim([0, 1.3] )

ax.set_xticks([p + 3 * width for p in pos])
ax.set_xticklabels(dfnew['init_cond'],fontsize=22)
ax.set_xlabel('initial condition',fontsize=22)
ax.set_ylabel('relative total amount',fontsize=22)
ax.set_yticklabels([0,0.2,0.4,0.6,0.8,1.0],fontsize=20)


ax.annotate('', xy=(-0.05, 0.75), xycoords='data',
             xytext=(0.25, 0.75), textcoords='data',
             arrowprops=dict(arrowstyle='|-|,widthA=0.2,widthB=0.2',
                            color='C3'))
ax.text(.045,0.95,r'TGF-$\beta$ Inh',color='C3',fontsize=16,rotation=45)

ax.annotate('', xy=(0.32, 0.88), xycoords='data',
             xytext=(0.62, 0.88), textcoords='data',
             arrowprops=dict(arrowstyle='|-|,widthA=0.2,widthB=0.2', 
                            color='C2'))
ax.text(.4,0.99,r'IFN-$\beta$',color='C2',fontsize=16,rotation=45)

ax.annotate('', xy=(0.94, 0.65), xycoords='data',
             xytext=(1.26, 0.65), textcoords='data',
             arrowprops=dict(arrowstyle='|-|,widthA=0.2,widthB=0.2',
                            color='C3'))
ax.text(1.05,0.85,r'TGF-$\beta$ Inh',color='C3',fontsize=16,rotation=45)

ax.annotate('', xy=(1.32, 0.6), xycoords='data',
             xytext=(1.62, 0.6), textcoords='data',
             arrowprops=dict(arrowstyle='|-|,widthA=0.2,widthB=0.2', 
                            color='C2'))
ax.text(1.4,0.71,r'IFN-$\beta$',color='C2',fontsize=16,rotation=45)

ax.annotate('', xy=(1.94, 0.76), xycoords='data',
             xytext=(2.26, 0.76), textcoords='data',
             arrowprops=dict(arrowstyle='|-|,widthA=0.2,widthB=0.2',
                            color='C3'))
ax.text(2.05,0.96,r'TGF-$\beta$ Inh',color='C3',fontsize=16,rotation=45)

ax.annotate('', xy=(2.32, 0.87), xycoords='data',
             xytext=(2.62, 0.87), textcoords='data',
             arrowprops=dict(arrowstyle='|-|,widthA=0.2,widthB=0.2', 
                            color='C2'))
ax.text(2.4,0.98,r'IFN-$\beta$',color='C2',fontsize=16,rotation=45)

ax.annotate('', xy=(2.94, 1.05), xycoords='data',
             xytext=(3.26, 1.05), textcoords='data',
             arrowprops=dict(arrowstyle='|-|,widthA=0.2,widthB=0.2',
                            color='C3'))
ax.text(3.05,1.25,r'TGF-$\beta$ Inh',color='C3',fontsize=16,rotation=45)

ax.annotate('', xy=(3.32, 1.05), xycoords='data',
             xytext=(3.62, 1.05), textcoords='data',
             arrowprops=dict(arrowstyle='|-|,widthA=0.2,widthB=0.2', 
                            color='C2'))
ax.text(3.4,1.16,r'IFN-$\beta$',color='C2',fontsize=16,rotation=45)

ax.text(0.28,1.46,r'TGF-$\beta$ Inh amount',color='C3',fontsize=18)
ax.text(2.4,1.46,r'IFN-$\beta$ amount',color='C2',fontsize=18)

plt.tight_layout()
plt.show()

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 13.tif',\
                 dpi=300, bbox_inches='tight', pad_inches=0)

os.chdir(workingdir)



