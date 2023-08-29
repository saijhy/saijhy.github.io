#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This code generates the optimal TGF-beta inhibitor  control only. It creates 
a "TGFInhdataSummary" containing a CSV file summarizing the information
obtained from PeriodAmtTGFInh.py after running the following cases 
    (1) tumorigenic region (Q1) with [C0,I0]=[4,0.3]
    (2) safe region (Q2) with [C0,I0]=[3,2]
    (3) anti-tumorigenic region (Q3) with [C0,I0]=[1,2]
    (4) safe region (Q4)with [C0,I0]=[1,0.5]
In addition, a folder "plots" is generated containing the plots obtained in
Figures 3(A), 3(B), 4(A) and 4(B).

@authors: Aurelio A. de los Reyes V and Yangjin Kim
"""

#==============================================================================

from __future__ import division     #floating point division
import numpy as np                  #library supporting arrays and matrices 
import matplotlib.pyplot as plt     #plotting library
import os as os
import pandas as pd
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#==============================================================================

def read_data(path,filename):
    workingdir = os.getcwd()

    os.chdir(path)
    df = pd.read_csv(filename)
    df.columns = ['TreatmentPeriod0','TreatmentPeriod1','TGFAmountTreatment',\
                  'RelapsePeriod0','RelapsePeriod1','TGFAmountRelapse',\
                  'RelapsePeriod','frequency','totTGFAmount'] 
    os.chdir(workingdir)

    return df

#==============================================================================

def read_profile(path,filename):
    workingdir = os.getcwd()

    os.chdir(path)
    dfprof = pd.read_csv(filename)
    dfprof.columns = ['time','Lsol','Gsol','Csol','Isol','Tsol','control'] 
    os.chdir(workingdir)

    return dfprof

#==============================================================================

def read_injection(path,filename):
    workingdir = os.getcwd()

    os.chdir(path)
    dfinj = pd.read_csv(filename)
    dfinj.columns = ['TGFinj_t1','TGFinj_t2'] 
    os.chdir(workingdir)

    return dfinj

#==============================================================================


def func(pct, amount):
    absamount = pct/100.*np.sum(amount)
    return "{:.1f}%\n({:0.2f} u)".format(pct, absamount)

#==============================================================================
#==============================================================================

#accessing the folder with the csv files
path = os.getcwd()+'/'+'TGFInhdatafile'

filename = ['TGFInhsummary_[C0,I0]=[3,2].csv',\
            'TGFInhsummary_[C0,I0]=[1,2].csv',\
            'TGFInhsummary_[C0,I0]=[1,0.5].csv',\
            'TGFInhsummary_[C0,I0]=[4,0.3].csv']

C0 = [3,1,1,4]
I0 = [2,2,0.5,0.3]

TreatmentPeriod0 = []
TreatmentPeriod1 = [] 
TGFAmountTreatment = []
RelapsePeriod0 = [] 
RelapsePeriod1 = [] 
TGFAmountRelapse = []
RelapsePeriod = []
frequency = []
totTGFAmount = []

for j in range(len(filename)):
    df = read_data(path,filename[j])
    TreatmentPeriod0 = np.append(TreatmentPeriod0,df['TreatmentPeriod0'])
    TreatmentPeriod1 = np.append(TreatmentPeriod1,df['TreatmentPeriod1'])
    TGFAmountTreatment = np.append(TGFAmountTreatment,df['TGFAmountTreatment'])
    RelapsePeriod0 = np.append(RelapsePeriod0,df['RelapsePeriod0'])
    RelapsePeriod1 = np.append(RelapsePeriod1,df['RelapsePeriod1'])
    TGFAmountRelapse = np.append(TGFAmountRelapse,df['TGFAmountRelapse'])
    RelapsePeriod = np.append(RelapsePeriod,df['RelapsePeriod'])
    frequency = np.append(frequency,df['frequency'])
    totTGFAmount = np.append(totTGFAmount,df['totTGFAmount'])

dfnew = pd.DataFrame({'C0':C0})
dfnew['I0'] = I0
dfnew['TreatmentPeriod0'] = TreatmentPeriod0
dfnew['TreatmentPeriod1'] = TreatmentPeriod1
dfnew['TGFAmountTreatment'] = TGFAmountTreatment
dfnew['RelapsePeriod0'] = RelapsePeriod0 
dfnew['RelapsePeriod1'] = RelapsePeriod1
dfnew['TGFAmountRelapse'] = TGFAmountRelapse
dfnew['RelapsePeriod'] = RelapsePeriod
dfnew['frequency'] = frequency
dfnew['totTGFAmount'] = totTGFAmount

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'TGFInhdataSummary'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

dfnew.to_csv('TGFInhsummary.csv',index=False)
 
os.chdir(workingdir)

#==============================================================================
#==============================================================================


plt.close('all')

plt.figure(figsize=(10,7))
mpl.rcParams['font.size'] = 24.0
myexplodes = (0.01, 0.01, 0.01, 0.01)
mycolors = ['lightcoral']*4
mylabels = [r'$Q_1$',r'$Q_2$',r'$Q_3$',r'$Q_4$']

bbox1 = dict(boxstyle="round,pad=0.3", fc="C1", ec="k", lw=2)
bbox2 = dict(boxstyle="round,pad=0.3", fc="C0", ec="k", lw=2)
bbox3 = dict(boxstyle="round,pad=0.3", fc="C2", ec="k", lw=2)
bbox4 = dict(boxstyle="round,pad=0.3", fc="C3", ec="k", lw=2)

arrowprops=dict(arrowstyle="-",connectionstyle="angle,angleA=0,angleB=90")
kw1 = dict(xycoords='data',textcoords='data',
          arrowprops=arrowprops, bbox=bbox1, zorder=0)
kw2 = dict(xycoords='data',textcoords='data',
          arrowprops=arrowprops, bbox=bbox2, zorder=0)
kw3 = dict(xycoords='data',textcoords='data',
          arrowprops=arrowprops, bbox=bbox3, zorder=0)
kw4 = dict(xycoords='data',textcoords='data',
          arrowprops=arrowprops, bbox=bbox4, zorder=0)

plt.pie(dfnew['totTGFAmount'],explode=myexplodes,\
        colors=mycolors,autopct=lambda pct: func(pct, totTGFAmount))
plt.axis('equal')
plt.gca().annotate(r"$Q_1$", xy=(0.5, 0), xytext=(1.2,  0.4), **kw1)
plt.gca().annotate(r"$Q_2$", xy=(-0.1, 0.9), xytext=(-1.5,  0.9), **kw2)
plt.gca().annotate(r"$Q_3$", xy=(-0.5, 0), xytext=( -1.5, -0.1), **kw3)
plt.gca().annotate(r"$Q_4$", xy=(0, -0.4), xytext=( 1.2, -0.9), **kw4)
plt.title(r'Total amount of TGF-$\beta$ inhibitor used for entire duration',\
          fontsize=22)

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 4A.tif',\
                 dpi=300, bbox_inches='tight', pad_inches=0)

os.chdir(workingdir)

#==============================================================================
#==============================================================================

#accessing the folder with the csv files
pathprof = os.getcwd()+'/'+'TGFInhdatafile'

fileprof = ['TGFInhdata_[C0,I0]=[3,2].csv',\
            'TGFInhdata_[C0,I0]=[1,2].csv',\
            'TGFInhdata_[C0,I0]=[1,0.5].csv',\
            'TGFInhdata_[C0,I0]=[4,0.3].csv']

dfprof1 = read_profile(pathprof,fileprof[0])
time1 = dfprof1['time']
Csol1 = dfprof1['Csol']
Isol1 = dfprof1['Isol']
Lsol1 = dfprof1['Lsol']
Gsol1 = dfprof1['Gsol']
uLcontrol1 = dfprof1['control']

dfprof2 = read_profile(pathprof,fileprof[1])
time2 = dfprof2['time']
Csol2 = dfprof2['Csol']
Isol2 = dfprof2['Isol']
Lsol2 = dfprof2['Lsol']
Gsol2 = dfprof2['Gsol']
uLcontrol2 = dfprof2['control']

dfprof3 = read_profile(pathprof,fileprof[2])
time3 = dfprof3['time']
Csol3 = dfprof3['Csol']
Isol3 = dfprof3['Isol']
Lsol3 = dfprof3['Lsol']
Gsol3 = dfprof3['Gsol']
uLcontrol3 = dfprof3['control']

dfprof4 = read_profile(pathprof,fileprof[3])
time4 = dfprof4['time']
Csol4 = dfprof4['Csol']
Isol4 = dfprof4['Isol']
Lsol4 = dfprof4['Lsol']
Gsol4 = dfprof4['Gsol']
uLcontrol4 = dfprof4['control']

#==============================================================================
#==============================================================================

#accessing the folder with the csv injection  files
pathprof = os.getcwd()+'/'+'TGFInhdatafile'

fileprof = ['TGFInhinjsummary_[C0,I0]=[3,2].csv',\
            'TGFInhinjsummary_[C0,I0]=[1,2].csv',\
            'TGFInhinjsummary_[C0,I0]=[1,0.5].csv',\
            'TGFInhinjsummary_[C0,I0]=[4,0.3].csv']

dfinj1 = read_injection(pathprof,fileprof[0])
TGFinj_t11 = dfinj1['TGFinj_t1']
TGFinj_t21 = dfinj1['TGFinj_t2']

dfinj2 = read_injection(pathprof,fileprof[1])
TGFinj_t12 = dfinj2['TGFinj_t1']
TGFinj_t22 = dfinj2['TGFinj_t2']

dfinj3 = read_injection(pathprof,fileprof[2])
TGFinj_t13 = dfinj3['TGFinj_t1']
TGFinj_t23 = dfinj3['TGFinj_t2']

dfinj4 = read_injection(pathprof,fileprof[3])
TGFinj_t14 = dfinj4['TGFinj_t1']
TGFinj_t24 = dfinj4['TGFinj_t2']

Cth = 1.81
Ith = 1.29

x=np.linspace(0,4.25,1000)
plt.figure(figsize=(10,7))
plt.plot(Csol1,Isol1,'C1',linestyle='-.',linewidth=3,label='$(C_0,I_0)=(3,2)$')
plt.plot(Csol2,Isol2,'C0',linewidth=3,label='$(C_0,I_0)=(1,2)$')
plt.plot(Csol3,Isol3,'C2',linestyle='--',linewidth=3,label='$(C_0,I_0)=(1,0.5)$')
plt.plot(Csol4,Isol4,'C3',linestyle=':',linewidth=3,label='$(C_0,I_0)=(4,0.3)$')
plt.scatter(Csol1[0],Isol1[0],s=130,c='C1',marker='*')
plt.scatter(Csol2[0],Isol2[0],s=130,c='C0',marker='*')
plt.scatter(Csol3[0],Isol3[0],s=130,c='C2',marker='*')
plt.scatter(Csol4[0],Isol4[0],s=130,c='C3',marker='*')
plt.axvline(Cth,color='k',linestyle='--')
plt.axhline(Ith,color='k',linestyle='--')
ticks=np.arange(0,5,1)
listticks=[int(ticks[i]) for i in range(len(ticks))]
plt.yticks(ticks,fontsize=24)
plt.xticks(ticks,fontsize=24)
plt.xlabel(r'N2 TANs $(C)$',fontsize=26)
plt.ylabel(r'N1 TANs $(I)$',fontsize=26)
plt.xlim((0,4.25))
plt.ylim((0,3.25))
plt.fill_between(x, 0, Ith, where=x > Cth, facecolor='C3', alpha=0.25)
plt.fill_between(x, Ith, 4.25, where=x < Cth, facecolor='C0', alpha=0.25)
plt.fill_between(x, Ith, 4.25, where=x > Cth, facecolor='C1', alpha=0.25)
plt.fill_between(x, 0, Ith, where=x < Cth, facecolor='C2', alpha=0.25)
plt.text(0.2,1.5,'anti-tumorigenic\n       region',fontsize=24,rotation=90)
plt.text(2.75,0.75,'tumorigenic\n     region',fontsize=24)
plt.text(C0[0]+0.1,I0[0]-0.05,r'$Q_1$',fontsize=20, bbox=bbox1)
plt.text(C0[1]+0.09,I0[1]-0.15,r'$Q_2$',fontsize=20, bbox=bbox2)
plt.text(C0[2]-0.25,I0[2]-0.15,r'$Q_3$',fontsize=20, bbox=bbox3)
plt.text(C0[3],I0[3]+0.17,r'$Q_4$',fontsize=20, bbox=bbox4)
plt.tight_layout()

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()

plt.savefig('Figure 4B.tif',\
                 dpi=300, bbox_inches='tight', pad_inches=0)

os.chdir(workingdir)


#==============================================================================
#==============================================================================

Cth = 1.81
Ith = 1.29
Tf = 30

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True,\
                       figsize=(10, 7))
ax1=plt.subplot(211)
ax1.plot(time1,uLcontrol1,'-.',color='gray')
ax1.fill_between(time1,uLcontrol1,0,\
        where=uLcontrol1>=0,color='lightcoral',interpolate=True)
ax1.set_ylim([0,1.3])
ax1.yaxis.set_tick_params(labelsize=20)
ax1.set_ylabel(r'TGF-$\beta$ Inh',color='C3',fontsize=20)
ax1.tick_params(axis='y', labelcolor='C3')
ax1.text(-4,1.2,r'$Q_1$',fontsize=24, bbox=bbox1)

ax1b = ax1.twinx()
ax1b.plot(time1,Csol1,color='C0',linewidth=3,label=r'$C$')
ax1b.plot(time1,Isol1,color='m',linewidth=3,label=r'$I$')
ax1b.plot(time1,Lsol1,color='C3',linewidth=3,label=r'$L$')
ax1b.plot(time1,Gsol1,color='C5',linewidth=3,label=r'$G$')
ax1b.axhline(Cth,color='C0',linestyle=':',linewidth=3,label=r'$C_{\rm th}$')
ax1b.axhline(Ith,color='m',linestyle=':',linewidth=3,label=r'$I_{\rm th}$')
ax1b.axvline(dfnew['TreatmentPeriod1'][0],color='k',linestyle='--',linewidth=2)
ax1b.axvline(dfnew['RelapsePeriod0'][0],color='k',linestyle='--',linewidth=2)
ax1b.annotate('', xy=(dfnew['TreatmentPeriod0'][0], 3.25), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][0], 3.25), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax1b.text(0,3.4,r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax1b.annotate('', xy=(dfnew['TreatmentPeriod1'][0], 3.25), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][0], 3.25), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax1b.text((dfnew['TreatmentPeriod1'][0]+dfnew['RelapsePeriod0'][0])/2.1,\
         3.4,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax1b.annotate('', xy=(dfnew['RelapsePeriod0'][0], 3.25), xycoords='data',
             xytext=(Tf, 3.25), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax1b.text((dfnew['RelapsePeriod0'][0]+Tf)/2.05,\
         3.4,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax1b.yaxis.set_tick_params(labelsize=20)
ax1b.set_xlim([0,Tf])
ax1b.set_ylim([0,3.9])
ax1b.set_ylabel('concentration',fontsize=20,rotation=270, labelpad=25)
ax1b.legend(bbox_to_anchor=(0, 1, 1, .225),loc=1, ncol=10, borderaxespad=0,\
           prop={'size': 16})


ax2=plt.subplot(212,sharex=ax1)
ax2.plot(time2,uLcontrol2,'-.',color='gray')
ax2.fill_between(time2,uLcontrol2,0,\
        where=uLcontrol2>=0,color='lightcoral',interpolate=True)
ax2.set_ylim([0,1.3])
ax2.yaxis.set_tick_params(labelsize=20)
ax2.set_ylabel(r'TGF-$\beta$ Inh',color='C3',fontsize=20)
ax2.tick_params(axis='y', labelcolor='C3')
ax2.text(-4,1.2,r'$Q_2$',fontsize=24, bbox=bbox2)

ax2b = ax2.twinx()
ax2b.plot(time2,Csol2,color='C0',linewidth=3,label=r'$C$')
ax2b.plot(time2,Isol2,color='m',linewidth=3,label=r'$I$')
ax2b.plot(time2,Lsol2,color='C3',linewidth=3,label=r'$L$')
ax2b.plot(time2,Gsol2,color='C5',linewidth=3,label=r'$G$')
ax2b.axhline(Cth,color='C0',linestyle=':',linewidth=3,label=r'$C_{\rm th}$')
ax2b.axhline(Ith,color='m',linestyle=':',linewidth=3,label=r'$I_{\rm th}$')
ax2b.axvline(dfnew['TreatmentPeriod1'][1],color='k',linestyle='--',linewidth=2)
ax2b.axvline(dfnew['RelapsePeriod0'][1],color='k',linestyle='--',linewidth=2)
ax2b.annotate('', xy=(0, 3.25), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][1], 3.25), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax2b.text(dfnew['TreatmentPeriod1'][1]/2.2,3.4,\
         r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax2b.annotate('', xy=(dfnew['TreatmentPeriod1'][1], 3.25), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][1], 3.25), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax2b.text((dfnew['TreatmentPeriod1'][1]+dfnew['RelapsePeriod0'][1])/2.1,\
         3.4,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax2b.annotate('', xy=(dfnew['RelapsePeriod0'][1], 3.25), xycoords='data',
             xytext=(Tf, 3.25), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax2b.text((dfnew['RelapsePeriod0'][1]+Tf)/2.05,\
         3.4,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax2b.yaxis.set_tick_params(labelsize=20)
ax2b.set_xlim([0,Tf])
ax2b.set_ylim([0,3.9])
ax2b.set_ylabel('concentration',fontsize=20,rotation=270, labelpad=25)
ax2.xaxis.set_tick_params(labelsize=20)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.subplots_adjust(bottom=0.13)

#fig.text(0.04, 0.5, 'concentration', va='center', rotation='vertical',\
#         fontsize=22)

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 3A-1.tif',dpi=300, pad_inches=0)

os.chdir(workingdir)

#==============================================================================
#==============================================================================

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True,\
                       figsize=(10, 7))
ax1=plt.subplot(211)
ax1.plot(time1,uLcontrol1,'-.',color='gray')
ax1.fill_between(time1,uLcontrol1,0,\
        where=uLcontrol1>=0,color='lightcoral',interpolate=True)
ax1.set_ylim([0,1.3])
ax1.yaxis.set_tick_params(labelsize=20)
ax1.set_ylabel(r'TGF-$\beta$ Inh',color='C3',fontsize=20)
ax1.tick_params(axis='y', labelcolor='C3')
ax1.text(-4,1.2,r'$Q_1$',fontsize=24, bbox=bbox1)

ax1b = ax1.twinx()
ax1b.plot(time1,Csol1/Isol1,color='k',linewidth=3,linestyle='-.',label='N2-N1 ratio')
ax1b.axvline(dfnew['TreatmentPeriod1'][0],color='k',linestyle='--',linewidth=2)
ax1b.axvline(dfnew['RelapsePeriod0'][0],color='k',linestyle='--',linewidth=2)
ax1b.annotate('', xy=(dfnew['TreatmentPeriod0'][0], 1.35), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][0], 1.35), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax1b.text(0,1.4,r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax1b.annotate('', xy=(dfnew['TreatmentPeriod1'][0], 1.35), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][0], 1.35), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax1b.text((dfnew['TreatmentPeriod1'][0]+dfnew['RelapsePeriod0'][0])/2.1,\
         1.4,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax1b.annotate('', xy=(dfnew['RelapsePeriod0'][0], 1.35), xycoords='data',
             xytext=(Tf, 1.35), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax1b.text((dfnew['RelapsePeriod0'][0]+Tf)/2.05,\
         1.4,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax1b.yaxis.set_tick_params(labelsize=20)
ax1b.set_xlim([0,Tf])
ax1b.set_ylim([0,1.6])
ax1b.set_ylabel('N2/N1',fontsize=20,rotation=270, labelpad=25)
ax1b.legend(bbox_to_anchor=(0, 1, 1, .225),loc=1, ncol=10, borderaxespad=0,\
           prop={'size': 16})


ax2=plt.subplot(212,sharex=ax1)
ax2.plot(time2,uLcontrol2,'-.',color='gray')
ax2.fill_between(time2,uLcontrol2,0,\
        where=uLcontrol2>=0,color='lightcoral',interpolate=True)
ax2.set_ylim([0,1.3])
ax2.yaxis.set_tick_params(labelsize=20)
ax2.set_ylabel(r'TGF-$\beta$ Inh',color='C3',fontsize=20)
ax2.tick_params(axis='y', labelcolor='C3')
ax2.text(-4,1.2,r'$Q_2$',fontsize=24, bbox=bbox2)

ax2b = ax2.twinx()
ax2b.plot(time2,Csol2/Isol2,color='k',linewidth=3,linestyle='-.')
ax2b.axvline(dfnew['TreatmentPeriod1'][1],color='k',linestyle='--',linewidth=2)
ax2b.axvline(dfnew['RelapsePeriod0'][1],color='k',linestyle='--',linewidth=2)
ax2b.annotate('', xy=(0, 1.35), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][1], 1.35), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax2b.text(dfnew['TreatmentPeriod1'][1]/2.2,1.4,\
         r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax2b.annotate('', xy=(dfnew['TreatmentPeriod1'][1], 1.35), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][1], 1.35), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax2b.text((dfnew['TreatmentPeriod1'][1]+dfnew['RelapsePeriod0'][1])/2.1,\
         1.4,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax2b.annotate('', xy=(dfnew['RelapsePeriod0'][1], 1.35), xycoords='data',
             xytext=(Tf, 1.35), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax2b.text((dfnew['RelapsePeriod0'][1]+Tf)/2.05,\
         1.4,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax2b.yaxis.set_tick_params(labelsize=20)
ax2b.set_xlim([0,Tf])
ax2b.set_ylim([0,1.6])
#ax2b.yaxis.set_ticklabels(np.arange(0,1.6,0.5),[0,0.5,1.0,1.5])
ax2b.set_ylabel('N2/N1',fontsize=20,rotation=270, labelpad=25)
ax2.xaxis.set_tick_params(labelsize=20)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.subplots_adjust(bottom=0.13)

#==============================================================================
#==============================================================================

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 3B-1.tif',dpi=300, pad_inches=0)

os.chdir(workingdir)

#==============================================================================
#==============================================================================

fig2, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True,\
                       figsize=(10, 7))

ax3=plt.subplot(211)
ax3.plot(time3,uLcontrol3,'-.',color='gray')
ax3.fill_between(time3,uLcontrol3,0,\
        where=uLcontrol3>=0,color='lightcoral',interpolate=True)
ax3.set_ylim([0,1.3])
ax3.yaxis.set_tick_params(labelsize=20)
ax3.set_ylabel(r'TGF-$\beta$ Inh',color='C3',fontsize=20)
ax3.tick_params(axis='y', labelcolor='C3')
ax3.text(-4,1.2,r'$Q_3$',fontsize=24, bbox=bbox3)

ax3b = ax3.twinx()
ax3b.plot(time3,Csol3,color='C0',linewidth=3,label=r'$C$')
ax3b.plot(time3,Isol3,color='m',linewidth=3,label=r'$I$')
ax3b.plot(time3,Lsol3,color='C3',linewidth=3,label=r'$L$')
ax3b.plot(time3,Gsol3,color='C5',linewidth=3,label=r'$G$')
ax3b.axhline(Cth,color='C0',linestyle=':',linewidth=3,label=r'$C_{\rm th}$')
ax3b.axhline(Ith,color='m',linestyle=':',linewidth=3,label=r'$I_{\rm th}$')
ax3b.axvline(dfnew['TreatmentPeriod1'][2],color='k',linestyle='--',linewidth=2)
ax3b.axvline(dfnew['RelapsePeriod0'][2],color='k',linestyle='--',linewidth=2)
ax3b.annotate('', xy=(0, 3.25), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][2], 3.25), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax3b.text(0,3.4,r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax3b.annotate('', xy=(dfnew['TreatmentPeriod1'][2], 3.25), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][2], 3.25), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax3b.text((dfnew['TreatmentPeriod1'][2]+dfnew['RelapsePeriod0'][2])/2.1,\
         3.4,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax3b.annotate('', xy=(dfnew['RelapsePeriod0'][2], 3.25), xycoords='data',
             xytext=(Tf, 3.25), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax3b.text((dfnew['RelapsePeriod0'][2]+Tf)/2.05,\
         3.4,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax3b.yaxis.set_tick_params(labelsize=20)
ax3b.set_xlim([0,Tf])
ax3b.set_ylim([0,3.9])
ax3b.set_ylabel('concentration',fontsize=20,rotation=270, labelpad=25)


ax4=plt.subplot(212,sharex=ax3)
ax4.plot(time4,uLcontrol4,'-.',color='gray')
ax4.fill_between(time4,uLcontrol4,0,\
        where=uLcontrol4>=0,color='lightcoral',interpolate=True)
ax4.set_ylim([0,1.3])
ax4.yaxis.set_tick_params(labelsize=20)
ax4.xaxis.set_tick_params(labelsize=20)
ax4.set_ylabel(r'TGF-$\beta$ Inh',color='C3',fontsize=20)
ax4.tick_params(axis='y', labelcolor='C3')
ax4.text(-4,1.2,r'$Q_4$',fontsize=24, bbox=bbox4)
ax4.set_xlabel('Time (days)',fontsize=20)

ax4b = ax4.twinx()
ax4b.plot(time4,Csol4,color='C0',linewidth=3,label=r'$C$')
ax4b.plot(time4,Isol4,color='m',linewidth=3,label=r'$I$')
ax4b.plot(time4,Lsol4,color='C3',linewidth=3,label=r'$L$')
ax4b.plot(time4,Gsol4,color='C5',linewidth=3,label=r'$G$')
ax4b.axhline(Cth,color='C0',linestyle=':',linewidth=3,label=r'$C_{\rm th}$')
ax4b.axhline(Ith,color='m',linestyle=':',linewidth=3,label=r'$I_{\rm th}$')
ax4b.axvline(dfnew['TreatmentPeriod1'][3],color='k',linestyle='--',linewidth=2)
ax4b.axvline(dfnew['RelapsePeriod0'][3],color='k',linestyle='--',linewidth=2)
ax4b.annotate('', xy=(dfnew['TreatmentPeriod0'][3], 4), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][3], 4), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax4b.text((dfnew['TreatmentPeriod0'][3]+dfnew['TreatmentPeriod1'][3])/2.4,\
         4.2,r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax4b.annotate('', xy=(dfnew['TreatmentPeriod1'][3], 4), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][3], 4), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax4b.text((dfnew['TreatmentPeriod1'][3]+dfnew['RelapsePeriod0'][3])/2.1,\
         4.2,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax4b.annotate('', xy=(dfnew['RelapsePeriod0'][3], 4), xycoords='data',
             xytext=(Tf, 4), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax4b.text((dfnew['RelapsePeriod0'][3]+Tf)/2.05,\
         4.2,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax4b.yaxis.set_tick_params(labelsize=20)
ax4b.set_xlim([0,Tf])
ax4b.set_ylim([0,4.75])
ax4b.yaxis.set_tick_params(labelsize=20)
ax4b.set_ylabel('concentration',fontsize=20,rotation=270, labelpad=25)

plt.setp(ax3.get_xticklabels(), visible=False)
plt.subplots_adjust(bottom=0.13)

#==============================================================================
#==============================================================================

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 3A-2.tif',dpi=300, pad_inches=0)

os.chdir(workingdir)

#==============================================================================
#==============================================================================

fig2, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True,\
                       figsize=(10, 7))

ax3=plt.subplot(211)
ax3.plot(time3,uLcontrol3,'-.',color='gray')
ax3.fill_between(time3,uLcontrol3,0,\
        where=uLcontrol3>=0,color='lightcoral',interpolate=True)
ax3.set_ylim([0,1.3])
ax3.yaxis.set_tick_params(labelsize=20)
ax3.set_ylabel(r'TGF-$\beta$ Inh',color='C3',fontsize=20)
ax3.tick_params(axis='y', labelcolor='C3')
ax3.text(-4,1.2,r'$Q_3$',fontsize=24, bbox=bbox3)

ax3b = ax3.twinx()
ax3b.plot(time3,Csol3/Isol3,color='k',linewidth=3,linestyle='-.')
ax3b.axvline(dfnew['TreatmentPeriod1'][2],color='k',linestyle='--',linewidth=2)
ax3b.axvline(dfnew['RelapsePeriod0'][2],color='k',linestyle='--',linewidth=2)
ax3b.annotate('', xy=(0, 1.35), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][2], 1.35), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax3b.text(0,1.4,r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax3b.annotate('', xy=(dfnew['TreatmentPeriod1'][2], 1.35), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][2], 1.35), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax3b.text((dfnew['TreatmentPeriod1'][2]+dfnew['RelapsePeriod0'][2])/2.1,\
         1.4,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax3b.annotate('', xy=(dfnew['RelapsePeriod0'][2], 1.35), xycoords='data',
             xytext=(Tf, 1.35), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax3b.text((dfnew['RelapsePeriod0'][2]+Tf)/2.05,\
         1.4,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax3b.yaxis.set_tick_params(labelsize=20)
ax3b.set_xlim([0,Tf])
ax3b.set_ylim([0,1.6])
ax3b.set_ylabel('N2/N1',fontsize=20,rotation=270, labelpad=25)


ax4=plt.subplot(212,sharex=ax3)
ax4.plot(time4,uLcontrol4,'-.',color='gray')
ax4.fill_between(time4,uLcontrol4,0,\
        where=uLcontrol4>=0,color='lightcoral',interpolate=True)
ax4.set_ylim([0,1.3])
ax4.yaxis.set_tick_params(labelsize=20)
ax4.xaxis.set_tick_params(labelsize=20)
ax4.set_ylabel(r'TGF-$\beta$ Inh',color='C3',fontsize=20)
ax4.tick_params(axis='y', labelcolor='C3')
ax4.text(-4,1.2,r'$Q_4$',fontsize=24, bbox=bbox4)
ax4.set_xlabel('Time (days)',fontsize=20)

ax4b = ax4.twinx()
ax4b.plot(time4,Csol4/Isol4,color='k',linewidth=3,linestyle='-.')
ax4b.axvline(dfnew['TreatmentPeriod1'][3],color='k',linestyle='--',linewidth=2)
ax4b.axvline(dfnew['RelapsePeriod0'][3],color='k',linestyle='--',linewidth=2)
ax4b.annotate('', xy=(dfnew['TreatmentPeriod0'][3], 11), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][3], 11), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax4b.text((dfnew['TreatmentPeriod0'][3]+dfnew['TreatmentPeriod1'][3])/2.4,\
         11.4,r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax4b.annotate('', xy=(dfnew['TreatmentPeriod1'][3], 11), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][3], 11), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax4b.text((dfnew['TreatmentPeriod1'][3]+dfnew['RelapsePeriod0'][3])/2.1,\
         11.4,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax4b.annotate('', xy=(dfnew['RelapsePeriod0'][3], 11), xycoords='data',
             xytext=(Tf, 11), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax4b.text((dfnew['RelapsePeriod0'][3]+Tf)/2.05,\
         11.4,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax4b.yaxis.set_tick_params(labelsize=20)
ax4b.set_xlim([0,Tf])
ax4b.set_ylim([0,13.3])
ax4b.yaxis.set_tick_params(labelsize=20)
ax4b.set_ylabel('N2/N1',fontsize=20,rotation=270, labelpad=25)

axins = zoomed_inset_axes(ax4b, 3.0, loc=7)
axins.plot(time4,Csol4/Isol4,color='k',linewidth=3,linestyle='-.')
axins.axvline(dfnew['RelapsePeriod0'][3],color='gray',linestyle='-.')
axins.axvline(dfnew['RelapsePeriod0'][3]+1,color='gray',linestyle='-.')
axins.axvspan(dfnew['RelapsePeriod0'][3],dfnew['RelapsePeriod0'][3]+1,\
              color='lightcoral')
axins.yaxis.set_tick_params(labelsize=16)
axins.xaxis.set_tick_params(labelsize=16)
axins.set_ylabel('N2/N1',fontsize=12)
x1, x2, y1, y2 = 19.8, 21.4, 0, 1.7 # specify the limits
axins.set_xlim(x1, x2) # apply the x-limits
axins.set_ylim(y1, y2) # apply the y-limits
mark_inset(ax4b, axins, loc1=2, loc2=4, fc="none", ec="0.5",ls='--')


plt.setp(ax3.get_xticklabels(), visible=False)
plt.subplots_adjust(bottom=0.13)

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 3B-2.tif',dpi=300, pad_inches=0)

os.chdir(workingdir)