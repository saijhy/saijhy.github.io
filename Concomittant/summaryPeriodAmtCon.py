#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This code generates the optimal concomitant administration of TGF-beta 
inhibitor and IFN-beta control.  It creates a folder "CondataSummary" 
containing a CSV file summarizing the information obtained from 
PeriodAmtCon.py after running the following cases 
    (1) tumorigenic region (Q1) with [C0,I0]=[4,0.3]
    (2) safe region (Q2) with [C0,I0]=[3,2]
    (3) anti-tumorigenic region (Q3) with [C0,I0]=[1,2]
    (4) safe region (Q4)with [C0,I0]=[1,0.5]
In addition, a folder "plots" is generated containing the plots obtained in
Figures 7(A), 7(B), 8(A), 8(B), 9(A), 9(B) and 9(C).

@authors: Aurelio A. de los Reyes V and Yangjin Kim
"""

from __future__ import division     #floating point division
import numpy as np                  #library supporting arrays and matrices 
import matplotlib.pyplot as plt     #plotting library
import os as os
import pandas as pd
import matplotlib as mpl
import matplotlib.patches as mpatches

#==============================================================================

def read_data(path,filename):
    workingdir = os.getcwd()

    os.chdir(path)
    df = pd.read_csv(filename)
    df.columns = ['TreatmentPeriod0','TreatmentPeriod1', 'TGFAmountTreatment',\
                  'IFNAmountTreatment', 'CancerFreePeriod0',\
                  'CancerFreePeriod1','RelapsePeriod0','TGFAmountRelapse',\
                  'IFNAmountRelapse','totTGFAmount','totIFNAmount',\
                  'RelapsePeriod','frequency'] 
    os.chdir(workingdir)

    return df

#==============================================================================

def read_profile(path,filename):
    workingdir = os.getcwd()

    os.chdir(path)
    dfprof = pd.read_csv(filename)
    dfprof.columns = ['time','Lsol','Gsol','Csol','Ssol','Isol','Tsol',\
                      'uLcontrol','uScontrol'] 
    os.chdir(workingdir)

    return dfprof


#==============================================================================

def read_injection(path,filename):
    workingdir = os.getcwd()

    os.chdir(path)
    dfinj = pd.read_csv(filename)
    dfinj.columns = ['Coninj_t1','Coninj_t2'] 
    os.chdir(workingdir)

    return dfinj

#==============================================================================

def func(pct, amount):
    absamount = pct/100.*np.sum(amount)
    return "{:.1f}%\n({:0.2f} u)".format(pct, absamount)

#==============================================================================

#accessing the folder with the csv files
path = os.getcwd()+'/'+'Condatafile'

filename = ['Consummary_[C0,I0]=[3,2].csv',\
            'Consummary_[C0,I0]=[1,2].csv',\
            'Consummary_[C0,I0]=[1,0.5].csv',\
            'Consummary_[C0,I0]=[4,0.3].csv']

C0 = [3,1,1,4]
I0 = [2,2,0.5,0.3]

TreatmentPeriod0 = [] 
TreatmentPeriod1 = []
TGFAmountTreatment = []
IFNAmountTreatment = []
CancerFreePeriod0 = []
CancerFreePeriod1 = []
RelapsePeriod0 = []
TGFAmountRelapse = []
IFNAmountRelapse = []
totTGFAmount = []
totIFNAmount = []
RelapsePeriod = []
frequency = []

for j in range(len(filename)):
    df = read_data(path,filename[j])
    TreatmentPeriod0 = np.append(TreatmentPeriod0,df['TreatmentPeriod0']) 
    TreatmentPeriod1 = np.append(TreatmentPeriod1,df['TreatmentPeriod1'])
    TGFAmountTreatment = np.append(TGFAmountTreatment,df['TGFAmountTreatment'])
    IFNAmountTreatment = np.append(IFNAmountTreatment,df['IFNAmountTreatment'])
    CancerFreePeriod0 = np.append(CancerFreePeriod0,df['CancerFreePeriod0'])
    CancerFreePeriod1 = np.append(CancerFreePeriod1,df['CancerFreePeriod1'])
    RelapsePeriod0 = np.append(RelapsePeriod0,df['RelapsePeriod0'])
    TGFAmountRelapse = np.append(TGFAmountRelapse,df['TGFAmountRelapse'])
    IFNAmountRelapse = np.append(IFNAmountRelapse,df['IFNAmountRelapse'])
    totTGFAmount = np.append(totTGFAmount,df['totTGFAmount'])
    totIFNAmount = np.append(totIFNAmount,df['totIFNAmount'])
    RelapsePeriod = np.append(RelapsePeriod,df['RelapsePeriod'])
    frequency = np.append(frequency,df['frequency'])
    
dfnew = pd.DataFrame({'C0':C0})
dfnew['I0'] = I0
dfnew['TreatmentPeriod0'] = TreatmentPeriod0
dfnew['TreatmentPeriod1'] = TreatmentPeriod1
dfnew['TGFAmountTreatment'] = TGFAmountTreatment
dfnew['IFNAmountTreatment'] = IFNAmountTreatment
dfnew['CancerFreePeriod0'] = CancerFreePeriod0
dfnew['CancerFreePeriod1'] = CancerFreePeriod1
dfnew['RelapsePeriod0'] = RelapsePeriod0
dfnew['TGFAmountRelapse'] = TGFAmountRelapse
dfnew['IFNAmountRelapse'] = IFNAmountRelapse
dfnew['totTGFAmount']= totTGFAmount
dfnew['totIFNAmount'] = totIFNAmount
dfnew['RelapsePeriod'] = RelapsePeriod
dfnew['frequency'] = frequency

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'CondataSummary'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

dfnew.to_csv('Consummary.csv',index=False)
 
os.chdir(workingdir)

#==============================================================================
#==============================================================================

plt.close('all')

plt.figure(figsize=(10,7))
mpl.rcParams['font.size'] = 24.0
myexplodes = (0.02, 0.02, 0.02, 0.02)
#mycolors = ['C3','C1','C0','C2']
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
plt.gca().annotate(r"$Q_1$", xy=(0.5, 0.5), xytext=(1.2,  0.6), **kw1)
plt.gca().annotate(r"$Q_2$", xy=(-0.1, 0.9), xytext=(-1.5,  0.7), **kw2)
plt.gca().annotate(r"$Q_3$", xy=(-0.5, 0), xytext=( -1.5, -0.6), **kw3)
plt.gca().annotate(r"$Q_4$", xy=(0, -0.4), xytext=( 1.2, -0.7), **kw4)
plt.title(r'Total TGF-$\beta$ inhibitor used during entire duration',\
          fontsize=22)

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 9A.tif',\
                 dpi=300, bbox_inches='tight', pad_inches=0)

os.chdir(workingdir)

#==============================================================================
#==============================================================================

plt.figure(figsize=(10,7))
mpl.rcParams['font.size'] = 24.0
myexplodes = (0.02, 0.02, 0.02, 0.02)
#mycolors = ['C3','C1','C0','C2']
mycolors = ['mediumseagreen']*4
mylabels = [r'$Q_1$',r'$Q_2$',r'$Q_3$',r'$Q_4$']

plt.pie(dfnew['totIFNAmount'],explode=myexplodes,\
        colors=mycolors,autopct=lambda pct: func(pct, totIFNAmount))
plt.axis('equal')
plt.gca().annotate(r"$Q_1$", xy=(0.5, 0.5), xytext=(1.2,  0.6), **kw1)
plt.gca().annotate(r"$Q_2$", xy=(-0.1, 0.9), xytext=(-1.5,  0.9), **kw2)
plt.gca().annotate(r"$Q_3$", xy=(-0.5, 0), xytext=( -1.5, -0.3), **kw3)
plt.gca().annotate(r"$Q_4$", xy=(0, -0.4), xytext=( 1.2, -0.9), **kw4)
plt.title(r'Total IFN-$\beta$ used during entire duration',\
          fontsize=22)

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()

plt.savefig('Figure 9B.tif',\
                 dpi=300, bbox_inches='tight', pad_inches=0)

os.chdir(workingdir)

#==============================================================================
#==============================================================================

#accessing the folder with the csv data files
pathprof = os.getcwd()+'/'+'Condatafile'

fileprof = ['Condata_[C0,I0]=[3,2].csv',\
            'Condata_[C0,I0]=[1,2].csv',\
            'Condata_[C0,I0]=[1,0.5].csv',\
            'Condata_[C0,I0]=[4,0.3].csv']

dfprof1 = read_profile(pathprof,fileprof[0])
time1 = dfprof1['time']
Lsol1 = dfprof1['Lsol']
Gsol1 = dfprof1['Gsol']
Csol1 = dfprof1['Csol']
Ssol1 = dfprof1['Ssol']
Isol1 = dfprof1['Isol']
uLcontrol1 = dfprof1['uLcontrol']
uScontrol1 = dfprof1['uScontrol']

dfprof2 = read_profile(pathprof,fileprof[1])
time2 = dfprof2['time']
Lsol2 = dfprof2['Lsol']
Gsol2 = dfprof2['Gsol']
Csol2 = dfprof2['Csol']
Ssol2 = dfprof2['Ssol']
Isol2 = dfprof2['Isol']
uLcontrol2 = dfprof2['uLcontrol']
uScontrol2 = dfprof2['uScontrol']

dfprof3 = read_profile(pathprof,fileprof[2])
time3 = dfprof3['time']
Lsol3 = dfprof3['Lsol']
Gsol3 = dfprof3['Gsol']
Csol3 = dfprof3['Csol']
Ssol3 = dfprof3['Ssol']
Isol3 = dfprof3['Isol']
uLcontrol3 = dfprof3['uLcontrol']
uScontrol3 = dfprof3['uScontrol']

dfprof4 = read_profile(pathprof,fileprof[3])
time4 = dfprof4['time']
Lsol4 = dfprof4['Lsol']
Gsol4 = dfprof4['Gsol']
Csol4 = dfprof4['Csol']
Ssol4 = dfprof4['Ssol']
Isol4 = dfprof4['Isol']
uLcontrol4 = dfprof4['uLcontrol']
uScontrol4 = dfprof4['uScontrol']

#==============================================================================
#==============================================================================

#accessing the folder with the csv injection  files
pathprof = os.getcwd()+'/'+'Condatafile'

fileprof = ['Coninjsummary_[C0,I0]=[3,2].csv',\
            'Coninjsummary_[C0,I0]=[1,2].csv',\
            'Coninjsummary_[C0,I0]=[1,0.5].csv',\
            'Coninjsummary_[C0,I0]=[4,0.3].csv']

dfinj1 = read_injection(pathprof,fileprof[0])
Coninj_t11 = dfinj1['Coninj_t1']
Coninj_t21 = dfinj1['Coninj_t2']

dfinj2 = read_injection(pathprof,fileprof[1])
Coninj_t12 = dfinj2['Coninj_t1']
Coninj_t22 = dfinj2['Coninj_t2']

dfinj3 = read_injection(pathprof,fileprof[2])
Coninj_t13 = dfinj3['Coninj_t1']
Coninj_t23 = dfinj3['Coninj_t2']

dfinj4 = read_injection(pathprof,fileprof[3])
Coninj_t14 = dfinj4['Coninj_t1']
Coninj_t24 = dfinj4['Coninj_t2']


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
plt.text(0.2,1.4,'anti-tumorigenic\n       region',fontsize=24,rotation=90)
plt.text(2.75,0.75,'tumorigenic\n    region',fontsize=24)
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

plt.savefig('Figure 9C.tif',\
                 dpi=300, bbox_inches='tight', pad_inches=0)

os.chdir(workingdir)

#==============================================================================
Tsim=30
TGF_patch = mpatches.Patch(color='lightcoral', label=r'TGF-$\beta$ Inh')
IFN_patch = mpatches.Patch(color='mediumseagreen', label=r'IFN-$\beta$')

#==============================================================================

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True,figsize=(10, 7))

ax1=plt.subplot(211)
ax1.plot(time1,uScontrol1,'-.',color='gray')
ax1.fill_between(time1,uScontrol1,0,\
        where=uScontrol1>=0,color='mediumseagreen',alpha=0.5,interpolate=True)
ax1.plot(time1,uLcontrol1,'-.',color='gray')
ax1.fill_between(time1,uLcontrol1,0,\
        where=uLcontrol1>=0,color='lightcoral',interpolate=True)
ax1.set_ylim([0,2.5])
ax1.yaxis.set_tick_params(labelsize=20)
ax1.set_ylabel(r'control',fontsize=20)
ax1.text(-3.5,2.2,r'$Q_1$',fontsize=24, bbox=bbox1)
ax1.legend(handles=[TGF_patch,IFN_patch],bbox_to_anchor=(0,0.74, 1, .6),\
           loc=2, ncol=1,borderaxespad=0,prop={'size': 15})

ax1b = ax1.twinx()
ax1b.plot(time1,Csol1,color='C0',linewidth=3,label=r'$C$')
ax1b.plot(time1,Isol1,color='m',linewidth=3,label=r'$I$')
ax1b.plot(time1,Lsol1,color='C3',linewidth=3,label=r'$L$')
ax1b.plot(time1,Gsol1,color='C5',linewidth=3,label=r'$G$')
ax1b.plot(time1,Ssol1,color='C2',linewidth=3,label=r'$S$')
ax1b.axhline(Cth,color='C0',linestyle=':',linewidth=3,label=r'$C_{\rm th}$')
ax1b.axhline(Ith,color='m',linestyle=':',linewidth=3,label=r'$I_{\rm th}$')
ax1b.axvline(dfnew['TreatmentPeriod1'][0],color='k',linestyle='--',linewidth=2)
ax1b.axvline(dfnew['RelapsePeriod0'][0],color='k',linestyle='--',linewidth=2)
ax1b.annotate('', xy=(0, 3), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][0], 3), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax1b.text(0,3.15,r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax1b.annotate('', xy=(dfnew['TreatmentPeriod1'][0], 3), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][0], 3), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax1b.text((dfnew['TreatmentPeriod1'][0]+dfnew['RelapsePeriod0'][0])/2.1,\
         3.15,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax1b.annotate('', xy=(dfnew['RelapsePeriod0'][0], 3), xycoords='data',
             xytext=(Tsim, 3), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax1b.text((dfnew['RelapsePeriod0'][0]+Tsim)/2.05,\
         3.15,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax1b.set_xlabel('Time (days)',fontsize=24)
ax1b.set_ylabel('Concentration',fontsize=24)
ax1b.xaxis.set_tick_params(labelsize=20)
ax1b.yaxis.set_tick_params(labelsize=20)
ax1b.set_xlim([0,Tsim])
ax1b.set_ylim([0,3.55])
ax1b.set_ylabel('concentration',fontsize=20,rotation=270, labelpad=25)
ax1b.legend(bbox_to_anchor=(0, 1.11, 1, .225), loc=1, ncol=4,borderaxespad=0,\
           prop={'size': 15})

ax2=plt.subplot(212,sharex=ax1)
ax2.plot(time2,uScontrol2,'-.',color='gray')
ax2.fill_between(time2,uScontrol2,0,\
        where=uScontrol2>=0,color='mediumseagreen',alpha=0.5,interpolate=True)
ax2.plot(time2,uLcontrol2,'-.',color='gray')
ax2.fill_between(time2,uLcontrol2,0,\
        where=uLcontrol2>=0,color='lightcoral',interpolate=True)
ax2.set_ylim([0,2.5])
ax2.yaxis.set_tick_params(labelsize=20)
ax2.set_ylabel(r'control',fontsize=20)
ax2.text(-3.5,2.2,r'$Q_2$',fontsize=24, bbox=bbox2)

ax2b = ax2.twinx()
ax2b.plot(time2,Csol2,color='C0',linewidth=3,label=r'$C$')
ax2b.plot(time2,Isol2,color='m',linewidth=3,label=r'$I$')
ax2b.plot(time2,Lsol2,color='C3',linewidth=3,label=r'$L$')
ax2b.plot(time2,Gsol2,color='C5',linewidth=3,label=r'$G$')
ax2b.plot(time2,Ssol2,color='C2',linewidth=3,label=r'$S$')
ax2b.axhline(Cth,color='C0',linestyle=':',linewidth=3,label=r'$C_{\rm th}$')
ax2b.axhline(Ith,color='m',linestyle=':',linewidth=3,label=r'$I_{\rm th}$')
ax2b.axvline(dfnew['TreatmentPeriod1'][1],color='k',linestyle='--',linewidth=2)
ax2b.axvline(dfnew['RelapsePeriod0'][1],color='k',linestyle='--',linewidth=2)
ax2b.annotate('', xy=(0, 3), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][1], 3), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax2b.text((dfnew['TreatmentPeriod1'][1])/2.3,\
          3.15,r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax2b.annotate('', xy=(dfnew['TreatmentPeriod1'][1], 3), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][1], 3), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax2b.text((dfnew['TreatmentPeriod1'][1]+dfnew['RelapsePeriod0'][1])/2.1,\
         3.15,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax2b.annotate('', xy=(dfnew['RelapsePeriod0'][1], 3), xycoords='data',
             xytext=(Tsim, 3), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax2b.text((dfnew['RelapsePeriod0'][1]+Tsim)/2.05,\
         3.15,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax2b.set_xlabel('Time (days)',fontsize=24)
ax2b.set_ylabel('Concentration',fontsize=24)
ax2b.xaxis.set_tick_params(labelsize=20)
ax2b.yaxis.set_tick_params(labelsize=20)
ax2b.set_xlim([0,Tsim])
ax2b.set_ylim([0,3.55])
ax2b.set_ylabel('concentration',fontsize=20,rotation=270, labelpad=25)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.subplots_adjust(bottom=0.13)

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 7A-1.tif', dpi=300,pad_inches=0)

os.chdir(workingdir)

#==============================================================================
#==============================================================================

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True,figsize=(10, 7))

ax3=plt.subplot(211)
ax3.plot(time3,uScontrol3,'-.',color='gray')
ax3.fill_between(time3,uScontrol3,0,\
        where=uScontrol3>=0,color='mediumseagreen',alpha=0.5,interpolate=True)
ax3.plot(time3,uLcontrol3,'-.',color='gray')
ax3.fill_between(time3,uLcontrol3,0,\
        where=uLcontrol3>=0,color='lightcoral',interpolate=True)
ax3.set_ylim([0,2.5])
ax3.yaxis.set_tick_params(labelsize=20)
ax3.set_ylabel(r'control',fontsize=20)
ax3.text(-3.5,2.2,r'$Q_3$',fontsize=24, bbox=bbox3)

ax3b = ax3.twinx()
ax3b.plot(time3,Csol3,color='C0',linewidth=3,label=r'$C$')
ax3b.plot(time3,Isol3,color='m',linewidth=3,label=r'$I$')
ax3b.plot(time3,Lsol3,color='C3',linewidth=3,label=r'$L$')
ax3b.plot(time3,Gsol3,color='C5',linewidth=3,label=r'$G$')
ax3b.plot(time3,Ssol3,color='C2',linewidth=3,label=r'$S$')
ax3b.axhline(Cth,color='C0',linestyle=':',linewidth=3,label=r'$C_{\rm th}$')
ax3b.axhline(Ith,color='m',linestyle=':',linewidth=3,label=r'$I_{\rm th}$')
ax3b.axvline(dfnew['TreatmentPeriod1'][2],color='k',linestyle='--',linewidth=2)
ax3b.axvline(dfnew['RelapsePeriod0'][2],color='k',linestyle='--',linewidth=2)
ax3b.annotate('', xy=(0, 3), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][2], 3), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax3b.text(0,3.15,r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax3b.annotate('', xy=(dfnew['TreatmentPeriod1'][2], 3), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][2], 3), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax3b.text((dfnew['TreatmentPeriod1'][2]+dfnew['RelapsePeriod0'][2])/2.1,\
         3.15,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax3b.annotate('', xy=(dfnew['RelapsePeriod0'][2], 3), xycoords='data',
             xytext=(Tsim, 3), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax3b.text((dfnew['RelapsePeriod0'][2]+Tsim)/2.05,\
         3.15,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax3b.set_xlabel('Time (days)',fontsize=24)
ax3b.set_ylabel('Concentration',fontsize=24)
ax3b.xaxis.set_tick_params(labelsize=20)
ax3b.yaxis.set_tick_params(labelsize=20)
ax3b.set_xlim([0,Tsim])
ax3b.set_ylim([0,3.55])
ax3b.set_ylabel('concentration',fontsize=20,rotation=270, labelpad=25)

ax4=plt.subplot(212,sharex=ax2)
ax4.plot(time4,uScontrol4,'-.',color='gray')
ax4.fill_between(time4,uScontrol4,0,\
        where=uScontrol4>=0,color='mediumseagreen',alpha=0.5,interpolate=True)
ax4.plot(time4,uLcontrol4,'-.',color='gray')
ax4.fill_between(time4,uLcontrol4,0,\
        where=uLcontrol4>=0,color='lightcoral',interpolate=True)
ax4.set_ylim([0,7.6])
ax4.xaxis.set_tick_params(labelsize=20)
ax4.yaxis.set_tick_params(labelsize=20)
ax4.set_ylabel(r'control',fontsize=20)
ax4.text(-3.5,7,r'$Q_4$',fontsize=24, bbox=bbox4)
ax4.set_xlabel('Time (days)',fontsize=20)

ax4b = ax4.twinx()
ax4b.plot(time4,Csol4,color='C0',linewidth=3,label=r'$C$')
ax4b.plot(time4,Isol4,color='m',linewidth=3,label=r'$I$')
ax4b.plot(time4,Lsol4,color='C3',linewidth=3,label=r'$L$')
ax4b.plot(time4,Gsol4,color='C5',linewidth=3,label=r'$G$')
ax4b.plot(time4,Ssol4,color='C2',linewidth=3,label=r'$S$')
ax4b.axhline(Cth,color='C0',linestyle=':',linewidth=3,label=r'$C_{\rm th}$')
ax4b.axhline(Ith,color='m',linestyle=':',linewidth=3,label=r'$I_{\rm th}$')
ax4b.axvline(dfnew['TreatmentPeriod1'][3],color='k',linestyle='--',linewidth=2)
ax4b.axvline(dfnew['RelapsePeriod0'][3],color='k',linestyle='--',linewidth=2)
ax4b.annotate('', xy=(0, 4), xycoords='data',
             xytext=(dfnew['TreatmentPeriod1'][3], 4), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='red'))
ax4b.text((dfnew['TreatmentPeriod1'][3])/2.3,\
          4.2,r'$\mathrm{\mathbb{P}}_{\rm t}$',color='red',fontsize=20)
ax4b.annotate('', xy=(dfnew['TreatmentPeriod1'][3], 4), xycoords='data',
             xytext=(dfnew['RelapsePeriod0'][3], 4), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='green'))
ax4b.text((dfnew['TreatmentPeriod1'][3]+dfnew['RelapsePeriod0'][3])/2.1,\
         4.2,r'$\mathrm{\mathbb{P}}_{\rm cf}$',color='green',fontsize=20)
ax4b.annotate('', xy=(dfnew['RelapsePeriod0'][3], 4), xycoords='data',
             xytext=(Tsim, 4), textcoords='data',
             arrowprops=dict(arrowstyle='<->', connectionstyle='arc3', 
                            color='blue'))
ax4b.text((dfnew['RelapsePeriod0'][3]+Tsim)/2.05,\
         4.2,r'$\mathrm{\mathbb{P}}_{\rm r}$',color='blue',fontsize=20)
ax4b.set_xlabel('Time (days)',fontsize=24)
ax4b.set_ylabel('Concentration',fontsize=24)
ax4b.xaxis.set_tick_params(labelsize=20)
ax4b.yaxis.set_tick_params(labelsize=20)
ax4b.set_xlim([0,Tsim])
ax4b.set_ylim([0,4.75])
ax4b.set_ylabel('concentration',fontsize=20,rotation=270, labelpad=25)

plt.setp(ax3.get_xticklabels(), visible=False)
plt.subplots_adjust(bottom=0.13)

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 7A-2.tif', dpi=300, pad_inches=0)

os.chdir(workingdir)

#==============================================================================
#==============================================================================

Tind = np.arange(3)
width = 0.6
#==============================================================================

t1idx_Q1 = np.where(time1>=dfnew['TreatmentPeriod1'][0])[0][1]
t2idx_Q1 = np.where(time1>=dfnew['RelapsePeriod0'][0])[0][1]
t3idx_Q1 = np.where(time1>=Tsim)[0][0]

aveN2N1ratio_Q1t1 = np.mean(Csol1[0:t1idx_Q1]/Isol1[0:t1idx_Q1])
aveN2N1ratio_Q1t2 = np.mean(Csol1[t1idx_Q1:t2idx_Q1]/Isol1[t1idx_Q1:t2idx_Q1])
aveN2N1ratio_Q1t3 = np.mean(Csol1[t2idx_Q1:t3idx_Q1]/Isol1[t2idx_Q1:t3idx_Q1])

aveN2N1ratio_Q1 = (aveN2N1ratio_Q1t1,aveN2N1ratio_Q1t2,aveN2N1ratio_Q1t3)
max_aveN2N1ratio_Q1 = np.max(aveN2N1ratio_Q1)
min_aveN2N1ratio_Q1 = np.min(aveN2N1ratio_Q1)
cost_aveN2N1ratio_Q1 = max_aveN2N1ratio_Q1/min_aveN2N1ratio_Q1
#==============================================================================

t1idx_Q2 = np.where(time2>=dfnew['TreatmentPeriod1'][1])[0][1]
t2idx_Q2 = np.where(time2>=dfnew['RelapsePeriod0'][1])[0][1]
t3idx_Q2 = np.where(time2>=Tsim)[0][0]

aveN2N1ratio_Q2t1 = np.mean(Csol2[0:t1idx_Q2]/Isol2[0:t1idx_Q2])
aveN2N1ratio_Q2t2 = np.mean(Csol2[t1idx_Q2:t2idx_Q2]/Isol2[t1idx_Q2:t2idx_Q2])
aveN2N1ratio_Q2t3 = np.mean(Csol2[t2idx_Q2:t3idx_Q2]/Isol2[t2idx_Q2:t3idx_Q2])

aveN2N1ratio_Q2 = (aveN2N1ratio_Q2t1,aveN2N1ratio_Q2t2,aveN2N1ratio_Q2t3)
max_aveN2N1ratio_Q2 = np.max(aveN2N1ratio_Q2)
min_aveN2N1ratio_Q2 = np.min(aveN2N1ratio_Q2)
cost_aveN2N1ratio_Q2 = max_aveN2N1ratio_Q2/min_aveN2N1ratio_Q2
#==============================================================================

t1idx_Q3 = np.where(time3>=dfnew['TreatmentPeriod1'][2])[0][1]
t2idx_Q3 = np.where(time3>=dfnew['RelapsePeriod0'][2])[0][1]
t3idx_Q3 = np.where(time3>=Tsim)[0][0]

aveN2N1ratio_Q3t1 = np.mean(Csol3[0:t1idx_Q3]/Isol3[0:t1idx_Q3])
aveN2N1ratio_Q3t2 = np.mean(Csol3[t1idx_Q3:t2idx_Q3]/Isol3[t1idx_Q3:t2idx_Q3])
aveN2N1ratio_Q3t3 = np.mean(Csol3[t2idx_Q3:t3idx_Q3]/Isol3[t2idx_Q3:t3idx_Q3])

aveN2N1ratio_Q3 = (aveN2N1ratio_Q3t1,aveN2N1ratio_Q3t2,aveN2N1ratio_Q3t3)
max_aveN2N1ratio_Q3 = np.max(aveN2N1ratio_Q3)
min_aveN2N1ratio_Q3 = np.min(aveN2N1ratio_Q3)
cost_aveN2N1ratio_Q3 = max_aveN2N1ratio_Q3/min_aveN2N1ratio_Q3
#==============================================================================

t1idx_Q4 = np.where(time4>=dfnew['TreatmentPeriod1'][3])[0][1]
t2idx_Q4 = np.where(time4>=dfnew['RelapsePeriod0'][3])[0][1]
t3idx_Q4 = np.where(time4>=Tsim)[0][0]

aveN2N1ratio_Q4t1 = np.mean(Csol4[0:t1idx_Q4]/Isol4[0:t1idx_Q4])
aveN2N1ratio_Q4t2 = np.mean(Csol4[t1idx_Q4:t2idx_Q4]/Isol4[t1idx_Q4:t2idx_Q4])
aveN2N1ratio_Q4t3 = np.mean(Csol4[t2idx_Q4:t3idx_Q4]/Isol4[t2idx_Q4:t3idx_Q4])

aveN2N1ratio_Q4 = (aveN2N1ratio_Q4t1,aveN2N1ratio_Q4t2,aveN2N1ratio_Q4t3)
max_aveN2N1ratio_Q4 = np.max(aveN2N1ratio_Q4)
min_aveN2N1ratio_Q4 = np.min(aveN2N1ratio_Q4)
cost_aveN2N1ratio_Q4 = max_aveN2N1ratio_Q4/min_aveN2N1ratio_Q4

cost_aveN2N1ratio = (cost_aveN2N1ratio_Q1,cost_aveN2N1ratio_Q2,\
                     cost_aveN2N1ratio_Q3,cost_aveN2N1ratio_Q4)
#==============================================================================


TGF_patch = mpatches.Patch(color='lightcoral', label=r'TGF-$\beta$ Inh')
IFN_patch = mpatches.Patch(color='mediumseagreen', label=r'IFN-$\beta$')

#==============================================================================

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True,figsize=(10, 7))

ax1=plt.subplot(211)
ax1.bar(Tind,aveN2N1ratio_Q1,width,color='C1', alpha=0.75, )
ax1.set_xticks(Tind)
ax1.set_xticklabels((r'$\tau_{\mathrm{\mathbb{P}}_{\rm t}}$', \
                     r'$\tau_{\mathrm{\mathbb{P}}_{\rm cf}}$', \
                     r'$\tau_{\mathrm{\mathbb{P}}_{\rm r}}$'))
ax1.xaxis.set_tick_params(labelsize=20)
ax1.yaxis.set_tick_params(labelsize=20)
ax1.set_ylim([0,1.6])
ax1.set_ylabel('average\nN2/N1',multialignment='center',fontsize=20)
ax1.text(-0.8,1.5,r'$Q_1$',fontsize=24, bbox=bbox1)

ax2=plt.subplot(212,sharex=ax1)
ax2.bar(Tind,aveN2N1ratio_Q2,width,color='C0', alpha=0.75, )
ax2.set_xticks(Tind)
ax2.set_xticklabels((r'$\tau_{\mathrm{\mathbb{P}}_{\rm t}}$', \
                     r'$\tau_{\mathrm{\mathbb{P}}_{\rm cf}}$', \
                     r'$\tau_{\mathrm{\mathbb{P}}_{\rm r}}$'))
ax2.xaxis.set_tick_params(labelsize=20)
ax2.yaxis.set_tick_params(labelsize=20)
ax2.set_ylim([0,1.6])
ax2.set_ylabel('average\nN2/N1',multialignment='center',fontsize=20)
ax2.text(-0.8,1.5,r'$Q_2$',fontsize=24, bbox=bbox2)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.subplots_adjust(bottom=0.13)

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 7B-1.tif', dpi=300,pad_inches=0)

os.chdir(workingdir)
#==============================================================================

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True,figsize=(10, 7))

ax3=plt.subplot(211)
ax3.bar(Tind,aveN2N1ratio_Q3,width,color='C2', alpha=0.75, )
ax3.set_xticks(Tind)
ax3.set_xticklabels((r'$\tau_{\mathrm{\mathbb{P}}_{\rm t}}$', \
                     r'$\tau_{\mathrm{\mathbb{P}}_{\rm cf}}$', \
                     r'$\tau_{\mathrm{\mathbb{P}}_{\rm r}}$'))
ax3.xaxis.set_tick_params(labelsize=20)
ax3.yaxis.set_tick_params(labelsize=20)
ax3.set_ylim([0,1.6])
ax3.set_ylabel('average\nN2/N1',multialignment='center',fontsize=20)
ax3.text(-0.8,1.5,r'$Q_3$',fontsize=24, bbox=bbox3)

ax4=plt.subplot(212,sharex=ax3)
ax4.bar(Tind,aveN2N1ratio_Q4,width,color='C3', alpha=0.75, )
ax4.set_xticks(Tind)
ax4.set_xticklabels((r'$\tau_{\mathrm{\mathbb{P}}_{\rm t}}$', \
                     r'$\tau_{\mathrm{\mathbb{P}}_{\rm cf}}$', \
                     r'$\tau_{\mathrm{\mathbb{P}}_{\rm r}}$'))
ax4.xaxis.set_tick_params(labelsize=26)
ax4.yaxis.set_tick_params(labelsize=20)
ax4.set_ylabel('average\nN2/N1',multialignment='center',fontsize=20)
ax4.text(-0.8,5.25,r'$Q_4$',fontsize=24, bbox=bbox4)

plt.setp(ax3.get_xticklabels(), visible=False)
plt.subplots_adjust(bottom=0.13)

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 7B-2.tif', dpi=300,pad_inches=0)

os.chdir(workingdir)

#==============================================================================

width=0.15
fig, ax = plt.subplots(figsize=(10, 7))
Q1 = ax.bar(Tind - 1.5*width, aveN2N1ratio_Q1, width, color='C1', alpha=0.75, label='Q1')
Q2 = ax.bar(Tind - width/2, aveN2N1ratio_Q2, width, color='C0', alpha=0.75, label='Q2')
Q3 = ax.bar(Tind + width/2, aveN2N1ratio_Q3, width, color='C2', alpha=0.75, label='Q3')
Q4 = ax.bar(Tind + 1.5*width, aveN2N1ratio_Q4, width, color='C3', alpha=0.75, label='Q4')
ax.set_xticks(Tind)
ax.set_xticklabels((r'$\tau_{\mathrm{\mathbb{P}}_{\rm t}}$', \
                     r'$\tau_{\mathrm{\mathbb{P}}_{\rm cf}}$', \
                     r'$\tau_{\mathrm{\mathbb{P}}_{\rm r}}$'))
ax.xaxis.set_tick_params(labelsize=30)
ax.yaxis.set_tick_params(labelsize=30)
ax.set_ylabel(r'average N2/N1',fontsize=32)
ax.legend(ncol=2,prop={'size': 24})

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 8A.tif', dpi=300,pad_inches=0)

os.chdir(workingdir)

#==============================================================================

mycolors = ('C1','C0','C2','C3')
Tindnew=np.arange(4)
fig, ax = plt.subplots(figsize=(10, 7))
cost = ax.bar(Tindnew, cost_aveN2N1ratio, 0.65, color=mycolors, alpha=0.75)
ax.set_xticks(Tindnew)
ax.set_xticklabels(('Q1','Q2','Q3','Q4'))
ax.xaxis.set_tick_params(labelsize=30)
ax.yaxis.set_tick_params(labelsize=30)
ax.set_ylabel(r'$\tau_{\mathrm{\max}}/\tau_{\mathrm{\min}}$',fontsize=32)

#save figure in "plots" folder
workingdir = os.getcwd() # accessing current directory
plotsdir = 'plots' 

if not os.path.exists(plotsdir): # making a folder named 'plots'
    os.makedirs(plotsdir) 

os.chdir(plotsdir)
    
plotsdir = os.getcwd()
    
plt.savefig('Figure 8B.tif', dpi=300,pad_inches=0)

os.chdir(workingdir)


