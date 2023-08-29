#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This code generates the optimal IFN-beta control only. It creates 
a "TGFInhdatafile folder" containing CSV files of all the information on (1)
state variables and control solution; (2) injection times; (3) period, 
frequency, amount, etc. The following cases are considered
    (1) tumorigenic region (Q1) with [C0,I0]=[4,0.3]
    (2) safe region (Q2) with [C0,I0]=[3,2]
    (3) anti-tumorigenic region (Q3) with [C0,I0]=[1,2]
    (4) safe region (Q4)with [C0,I0]=[1,0.5]
The above-mentioned cases can be run by simply choosing and uncommenting 
lines 290-293.

@authors: Aurelio A. de los Reyes V and Yangjin Kim
"""

from __future__ import division     #floating point division
import numpy as np                  #library supporting arrays and matrices 
import os as os
import pandas as pd

#==============================================================================

def load_parameters():
    #model parameters stored in a dictionary
    #N1/N2 modules
    params = {'lamda':0.01,'k1':4.0,'k3':1.0,'alpha':1.5,'k2':4.0,'k4':1.0,\
            'beta':1.0,'mu':1.0,'lamdaG':1,'lamdaS':1,'G':0.45}
    #tumor module
    p1 = {'r':0.05,'K':1.0,'gamma1':0.1,'T0':100}
    #therapeutics
    p2 = {'Gs':0.826,'muG':0.826,'gammaL':100,'LS':13,'muL':6.6,\
                    'muS':3.96}
    #tumor decay
    p3 = {'decayT':0.005}
    
    params.update(p1)
    params.update(p2)
    params.update(p3)

    return params

#==============================================================================

def LungCancerModel(y,t,uS,params):
    #model equations
    C = y[0]
    S = y[1]
    I = y[2]
    T = y[3]
    z = y[4]
      
    #dimensionless ODE model for C(N2 complex), I(N1 complex) and tumor (T)
    dCdt = params['lamda'] + params['lamdaG']*params['G'] +\
        params['k1']/(params['k3']**2 + params['alpha']*I**2)-C
    dSdt = uS - params['muS']*S
    dIdt = params['lamdaS']*S +\
        params['k2']/(params['k4']**2 + params['beta']*C**2) - params['mu']*I
    dTdt = params['r']*(1 + C/(params['K']+params['gamma1']*I))\
        *T*(1-T/params['T0']) - params['decayT']*I*T
    dzdt = uS
        
    dy = np.zeros(len(y))
    dy[0] = dCdt
    dy[1] = dSdt
    dy[2] = dIdt
    dy[3] = dTdt
    dy[4] = dzdt
   
    return dy

#==============================================================================

def RK4ForwardState(y,t,uS,dt,N,params):
    #discretizing the system via RK4
    for i in range(1,N):
        k1 = LungCancerModel(y[i-1,:],t[i-1],uS[i-1],params)   
        k2 = LungCancerModel(y[i-1,:]+(dt/2)*k1,t[i-1]+(dt/2),\
                             0.5*(uS[i-1]+uS[i]),params)
        k3 = LungCancerModel(y[i-1,:]+(dt/2)*k2,t[i-1]+(dt/2),\
                             0.5*(uS[i-1]+uS[i]),params)
        k4 = LungCancerModel(y[i-1,:]+dt*k3,t[i-1]+dt,uS[i],params)
    
        y[i,:] = y[i-1,:] + (dt/6)*(k1 + 2*k2 + 2*k3 +k4)
               
    return y

#==============================================================================

def AdjointFunc(y,t,uS,adjoint,params):
    C = y[0]
    S = y[1]
    I = y[2]
    T = y[3]
    z = y[4]
    
    adjoint1 = adjoint[0]
    adjoint2 = adjoint[1]
    adjoint3 = adjoint[2]
    adjoint4 = adjoint[3]
    adjoint5 = adjoint[4]
    
    
    adjointprime1 = adjoint1 + adjoint3*((2*params['k2']*params['beta']*C)/\
        (params['k4']**2+params['beta']*C**2)**2) - adjoint4*params['r']*\
        T*(1-T/params['T0'])*(1/(params['K'] + params['gamma1']*I))
    adjointprime2 = adjoint2*params['muS'] - adjoint3*params['lamdaS']
    adjointprime3 = adjoint1*((2*params['k1']*params['alpha']*I)/\
        (params['k3']**2 + params['alpha']*I**2)**2) + adjoint3*params['mu'] \
        + adjoint4*(params['r']*T*(1-T/params['T0'])*((params['gamma1']*C)/\
        (params['K']+params['gamma1']*I)**2) - params['decayT']*T)
    adjointprime4 = -1 - adjoint4*(params['r']*(1+C/(params['K'] + \
        params['gamma1']*I))*(1-2*T/params['T0']) - params['decayT']*I) 
    adjointprime5 = 0
     
    adjointprime = np.zeros(len(adjoint))
    adjointprime[0] = adjointprime1
    adjointprime[1] = adjointprime2
    adjointprime[2] = adjointprime3
    adjointprime[3] = adjointprime4
    adjointprime[4] = adjointprime5
    
    return adjointprime

#==============================================================================

def RK4BackwardAdjoint(y,t,uS,adjoint,params,dt,N):
    for i in range(1,N):
        j = N + 1 - i
        k1 = AdjointFunc(y[j-1,:],t[j-1],uS[j-1],adjoint[j-1,:],params)
        
        k2 = AdjointFunc(0.5*(y[j-1,:]+y[j-2,:]),t[j-1]-(dt/2),\
                         0.5*(uS[j-1]+uS[j-2]),adjoint[j-1,:]-(dt/2)*k1,params)
        
        k3 = AdjointFunc(0.5*(y[j-1,:]+y[j-2,:]),t[j-1]-(dt/2),\
                         0.5*(uS[j-1]+uS[j-2]),adjoint[j-1,:]-(dt/2)*k2,params)

        k4 = AdjointFunc(y[j-2,:],t[j-2],uS[j-2],adjoint[j-1,:]-dt*k3,params)
        
        adjoint[j-2,:] = adjoint[j-1,:] - (dt/6)*(k1 + 2*k2 + 2*k3 +k4)
        
    
    return adjoint

#==============================================================================

def OptiControl(y,T,N,t,dt,params,A,B,theta):
    
    test = -1                       #convergence test variable
    
    delta = 0.001                   #tolerance value
    
    uS = np.zeros(N) 
    #y[0,:] = [C0,I0,T0]          #value of S at t=0 (initial value)
    
    adjoint = np.zeros((N,5))       #declaration for adjoint
    
    adjoint1 = adjoint[:,0]
    adjoint2 = adjoint[:,1]
    adjoint3 = adjoint[:,2]
    adjoint4 = adjoint[:,3]
    adjoint5 = adjoint[:,4]
    
    adjoint5[-1] = theta
    
    #note:values of x and adjoint will be overwritten in the sweep process
    
    C=y[:,0]
    S=y[:,1]
    I=y[:,2]
    T=y[:,3]
    z=y[:,4]
    
#    c=0.9
#    iter=1
    #Forward-Backward Sweep Method (FBSM)
    while (test<0):
        #store previous values of u, x, and adjoint as oldu, oldy, oldadjoint, 
        #respectively
        olduS = uS
        oldC = C
        oldS = S
        oldI = I
        oldT = T
        oldz = z
        oldadjoint1 = adjoint1
        oldadjoint2 = adjoint2
        oldadjoint3 = adjoint3
        oldadjoint4 = adjoint4
        oldadjoint5 = adjoint5
            
        #solve for x forward in time using Runge-Kutta method
        y = RK4ForwardState(y,t,uS,dt,N,params)
        
        #solve for adjoint backward in time using Runge-Kutta method
        adjoint = RK4BackwardAdjoint(y,t,uS,adjoint,params,dt,N)
    
        adjoint1 = adjoint[:,0]
        adjoint2 = adjoint[:,1]
        adjoint3 = adjoint[:,2]
        adjoint4 = adjoint[:,3]
        adjoint5 = adjoint[:,4]
    
        uS1 = max(A,min(-(adjoint2+adjoint5)/B),0)
        uS = 0.5*(uS1 + olduS) 
        
        
        #convergence test parameters for the variables
        temp1 = delta*np.sum(np.abs(uS)) - np.sum(np.abs(olduS-uS))
        temp2 = delta*np.sum(np.abs(C)) - np.sum(np.abs(oldC-C))
        temp3 = delta*np.sum(np.abs(S)) - np.sum(np.abs(oldS-S))
        temp4 = delta*np.sum(np.abs(I)) - np.sum(np.abs(oldI-I))
        temp5 = delta*np.sum(np.abs(T)) - np.sum(np.abs(oldT-T))
        temp6 = delta*np.sum(np.abs(z)) - np.sum(np.abs(oldz-z))
        temp7 = delta*np.sum(np.abs(adjoint1)) - \
            np.sum(np.abs(oldadjoint1-adjoint1))
        temp8 = delta*np.sum(np.abs(adjoint2)) - \
            np.sum(np.abs(oldadjoint2-adjoint2))
        temp9 = delta*np.sum(np.abs(adjoint3)) - \
            np.sum(np.abs(oldadjoint3-adjoint3))
        temp10 = delta*np.sum(np.abs(adjoint4)) - \
            np.sum(np.abs(oldadjoint4-adjoint4))
        temp11 = delta*np.sum(np.abs(adjoint5)) - \
            np.sum(np.abs(oldadjoint5-adjoint5))
            
        #minimum among the test values
        test = np.minimum(temp1, np.minimum(temp2, np.minimum(temp3, \
                np.minimum(temp4, np.minimum(temp5,np.minimum(temp6,\
                np.minimum(temp7,np.minimum(temp8,np.minimum(temp9,\
                np.minimum(temp10,temp11))))))))))
        
    return [t,y,uS]

#==============================================================================

def secantmethod(y,Tend,N,t,dt,params,A,B,guess1,guess2):
    flag = -1
    
    [ta,ya,uSa] = OptiControl(y,Tend,N,t,dt,params,A,B,guess1)
    za = ya[-1,4] - A
    [tb,yb,uSb] = OptiControl(y,Tend,N,t,dt,params,A,B,guess2)
    zb = yb[-1,4] - A
    
    while(flag<0):
        if (np.abs(za) > np.abs(zb)):
            k = guess1
            guess1 = guess2
            guess2 = k
            za = zb
            zb = k;
            
        d = za*(guess2-guess1)/(zb-za)
        guess2 = guess1
        zb = za
        guess1 = guess1 - d
        [t,y,uS] = OptiControl(y,Tend,N,t,dt,params,A,B,guess1)
        za = y[-1,4] - A
        
        if (np.abs(za) < 1e-10):
            flag=1
            
            
    return [t,y,uS]

#==============================================================================

params = load_parameters()
Cth = 1.81
Ith = 1.29

Tcontrol = 1
Tnocontrol = 1
Tnocontrol2 = 20
Tfinal1 = 15
Tfinal2 = 60

N = 1000
Tinit=0
Tend = 1
Tendnew = 1

S0 = 0
T0 = 0.2
z0 = 0

guess1 = 0.5
guess2 = 1.5

[C0,I0]=[4,0.3] #tumorigenic region (Q1)
#[C0,I0]=[3,2] #safe region (Q2)
#[C0,I0]=[1,2] #anti-tumorigenic region (Q3)
#[C0,I0]=[1,0.5] #safe region (Q4)

#idxset = [0]
Amin=2
Ainc=2
A=Amin
Amax=10
B=1
#umax=10
y = np.zeros((N,5))             #declaration for S 
y[0,:] = [C0,S0,I0,T0,z0] 
t = np.linspace(Tinit,Tend,N)
dt = (Tend-Tinit)/(N-1) 
[t1,y1,uS1] = OptiControl(y,Tend,N,t,dt,params,A,B,-A)#secantmethod(y,Tend,N,t,dt,params,A,B,guess1,guess2)
C = y1[:,0]
S = y1[:,1]
I = y1[:,2]
T = y1[:,3]
z = y1[:,4]

Tendnew = Tend+Tnocontrol
tnew = np.linspace(t1[-1],Tendnew,N)
dtnew = (Tendnew-t1[-1])/(N-1)    
ynew = np.zeros((N,5))           
ynew[0,:] = [C[-1],S[-1],I[-1],T[-1],0]
uSnew = np.zeros(N)
ynew = RK4ForwardState(ynew,tnew,uSnew,dtnew,N,params)

tsol = np.append(t,tnew)
Csol = np.append(C,ynew[:,0])
Ssol = np.append(S,ynew[:,1])
Isol = np.append(I,ynew[:,2])
Tsol = np.append(T,ynew[:,3])
zsol = np.append(z,ynew[:,4])
control = np.append(uS1,uSnew)

Tinit = tsol[-1]
Tend = tsol[-1]+Tcontrol


IFNinj_init = t1[0] #first IFN beta injection
IFNinj_end = t1[-1] #end of first IFN beta injection  

IFNamount = (IFNinj_end-IFNinj_init)*A
usedA = A #amount of IFN beta used

while Tend<Tfinal1:# or (Csol[-1]>Cth) and (Isol[-1]<Ith):
    if ((Csol[-1]>Cth) or (Isol[-1]<Ith) and (A<Amax)):
        A+=Ainc
        if A>Amax:
            A=Amax
    else:
        A=Amin

    y[0,:] = [Csol[-1],Ssol[-1],Isol[-1],Tsol[-1],0] 
    t = np.linspace(Tinit,Tend,N)
    dt = (Tend-Tinit)/(N-1) 
    [t1,y1,uS1] = OptiControl(y,Tend,N,t,dt,params,A,B,-A)#secantmethod(y,Tend,N,t,dt,params,A,B,guess1,guess2)
    Csol = np.append(Csol,y1[:,0])
    Ssol = np.append(Ssol,y1[:,1])
    Isol = np.append(Isol,y1[:,2])
    Tsol = np.append(Tsol,y1[:,3])
    zsol = np.append(zsol,y1[:,4])
    tsol = np.append(tsol,t1)
    control = np.append(control,uS1)
    
    IFNinj_init = np.append(IFNinj_init,t1[0])
    IFNinj_end = np.append(IFNinj_end,t1[-1])
    
    usedA = np.append(usedA,A)
    
    Tendnew = Tend+Tnocontrol
    tnew = np.linspace(t1[-1],Tendnew,N)
    dtnew = (Tendnew-t1[-1])/(N-1)    
    ynew = np.zeros((N,5))           
    ynew[0,:] = [Csol[-1],Ssol[-1],Isol[-1],Tsol[-1],0]
    uSnew = np.zeros(N)
    ynew = RK4ForwardState(ynew,tnew,uSnew,dtnew,N,params)
    
    tsol = np.append(tsol,tnew)
    Csol = np.append(Csol,ynew[:,0])
    Ssol = np.append(Ssol,ynew[:,1])
    Isol = np.append(Isol,ynew[:,2])
    Tsol = np.append(Tsol,ynew[:,3])
    zsol = np.append(zsol,ynew[:,4])
    control = np.append(control,uSnew)
    
    Tinit = tnew[-1]
    Tend = Tendnew+Tcontrol

idxCsol = np.where(Csol<Cth)[-1][0]
idxIsol=np.where(Isol>Ith)[-1][0]
if [C0,I0]==[4,0.3]:
    idx = np.minimum(idxCsol,idxIsol)
else:
    idx = np.maximum(idxCsol,idxIsol)


T1inj2 = tsol[idx] #last IFN beta injection time
IFNinj_init = IFNinj_init[IFNinj_init<T1inj2] 
if [C0,I0]==[3,2]:
    IFNinj_end = IFNinj_end[IFNinj_end<T1inj2]
else:
    IFNinj_end = np.append(IFNinj_end[IFNinj_end<T1inj2],T1inj2)

IFNamount = (IFNinj_end-IFNinj_init)*usedA[:len(IFNinj_init)]
totIFNamount = np.sum(IFNamount)

IFNinj2_init = []
IFNinj2_end = []

A=Amin
usedA2 = A
while Tend<Tfinal2: 
    Tend = tsol[idx]
    Tendnew = Tend+Tnocontrol2
    tnew = np.linspace(Tend,Tendnew,N)
    dtnew = (Tendnew-Tend)/(N-1)    
    ynew = np.zeros((N,5))           
    ynew[0,:] = [Csol[idx],Ssol[idx],Isol[idx],Tsol[idx],0]
    uSnew = np.zeros(N)
    ynew = RK4ForwardState(ynew,tnew,uSnew,dtnew,N,params)
    
    tsol = np.append(tsol[:idx],tnew)
    Csol = np.append(Csol[:idx],ynew[:,0])
    Ssol = np.append(Ssol[:idx],ynew[:,1])
    Isol = np.append(Isol[:idx],ynew[:,2])
    Tsol = np.append(Tsol[:idx],ynew[:,3])
    zsol = np.append(zsol[:idx],ynew[:,4])
    control = np.append(control[:idx],uSnew)
    
    idx = np.where(tsol==Tend)[0][0]
    idxCsol = np.where(Csol[idx:]<Cth+0.09)[-1][-1]
    idxIsol = np.where(Isol[idx:]>Ith+0.03)[-1][-1]
    idx = np.minimum(idxCsol,idxIsol)
    idx = np.where(tsol==tnew[idx])[0][0]

    tsol = tsol[:idx]
    Csol = Csol[:idx]
    Ssol = Ssol[:idx]
    Isol = Isol[:idx]
    Tsol = Tsol[:idx]
    zsol = zsol[:idx]
    control = control[:idx]
    
    Tinit = tsol[-1]
    Tend = Tinit+Tcontrol
    tidx = np.where(tsol==Tinit)[0][0]
    y[0,:] = [Csol[-1],Ssol[-1],Isol[-1],Tsol[-1],0] 
    t = np.linspace(Tinit,Tend,N)
    dt = (Tend-Tinit)/(N-1) 
    [t1,y1,uS1] = OptiControl(y,Tend,N,t,dt,params,A,B,-A)#secantmethod(y,Tend,N,t,dt,params,A,B,guess1,guess2)
    Csol = np.append(Csol,y1[:,0])
    Ssol = np.append(Ssol,y1[:,1])
    Isol = np.append(Isol,y1[:,2])
    Tsol = np.append(Tsol,y1[:,3])
    zsol = np.append(zsol,y1[:,4])
    tsol = np.append(tsol,t1)
    control = np.append(control,uS1)
    
    IFNinj2_init = np.append(IFNinj2_init,Tinit) #injection time
    IFNinj2_end = np.append(IFNinj2_end,Tend) #end injection time
    
    usedA2 = np.append(usedA2,A)
    
    idx = np.where(tsol==tsol[-1])[0][0]
    
IFNinj2_init = IFNinj2_init[IFNinj2_init<Tfinal2-Tcontrol]
IFNinj2_end = IFNinj2_end[IFNinj2_end<Tfinal2]
IFNamount2 = (IFNinj2_end-IFNinj2_init)*usedA2[:len(IFNinj2_init)]

totIFNamount2 = np.sum(IFNamount2)

relapsePeriod = IFNinj2_init[-1]-IFNinj2_init[-2]
freq = 1/relapsePeriod

IFNinj_t1 = np.append(IFNinj_init,IFNinj2_init)
IFNinj_t2 = np.append(IFNinj_end,IFNinj2_end)

#storing data summary in a data frame
if [C0,I0]==[4,0.3] or [C0,I0]==[3,2] or [C0,I0]==[1,0.5]:
    df = pd.DataFrame({'TreatmentPeriod0':IFNinj_init[0]},index=[0])
    df['TreatmentPeriod1'] = IFNinj_end[-1]
    df['IFNAmountTreatment'] = totIFNamount
    df['CritPeriod0'] = IFNinj_end[-1]
    df['CritPeriod1'] = IFNinj2_init[2]
    df['IFNAmountCrit'] = np.sum(IFNamount2[:2])
    df['MainPeriod0'] = IFNinj2_init[2]
    df['MainPeriod1'] = IFNinj2_end[-1]
    df['IFNAmountMain'] = np.sum(IFNamount2[2:])
    df['totIFNAmount'] = totIFNamount+totIFNamount2
    df['RelapsePeriod'] = relapsePeriod
    df['frequency'] = freq
    
    df1 = pd.DataFrame({'IFNinj_t1':IFNinj_t1},index=range(len(IFNinj_t1)))
    df1['IFNinj_t2'] = IFNinj_t2
    
else:# [C0,I0]==[1,2]:
    df = pd.DataFrame({'TreatmentPeriod0':IFNinj2_init[0]},index=[0])
    df['TreatmentPeriod1'] = IFNinj2_init[1]
    df['IFNAmountTreatment'] = IFNamount2[0]
    df['CritPeriod0'] = IFNinj2_init[1]
    df['CritPeriod1'] = IFNinj2_init[2]
    df['IFNAmountCrit'] = IFNamount2[1]
    df['MainPeriod0'] = IFNinj2_init[2]
    df['MainPeriod1'] = IFNinj2_end[-1]
    df['IFNAmountMain'] = np.sum(IFNamount2[1:])
    df['totIFNAmount'] = totIFNamount+totIFNamount2
    df['RelapsePeriod'] = relapsePeriod
    df['frequency'] = freq
    
    df1 = pd.DataFrame({'IFNinj_t1':IFNinj2_init},index=range(len(IFNinj2_init)))
    df1['IFNinj_t2'] = IFNinj2_end
    
#==============================================================================
#==============================================================================

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'IFNdatafile'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

df.to_csv('IFNsummary_[C0,I0]=['+str(C0)+','+str(I0)+'].csv',index=False)
 
os.chdir(workingdir)

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'IFNdatafile'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

df1.to_csv('IFNinjsummary_[C0,I0]=['+str(C0)+','+str(I0)+'].csv',index=False)
 
os.chdir(workingdir)

#==============================================================================
#==============================================================================

#storing all data in a data frame
dfdata = pd.DataFrame({'time':tsol})
dfdata['Csol'] = Csol
dfdata['Ssol'] = Ssol
dfdata['Isol'] = Isol
dfdata['Tsol'] = Tsol
dfdata['control'] = control

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'IFNdatafile'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

dfdata.to_csv('IFNdata_[C0,I0]=['+str(C0)+','+str(I0)+'].csv',index=False)
 
os.chdir(workingdir)
        
