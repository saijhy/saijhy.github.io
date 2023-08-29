#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This code generates the optimal TGF-beta inhibitor  control only. It creates 
a "TGFInhdatafile folder" containing CSV files of all the information on (1)
state variables and control solution; (2) injection times; (3) period, 
frequency, amount, etc. The following cases are considered
    (1) tumorigenic region (Q1) with [C0,I0]=[4,0.3]
    (2) safe region (Q2) with [C0,I0]=[3,2]
    (3) anti-tumorigenic region (Q3) with [C0,I0]=[1,2]
    (4) safe region (Q4)with [C0,I0]=[1,0.5]
The above-mentioned cases can be run by simply choosing and uncommenting 
lines 307-310.

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
            'beta':1.0,'mu':1.0,'lamdaG':1,'lamdaS':1,'S':0.2}
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

def LungCancerModel(y,t,uL,params):
    #model equations
    L = y[0]
    G = y[1]
    C = y[2]
    I = y[3]
    T = y[4]
    z = y[5]
      
    #dimensionless ODE model for C(N2 complex), I(N1 complex) and tumor (T)
    dLdt = uL - params['muG']*L
    dGdt = params['Gs'] - params['muG']*G - params['gammaL']*L*G
    dCdt = params['lamda'] + params['lamdaG']*G +\
        params['k1']/(params['k3']**2 + params['alpha']*I**2)-C
    dIdt = params['lamdaS']*params['S'] +\
        params['k2']/(params['k4']**2 + params['beta']*C**2) - params['mu']*I
    dTdt = params['r']*(1 + C/(params['K']+params['gamma1']*I))\
        *T*(1-T/params['T0']) - params['decayT']*I*T
    dzdt = uL
        
    dy = np.zeros(len(y))
    dy[0] = dLdt
    dy[1] = dGdt
    dy[2] = dCdt
    dy[3] = dIdt
    dy[4] = dTdt
    dy[5] = dzdt
   
    return dy

#==============================================================================

def RK4ForwardState(y,t,uL,dt,N,params):
    #discretizing the system via RK4
    for i in range(1,N):
        k1 = LungCancerModel(y[i-1,:],t[i-1],uL[i-1],params)   
        k2 = LungCancerModel(y[i-1,:]+(dt/2)*k1,t[i-1]+(dt/2),\
                             0.5*(uL[i-1]+uL[i]),params)
        k3 = LungCancerModel(y[i-1,:]+(dt/2)*k2,t[i-1]+(dt/2),\
                             0.5*(uL[i-1]+uL[i]),params)
        k4 = LungCancerModel(y[i-1,:]+dt*k3,t[i-1]+dt,uL[i],params)
    
        y[i,:] = y[i-1,:] + (dt/6)*(k1 + 2*k2 + 2*k3 +k4)
               
    return y

#==============================================================================

def AdjointFunc(y,t,uL,adjoint,params):
    L = y[0]
    G = y[1]
    C = y[2]
    I = y[3]
    T = y[4]
    z = y[5]
    
    adjoint1 = adjoint[0]
    adjoint2 = adjoint[1]
    adjoint3 = adjoint[2]
    adjoint4 = adjoint[3]
    adjoint5 = adjoint[4]
    adjoint6 = adjoint[5]
    
    adjointprime1 = adjoint1*params['muL'] + adjoint2*params['gammaL']*G
    adjointprime2 = adjoint2*(params['muG']+params['gammaL']*L)\
        - adjoint3*params['lamdaG']
    adjointprime3 = adjoint3 + adjoint4*((2*params['k2']*params['beta']*C)/\
        (params['k4']**2+params['beta']*C**2)**2) - adjoint5*params['r']*\
        T*(1-T/params['T0'])*(1/(params['K'] + params['gamma1']*I))
    adjointprime4 = adjoint3*((2*params['k1']*params['alpha']*I)/\
        (params['k3']**2 + params['alpha']*I**2)**2) + adjoint4*params['mu'] \
        + adjoint5*(params['r']*T*(1-T/params['T0'])*((params['gamma1']*C)/\
        (params['K']+params['gamma1']*I)**2) - params['decayT']*T)
    adjointprime5 = -1 - adjoint5*(params['r']*(1+C/(params['K'] + \
        params['gamma1']*I))*(1-2*T/params['T0']) - params['decayT']*I)
    adjointprime6 = 0
     
    adjointprime = np.zeros(len(adjoint))
    adjointprime[0] = adjointprime1
    adjointprime[1] = adjointprime2
    adjointprime[2] = adjointprime3
    adjointprime[3] = adjointprime4
    adjointprime[4] = adjointprime5
    adjointprime[5] = adjointprime6
    
    return adjointprime

#==============================================================================

def RK4BackwardAdjoint(y,t,uL,adjoint,params,dt,N):
    for i in range(1,N):
        j = N + 1 - i
        k1 = AdjointFunc(y[j-1,:],t[j-1],uL[j-1],adjoint[j-1,:],params)
        
        k2 = AdjointFunc(0.5*(y[j-1,:]+y[j-2,:]),t[j-1]-(dt/2),\
                         0.5*(uL[j-1]+uL[j-2]),adjoint[j-1,:]-(dt/2)*k1,params)
        
        k3 = AdjointFunc(0.5*(y[j-1,:]+y[j-2,:]),t[j-1]-(dt/2),\
                         0.5*(uL[j-1]+uL[j-2]),adjoint[j-1,:]-(dt/2)*k2,params)

        k4 = AdjointFunc(y[j-2,:],t[j-2],uL[j-2],adjoint[j-1,:]-dt*k3,params)
        
        adjoint[j-2,:] = adjoint[j-1,:] - (dt/6)*(k1 + 2*k2 + 2*k3 +k4)
        
    
    return adjoint

#==============================================================================

def OptiControl(y,Tend,N,t,dt,params,A2,B2,theta):
    
    test = -1                       #convergence test variable
    
    delta = 0.001                   #tolerance value
    
    uL = np.zeros(N) 
    #y[0,:] = [C0,I0,T0]          #value of S at t=0 (initial value)
    
    adjoint = np.zeros((N,6))       #declaration for adjoint
    
    adjoint1 = adjoint[:,0]
    adjoint2 = adjoint[:,1]
    adjoint3 = adjoint[:,2]
    adjoint4 = adjoint[:,3]
    adjoint5 = adjoint[:,4]
    adjoint6 = adjoint[:,5]
    
    adjoint6[-1] = -theta
    
    #note:values of x and adjoint will be overwritten in the sweep process
    
    L=y[:,0]
    G=y[:,1]
    C=y[:,2]
    I=y[:,3]
    T=y[:,4]
    z=y[:,5]
    
#    c=0.9
#    iter=1
    #Forward-Backward Sweep Method (FBSM)
    while (test<0):
        #store previous values of u, x, and adjoint as oldu, oldy, oldadjoint, 
        #respectively
        olduL = uL
        oldL = L
        oldG = G
        oldC = C
        oldI = I
        oldT = T
        oldz = z
        oldadjoint1 = adjoint1
        oldadjoint2 = adjoint2
        oldadjoint3 = adjoint3
        oldadjoint4 = adjoint4
        oldadjoint5 = adjoint5
        oldadjoint6 = adjoint6
            
        #solve for x forward in time using Runge-Kutta method
        y = RK4ForwardState(y,t,uL,dt,N,params)
        
        #solve for adjoint backward in time using Runge-Kutta method
        adjoint = RK4BackwardAdjoint(y,t,uL,adjoint,params,dt,N)
    
        adjoint1 = adjoint[:,0]
        adjoint2 = adjoint[:,1]
        adjoint3 = adjoint[:,2]
        adjoint4 = adjoint[:,3]
        adjoint5 = adjoint[:,4]
        adjoint6 = adjoint[:,5]
    
        uL1 = -(adjoint1+adjoint6)/B2
        uL = 0.5*(uL1 + olduL) 
        
        
        #convergence test parameters for the variables
        temp1 = delta*np.sum(np.abs(uL)) - np.sum(np.abs(olduL-uL))
        temp2 = delta*np.sum(np.abs(L)) - np.sum(np.abs(oldL-L))
        temp3 = delta*np.sum(np.abs(G)) - np.sum(np.abs(oldG-G))
        temp4 = delta*np.sum(np.abs(C)) - np.sum(np.abs(oldC-C))       
        temp5 = delta*np.sum(np.abs(I)) - np.sum(np.abs(oldI-I))
        temp6 = delta*np.sum(np.abs(T)) - np.sum(np.abs(oldT-T))
        temp7 = delta*np.sum(np.abs(z)) - np.sum(np.abs(oldz-z))
        temp8 = delta*np.sum(np.abs(adjoint1)) - \
            np.sum(np.abs(oldadjoint1-adjoint1))
        temp9 = delta*np.sum(np.abs(adjoint2)) - \
            np.sum(np.abs(oldadjoint2-adjoint2))
        temp10 = delta*np.sum(np.abs(adjoint3)) - \
            np.sum(np.abs(oldadjoint3-adjoint3))
        temp11 = delta*np.sum(np.abs(adjoint4)) - \
            np.sum(np.abs(oldadjoint4-adjoint4))
        temp12 = delta*np.sum(np.abs(adjoint5)) - \
            np.sum(np.abs(oldadjoint5-adjoint5))
        temp13 = delta*np.sum(np.abs(adjoint6)) - \
            np.sum(np.abs(oldadjoint6-adjoint6))
            
        #minimum among the test values
        test = np.minimum(temp1, np.minimum(temp2, np.minimum(temp3, \
                np.minimum(temp4, np.minimum(temp5,np.minimum(temp6,\
                np.minimum(temp7,np.minimum(temp8,np.minimum(temp9,\
                np.minimum(temp10,np.minimum(temp11,np.minimum(temp12,\
                temp13))))))))))))
        
    return [t,y,uL]

#==============================================================================

def secantmethod(y,Tend,N,t,dt,params,A,B,guess1,guess2):
    flag = -1
    
    [ta,ya,uSa] = OptiControl(y,Tend,N,t,dt,params,A,B,guess1)
    za = ya[-1,5] - A
    [tb,yb,uSb] = OptiControl(y,Tend,N,t,dt,params,A,B,guess2)
    zb = yb[-1,5] - A
    
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
        za = y[-1,5] - A
        
        if (np.abs(za) < 1e-10):
            flag=1
            
            
    return [t,y,uS]

#==============================================================================

params = load_parameters()
Cth = 1.81
Ith = 1.29

Tcontrol = 1
Tnocontrol = 1
Tnocontrol2 = 15
Tfinal1 = 20
Tfinal2 = 60 #[40,60]

N = 1000
Tinit=0
Tend = 1
Tendnew = 1

L0 = 0
G0 = 0
T0 = 0.2
z0 = 0

guess1 = 0.5
guess2 = 1.5

[C0,I0]=[4,0.3] #tumorigenic region (Q1)
#[C0,I0]=[3,2] #safe region (Q2)
#[C0,I0]=[1,2] #anti-tumorigenic region (Q3)
#[C0,I0]=[1,0.5] #safe region (Q4)

#idxset = [0]
A2min=1
A2=A2min
B2=1
#umax=10
y = np.zeros((N,6))             #declaration for S 
y[0,:] = [L0,G0,C0,I0,T0,z0] 
t = np.linspace(Tinit,Tend,N)
dt = (Tend-Tinit)/(N-1) 
[t1,y1,uL1] = OptiControl(y,Tend,N,t,dt,params,A2,B2,A2)#secantmethod(y,Tend,N,t,dt,params,A2,B2,guess1,guess2)#OptiControl(y,Tend,N,t,dt,params,A2,B2)
L = y1[:,0]
G = y1[:,1]
C = y1[:,2]
I = y1[:,3]
T = y1[:,4]
z = y1[:,5]

Tendnew = Tend+Tnocontrol
tnew = np.linspace(t1[-1],Tendnew,N)
dtnew = (Tendnew-t1[-1])/(N-1)    
ynew = np.zeros((N,6))           
ynew[0,:] = [L[-1],G[-1],C[-1],I[-1],T[-1],0]
uLnew = np.zeros(N)
ynew = RK4ForwardState(ynew,tnew,uLnew,dtnew,N,params)

tsol = np.append(t,tnew)
Lsol = np.append(L,ynew[:,0])
Gsol = np.append(G,ynew[:,1])
Csol = np.append(C,ynew[:,2])
Isol = np.append(I,ynew[:,3])
Tsol = np.append(T,ynew[:,4])
zsol = np.append(z,ynew[:,5])
control = np.append(uL1,uLnew)

Tinit = tsol[-1]
Tend = tsol[-1]+1

TGFinj_init = t1[0] #first TGF beta inhibitor injection
TGFinj_end = t1[-1]

usedA2 = A2
while Tend<Tfinal1:
    if ((Csol[-1]>Cth) or (Isol[-1]<Ith)):
        A2=A2min
    else: 
        A2=0
    y[0,:] = [Lsol[-1],Gsol[-1],Csol[-1],Isol[-1],Tsol[-1],0] 
    t = np.linspace(Tinit,Tend,N)
    dt = (Tend-Tinit)/(N-1) 
    [t1,y1,uL1] = OptiControl(y,Tend,N,t,dt,params,A2,B2,A2)#secantmethod(y,Tend,N,t,dt,params,A2,B2,guess1,guess2)#
    Lsol = np.append(Lsol,y1[:,0])
    Gsol = np.append(Gsol,y1[:,1])
    Csol = np.append(Csol,y1[:,2])
    Isol = np.append(Isol,y1[:,3])
    Tsol = np.append(Tsol,y1[:,4])
    zsol = np.append(zsol,y1[:,5])
    tsol = np.append(tsol,t1)
    control = np.append(control,uL1)
    
    TGFinj_init = np.append(TGFinj_init,t1[0]) #initial injection time
    TGFinj_end = np.append(TGFinj_end,t1[-1]) #end injection time
    
    usedA2 = np.append(usedA2,A2)
        
    Tendnew = Tend+Tnocontrol
    tnew = np.linspace(t1[-1],Tendnew,N)
    dtnew = (Tendnew-t1[-1])/(N-1)    
    ynew = np.zeros((N,6))           
    ynew[0,:] = [Lsol[-1],Gsol[-1],Csol[-1],Isol[-1],Tsol[-1],0]
    uLnew = np.zeros(N)
    ynew = RK4ForwardState(ynew,tnew,uLnew,dtnew,N,params)
    
    tsol = np.append(tsol,tnew)
    Lsol = np.append(Lsol,ynew[:,0])
    Gsol = np.append(Gsol,ynew[:,1])
    Csol = np.append(Csol,ynew[:,2])    
    Isol = np.append(Isol,ynew[:,3])
    Tsol = np.append(Tsol,ynew[:,4])
    zsol = np.append(zsol,ynew[:,5])
    control = np.append(control,uLnew)
    
    Tinit = tnew[-1]
    Tend = Tendnew+1


idxCsol = np.where(Csol<Cth)[-1][0]
idxIsol=np.where(Isol>Ith)[-1][0]
if [C0,I0]==[4,0.3]:
    idx = np.minimum(idxCsol,idxIsol)
else:
    idx = np.maximum(idxCsol,idxIsol)


T1inj2 = tsol[idx] #last injection time

TGFinj_init = TGFinj_init[TGFinj_init<T1inj2] 
if [C0,I0]==[1,0.5]:
    TGFinj_end = TGFinj_end[TGFinj_end<T1inj2]
else:
    TGFinj_end = np.append(TGFinj_end[TGFinj_end<T1inj2],T1inj2)
TGFamount = (TGFinj_end-TGFinj_init)*usedA2[:len(TGFinj_init)]

totTGFamount = np.sum(TGFamount)

TGFinj2_init = []
TGFinj2_end = []
A2=A2min
usedA2R=A2
while Tend<Tfinal2: 
    Tend = tsol[idx]
    Tendnew = Tend+Tnocontrol2
    tnew = np.linspace(Tend,Tendnew,N)
    dtnew = (Tendnew-Tend)/(N-1)    
    ynew = np.zeros((N,6))           
    ynew[0,:] = [Lsol[idx],Gsol[idx],Csol[idx],Isol[idx],Tsol[idx],0]
    uLnew = np.zeros(N)
    ynew = RK4ForwardState(ynew,tnew,uLnew,dtnew,N,params)
    
    tsol = np.append(tsol[:idx],tnew)
    Lsol = np.append(Lsol[:idx],ynew[:,0])
    Gsol = np.append(Gsol[:idx],ynew[:,1])
    Csol = np.append(Csol[:idx],ynew[:,2])
    Isol = np.append(Isol[:idx],ynew[:,3])
    Tsol = np.append(Tsol[:idx],ynew[:,4])
    zsol = np.append(zsol[:idx],ynew[:,5])
    control = np.append(control[:idx],uLnew)
    
    
    idx = np.where(tsol==Tend)[0][0]
    idxCsol = np.where(Csol[idx:]<Cth-0.04)[-1][-1]
    idxIsol = np.where(Isol[idx:]>Ith+0.05)[-1][-1]
    idx = np.minimum(idxCsol,idxIsol)
    idx = np.where(tsol==tnew[idx])[0][0]
          
    tsol = tsol[:idx]
    Lsol = Lsol[:idx]
    Gsol = Gsol[:idx]
    Csol = Csol[:idx]
    Isol = Isol[:idx]
    Tsol = Tsol[:idx]
    zsol = zsol[:idx]
    control = control[:idx]
    
    Tinit = tsol[-1]
    Tend = Tinit+Tcontrol
    y[0,:] = [Lsol[-1],Gsol[-1],Csol[-1],Isol[-1],Tsol[-1],0]
    t = np.linspace(Tinit,Tend,N)
    dt = (Tend-Tinit)/(N-1) 
    [t1,y1,uL1] = OptiControl(y,Tend,N,t,dt,params,A2,B2,A2)#secantmethod(y,Tend,N,t,dt,params,A2,B2,guess1,guess2)#
    Lsol = np.append(Lsol,y1[:,0])
    Gsol = np.append(Gsol,y1[:,1])
    Csol = np.append(Csol,y1[:,2])
    Isol = np.append(Isol,y1[:,3])
    Tsol = np.append(Tsol,y1[:,4])
    zsol = np.append(zsol,y1[:,5])
    tsol = np.append(tsol,t1)
    control = np.append(control,uL1)

    
    TGFinj2_init = np.append(TGFinj2_init,Tinit) #injection time
    TGFinj2_end = np.append(TGFinj2_end,Tend) #end injection time
    
    usedA2R = np.append(usedA2R,A2)
    
    idx = np.where(tsol==tsol[-1])[0][0]


TGFinj2_init = TGFinj2_init[TGFinj2_init<Tfinal2-Tcontrol]
TGFinj2_end = TGFinj2_end[TGFinj2_end<Tfinal2]
TGFamount2 = (TGFinj2_end-TGFinj2_init)*usedA2R[:len(TGFinj2_init)]

totTGFamount2 = np.sum(TGFamount2)

relapsePeriod = TGFinj2_init[-1]-TGFinj2_init[-2]
freq = 1/relapsePeriod

TGFinj_t1 = np.append(TGFinj_init,TGFinj2_init)
TGFinj_t2 = np.append(TGFinj_end,TGFinj2_end)

#storing data summary in a data frame
if [C0,I0]==[4,0.3] or [C0,I0]==[3,2] or [C0,I0]==[1,0.5]:
    df = pd.DataFrame({'TreatmentPeriod0':TGFinj_init[0]},index=[0])
    df['TreatmentPeriod1'] = TGFinj_end[-1]
    df['TGFAmountTreatment'] = totTGFamount
    df['RelapsePeriod0'] = TGFinj2_init[0]
    df['RelapsePeriod1'] = TGFinj2_init[-1]
    df['TGFAmountRelapse'] = totTGFamount2
    df['RelapsePeriod'] = relapsePeriod
    df['frequency'] = freq
    df['totTGFAmount'] = totTGFamount+totTGFamount2
    
    df1 = pd.DataFrame({'TGFinj_t1':TGFinj_t1},index=range(len(TGFinj_t1)))
    df1['TGFinj_t2'] = TGFinj_t2
else:
    df = pd.DataFrame({'TreatmentPeriod0':TGFinj2_init[0]},index=[0])
    df['TreatmentPeriod1'] = TGFinj2_end[0]
    df['TGFAmountTreatment'] = (TGFinj2_end[0]-TGFinj2_init[0])*A2min
    df['RelapsePeriod0'] = TGFinj2_init[1]
    df['RelapsePeriod1'] = TGFinj2_init[-1]
    df['TGFAmountRelapse'] = np.sum([(TGFinj2_end[i]-TGFinj2_init[i])*A2min\
      for i in range(1,len(TGFinj2_init))])
    df['RelapsePeriod'] = relapsePeriod
    df['frequency'] = freq
    df['totTGFAmount'] = np.sum([(TGFinj2_end[i]-TGFinj2_init[i])*A2min\
      for i in range(len(TGFinj2_init))])
        
    df1 = pd.DataFrame({'TGFinj_t1':TGFinj2_init},index=range(len(TGFinj2_init)))
    df1['TGFinj_t2'] = TGFinj2_end
    
#==============================================================================
#==============================================================================

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'TGFInhdatafile'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

df.to_csv('TGFInhsummary_[C0,I0]=['+str(C0)+','+str(I0)+'].csv',index=False)
 
os.chdir(workingdir)

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'TGFInhdatafile'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

df1.to_csv('TGFInhinjsummary_[C0,I0]=['+str(C0)+','+str(I0)+'].csv',index=False)
 
os.chdir(workingdir)


#==============================================================================
#==============================================================================

#storing all data in a data frame
dfdata = pd.DataFrame({'time':tsol})
dfdata['Lsol'] = Lsol
dfdata['Gsol'] = Gsol
dfdata['Csol'] = Csol
dfdata['Isol'] = Isol
dfdata['Tsol'] = Tsol
dfdata['control'] = control

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'TGFInhdatafile'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

dfdata.to_csv('TGFInhdata_[C0,I0]=['+str(C0)+','+str(I0)+'].csv',index=False)
 
os.chdir(workingdir)
        
