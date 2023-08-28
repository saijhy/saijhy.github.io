#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This code generates the optimal concomitant administration of TGF-beta 
inhibitor and IFN-beta control. It creates a "Condatafile folder" containing 
CSV files of all the information on (1) state variables and control solution; 
(2) injection times; (3) period, frequency, amount, etc. 
The following cases are considered
    (1) tumorigenic region (Q1) with [C0,I0]=[4,0.3]
    (2) safe region (Q2) with [C0,I0]=[3,2]
    (3) anti-tumorigenic region (Q3) with [C0,I0]=[1,2]
    (4) safe region (Q4)with [C0,I0]=[1,0.5]
The above-mentioned cases can be run by simply choosing and uncommenting 
lines 321-324.

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
            'beta':1.0,'mu':1.0,'thC':1.5,'lamdaG':1,'lamdaS':1,'G':0.45}
    #tumor module
    p1 = {'r':0.001,'K':1.0,'gamma1':0.1,'T0':100}
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

def LungCancerModel(y,t,uL,uS,params):
    #model equations
    L = y[0]
    G = y[1]
    C = y[2]
    S = y[3]
    I = y[4]
    T = y[5]
    z1 = y[6]
    z2 = y[7]
      
    #dimensionless ODE model for C(N2 complex), I(N1 complex) and tumor (T)
    dLdt = uL - params['muG']*L
    dGdt = params['Gs'] - params['muG']*G - params['gammaL']*L*G
    dCdt = params['lamda'] + params['lamdaG']*G +\
        params['k1']/(params['k3']**2 + params['alpha']*I**2)-C
    dSdt = uS - params['muS']*S
    dIdt = params['lamdaS']*S +\
        params['k2']/(params['k4']**2 + params['beta']*C**2) - params['mu']*I
    dTdt = params['r']*(1 + C/(params['K']+params['gamma1']*I))\
        *T*(1-T/params['T0']) - params['decayT']*I*T
    dz1dt = uL
    dz2dt = uS
        
    dy = np.zeros(len(y))
    dy[0] = dLdt
    dy[1] = dGdt
    dy[2] = dCdt
    dy[3] = dSdt
    dy[4] = dIdt
    dy[5] = dTdt
    dy[6] = dz1dt
    dy[7] = dz2dt
   
    return dy

#==============================================================================

def RK4ForwardState(y,t,uL,uS,dt,N,params):
    #discretizing the system via RK4
    for i in range(1,N):
        k1 = LungCancerModel(y[i-1,:],t[i-1],uL[i-1],uS[i-1],params)   
        k2 = LungCancerModel(y[i-1,:]+(dt/2)*k1,t[i-1]+(dt/2),\
                             0.5*(uL[i-1]+uL[i]),0.5*(uS[i-1]+uS[i]),params)
        k3 = LungCancerModel(y[i-1,:]+(dt/2)*k2,t[i-1]+(dt/2),\
                             0.5*(uL[i-1]+uL[i]),0.5*(uS[i-1]+uS[i]),params)
        k4 = LungCancerModel(y[i-1,:]+dt*k3,t[i-1]+dt,uL[i],uS[i],params)
    
        y[i,:] = y[i-1,:] + (dt/6)*(k1 + 2*k2 + 2*k3 +k4)
               
    return y

#==============================================================================

def AdjointFunc(y,t,uL,uS,adjoint,params):
    L = y[0]
    G = y[1]
    C = y[2]
    S = y[3]
    I = y[4]
    T = y[5]
    z1 = y[6]
    z2 = y[7]
    
    adjoint1 = adjoint[0]
    adjoint2 = adjoint[1]
    adjoint3 = adjoint[2]
    adjoint4 = adjoint[3]
    adjoint5 = adjoint[4]
    adjoint6 = adjoint[5]
    adjoint7 = adjoint[6]
    adjoint8 = adjoint[7]
    
    
    adjointprime1 = adjoint1*params['muL'] + adjoint2*params['gammaL']*G
    adjointprime2 = adjoint2*(params['muG']+params['gammaL']*L)\
        - adjoint3*params['lamdaG']
    adjointprime3 = adjoint3 + adjoint5*((2*params['k2']*params['beta']*C)/\
        (params['k4']**2+params['beta']*C**2)**2) - adjoint6*params['r']*\
        T*(1-T/params['T0'])*(1/(params['K'] + params['gamma1']*I))
    adjointprime4 = adjoint4*params['muS'] - adjoint5*params['lamdaS']
    adjointprime5 = adjoint3*((2*params['k1']*params['alpha']*I)/\
        (params['k3']**2 + params['alpha']*I**2)**2) + adjoint5*params['mu'] \
        + adjoint6*params['r']*T*(1-T/params['T0'])*((params['gamma1']*C)/\
        (params['K']+params['gamma1']*I)**2) - params['decayT']*T
    adjointprime6 = -1 - adjoint6*(params['r']*(1+C/(params['K'] + \
        params['gamma1']*I))*(1-2*T/params['T0']) - params['decayT']*I)
    adjointprime7 = 0
    adjointprime8 = 0
     
    adjointprime = np.zeros(len(adjoint))
    adjointprime[0] = adjointprime1
    adjointprime[1] = adjointprime2
    adjointprime[2] = adjointprime3
    adjointprime[3] = adjointprime4
    adjointprime[4] = adjointprime5
    adjointprime[5] = adjointprime6
    adjointprime[6] = adjointprime7
    adjointprime[7] = adjointprime8
    
    return adjointprime

#==============================================================================

def RK4BackwardAdjoint(y,t,uL,uS,adjoint,params,dt,N):
    for i in range(1,N):
        j = N + 1 - i
        k1 = AdjointFunc(y[j-1,:],t[j-1],uL[j-1],uS[j-1],adjoint[j-1,:],params)
        
        k2 = AdjointFunc(0.5*(y[j-1,:]+y[j-2,:]),t[j-1]-(dt/2),\
                         0.5*(uL[j-1]+uL[j-2]),0.5*(uS[j-1]+uS[j-2]),\
                         adjoint[j-1,:]-(dt/2)*k1,params)
        
        k3 = AdjointFunc(0.5*(y[j-1,:]+y[j-2,:]),t[j-1]-(dt/2),\
                         0.5*(uL[j-1]+uL[j-2]),0.5*(uS[j-1]+uS[j-2]),\
                         adjoint[j-1,:]-(dt/2)*k2,params)

        k4 = AdjointFunc(y[j-2,:],t[j-2],uL[j-2],uS[j-2],\
                         adjoint[j-1,:]-dt*k3,params)
        
        adjoint[j-2,:] = adjoint[j-1,:] - (dt/6)*(k1 + 2*k2 + 2*k3 +k4)
        
    
    return adjoint

#==============================================================================

def OptiControl(y,T,N,t,dt,params,A1,A2,B1,B2):
    
    test = -1                       #convergence test variable
    
    delta = 0.001                   #tolerance value
    
    uL = np.zeros(N)
    uS = np.zeros(N) 
    #y[0,:] = [C0,I0,T0]          #value of S at t=0 (initial value)
    
    adjoint = np.zeros((N,8))       #declaration for adjoint
    
    adjoint1 = adjoint[:,0]
    adjoint2 = adjoint[:,1]
    adjoint3 = adjoint[:,2]
    adjoint4 = adjoint[:,3]
    adjoint5 = adjoint[:,4]
    adjoint6 = adjoint[:,5]
    adjoint7 = adjoint[:,6]
    adjoint8 = adjoint[:,7]
    
    adjoint7[-1] = -A1
    adjoint8[-1] = -A2
    
    #note:values of x and adjoint will be overwritten in the sweep process
    
    L=y[:,0]
    G=y[:,1]
    C=y[:,2]
    S=y[:,3]
    I=y[:,4]
    T=y[:,5]
    z1=y[:,6]
    z2=y[:,7]
    
#    c=0.9
#    iter=1
    #Forward-Backward Sweep Method (FBSM)
    while (test<0):
        #store previous values of u, x, and adjoint as oldu, oldy, oldadjoint, 
        #respectively
        olduL = uL
        olduS = uS
        oldL = L
        oldG = G
        oldC = C
        oldS = S
        oldI = I
        oldT = T
        oldz1 = z1
        oldz2 = z2
        oldadjoint1 = adjoint1
        oldadjoint2 = adjoint2
        oldadjoint3 = adjoint3
        oldadjoint4 = adjoint4
        oldadjoint5 = adjoint5
        oldadjoint6 = adjoint6
        oldadjoint7 = adjoint7
        oldadjoint8 = adjoint8
            
        #solve for x forward in time using Runge-Kutta method
        y = RK4ForwardState(y,t,uL,uS,dt,N,params)
        
        #solve for adjoint backward in time using Runge-Kutta method
        adjoint = RK4BackwardAdjoint(y,t,uL,uS,adjoint,params,dt,N)
    
        adjoint1 = adjoint[:,0]
        adjoint2 = adjoint[:,1]
        adjoint3 = adjoint[:,2]
        adjoint4 = adjoint[:,3]
        adjoint5 = adjoint[:,4]
        adjoint6 = adjoint[:,5]
        adjoint7 = adjoint[:,6]
        adjoint8 = adjoint[:,7]
    
        uL1 = -(adjoint1+adjoint7)/B1
        uL = 0.5*(uL1 + olduL) 
        
        uS1 = -(adjoint4+adjoint8)/B2
        uS = 0.5*(uS1 + olduS) 
        
        
        #convergence test parameters for the variables
        temp1 = delta*np.sum(np.abs(uL)) - np.sum(np.abs(olduL-uL))
        temp2 = delta*np.sum(np.abs(uS)) - np.sum(np.abs(olduS-uS))
        temp3 = delta*np.sum(np.abs(L)) - np.sum(np.abs(oldL-L))
        temp4 = delta*np.sum(np.abs(G)) - np.sum(np.abs(oldG-G))
        temp5 = delta*np.sum(np.abs(C)) - np.sum(np.abs(oldC-C))
        temp6 = delta*np.sum(np.abs(S)) - np.sum(np.abs(oldS-S))
        temp7 = delta*np.sum(np.abs(I)) - np.sum(np.abs(oldI-I))
        temp8 = delta*np.sum(np.abs(T)) - np.sum(np.abs(oldT-T))
        temp9 = delta*np.sum(np.abs(z1)) - np.sum(np.abs(oldz1-z1))
        temp10 = delta*np.sum(np.abs(z2)) - np.sum(np.abs(oldz2-z2))
        temp11 = delta*np.sum(np.abs(adjoint1)) - \
            np.sum(np.abs(oldadjoint1-adjoint1))
        temp12 = delta*np.sum(np.abs(adjoint2)) - \
            np.sum(np.abs(oldadjoint2-adjoint2))
        temp13 = delta*np.sum(np.abs(adjoint3)) - \
            np.sum(np.abs(oldadjoint3-adjoint3))
        temp14 = delta*np.sum(np.abs(adjoint4)) - \
            np.sum(np.abs(oldadjoint4-adjoint4))
        temp15 = delta*np.sum(np.abs(adjoint5)) - \
            np.sum(np.abs(oldadjoint5-adjoint5))
        temp16 = delta*np.sum(np.abs(adjoint6)) - \
            np.sum(np.abs(oldadjoint6-adjoint6))
        temp17 = delta*np.sum(np.abs(adjoint7)) - \
            np.sum(np.abs(oldadjoint7-adjoint7))
        temp18 = delta*np.sum(np.abs(adjoint8)) - \
            np.sum(np.abs(oldadjoint8-adjoint8))
            
        #minimum among the test values
        test = np.minimum(temp1, np.minimum(temp2, np.minimum(temp3, \
                np.minimum(temp4, np.minimum(temp5,np.minimum(temp6,\
                np.minimum(temp7,np.minimum(temp8,np.minimum(temp9,\
                np.minimum(temp10,np.minimum(temp11,np.minimum(temp12,\
                np.minimum(temp13,np.minimum(temp14,np.minimum(temp15,\
                np.minimum(temp16,np.minimum(temp17,temp18)))))))))))))))))
        
    return [t,y,uL,uS]

#==============================================================================

params = load_parameters()
Cth = 1.81
Ith = 1.29

Tnocontrol = 1
Tnocontrol2 = 15
Tfinal1 = 14#18
Tfinal2 = 60 #[32,60]
Tcontrol = 1

N = 1000
Tinit=0
Tend = 1
Tendnew = 1
#t = np.linspace(Tinit,Tend,N)          #creates N equally spaced nodes bet. 0 & 1
#dt = Tend/N                        #spaces between nodes

L0 = 0
G0 = 0
S0 = 0
T0 = 0.2
z10 = 0
z20 = 0

[C0,I0]=[4,0.3] #tumorigenic region (Q1)
#[C0,I0]=[3,2] #safe region (Q2)
#[C0,I0]=[1,2] #anti-tumorigenic region (Q3)
#[C0,I0]=[1,0.5] #safe region (Q4)

A1min=1
A1= A1min
B1=1

A2min=2
A2inc=2
A2=A2min
A2max=10
B2=1

y = np.zeros((N,8))             
y[0,:] = [L0,G0,C0,S0,I0,T0,z10,z20] 
t = np.linspace(Tinit,Tend,N)
dt = (Tend-Tinit)/(N-1) 
[t1,y1,uL1,uS1] = OptiControl(y,Tend,N,t,dt,params,A1,A2,B1,B2)
L = y1[:,0]
G = y1[:,1]
C = y1[:,2]
S = y1[:,3]
I = y1[:,4]
T = y1[:,5]
z1 = y1[:,6]
z2 = y1[:,7]

Tendnew = Tend+Tnocontrol
tnew = np.linspace(t1[-1],Tendnew,N)
dtnew = (Tendnew-t1[-1])/(N-1)    
ynew = np.zeros((N,8))           
ynew[0,:] = [L[-1],G[-1],C[-1],S[-1],I[-1],T[-1],0,0]
uLnew = np.zeros(N)
uSnew = np.zeros(N)
ynew = RK4ForwardState(ynew,tnew,uLnew,uSnew,dtnew,N,params)

tsol = np.append(t,tnew)
Lsol = np.append(L,ynew[:,0])
Gsol = np.append(G,ynew[:,1])
Csol = np.append(C,ynew[:,2])
Ssol = np.append(S,ynew[:,3])
Isol = np.append(I,ynew[:,4])
Tsol = np.append(T,ynew[:,5])
z1sol = np.append(z1,ynew[:,6])
z2sol = np.append(z2,ynew[:,7])
uLcontrol = np.append(uL1,uLnew)
uScontrol = np.append(uS1,uSnew)

Tinit = tsol[-1]
Tend = tsol[-1]+Tnocontrol

Coninj_init = t1[0] #first IFN and TGF beta inhibitor injection
Coninj_end = t1[-1]

IFNamount = (Coninj_end-Coninj_init)*A2
usedA2 = A2 #amount of IFN beta used

idx=-1
while Tend<Tfinal1:
    if ((Csol[idx]>Cth) or (Isol[idx]<Ith) and (A2<A2max)):
        A2+=A2inc
        A1=A1min
        if A2>A2max:
            A2=A2max
    elif ((Csol[-1]<Cth) and (Isol[-1]>Ith)):
        A1=0
        A2=0
    else:
        A1=A1min
        A2=A2min

    y[0,:] = [Lsol[-1],Gsol[-1],Csol[-1],Ssol[-1],Isol[-1],Tsol[-1],0,0] 
    t = np.linspace(Tinit,Tend,N)
    dt = (Tend-Tinit)/(N-1) 
    [t1,y1,uL1,uS1] = OptiControl(y,Tend,N,t,dt,params,A1,A2,B1,B2)
    Lsol = np.append(Lsol,y1[:,0])
    Gsol = np.append(Gsol,y1[:,1])
    Csol = np.append(Csol,y1[:,2])
    Ssol = np.append(Ssol,y1[:,3])
    Isol = np.append(Isol,y1[:,4])
    Tsol = np.append(Tsol,y1[:,5])
    z1sol = np.append(z1sol,y1[:,6])
    z2sol = np.append(z2sol,y1[:,7])
    tsol = np.append(tsol,t1)
    uLcontrol = np.append(uLcontrol,uL1)
    uScontrol = np.append(uScontrol,uS1)
    
    Coninj_init = np.append(Coninj_init,t1[0]) #injection time
    Coninj_end = np.append(Coninj_end,t1[-1])
    
    usedA2 = np.append(usedA2,A2)
    
    Tendnew = Tend+Tnocontrol
    tnew = np.linspace(t1[-1],Tendnew,N)
    dtnew = (Tendnew-t1[-1])/(N-1)    
    ynew = np.zeros((N,8))           
    ynew[0,:] = [Lsol[-1],Gsol[-1],Csol[-1],Ssol[-1],Isol[-1],Tsol[-1],0,0]
    uLnew = np.zeros(N)
    uSnew = np.zeros(N)
    ynew = RK4ForwardState(ynew,tnew,uLnew,uSnew,dtnew,N,params)
    
    tsol = np.append(tsol,tnew)
    Lsol = np.append(Lsol,ynew[:,0])
    Gsol = np.append(Gsol,ynew[:,1])
    Csol = np.append(Csol,ynew[:,2])
    Ssol = np.append(Ssol,ynew[:,3])
    Isol = np.append(Isol,ynew[:,4])
    Tsol = np.append(Tsol,ynew[:,5])
    z1sol = np.append(z1sol,ynew[:,6])
    z2sol = np.append(z2sol,ynew[:,7])
    uLcontrol = np.append(uLcontrol,uLnew)
    uScontrol = np.append(uScontrol,uSnew)
    
    Tinit = tnew[-1]
    Tend = Tendnew+1

idxCsol = np.where(Csol<Cth)[-1][0]
idxIsol=np.where(Isol>Ith)[-1][0]
if [C0,I0]==[4,0.3]:
    idx = np.minimum(idxCsol,idxIsol)
else:
    idx = np.maximum(idxCsol,idxIsol)

T1inj2 = tsol[idx] #last injection time
Coninj_init = Coninj_init[Coninj_init<T1inj2] 
if [C0,I0]==[4,0.3] or [C0,I0]==[3,2] or [C0,I0]==[1,0.5]:
    Coninj_end = np.append(Coninj_end[Coninj_end<T1inj2],T1inj2)
else:
    Coninj_end = Coninj_end[Coninj_end<T1inj2]

IFNTreatmentAmount = (Coninj_end-Coninj_init)*usedA2[:len(Coninj_init)]
totIFNTreatmentAmount = np.sum(IFNTreatmentAmount)

TGFTreatmentAmount = (Coninj_end-Coninj_init)*A1min
totTGFTreatmentAmount = np.sum(TGFTreatmentAmount)

Coninj2_init = []
Coninj2_end = []


A1=A1min
A2=A2min
usedA2R=A2

while Tend<Tfinal2: 
    Tend = tsol[idx]
    Tendnew = Tend+Tnocontrol2
    tnew = np.linspace(Tend,Tendnew,N)
    dtnew = (Tendnew-Tend)/(N-1)    
    ynew = np.zeros((N,8))           
    ynew[0,:] = [Lsol[idx],Gsol[idx],Csol[idx],Ssol[idx],Isol[idx],Tsol[idx],0,0]
    uLnew = np.zeros(N)
    uSnew = np.zeros(N)
    ynew = RK4ForwardState(ynew,tnew,uLnew,uSnew,dtnew,N,params)
    
    tsol = np.append(tsol[:idx],tnew)
    Lsol = np.append(Lsol[:idx],ynew[:,0])
    Gsol = np.append(Gsol[:idx],ynew[:,1])
    Csol = np.append(Csol[:idx],ynew[:,2])
    Ssol = np.append(Ssol[:idx],ynew[:,3])
    Isol = np.append(Isol[:idx],ynew[:,4])
    Tsol = np.append(Tsol[:idx],ynew[:,5])
    z1sol = np.append(z1sol[:idx],ynew[:,6])
    z2sol = np.append(z2sol[:idx],ynew[:,7])
    uLcontrol = np.append(uLcontrol[:idx],uLnew)
    uScontrol = np.append(uScontrol[:idx],uSnew)
    
    idx = np.where(tsol==Tend)[0][0]
    idxCsol = np.where(Csol[idx:]<Cth)[-1][-1]
    idxIsol = np.where(Isol[idx:]>Ith+0.05)[-1][-1]
    idx = np.minimum(idxCsol,idxIsol)
    idx = np.where(tsol==tnew[idx])[0][0]
          
    tsol = tsol[:idx]
    Lsol = Lsol[:idx]
    Gsol = Gsol[:idx]
    Csol = Csol[:idx]
    Ssol = Ssol[:idx]
    Isol = Isol[:idx]
    Tsol = Tsol[:idx]
    z1sol = z1sol[:idx]
    z2sol = z2sol[:idx]
    uLcontrol = uLcontrol[:idx]
    uScontrol = uScontrol[:idx]
       
    Tinit = tsol[-1]
    Tend = Tinit+Tcontrol
    y[0,:] = [Lsol[-1],Gsol[-1],Csol[-1],Ssol[-1],Isol[-1],Tsol[-1],0,0]
    t = np.linspace(Tinit,Tend,N)
    dt = (Tend-Tinit)/(N-1) 
    [t1,y1,uL1,uS1] = OptiControl(y,Tend,N,t,dt,params,A1,A2,B1,B2)
    Lsol = np.append(Lsol,y1[:,0])
    Gsol = np.append(Gsol,y1[:,1])
    Csol = np.append(Csol,y1[:,2])
    Ssol = np.append(Ssol,y1[:,3])
    Isol = np.append(Isol,y1[:,4])
    Tsol = np.append(Tsol,y1[:,5])
    z1sol = np.append(z1sol,y1[:,6])
    z2sol = np.append(z2sol,y1[:,7])
    tsol = np.append(tsol,t1)
    uLcontrol = np.append(uLcontrol,uL1)
    uScontrol = np.append(uScontrol,uS1)
    
    Coninj2_init = np.append(Coninj2_init,Tinit) #injection time
    Coninj2_end = np.append(Coninj2_end,Tend) #end injection time
    
    usedA2R = np.append(usedA2R,A2)
    
    idx = np.where(tsol==tsol[-1])[0][0]


Coninj2_init = Coninj2_init[Coninj2_init<Tfinal2-Tcontrol]
Coninj2_end = Coninj2_end[Coninj2_end<Tfinal2]

IFNRelapseAmount = (Coninj2_end-Coninj2_init)*usedA2R[:len(Coninj2_init)]
totIFNRelapseAmount = np.sum(IFNRelapseAmount)

TGFRelapseAmount = (Coninj2_end-Coninj2_init)*A1min
totTGFRelapseAmount = np.sum(TGFRelapseAmount)

relapsePeriod = Coninj2_init[-1]-Coninj2_init[-2]
freq = 1/relapsePeriod

Coninj_t1 = np.append(Coninj_init,Coninj2_init)
Coninj_t2 = np.append(Coninj_end,Coninj2_end)

#==============================================================================

#storing data summary in a data frame
if [C0,I0]==[4,0.3] or [C0,I0]==[3,2] or [C0,I0]==[1,0.5]:
    df = pd.DataFrame({'TreatmentPeriod0':Coninj_init[0]},index=[0])
    df['TreatmentPeriod1'] = Coninj_end[-1]
    df['TGFAmountTreatment'] = totTGFTreatmentAmount
    df['IFNAmountTreatment'] = totIFNTreatmentAmount
    df['CancerFreePeriod0'] = Coninj_end[-1]
    df['CancerFreePeriod1'] = Coninj2_init[0]
    df['RelapsePeriod0'] = Coninj2_init[0]
    df['TGFAmountRelapse'] = totTGFRelapseAmount
    df['IFNAmountRelapse'] = totIFNRelapseAmount
    df['totTGFAmount'] = totTGFTreatmentAmount+totTGFRelapseAmount
    df['totIFNAmount'] = totIFNTreatmentAmount+totIFNRelapseAmount
    df['RelapsePeriod'] = relapsePeriod
    df['frequency'] = freq
    
    df1 = pd.DataFrame({'Coninj_t1':Coninj_t1},index=range(len(Coninj_t1)))
    df1['Coninj_t2'] = Coninj_t2
    
else:
    df = pd.DataFrame({'TreatmentPeriod0':Coninj2_init[0]},index=[0])
    df['TreatmentPeriod1'] = Coninj2_end[0]
    df['TGFAmountTreatment'] = TGFRelapseAmount[0]
    df['IFNAmountTreatment'] = IFNRelapseAmount[0]
    df['CancerFreePeriod0'] = Coninj2_end[0]
    df['CancerFreePeriod1'] = Coninj2_init[1]
    df['RelapsePeriod0'] = Coninj2_init[1]
    df['TGFAmountRelapse'] = np.sum(TGFRelapseAmount[1:])
    df['IFNAmountRelapse'] = np.sum(IFNRelapseAmount[1:])
    df['totTGFAmount'] = totTGFRelapseAmount
    df['totIFNAmount'] = totIFNRelapseAmount
    df['RelapsePeriod'] = relapsePeriod
    df['frequency'] = freq
    
    df1 = pd.DataFrame({'Coninj_t1':Coninj2_init},index=range(len(Coninj_t1)))
    df1['Coninj_t2'] = Coninj2_end
    
#==============================================================================
#==============================================================================

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'Condatafile'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

df.to_csv('Consummary_[C0,I0]=['+str(C0)+','+str(I0)+'].csv',index=False)
 
os.chdir(workingdir)

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'Condatafile'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

df1.to_csv('Coninjsummary_[C0,I0]=['+str(C0)+','+str(I0)+'].csv',index=False)
 
os.chdir(workingdir)

#==============================================================================
#==============================================================================

#storing all data in a data frame
dfdata = pd.DataFrame({'time':tsol})
dfdata['Lsol'] = Lsol
dfdata['Gsol'] = Gsol
dfdata['Csol'] = Csol
dfdata['Ssol'] = Ssol
dfdata['Isol'] = Isol
dfdata['Tsol'] = Tsol
dfdata['uLcontrol'] = uLcontrol
dfdata['uScontrol'] = uScontrol

# saving data results in a folder "datafile" as a csv file
workingdir = os.getcwd() # accessing current directory
datadir = 'Condatafile'

if not os.path.exists(datadir): # making a folder named 'plots'
    os.makedirs(datadir) 

os.chdir(datadir)

dfdata.to_csv('Condata_[C0,I0]=['+str(C0)+','+str(I0)+'].csv',index=False)
 
os.chdir(workingdir)