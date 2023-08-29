#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This program creates the dynamics in N1-N2 plane showing the 
nullclines, steady states and the direction fields. The function name to 
be called is mainN1N2dynamics(G) where G is the input such as 0.1, 0.4 or 1.0

@authors: Aurelio A. de Los Reyes V and Yangjin Kim
"""

from __future__ import division     #floating point division
#from scipy.integrate import odeint
import numpy as np                  #library supporting arrays and matrices 
import matplotlib.pyplot as plt     #plotting library
import os as os
import sympy as sm
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm


#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------

def dCdt(C,I,G,params):
    return params['lamda'] + params['lamdaG']*G +\
        params['k1']/(params['k3']**2 + params['alpha']*I**2)-C

def dIdt(C,I,G,params): 
    return params['lamdaS']*params['S'] +\
        params['k2']/(params['k4']**2 + params['beta']*C**2) - params['mu']*I
             
def dTdt(C,I,T,G,params): 
    return params['r']*(1 + C/(params['K']+params['gamma1']*I))\
        *T*(1-T/params['T0']) - params['decayT']*I*T

#------------------------------------------------------------------------------
def ODEsystem(state,t,G,params):
    C, I, T = state 
    dC = dCdt(C,I,G,params)   
    dI = dIdt(C,I,G,params)
    dT = dTdt(C,I,T,G,params)
   
    return [dC, dI, dT]

#------------------------------------------------------------------------------

def mainN1N2dynamics(G):

    if G==0.1:
        init = [1,2.5,0.2]
    elif G==1:
        init = [4,1.0,0.2]
    elif G==0.4:
        init = [1,2.5,0.2]
    else:
        pass
        
   
    params = load_parameters()  
    
    C_1 = np.linspace(0,5,20)
    I_1 = np.linspace(0,5,20)

    
    Cgrid,Igrid  = np.meshgrid(C_1, I_1) 
    
    dCquiv = dCdt(Cgrid,Igrid,G,params)
    dIquiv = dIdt(Cgrid,Igrid,G,params)
    MCI = (np.hypot(dCquiv, dIquiv))                          # norm growth rate 
    MCI[ MCI == 0] = 1.                                 # avoid zero divisio#n errors 
    dCquiv /= MCI                                         # normalize each arrows
    dIquiv /= MCI
    
    Cgrid2, Igrid2 = np.meshgrid(np.arange(0,5,0.01), np.arange(0,5,0.01))
    dCstream = dCdt(Cgrid2,Igrid2,G,params)
    dIstream = dIdt(Cgrid2,Igrid2,G,params)
    
    
    
    Csym,Isym,Tsym = sm.symbols('Csym,Isym,Tsym',negative=False)
    Ceq = dCdt(Csym,Isym,G,params)
    Ieq = dIdt(Csym,Isym,G,params)
    Teq = dTdt(Csym,Isym,Tsym,G,params)
    
    CEqual = sm.Eq(Ceq, 0)
    IEqual = sm.Eq(Ieq, 0)
    TEqual = sm.Eq(Teq, 0)
    
        
    neweq=[]
    equilibria=[]
    if G==0.1 or G==1:   
        equilibria = sm.nsolve([CEqual, IEqual, TEqual],\
                               [Csym, Isym, Tsym],init, maxsteps=50 )
    else:
        init1 = [1,2.5,0.2]
        init2 = [2.0,1.0,0.2]
        init3 = [4,1.0,0.2]
        init = [init1,init2,init3]
        for i in range(len(init)):
            neweq = sm.nsolve([CEqual, IEqual, TEqual],\
                               [Csym, Isym, Tsym],init[i], maxsteps=50 )
            equilibria.append(neweq)
    
    
    Cth = 1.81
    Ith = 1.29
#    Tf = 30

    
    plt.close('all')
    plt.figure(figsize=(10,7))
    plt.quiver(Cgrid, Igrid, dCquiv, dIquiv, MCI,cmap='gist_gray_r',pivot='mid')
    plt.streamplot(Cgrid2,Igrid2,dCstream,dIstream,color=Cgrid2,density=0.5,\
                   cmap='gist_earth')
    plt.contour(Cgrid2, Igrid2, dCstream, levels=[0], linewidths=3, colors='blue')
    plt.contour(Cgrid2, Igrid2, dIstream, levels=[0], linewidths=3, colors='magenta')
    plt.text(3.5, 1, r'tumorigenic region',
             {'color': 'k', 'fontsize': 16, 'ha': 'center', 'va': 'center',
              'bbox': dict(boxstyle="round", fc="w", ec="k", pad=0.2)})
    plt.text(0.9, 4.7, r'anti-tumorigenic region',
             {'color': 'k', 'fontsize': 16, 'ha': 'center', 'va': 'center',
              'bbox': dict(boxstyle="round", fc="w", ec="k", pad=0.2)}) 
    plt.text(4.37,4.5, r'$C$-nullcline', {'color': 'b', 'fontsize': 24,\
            'ha': 'center', 'va': 'center', 'bbox': dict(boxstyle="round",\
            fc="w", ec="w", pad=0.2)})
    plt.text(4.37,4.1, r'$I$-nullcline', {'color': 'm', 'fontsize': 24,\
            'ha': 'center', 'va': 'center', 'bbox': dict(boxstyle="round",\
            fc="w", ec="w", pad=0.2)})
    
    if G==0.1:
        plt.plot(equilibria[0],equilibria[1],'.', color='k', markersize=18)
        plt.text(0.4,3.85,r'SS$^{(s)}$',fontsize=20)
    elif G==0.4:
        plt.plot(equilibria[0][0],equilibria[0][1],'.', color='k', markersize=18)
        plt.plot(equilibria[1][0],equilibria[1][1],'.', color='k', markersize=18,markerfacecolor="None", markeredgewidth=2)
        plt.plot(equilibria[2][0],equilibria[2][1],'.', color='k', markersize=18)
        plt.text(0.8,2.85,r'SS$_1^{(s)}$',fontsize=20)
        plt.text(1.75,0.55,r'SS$^{(u)}$',fontsize=20)
        plt.text(2.8,0.1,r'SS$_2^{(s)}$',fontsize=20)
    else:
        plt.plot(equilibria[0],equilibria[1],'.', color='k', markersize=18)
        plt.text(4.27,0.5,r'SS$^{(s)}$',fontsize=20)
        
    #plt.plot(equilibria[0],equilibria[1],'.', color='r',  markerfacecolor="None",markersize=16)
    plt.axvline(Cth, color='k', linestyle='-.')
    plt.axhline(Ith, color='k', linestyle='-.')
    plt.xlim([0, 5])
    plt.ylim([0, 5])   
    plt.xlabel(r'N2 complex $(C)$',fontsize=28)
    plt.ylabel(r'N1 complex $(I)$',fontsize=28)
    plt.title(r"$G$ = %.1f" %G, fontsize=32)
    plt.tick_params(axis='both',labelsize=26)
    
    x=np.linspace(0,5,1000)
    plt.fill_between(x, 0, Ith, where=x > Cth, facecolor='C3', alpha=0.25)
    plt.fill_between(x, Ith, 5, where=x < Cth, facecolor='C0', alpha=0.25)
    
    plt.subplots_adjust(bottom=0.15)
    
    
    #save figure in "plots" folder
    workingdir = os.getcwd() # accessing current directory
    plotsdir = 'nullclines_varyG' 
    
    if not os.path.exists(plotsdir): # making a folder named 'plots'
        os.makedirs(plotsdir) 
    
    os.chdir(plotsdir)
        
    plotsdir = os.getcwd()
    
    plt.savefig('N1N2dynamics_G='+str(G)+'.tif',dpi=300, pad_inches=0,\
                bbox_inches='tight')
        
    os.chdir(workingdir)
#    