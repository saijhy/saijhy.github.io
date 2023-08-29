#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This program creates a bifurcation diagram varying both the parameters G and S.

@authors: Aurelio A. de Los Reyes V and Yangjin Kim
"""

import os as os
import sys
import numpy as np   

myPyDSToolfile='/Users/aurelio/Documents/python/robclewley-pydstool-9ecaefa/PyDSTool/setup.py'

path = os.path.abspath(os.path.join(os.path.dirname(myPyDSToolfile), '..'))
if not path in sys.path:
    sys.path.insert(1, path)
del path

import PyDSTool

from PyDSTool import *
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt


plt.close('all')
plt.figure(figsize=(10,7))

mylinestyle = ['--','-.','-']


Svalue =[0.05,0.2,0.356592] 
for i in range(len(Svalue)):   
    Sval = Svalue[i]
    #mylabel = r"$S$ = %.2f" %S
    
    lamda = Par(0.01, 'lamda')
    G = Par(0.45, 'G')
    k1 = Par(4.0, 'k1')
    k3 = Par(1.0, 'k3')
    alpha = Par(1.5, 'alpha')
    k2 = Par(4.0, 'k2')
    k4 = Par(1.0, 'k4')
    beta = Par(1.0, 'beta')
    S = Par(Sval, 'S')
    mu = Par(1.0, 'mu')
    r = Par(0.05, 'r')
    K = Par(1.0, 'K')
    gam1 = Par(0.1,'gam1')
    TK = Par(100,'TK')
    lamdaG = Par(1.0, 'lamdaG')
    lamdaS = Par(1.0, 'lamdaS')
    decayT = Par(0.005, 'decayT')
    
    
    ## Compute nontrivial boundary equilibrium initial condition from parameters (see reference)
    C0 = 1.0
    I0 = 2.5
    T0 = 0.2
    
    
    # Declare symbolic variables
    C = Var('C')
    I = Var('I')
    T = Var('T')
    
    
    # Create Symbolic Quantity objects for definitions
    Crhs = lamda + lamdaG*G + (k1/(k3**2 + alpha*I**2)) - C
    Irhs = lamdaS*S + (k2/(k4**2 + beta*C**2)) - mu*I
    Trhs = r*(1 + C/(K+gam1*I))*T*(1-T/TK) - decayT*I*T
    
    
    # Build Generator
    DSargs = args(name='N1N2dynamics')
    DSargs.pars = [lamda,G,k1,k3,alpha,k2,k4,beta,S,mu,r,K,gam1,TK,lamdaG,\
                   lamdaS,decayT]
    DSargs.varspecs = args(C=Crhs,I=Irhs,T=Trhs)
    
    
    # Use eval method to get a float value from the symbolic definitions given in
    # terms of parameter values
    DSargs.ics = args(C=C0, I=I0, T=T0)
    ode = Generator.Vode_ODEsystem(DSargs)
    
    DSargs.tdomain = [0,500]                             # set the range of integration.
    ode1  = PyDSTool.Generator.Vode_ODEsystem(DSargs)    # an instance of the 'Generator' class.
    
    traj = ode1.compute('polarization')                  # 
    pd   = traj.sample()    
    
    
    #plt.close('all')  
    
    # Prepare the system to start close to a steady state
    ode.set(pars = {'G': 0.45} )       # Lower bound of the control parameter 'i'
    ode.set(ics =  {'C': pd['C'][0]} )       # Close to one of the steady states present for G=0.45
    
    PyCont = PyDSTool.ContClass(ode)                 # Set up continuation class
    
    PCargs = PyDSTool.args(name='EQ1', type='EP-C')  # 'EP-C' stands for Equilibrium Point Curve. The branch will be labeled 'EQ1'.
    PCargs.freepars     = ['G']                      # control parameter(s) (it should be among those specified in DSargs.pars)
    
    PCargs.MaxNumPoints = 300                       # The following 3 parameters are set after trial-and-error
    PCargs.MaxStepSize  = 0.01
    PCargs.MinStepSize  = 1e-5
    PCargs.StepSize     = 2e-2
    #PCargs.LocBifPoints = 'LP'                       # detect limit points / saddle-node bifurcations
    PCargs.SaveEigen    = True                       # to tell unstable from stable branches
    
    PyCont.newCurve(PCargs)
    PyCont['EQ1'].forward()
    PyCont['EQ1'].backward()
    
    PyCont.display(['G','C'],color='blue',linewidth=3,linestyle=mylinestyle[i],\
                   label=r"$S$ = %.2f" %Svalue[i])        # stable and unstable branches as solid and dashed curves, resp.
    
Cth = 1.81
Ith = 1.29


plt.axhline(Cth,color='blue',linestyle=':',linewidth=2)
plt.text(0.85,1.9, r"$C_{\rm th}\approx 1.81$", color='blue',fontsize='16')

plt.xticks(np.arange(0,1.2,0.2),[0,0.2,0.4,0.6,0.8,1.0],fontsize=22)
plt.yticks(fontsize=22)
plt.xlim([0,1.0])
plt.ylim([0,5])
#plt.scatter(0.658443,1.808470,s=40,marker='h',c='red')
plt.title("")    
plt.ylabel(r'N2 TANs $(C)$',fontsize=24)
plt.xlabel(r'TGF-$\beta  (G)$',fontsize=24)
plt.legend(ncol=3,prop={'size': 16})
plt.tight_layout()

#save figure in "BifurcationPlots" folder
workingdir = os.getcwd() # accessing current directory
resultsdir = 'BifurcationPlots' 

if not os.path.exists(resultsdir): # making a folder named 'plots'
    os.makedirs(resultsdir) 

os.chdir(resultsdir)
    
plotdir = os.getcwd()


plt.savefig('N1N2bifurcation_varyGvaryS.tif',dpi=300, pad_inches=0)
    
os.chdir(workingdir)
    
