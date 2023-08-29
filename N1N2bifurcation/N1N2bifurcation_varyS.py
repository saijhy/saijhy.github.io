#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This program creates a bifurcation diagram varying the parameter S.

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


lamda = Par(0.01, 'lamda')
G = Par(0.45, 'G')
k1 = Par(4.0, 'k1')
k3 = Par(1.0, 'k3')
alpha = Par(1.5, 'alpha')
k2 = Par(4.0, 'k2')
k4 = Par(1.0, 'k4')
beta = Par(1.0, 'beta')
S = Par(0.2, 'S')
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

DSargs.tdomain = [0,100]                             # set the range of integration.
ode1  = PyDSTool.Generator.Vode_ODEsystem(DSargs)    # an instance of the 'Generator' class.

traj = ode1.compute('polarization')                  # 
pd   = traj.sample()    


#plt.close('all')  

# Prepare the system to start close to a steady state
ode.set(pars = {'S': 0.2} )       # Lower bound of the control parameter 'i'
ode.set(ics =  {'C': pd['C'][-1]} )       # Close to one of the steady states present for G=0.45

PyCont = PyDSTool.ContClass(ode)                 # Set up continuation class

PCargs = PyDSTool.args(name='EQ1', type='EP-C')  # 'EP-C' stands for Equilibrium Point Curve. The branch will be labeled 'EQ1'.
PCargs.freepars     = ['S']                      # control parameter(s) (it should be among those specified in DSargs.pars)

PCargs.MaxNumPoints = 300                       # The following 3 parameters are set after trial-and-error
PCargs.MaxStepSize  = 0.01
PCargs.MinStepSize  = 1e-5
PCargs.StepSize     = 2e-2
#PCargs.LocBifPoints = 'LP'                       # detect limit points / saddle-node bifurcations
PCargs.SaveEigen    = True                       # to tell unstable from stable branches

PyCont.newCurve(PCargs)
PyCont['EQ1'].forward()
PyCont['EQ1'].backward()

plt.close('all')
plt.figure(figsize=(10,7))

PyCont.display(['S','C'],color='blue',linewidth=3,label=r'$C$')        # stable and unstable branches as solid and dashed curves, resp.
PyCont.display(['S','I'],color='magenta',linewidth=3,label=r'$I$',linestyle='-.') 

Cth = 1.81
Ith = 1.29

#PyCont['EQ1'].getSpecialPoint('LP1')
LP1=0.05992734846244219
LP2=0.24965629260163663


plt.axvline(x=LP1,linestyle='--',color='k')
plt.axvline(x=LP2,linestyle='--',color='k')
plt.text(0.145,4.5, r"W$^b$", color='k',fontsize='22')
plt.arrow(0.143,4.58, -0.07, 0, head_width=0.05, head_length=0.01, fc='k', ec='k')
plt.arrow(0.168,4.58, 0.07, 0, head_width=0.05, head_length=0.01, fc='k', ec='k')

#arrows in (G-C) hysteresis curve
plt.annotate('', xy=(LP1+0.01,0.65), xytext=(LP2-0.01,0.55),
            arrowprops={'arrowstyle': '-|>', 'lw': 3, 'color': 'blue'},
            va='center')
plt.annotate('', xy=(LP2-0.01,3.2), xytext=(LP1+0.01,4.25),
            arrowprops={'arrowstyle': '-|>', 'lw': 3, 'color': 'blue'},
            va='center')
plt.annotate('', xy=(LP2+0.01,0.9), xytext=(LP2+0.01,2.4),
            arrowprops={'arrowstyle': '-|>', 'lw': 3, 'color': 'blue'},
            va='center')
plt.annotate('', xy=(LP1-0.01,3.9), xytext=(LP1-0.01,1.07),
            arrowprops={'arrowstyle': '-|>', 'lw': 3, 'color': 'blue'},
            va='center')

#arrows in (G-I) hysteresis curve
plt.annotate('', xy=(LP2-0.01,0.35), xytext=(LP1+0.01,0.05),
            arrowprops={'arrowstyle': '-|>', 'lw': 3, 'ls':'-.','color': 'magenta'},
            va='center')
plt.annotate('', xy=(LP1+0.01,2.42), xytext=(LP2-0.01,2.92),
            arrowprops={'arrowstyle': '-|>', 'lw': 3, 'ls':'-.','color': 'magenta'},
            va='center')
plt.annotate('', xy=(LP2+0.02,2.6), xytext=(LP2+0.02,0.9),
            arrowprops={'arrowstyle': '-|>', 'lw': 3, 'ls':'-.','color': 'magenta'},
            va='center')
plt.annotate('', xy=(LP1-0.02,0.35), xytext=(LP1-0.02,1.92),
            arrowprops={'arrowstyle': '-|>', 'lw': 3, 'ls':'-.','color': 'magenta'},
            va='center')


plt.axhline(Cth,color='blue',linestyle=':',linewidth=2)
plt.text(0.32,1.9, r"$C_{\rm th}\approx 1.81$", color='blue',fontsize='20')
plt.axhline(Ith,color='magenta',linestyle=':',linewidth=2)
plt.text(0.325,1.38, r"$I_{\rm th}\approx 1.29$", color='magenta',fontsize='20')

plt.xticks(np.arange(0,0.5,0.1),[0,0.1,0.2,0.3,0.4,0.5])
plt.xlim([0,0.4])
plt.ylim([0,5])
plt.title("")    
plt.ylabel('concentration',fontsize=28)
plt.xlabel(r'IFN-$\beta (S)$',fontsize=28)
plt.legend(prop={'size': 24})
plt.tick_params(axis='both',labelsize=26)
plt.subplots_adjust(bottom=0.15)

#save figure in "BifurcationPlots" folder
workingdir = os.getcwd() # accessing current directory
resultsdir = 'BifurcationPlots' 

if not os.path.exists(resultsdir): # making a folder named 'plots'
    os.makedirs(resultsdir) 

os.chdir(resultsdir)
    
plotdir = os.getcwd()


plt.savefig('N1N2bifurcation_varyS.tif',dpi=300, pad_inches=0)
    
os.chdir(workingdir)
    
