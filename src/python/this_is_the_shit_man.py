#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 18:31:08 2021

@author: os18o068
"""

import numpy as np
from matplotlib import pyplot as plt
from PIMPphysicalparams import m_earth, r_earth, mFe, mMg, mO, mSi

d1 = 10e3
d2 = 6e3


def R1(C):
    return ((C*m_earth*r_earth**2-r_earth**5*d2)/(d1-d2))**(1./5.)

def delta_R1(C, dC):
    return dC*C*m_earth*r_earth**2/(5*(d1-d2))*R1(C)**(-4.)

def delta_Mg(C, dC):
    r1 = R1(C)
    dr1 = delta_R1(C, dC)
    
    dMgFe = -3*d2/d1*r_earth**3/r1**4*r1/2
    
    MgFe = d2/d1*(r_earth**3-r1**3)/r1**3*mFe/(mMg+mSi+3*mO)
    print (MgFe/(1+MgFe))
    return dMgFe*(1/(1+MgFe) - MgFe/(1+MgFe)**2)

def plot():
    x = np.linspace(.3, .34)
    y = np.linspace(.0, .05)
    
    z1 = np.empty([len(x), len(y)])
    z2 = np.empty([len(x), len(y)])
    
    for i in range(len(x)):
        for j in range(len(y)):
            z1[i][j] = delta_Mg(x[i], y[j])
            z2[i][j] = delta_R1(x[i], y[j])
        
    plt.imshow(abs(z1.T), origin='lower')
    
    fig, ax = plt.subplots()
    
    ax.plot(x, abs(z1.T[0]))
    ax.set_xlabel(r'$C/MR^2$')
    ax.set_ylabel(r'$\delta_{\rm max}\rm Mg\#$')
    #plt.imshow(abs(z2.T), origin='lower')    
    
    