#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 12:36:25 2021

@author: os18o068
"""

import numpy as np
from matplotlib import pyplot as plt
from PIMPphysicalparams import M_solar, R_solar, MoI_solar, J2_solar, ecc_solar,\
    orbit_solar, axes_solar, hosts_solar, G, m_earth, rotation_solar
import matplotlib

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

def plot_orbits():
    
    objects = ['moon', 'titan', 'callisto', 'io', 'ganymede', 'europa']
    
    ecc = []
    host_mass = []
    distances = []
    fracs = []
    synch = []    
    
    for key in objects:
        ecc.append(ecc_solar[key])
        host_mass.append(M_solar[hosts_solar[key]])
        distances.append(axes_solar[key])
        fracs.append(M_solar[hosts_solar[key]]/(axes_solar[key])**2/M_solar['sun']*(axes_solar[hosts_solar[key]])**2)
        synch.append(rotation_solar[key]/orbit_solar[key])
        
    lists = sorted(M_solar.items())
    
    x, y = zip(*lists)
    
    print (x, y)
    
    fig, ax = plt.subplots(2,2)
    
    ax[0][0].set_ylabel(r'$\rm Eccentricity$')
    ax[0][1].set_ylabel(r'$g_{\rm host}/g_{\rm sun}$')
    ax[1][0].set_ylabel(r'$P_{\rm Rot}/P_{\rm Orb}$')
    
    ax[0][0].scatter(objects, ecc)
    ax[0][1].scatter(objects, fracs)
    ax[1][0].scatter(objects, synch)

    ax[0][1].plot([-1, len(objects)+1], [1., 1.])
    ax[0][1].set_xlim(-.5, len(objects)-.5)
    
    ax[0][0].set_yscale('log')
    ax[0][1].set_yscale('log')    
    
    