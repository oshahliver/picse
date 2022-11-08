#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 15:08:57 2019

@author: oshah
contact: oliver.shah@space.unibe.ch

This script reads in ascii data files generated with Planet.write(). The
simplest method for usage is shown in the example below for reading in
a data file yourFile located in yourPath:
    
    >import readPlanet as rd
    
    #loads data into multidimensional numpy array for further processing
    >data = rd.read(yourPath/yourFile)
    
    #do some stuff with the data
    #...

Alternatively you can directly plot the data:
    
    >import readPlanet as rd
    
    #visualize the data in 2x3 plots
    >rd.plot(yourPath/yourFile)


NOTE: The <astropy> package is required to use this routine. It can conveniently
be installed from the terminal via:
    
    > pip install astropy

"""

import numpy as np
from astropy.io import ascii
from matplotlib import pyplot as plt

def read(dataFile):
    """Read in standard data file generated with Planet.write() method.
    
        By default, the data is organized as follows:
                
    R           M           T       P       rho         v_esc       g
    (r_earth)   (m_earth)   (K)     (GPa)   (kg/m3)     (km/s)      (m/s2)
    -----------------------------------------------------------------------        
    R_center    M_center    ...
    ...         ...         ...
    R_surface   M_surface   ...
    
    Note:   r_earth = 6.371e6 m
            m_earth = 5.9722e24 kg

    """
    
    print ('\nReading in data file...', dataFile)
    dataTable = ascii.read(dataFile)
    print ('--> done')
    
    print ('\nPlanet consists of', len(dataTable), 'shells')
    
    #display contents
    print ('\nExtracted data table contents:')
    print ('\n', dataTable)
    
    dataArray = np.zeros([len(dataTable[0]), len(dataTable)])
    
    #gather data for plots
    print ('\nGathering data in numpy array...')
    #loop over all shells of the planet
    for p in range(len(dataArray[0])):
        for d in range(len(dataArray)):
            #note that the order of the indices must be inverted due to how
            #astropy saves the data into the ascii file
            dataArray[d][p] = dataTable[p][d]
            
    print ('--> done')

    print ('\nThe data is now ready for processing')
    print ('\nAccess parameter via the first index and the shell via the second')
    print ('\nRemember that python starts indexing at 0!')
    print ('''\nExample 1: dataArray[2][67] will give T (K) in the 68th shell
           counted from the center''')
    
    print ('''\nExample 2: dataArray[0][-1] will give R (R_earth) at the 
           surface of the planet (last shell)''')

    return  dataArray


def plot(dataFile):
    """Read in data file and plot the data
    """
    #read in data file
    plotData = read(dataFile)

    #generate plot labels for the different parameters
    labels = [r'$R \ [R_{\oplus}]$', 
              r'$M \ [M_{\oplus}]$', 
              r'$T \ [K]$', 
              r'$P \ [GPa]$',
              r'$\rho \ [kg/m^3]$', 
              r'$v_{esc} \ [km/s]$', 
              r'$g \ [m/s^2]$'
              ]
    
    #define plot colors
    pres_color = (.4, .5, .7)
    dens_color = (.8, .4, .0)
    mass_color = (.2, .8, .2)
    temp_color = (.8, .2, .2)
    vesc_color = (.5, .5, .5)
    grav_color = (.5, .5, .0)
    param_colors = [mass_color,
                    temp_color,
                    pres_color,
                    dens_color,
                    vesc_color, 
                    grav_color]
    
    #plot the data
    fig, ax = plt.subplots(2, 3)
    fig.tight_layout()
    fig.subplots_adjust(wspace=.5)
    
    for i in range(2):
        for j in range(3):
            #index goes from 1 to 6 to account for all parameters (y-axis)
            #index 0 corresponds to radius (x-axis)
            ind = (i+1)*(j+1) + (2-j)*i
            ax[i][j].plot(plotData[0], plotData[ind], color=param_colors[ind-1])
            ax[i][j].set_ylabel(labels[ind])
            ax[i][j].tick_params(top=True, right=True, direction='in')
            ax[i][j].grid(color = (.9, .9, .9))
            
            if i == 1:
                ax[i][j].set_xlabel(labels[0])
    
    return plotData