#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 21:02:11 2020

@author: os18o068
"""

import numpy as np
from PIMPphysicalparams import m_earth, r_earth
import pandas as pd
import PlanetLight
import os
from PIMPrunparams import color_list
import shutil
import sys
import time
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath, amssymb}"]

file_dir = '/home/os18o068/Documents/PHD/Abbildungen/Figures_paper_1/'

linestyles = ('-', '--', ':')
mark = 's'
fnts1=9
fnts2 = 6
fnts3 = 5
lwdth = .5

legendalpha = .5

'''
planet_dir = 'ocean_planets'


count = 0
#gather data in ocean directory
for root, dirs, files in os.walk(planet_dir):
    
    for file in files:
        count += 1

data = np.empty([7, count])

count = 0 
for root, dirs, files in os.walk(planet_dir):
    
    for file in files:
        pl = PlanetLight.Planet()
        pl.load(loc = planet_dir+'/', file_name='planet_'+str(count))
        data[0][count] = pl.finals['M_surface_is']
        data[1][count] = pl.finals['R_surface_is']
        data[2][count] = (pl.finals['layer_properties'][3]['indigenous_mass']+\
                            pl.finals['layer_properties'][2]['indigenous_mass'] +\
                            pl.finals['layer_properties'][1]['indigenous_mass'] +\
                            pl.finals['layer_properties'][0]['indigenous_mass'])/m_earth
                            
        data[3][count] = pl.finals['layer_properties'][3]['R_outer']/r_earth
        data[4][count] = pl.finals['layer_properties'][4]['indigenous_mass']/m_earth
        data[5][count] = pl.finals['layer_properties'][4]['R_outer']/r_earth
        data[6][count] = pl.finals['Mg_number_is']

        
        count += 1

np.save('ocean_planet_outputs.npy', data)

print ('count =', count)

'''

data = np.load('ocean_planet_outputs.npy')

fig, ax = plt.subplots(1,2)

ax[0].scatter(data[0], data[2]/(4/3*np.pi*data[3]**3)*m_earth/r_earth**3, c=data[-1], cmap='rainbow')
ax[1].scatter(data[0], data[4]/(4/3*np.pi*(data[5]**3-data[3]**3))*m_earth/r_earth**3, c=data[-1], cmap='rainbow')



ax[0].set_xlim(.1, 10.)
#ax[0].set_ylim(0., 1.)

ax[1].set_xlim(.1, 10.)
#ax[1].set_ylim(0., 1.)


plt.save('ocean_estimate.pdf')
plt.close('all')