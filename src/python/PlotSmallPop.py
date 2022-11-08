#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 21:02:11 2020

@author: os18o068
"""

import numpy as np
from PIMPphysicalparams import m_earth, r_earth
import pandas as pd
import PlanetFort
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

def delta_T_CMB(M):
    return 1400.*M**.75

linestyles = ('-', '--', ':')
mark = 's'
fnts1=9
fnts2 = 6
fnts3 = 5
lwdth = .5

legendalpha = .5

temp_label_color = 'grey'

widths = [4., 4., 4., .5]
heights = [1., 1., 1., 1.]

axes_ylims = [[.0, 20.], [.0, 1600.], [0., 20.]]

temp_labels = [r'${\Delta T_{\scriptscriptstyle \rm TBL} = 200 \ \rm K}$',
                r'${\Delta T_{\rm TBL} = 700 \ \rm K}$',
                r'${\Delta T_{\rm TBL} = 1200 \ \rm K}$',
                r'${\Delta T_{\rm TBL} = 1700 \ \rm K}$']

mass_labels = [r'$M/M_\oplus = 0.1$',
               r'$M/M_\oplus = 0.6$',
               r'$M/M_\oplus = 3$']

format = 'pdf'
file_dir = '/home/os18o068/Documents/PHD/Abbildungen/Figures_paper_1/'

Mg_numbers = np.linspace(.2, .7, 6)
masses = np.array([.1, .6, 3.])
temps = np.linspace(500., 2000., 4)

planets = [[], [], []]

axes_lims = np.zeros([2, 16, 2])

for k in range(len(masses)):
    m = masses[k]
    ax = []
    
    label_list1 = [r'$\rm Mg\# = \ $'+str(round(i*.1+.2, 1)) for i in range(6)]
    label_list2 = ['Type-1', 'Type-2', 'Type-3']
    
    #fig, ax = plt.subplots(len(temps), 4, sharex = True)
    fig = plt.figure()
    fig.subplots_adjust(hspace=.2, wspace=.5)
    
    #plt.title(r'$M/M_\oplus = $'+str(round(masses[k], 3)), pad=10)
    #plt.axis('off')
    
    gs = fig.add_gridspec(ncols=4, nrows=4, width_ratios=widths,
                          height_ratios=heights)

    for i in range(len(temps)):
        t = temps[i]
        ax.append([])
        
        for row in range(4):
            for col in range(4):
                ax[i].append(fig.add_subplot(gs[row,col]))     
        
        ax[i][3+4*i].text(-2., .5*heights[-1], 
                          temp_labels[i], 
                          rotation=90, 
                          fontsize=fnts2, 
                          verticalalignment='center',
                          color = temp_label_color)
          
        ax[i][3+4*i].axis('off')
        
        ax[i][0+4*i].set_ylabel(r'$\rm Temperature \ [10^{3} \ K]$', 
                                fontsize = fnts2)
        ax[i][1+4*i].set_ylabel(r'$\rm Pressure \ [GPa]$', 
                                fontsize = fnts2)
        ax[i][2+4*i].set_ylabel(r'$\rm Density \ \rm [gcc]$',
                                fontsize = fnts2)
        
    
        if i == len(temps)-1:
            for ii in range(3):
                ax[i][12+ii].set_xlabel(r'$R/R_\oplus$', fontsize = fnts2)
        
        plot_list1 = []
        plot_list2 = []
        
        for j in range(len(Mg_numbers)):
            mg = Mg_numbers[j]    
            
            planets = []
            loc = []
            dfs = []
            
            #loop over types
            for s in range(3):
                
                pl = PlanetFort.Planet()
                
                loc = './M_'+str(round(m, 2)).replace('.', '_')+'/Tsurf_'
                loc = loc + str(int(t))+'/Type_'+str(s+1)+'/'
    
                pl.load(loc=loc, file_name = 'planet_'+str(j))
                
                if s == 1:
                    pl.correct_density_profiles()
                
                df = pl.data_table.to_pandas()
                
                #For type 1 and type 2 planets, the outer temperature is
                #the temperature before the TBL. Assuming T_surf of 300
                #we must substract (T_TBl -300) from that value in order to
                #get 300 K at the surfaces of the planets
                if s < 2:
                    df['T (K)'][len(df['T (K)'])-1] -= t - 300.
                
                pl1, = ax[i][0+4*i].plot(df['R (r_earth)']/r_earth, df['T (K)']/1000,
                              linestyle = linestyles[s],
                              color=color_list[j],
                              linewidth=lwdth)
                
                
                pl2, = ax[i][1+4*i].plot(df['R (r_earth)']/r_earth, df['P (GPa)']*1.0e-9,
                              linestyle = linestyles[s],
                              color=color_list[j],
                              linewidth=lwdth)
                
                
                ax[i][2+4*i].plot(df['R (r_earth)']/r_earth, df['rho (kg/m3)']/1000,
                              linestyle = linestyles[s],
                              color=color_list[j],
                              linewidth=lwdth)
                
                if s == 0:
                    plot_list1.append(pl1)
                    
                if j == 0:
                    plot_list2.append(pl2)
                
                for ii in range(4):
                    ax[i][ii+4*i].tick_params(direction='in', top=True, right=True,
                                     which='both', labelsize=fnts3)
                    
    legend1 = ax[0][12].legend(plot_list1, label_list1, fontsize=fnts2,
                              loc=2, framealpha = legendalpha,
                              bbox_to_anchor=(0., -.5),
                              borderaxespad=0.,
                              edgecolor='None')
    
    legend2 = ax[0][12].legend(plot_list2, label_list2, fontsize=fnts2,
                              loc=2, framealpha = legendalpha,
                              bbox_to_anchor=(1., -.5),
                              edgecolor='None')
    
    ax[0][12].add_artist(legend1)
    ax[0][12].add_artist(legend2)
    
    bbox_props = dict(edgecolor='k', facecolor='white',
                      linewidth = .5)
    
    ax[0][3].text(-30, 1.1, mass_labels[k],
                  bbox=bbox_props, va='center')
    
    fig.align_labels()
    fig.savefig(file_dir+'profiles_'+str(round(m,2)).replace('.', '_')+'.'+format,
                bbox_inches='tight')
                
plt.close('all')