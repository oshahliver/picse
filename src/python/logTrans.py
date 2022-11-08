#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 13:20:35 2019

@author: oshah
"""

import numpy as np
import math
from matplotlib import pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable


def transform(data=None, param=2, def_zero=1.0e-3, dim=2):
    """Takes data as input and transforms it to logarithmic taking into account
    possible sign inversions which would make it impossible to just directly
    take log(data).
    """
    
    if dim == 1:
        nx = len(data)
        data_dummy = np.zeros([nx])
        data_dummy_dummy = np.zeros([nx])
        sign_inversion = False
        
        #gather data in 1d array
        for i in range(nx):
            data_dummy[i] = data[i][param]
            data_dummy_dummy[i] = data[i][param]
        
        dat_min = np.nanmin(data_dummy)
        dat_max = np.nanmax(data_dummy)
        
        print ('dat_min=', dat_min)
        print ('dat_max=', dat_max)
        
        #check for sign inversion
        for i in range(nx):
            if not math.isnan(data_dummy[i]):
                if data_dummy[0]/abs(data_dummy[0]) ==\
                data_dummy[i]/abs(data_dummy[i]):
                    pass
                
                else:
                    sign_inversion = True
            
            else:
                pass
            
        if sign_inversion:
            for i in range(nx):
                dat = data_dummy[i]
              
                #work on data set to determine min and max exponent
                if abs(dat) <= def_zero:
                    data_dummy_dummy[i] = np.log10(def_zero)
                    
                else:
                    data_dummy_dummy[i] = np.log10(abs(dat))
            
            exp_min = np.nanmin(data_dummy_dummy)
            exp_max = np.nanmax(data_dummy_dummy)
                        
        else:
            for i in range(nx):
                dat = data_dummy[i]
                #if value is below defzero keep it zero in the dummy array
                if abs(dat) <= def_zero:
                    data_dummy[i] = np.log10(def_zero)
                    
                else:
                    data_dummy[i] = np.log10(abs(dat))  
    
            exp_min = np.nanmin(data_dummy)
            exp_max = np.nanmax(data_dummy)
    
        print ('exp_max=', exp_max)
        print ('exp_min=', exp_min)
        
        #transform data to exp>=1 if
        if sign_inversion:
            for i in range(nx):
                dat = data_dummy[i]
                if not dat == 0.:
                    if abs(dat) <= def_zero:
                        data_dummy[i] = def_zero*dat/abs(dat)
                        dat = data_dummy[i]
                    
                    data_dummy[i] = (np.log10(abs(dat)) - \
                              int(exp_min))*dat/abs(dat)
            
                else:
                    pass        
    
    if dim == 2:
        nx = len(data)
        ny = len(data[0])
        data_dummy = np.zeros([nx, ny])
        data_dummy_dummy = np.zeros([nx, ny])
        sign_inversion = False
        
        #gather data in 2d array
        for i in range(nx):
            for j in range(ny):
                data_dummy[i][j] = data[i][j][param]
                data_dummy_dummy[i][j] = data[i][j][param]
        
        dat_min = np.nanmin(data_dummy)
        dat_max = np.nanmax(data_dummy)
        
        print ('dat_min=', dat_min)
        print ('dat_max=', dat_max)
        
        #check for sign inversion
        for i in range(nx):
            for j in range(ny):
                
                if not math.isnan(data_dummy[i][j]):
                    if data_dummy[0][0]/abs(data_dummy[0][0]) ==\
                    data_dummy[i][j]/abs(data_dummy[i][j]):
                        pass
                    
                    else:
                        sign_inversion = True
                
                else:
                    pass
            
        if sign_inversion:
            for i in range(nx):
                for j in range(ny):
                    dat = data_dummy[i][j]
                  
                    #work on data set to determine min and max exponent
                    if abs(dat) <= def_zero:
                        data_dummy_dummy[i][j] = np.log10(def_zero)
                        
                    else:
                        data_dummy_dummy[i][j] = np.log10(abs(dat))
            
            exp_min = np.nanmin(data_dummy_dummy)
            exp_max = np.nanmax(data_dummy_dummy)
                        
        else:
            for i in range(nx):
                for j in range(ny):
                    dat = data_dummy[i][j]
                    #if value is below defzero keep it zero in the dummy array
                    if abs(dat) <= def_zero:
                        data_dummy[i][j] = np.log10(def_zero)
                        
                    else:
                        data_dummy[i][j] = np.log10(abs(dat))  
    
            exp_min = np.nanmin(data_dummy)
            exp_max = np.nanmax(data_dummy)
    
        print ('exp_max=', exp_max)
        print ('exp_min=', exp_min)
        
        #transform data to exp>=1 if
        if sign_inversion:
            for i in range(nx):
                for j in range(ny):
                    dat = data_dummy[i][j]
                    if not dat == 0.:
                        if abs(dat) <= def_zero:
                            data_dummy[i][j] = def_zero*dat/abs(dat)
                            dat = data_dummy[i][j]
                        
                        data_dummy[i][j] = (np.log10(abs(dat)) - \
                                  int(exp_min))*dat/abs(dat)
                
                    else:
                        pass
                
    return data_dummy, exp_min, exp_max, sign_inversion


def plot(data, exp_min, exp_max, sign_inversion, cmap='rainbow', 
         symmetric_axis=False, major_ticks=[], major_tick_labels=[], 
         cbar_label='', material=None):
    
    #print ('plot data=', data)
    
    if sign_inversion:
        exp_max_plot = round(np.nanmax(data)+.5, 0)
        exp_min_plot = round(np.nanmin(data)-.5, 0)
        
        if symmetric_axis:
            if abs(exp_max_plot) > abs(exp_min_plot):
                exp_min_plot = -exp_max_plot
                    
            elif abs(exp_max_plot) < abs(exp_min_plot):
                exp_max_plot = -exp_min_plot
        
        n_ticks = int(exp_max_plot - exp_min_plot) + 1
  
    else:
        exp_max_plot = round(exp_max+.5, 0)
        exp_min_plot = round(exp_min-.5, 0)
        n_ticks = int(exp_max_plot - exp_min_plot)+1
    
    #print ('n_ticks=', n_ticks)
    #print ('exp_max_plot=', exp_max_plot)
    #print ('exp_min_plot=', exp_min_plot)
    ticks = np.linspace(exp_min_plot, exp_max_plot, n_ticks)

    norm = mpl.colors.Normalize(vmin=exp_min_plot, vmax=exp_max_plot)
    
    tick_labels = []
    for t in ticks:
        if sign_inversion:
            if not t == 0:
                if abs(t)/t < 0:
                    sign = '-'
                else:
                    sign = ''
                    
                tick_labels.append(sign+'10e'+str(int(abs(t))+int(exp_min)))
                
            else:
                tick_labels.append('0')

        else:
            tick_labels.append('10e'+str(t))

    fig, ax = plt.subplots()
    im = ax.imshow(data.T, cmap=cmap, extent=[0.,1.,0.,1.],
            norm=norm, origin='lower')              

    divider=make_axes_locatable(ax)
    cax=divider.append_axes('right', size='5%', pad=.5)
    cbar=fig.colorbar(im, cax=cax, cmap=cmap, norm=norm, ticks=ticks)
    
    cbar.ax.tick_params(size=0.)
    cbar.ax.set_yticklabels(tick_labels)
    cbar.ax.set_ylabel(cbar_label)
    
    #set major ticks and labels for x and y axis
    if not len(major_ticks) == 0 and not len(major_tick_labels) == 0:
        ax.set_xticks(major_ticks[0])
        ax.set_yticks(major_ticks[1])
        
        ax.set_xticklabels(major_tick_labels[0])
        ax.set_yticklabels(major_tick_labels[1])

        ax.set_xlabel(r'$log(T) \ [K]$')
        ax.set_ylabel(r'$log(P) \ [Pa]$')

        ax.set_title(material)
        
        ax.tick_params(right=True, top=True, labeltop=False, labelright=False)