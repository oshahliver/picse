#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 19:12:30 2021

@author: os18o068
"""


import pandas as pd
import numpy as np
from astropy.io import ascii
from matplotlib import pyplot as plt
import matplotlib as mpl
from tabulate import tabulate
from phase_transitions_water_Wagner2002 import T_critical, P_critical

#mpl.rc('text', usetex=True)
#mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath, amssymb}"]


def filter(M, T, **kwargs):
    a, b, c, d = .1468, .2694, .05546, -.03538
    return 10**(a + b * np.log10(M) + c * T / 1000 + d * np.log10(M) * T / 1000)


class archive():
    def __init__(self):
        self.data = ascii.read('/home/os18o068/Documents/PHD/Projects/Planets/Code/PS_2021.12.02_02.23.23.csv')
        self.df = self.data.to_pandas()
        self.planets = []
        self.names = []
        self.hostnames = []
        self.filtered_data = []
        self.distances = []
        self.T_max = T_critical
        
        i = 0
        while True:
            if np.isnan(self.df['pl_bmasse'].values[i]) or np.isnan(self.df['pl_rade'].values[i]) \
            or np.isnan(self.df['pl_radeerr1'][i]) or np.isnan(self.df['pl_radeerr2'][i]) \
            or pd.isnull(self.df['pl_eqt'][i]) or np.isnan(self.df['pl_bmasseerr1'][i]) \
            or np.isnan(self.df['pl_bmasseerr2'][i]):
                pass
            
            else:
                if self.df['pl_bmasse'].values[i] < 10.:
                    x1 = self.df['pl_bmasse'][i]
                    x2 = self.df['pl_bmasseerr1'][i]
                    x3 = self.df['pl_bmasseerr2'][i]
                    x4 = self.df['pl_rade'][i]
                    x5 = self.df['pl_radeerr1'][i]
                    x6 = self.df['pl_radeerr2'][i]
                    x7 = self.df['pl_eqt'][i]
                    x8 = self.df['pl_eqterr1'][i]
                    x9 = self.df['pl_eqterr2'][i]
                    x10 = self.df['pl_name'][i]
                    x11 = self.df['hostname'][i]
                    x12 = self.df['sy_dist'][i]
                    self.planets.append([x1, x2, x3, x4, x5, x6, x7, x8, x9])
                    self.names.append(x10)
                    self.hostnames.append(x11)
                    self.distances.append(x12)
            i += 1
            if i == self.df.shape[0] - 1:
                break        


    def plot_filter(self, ax = None, vmin = 0, vmax = 2500, cmap = 'jet'):
        x = np.linspace(.0001, 10, 100)
        t = np.linspace(vmin, vmax, 6)
        
        data = np.empty([len(t), len(x)])
        
        for i in range(len(t)):
            data[i] = filter(x, t[i])
        
        if ax == None:
            fig, ax = plt.subplots()
            
        cm = mpl.cm.get_cmap(cmap)
        for i in range(len(t)):
            c = (t[i] - vmin) / (vmax - vmin)
            ax.plot(x, data[i], color = cm(c), zorder = 0)
    
    
    def prt(self):
        dat = []
        for i in range(len(self.filtered_names)):
            dat.append([])
            dat[i].append(self.filtered_names[i])
            for j in range(len(self.filtered_data[i])):
                dat[i].append(self.filtered_data[i][j])
            
            dat[i].append(int(self.filtered_data[i][1] / self.filtered_data[i][0] * 100))
            dat[i].append(int(self.filtered_data[i][2] / self.filtered_data[i][0] * 100))
            
        tab = tabulate(dat, headers=['Name', 'Mass', '+err', '-err', 'Radius', 
                                     '+err', '-err', 'T_eq', '+err', '-err',
                                     '+%', '-%'])
           
        print()
        print (f"{tab}")
        print ()
    
    
    def filter(self, fct = filter, max_mass_err = .1, **kwargs):
        """Apply filter function to filter out unwanted data points
        """
        count = 0
        self.filtered_data = []
        self.filtered_names = []
        self.filtered_hosts = []
        self.filtered_distances = []
        for p in range(len(self.planets)):
            pl = self.planets[p]
            
            #extract possible mass and temp range for planet
            masses = np.array(pl[1:3]) + pl[0]
            temps = np.array(pl[7:9]) + pl[6]
            
            #check if pure water sphere would be possible in this range
            diffs = []
            for i in range(2):
                for j in range(2):
                    rad = fct(masses[i], temps[j], **kwargs)
                    for k in range(2):
                        diffs.append(min(pl[3] + pl[4 + k] - rad, 0))
            
            #check if mass error is below set limit
            errs = [abs(pl[2] / pl[0]), abs(pl[1] / pl[0])]
            
            if sum(diffs) < 0 and temps[0] < self.T_max and masses[0] < 10. \
                and errs[0] < max_mass_err and errs[1] < max_mass_err \
                and not self.names[p] == 'HD 219134 f':
                count += 1
                self.filtered_data.append(pl)
                self.filtered_names.append(self.names[p])
                self.filtered_hosts.append(self.hostnames[p])
                self.filtered_distances.append(self.distances[p])
                
        self.filtered_data = np.array(self.filtered_data)
        print ('count = ', count)
        
        
    def plot_filtered(self, cmap = 'jet', ax = None, annot = False,
                     vmin = 0, vmax = 2500, fnts = 16, n_ticks = 6):
        #self.filter()
        
        if ax == None:
            fig, ax = plt.subplots()
        self.plot_filter(ax = ax, vmin = vmin, vmax = vmax, cmap = cmap)    
        ax.set_xlim(.1, 10)
        ax.set_ylim(.5, 4.)
        
        yerr = [abs(self.filtered_data.T[4]), abs(self.filtered_data.T[5])]
        xerr = [abs(self.filtered_data.T[1]), abs(self.filtered_data.T[2])]

        ax.errorbar(self.filtered_data.T[0], self.filtered_data.T[3], 
                    yerr =  yerr, xerr = xerr, ecolor = 'k', 
                    marker = '', linestyle = '', zorder = 1,
                    linewidth = 1)
    
        sc = ax.scatter(self.filtered_data.T[0], self.filtered_data.T[3], 
                        c = self.filtered_data.T[6], cmap = cmap, zorder = 2,
                        marker = 's', s = 30, vmin = vmin, vmax = vmax,
                        edgecolor = 'k')
        
        if annot:
            for i in range(len(self.filtered_names)):
                ax.annotate(self.filtered_names[i], (self.filtered_data.T[0][i], 
                                                 self.filtered_data.T[3][i]))
        
        inset_ax = ax.inset_axes((.05, 0.6, .05, .35))
        cbar = plt.colorbar(sc, cax = inset_ax, ticks=np.linspace(vmin, vmax, n_ticks))
        cbar.ax.set_ylabel(r'$T_{\rm eq} \ \rm [K]$', fontsize = fnts-2)
        cbar.ax.tick_params(labelsize = fnts-2)
        
        ticklabels = np.linspace(vmin, vmax, n_ticks).astype(int)
        cbar.ax.set_yticklabels(ticklabels)
        
        x = np.linspace(.1, 10, 100)
        m = np.linspace(0, 100, 3)
        
        '''
        for i in range(len(m)):
            y = MR_Grasset_2009(x, m[i])
            ax.plot(x, y, color = cols[i])
            ax.text(2., y[-1]-.75, str(int(m[i]))+r'% $\rm H_2 O$', color = cols[i])
        '''
        #ax.vlines(data['pl_bmasse'], df['pl_rade']- df['pl_radeerr2'],   df['pl_rade']+df['pl_radeerr1'])
        #ax.hlines(data['pl_rade'], data['pl_bmasse']- data['pl_bmasseerr1'],   data['pl_bmasse']+data['pl_bmasseerr2'])
        #ax.scatter(data['pl_bmasse'], data['pl_rade'])
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlabel(r'${\rm Mass} \ [M_\oplus]$', fontsize = fnts)
        ax.set_ylabel(r'${\rm Radius} \ [R_\oplus]$', fontsize = fnts)        
        ax.tick_params(labelsize = fnts)
        
        
    def plot(self, cmap = 'jet', ax = None, vmin = 0, vmax = 2500,
             fnts = 16):
        
        yerr = [[], []]
        xerr = [[], []]
        x_data = []
        y_data = []
        z_data = []        
        #yerr_low = np.nan_to_num(df['pl_radeerr2'].values, nan= 0., copy = True)
        #yerr_up = np.nan_to_num(df['pl_radeerr1'].values, nan= 0., copy = True)
        i = 0
        while True:
            if np.isnan(self.df['pl_bmasse'].values[i]) or np.isnan(self.df['pl_rade'].values[i]) \
            or np.isnan(self.df['pl_radeerr1'][i]) or np.isnan(self.df['pl_radeerr2'][i]) \
            or pd.isnull(self.df['pl_eqt'][i]) or np.isnan(self.df['pl_bmasseerr1'][i]) \
            or np.isnan(self.df['pl_bmasseerr2'][i]):
                pass
            
            else:
                if self.df['pl_bmasse'].values[i] < 10.:
                    x_data.append(self.df['pl_bmasse'].values[i])
                    y_data.append(self.df['pl_rade'].values[i])
                    z_data.append(self.df['pl_eqt'].values[i])
                    yerr[0].append(abs(self.df['pl_radeerr2'][i]))
                    yerr[1].append(self.df['pl_radeerr1'][i])
                    xerr[0].append(abs(self.df['pl_bmasseerr2'][i]))
                    xerr[1].append(self.df['pl_bmasseerr1'][i])
                   
            i += 1
            if i == self.df.shape[0] - 1:
                break
        
        if ax == None:
            fig, ax = plt.subplots()
        self.plot_filter(ax = ax, vmin = vmin, vmax = vmax, cmap = cmap)
        ax.set_xlim(.1, 10)
        ax.set_ylim(.5, 4.)
        ax.errorbar(x_data, y_data, yerr =  yerr, xerr = xerr, ecolor = 'k', 
                    marker = '', linestyle = '', zorder = 1,
                    linewidth = 1)
         
        sc = ax.scatter(x_data, y_data, c = z_data, cmap = cmap, zorder = 2,
                        marker = 's', s = 30, vmin = vmin, vmax = vmax,
                        edgecolor = 'k')
        
        inset_ax = ax.inset_axes((.05, 0.6, .05, .35))
        cbar = plt.colorbar(sc, cax = inset_ax, ticks=np.linspace(vmin, vmax, 6))
        cbar.ax.set_ylabel(r'$T_{\rm eq} \ \rm [K]$', fontsize = fnts-2)
        cbar.ax.tick_params(labelsize = fnts-2)
        
        ticklabels = np.linspace(vmin, vmax, 6).astype(int)
        cbar.ax.set_yticklabels(ticklabels)
        
        #ax.vlines(data['pl_bmasse'], df['pl_rade']- df['pl_radeerr2'],   df['pl_rade']+df['pl_radeerr1'])
        #ax.hlines(data['pl_rade'], data['pl_bmasse']- data['pl_bmasseerr1'],   data['pl_bmasse']+data['pl_bmasseerr2'])
        #ax.scatter(data['pl_bmasse'], data['pl_rade'])
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        ax.set_xlabel(r'${\rm Mass} \ [M_\oplus]$', fontsize = fnts)
        ax.set_ylabel(r'${\rm Radius} \ [R_\oplus]$', fontsize = fnts)
        ax.tick_params(labelsize = fnts)
    
    
def MR_Grasset_2009(M, X_w):
    """Empirical M-R relation for ocean planets up to 100 M_E from 
    Grasset et al. 2009.
    
    M: Planet mass in M_earth
    X_w: Water content in wt%
    """
    
    all_xis = np.array([[1.010, 2.859e-1, -5.518e-4, 2.096e-6], 
          [5.714e-3, -2.116e-4, -4.221e-6, 3.085e-8], 
          [-1.704e-5, -9.015e-7, 5.397e-8, -3.166e-10]
          ])
    
    def xi(x, xis):
        dummy = [xis[i] * x**i for i in range(3)]
        return sum(dummy)
    
    
    coefs = [xi(X_w, all_xis.T[i]) for i in range(4)]
    
    a, b, c, d = coefs
    
    return 10**(np.log10(a) + (b + c * M + d * M**2) * np.log10(M))


def plot_MR(N_mass = 20, N_oc = 3, oc = [0., 100.]):
    masses = np.linspace(1., 10., N_mass)
    waters = np.linspace(oc[0], oc[1], N_oc)
    
    data = np.empty([len(waters), len(masses)])
    fig, ax = plt.subplots()
    ax.set_xlabel(r'Mass [$M_\oplus$]')
    ax.set_ylabel(r'Radius [$R_\oplus$]')
    ax.set_xlim(1., 10.)
    ax.set_ylim(.75, 2.5)
    
    for i in range(len(waters)):
        for j in range(len(masses)):
            data[i][j] = MR_Grasset_2009(masses[j], waters[i])
    
    ticklabels = np.linspace(0, 100, 6)
    cmap = plt.cm.get_cmap('viridis')
    sm = plt.cm.ScalarMappable(cmap = cmap)
    cbar = plt.colorbar(sm, ticks = np.linspace(0, 1, 6))
    cbar.set_label(r'Water mass fraction [%]')
    cbar.ax.set_yticklabels(ticklabels.astype(int))
    
    c1 = 0
    c2 = 100
    for i in range(len(waters)):
        c = (waters[i] - c1) / (c2 -c1)
        print ('c =', c)
        ax.plot(masses, data[i], color = cmap(c))
    