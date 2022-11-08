#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 17:06:08 2020

@author: os18o068
"""

import ternary
import numpy as np
from matplotlib import pyplot as plt
import Material
import abundances
import MyCmaps as mcm
import pyrolite.plot
import pandas as pd
from mpltern.ternary.datasets import get_scatter_points
from pyrolite.util.math import flattengrid
from pyrolite.util.plot.axes import axes_to_ternary, share_axes
from pyrolite.comp.codata import ILR, inverse_ILR
from pyrolite.plot.density.ternary import ternary_heatmap


np.random.seed(43)

df = pd.DataFrame(np.array([*get_scatter_points(n=80)]).T, columns=['A', 'B', 'C'])
df = df.loc[(df > .1).all(axis=1), :]

N_pnts = 100
Si_numbers = np.linspace(.4, .6, 5)
xi_Stv = np.linspace(0., 1., N_pnts)
points = np.empty([N_pnts, 3])
scale=25
fnts=14
off = .2
        
    
def get_SiMg(p, Fe_number, mantle='upper'):
    if mantle == 'upper':
        Mg = 2*(p[0]+p[1])*(1.-Fe_number)
        Si = p[2] + p[0] + 2*p[1]
        
    elif mantle == 'lower':
        Mg = (p[0]+p[1])*(1.-Fe_number)
        Si = p[2] + p[1]
    
    try:
        return Si/Mg
    
    except ZeroDivisionError:
        return scale

def f(p):
    
    SiMg = get_SiMg([p[i]/scale for i in range(len(p))], 0., 
                    mantle='upper')

    pp = abundances.abund.compute_abundance_vector(simg=SiMg, 
                                          femg=0., 
                                          n_mats=3, 
                                          contents=[6,7,3],
                                          additional=[p[2]/scale], 
                                          xifei=[0., 0., 0.])
    
    check_validity = 0
    for i in range(len(pp)):
        if pp[i] < 0.:
            check_validity += 1
            
    if check_validity > 0:
        SiMg = SiMg
    
    s=0.
    for i in range(len(p)):
        try:
            s += p[i]
            
        except ValueError:
            print ('no')
            continue
        
    try:
        Si_number = SiMg/(1.+SiMg)
    except TypeError:
        Si_number = 0.333
    
    return Si_number


cm = mcm.CMAP().cm


ticks = np.linspace(0., 2., 10)
ticklabels = [str(t) for t in ticks]

figure, tax = ternary.figure(scale=scale)

tax.heatmapf(f, boundary=True, style='hexagonal', 
             cmap=cm, vmin=.0, vmax=1., colorbar=False)

tax.gridlines(color='w', multiple = scale/10, linewidth=.5, 
              linestyle='-')

tax.ticks(axis='lbr', multiple=scale/10, offset=.025)

for l in range(len(Si_numbers)):
    for i in range(N_pnts):
        points[i] = abundances.abund.compute_abundance_vector(simg=Si_numbers[l]/(1.-Si_numbers[l]),
                                                  femg=.0, contents=[6,7,3],
                                                  additional=[xi_Stv[i]],
                                                  xifei=[.0, .0, .0])
    
    
    
        check_validity = 0
        for j in range(len(points[i])):
            if points[i][j] < 0.:
                check_validity += 1
                
        if check_validity > 0:
            points[i] = [None, None, None]
            
    t = tax.plot(points*scale, linestyle='-', linewidth=2., color='k')

ax = tax.get_axes()

im=ax.imshow(points, cmap=cm, vmin=.0, vmax=1.)

cb = figure.colorbar(im, cmap=cm, orientation='horizontal', shrink=.75,
                     label='Si#', pad=0.01)

cb.set_label(label='Si#', size=fnts)

ax.text(5., 1., 'Si# = 0.6', rotation = 105)
ax.text(7., 1., 'Si# = 0.55', rotation = 100)
ax.text(10., 1., 'Si# = 0.5', rotation = 90)
ax.text(12., 1., 'Si# = 0.45', rotation = 85)
ax.text(14., 1., 'Si# = 0.4', rotation = 80)

tax.boundary(linewidth=1.)
tax.get_axes().axis('off')
tax.set_title('Upper Mantle', pad=50., fontsize=2*fnts)
tax.right_corner_label('Ol', fontsize = fnts, offset=off)
tax.left_corner_label('Stv', fontsize = fnts, offset=off)
tax.top_corner_label('Pyr', fontsize = fnts, offset=off)

#tax.set_axis_limits({'b':[0,100], 'l':[0,100], 'r':[0,100]})

#tax.get_ticks_from_axis_limits()
#tax.set_custom_ticks(offset=.025, multiple=scale/10.)
tax.show()

fig, ax = plt.subplots()


ax = axes_to_ternary(ax)

