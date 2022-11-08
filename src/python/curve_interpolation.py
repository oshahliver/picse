#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 17:38:45 2021

@author: os18o068
"""

import matplotlib.pyplot as plt
import numpy as np

def interpolate(a1, a2, poly_deg=3, n_points=100, plot=True, method = 'poly'):

    min_a1_x, max_a1_x = min(a1[:,0]), max(a1[:,0])
    min_a2_x, max_a2_x = min(a2[:,0]), max(a2[:,0])    

    new_a1_x = np.linspace(min_a1_x, max_a1_x, n_points)
    new_a2_x = np.linspace(min_a2_x, max_a2_x, n_points)
    
    if method == 'poly':
        a1_coefs = np.polyfit(a1[:,0],a1[:,1], poly_deg)
        a2_coefs = np.polyfit(a2[:,0],a2[:,1], poly_deg)

        new_a1_y = np.polyval(a1_coefs, new_a1_x)
        new_a2_y = np.polyval(a2_coefs, new_a2_x)
        
        midx = [np.mean([new_a1_x[i], new_a2_x[i]]) for i in range(n_points)]
        midy = [np.mean([new_a1_y[i], new_a2_y[i]]) for i in range(n_points)]

    elif method == 'interp':
        new_a1_y = np.interp(new_a1_x, a1[:,0], a1[:,1])
        new_a2_y = np.interp(new_a2_x, a2[:,0], a2[:,1])
        
        midx = [np.mean([new_a1_x[i], new_a2_x[i]]) for i in range(100)]
        midy = [np.mean([new_a1_y[i], new_a2_y[i]]) for i in range(100)]       

    if plot:
        plt.plot(a1[:,0], a1[:,1],c='black')
        plt.plot(a2[:,0], a2[:,1],c='black')
        plt.plot(midx, midy, '--', c='black')
        plt.show()

    return np.array([[x, y] for x, y in zip(midx, midy)])