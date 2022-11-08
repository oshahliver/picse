#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 10:28:12 2020

@author: os18o068
"""
import numpy as np
import matplotlib as mpl

class CMAP():
    
    def __init__(self):
        border_colors = list(reversed([(.7, .2, .6),
                (.6, .4, 1.), 
              (.2, .4, 1.), 
              (.2,.6,.9), 
              (.2, .8, .8), 
              (.2, .8, .4), 
              (.6, .8, .4), 
              (.6, .6, .2), 
              (.8, .4, .2),
              (1., .2, .2),
              (1., .5, .5),
              (.7, .7, .7),
              (.5, .5, .5),
              (.2, .2, .2),
              (.0, .0, .0)
              ]))
        
        bin_refine_exponent = 4
        n_additional_bins = 2**bin_refine_exponent
        colors = []
        

        for i in range(len(border_colors)-1):
            colors.append(border_colors[i])
            for r in range(n_additional_bins):
                colors.append(np.asarray(border_colors[i])+\
                            (r+1)/(n_additional_bins+1)*(np.asarray(border_colors[i+1]) - \
                              np.asarray(border_colors[i])))
            
        colors.append(border_colors[-1])
            
        self.cm = mpl.colors.ListedColormap(colors)