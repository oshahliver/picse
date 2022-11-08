#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 11:39:37 2021

@author: os18o068
"""

import numpy as np
import numpy.random
import matplotlib.pyplot as plt
from matplotlib import cm

def selection(XY, limitXY=[[-2,+2],[-2,+2]]):
    XY_select = []
    for elt in XY:
        if elt[0] > limitXY[0][0] and elt[0] < limitXY[0][1] and elt[1] > limitXY[1][0] and elt[1] < limitXY[1][1]:
            XY_select.append(elt)

        return np.array(XY_select)

def plot(x, y, N = 8, limitXY = [[-1,+1], [-1,+1]], title = 'hist3d.png'):
    XY = np.stack((x,y),axis=-1)
    limitXY = [[min(x), max(x)], [min(y), max(y)]]
    print ('XY =', XY)
    XY_select = selection(XY, limitXY=limitXY)
    print ('XY_select =', XY_select)
    
    fig = plt.figure() #create a canvas, tell matplotlib it's 3d
    ax = fig.add_subplot(111, projection='3d')
    
    hist, xedges, yedges = np.histogram2d(x, y, bins=(N,N), range = limitXY) # you can change your bins, and the range on which to take data
    # hist is a 7X7 matrix, with the populations for each of the subspace parts.
    xpos, ypos = np.meshgrid(xedges[:-1]+xedges[1:], yedges[:-1]+yedges[1:]) -(xedges[1]-xedges[0])
    
    xpos = xpos.flatten()*1./2
    ypos = ypos.flatten()*1./2
    zpos = np.zeros_like (xpos)
    
    dx = xedges [1] - xedges [0]
    dy = yedges [1] - yedges [0]
    dz = hist.flatten()
    print ('xpos =', xpos)
    cmap = cm.get_cmap('viridis') # Get desired colormap - you can change this!
    max_height = np.max(dz)   # get range of colorbars so we can normalize
    min_height = np.min(dz)
    # scale each z to [0,1], and get their rgb values
    rgba = [cmap((k-min_height)/max_height) for k in dz] 
    
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=rgba, zsort='average')
    plt.xlabel(r'Total mass [$M_\oplus$]')
    plt.ylabel(r'Surface temperature [K]')
    plt.title(title)
    plt.savefig(title)
    plt.show()