#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 16:24:25 2019

@author: oshah
"""

import numpy as np
import random
from PIMPphysicalparams import kB, c_light, hPl, sigmaSB
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from PIMPrunparams import color_list
import functionTools as ftool

dl_min = .2
dl_max = 50.

def Wien(T):
    """Computes position of maximum in m in Planck Spectrum according to Wien's
    law
    """
    return 2897.8e-6/T

def B_Pl(T=None, nu=None):
    """Compute Planck spectrum at given temperature T in K and frequency
    in s-1
    """
    return 2.*hPl*nu**3/c_light**2 / (np.exp(hPl*nu/(kB*T)) -1.)


def f_Pl(T=None, nu=None, N=10, l_min=0.1e-6, l_max=100.0e-6,
         scale = 'linear'):
    """Approximate B_Pl as stepwise constant to sample from distribution.
    The Planck Spectrum is divided into N bins between l_min and l_max
    given in microns. Relevant wavelength range for atmospheric modelling is
    between .1 and 100 microns
    """
    
    #this array will contain the distribution values
    vals = np.zeros([N])
    
    #this array will contain the wavelength at the bin centers
    bins = np.zeros([N])
        
    if scale == 'linear':
        dl = (l_max - l_min)/N
        bins = np.linspace(l_min + .5 * dl, l_max - .5 * dl, N)

        for i in range(N):
            vals[i] = B_Pl(T=T, nu = c_light/bins[i])

    elif scale == 'log':
        dl = np.log10(l_max/l_min)/N
        bins = np.linspace(np.log10(l_min) + .5 * dl, 
                           np.log10(l_max) - .5 * dl, N)

    
        for i in range(N):
            vals[i] = B_Pl(T=T, nu = c_light/10**(bins[i]))
    
    #compute self consistent normalization of N-bin approximation
    norm = 0.0
    for i in range(N):
        if scale == 'linear':
            norm += dl*vals[i]
            
        elif scale == 'log':
            #convert constant dl in log space back into lin space
            dl_lin = 10**(bins[i] + .5 * dl) - 10**(bins[i] - .5 * dl)
            norm += dl_lin*vals[i]
            
    
    #return everything in units of microns
    return bins, vals/norm, norm, dl


def cdf_Pl(T=None, N=20, l_min=0.1e-6, l=100.0e-6, scale = 'linear'):
    """Computes cumualative distribution function by approximating the pdf
    with a finite number of N bins and computing the integral as a the sum of
    the bins.
    """
    bins, vals, norm, dl = f_Pl(T=T, N=N, l_min=l_min, l_max=l, scale = scale)
    cdf_list = np.zeros([N])
    
    for i in range(N):
        for j in range(i):
            if scale == 'linear':
                #add bin area to integral
                cdf_list[i] += vals[j]*dl
                
            elif scale == 'log':
                #convert constant dl in log space back into lin space
                dl_lin = 10**(bins[j] + .5 * dl) - 10**(bins[j] - .5 * dl)
                
                #add bin area to integral
                cdf_list[i] += vals[j]*dl_lin
               
    return bins, cdf_list


def fit_cdf_Pl(T=None, N=20, l_min=0.1e-6, l_max=100.0e-6, order=4,
               scale = 'linear'):
    
    bins, cdf_list = cdf_Pl(T=T, N=N, l_min=l_min, l=l_max, scale = scale)
    
    z = np.polyfit(bins, cdf_list, order)
    poly = np.poly1d(z)
    
    return poly


def cdf_Pl_spec(T=None, N=20, l_min=0.1e-6, l_max=None, l = None,
                scale = 'linear'):
    """Generates cdf_Pl between l_min and l_max in N bins and computes the
    value at the wavelength l in micrometers
    """

    #estimate maximum of spectrum using Wien's law
    if l_max == None:
        ll = Wien(T)
        l_min = ll * dl_min
        l_max = ll * dl_max
    
    bins, cdf_list = cdf_Pl(T=T, N=N, l_min=l_min, l=l_max, scale=scale)
    
    #find closest value to l in bins
    diff = np.zeros([N])
    for d in range(N):
        diff[d] = abs(bins[d] - l)
        
    ind = np.where(diff == np.amin(diff))[0][0]
    return cdf_list[ind]


def fitted_cdf(T=None, l=None, order=4, N=50, poly=[], scale = 'linear'):
    #estimate maximum of spectrum using Wien's law
    ll = Wien(T)
    l_min = ll * dl_min
    l_max = ll * dl_max

    #if no poly fit has been passed, generate it first
    if len(poly) == 0:
        poly = fit_cdf_Pl(T=T, N=50, l_min=l_min, l_max=l_max, order=4,
                          scale = scale)
        
    else:
        pass
    
    if scale == 'linear':
        return np.maximum(poly(l), 0.001)
    
    elif scale == 'log':
        return np.maximum(poly(np.log10(l)), 0.001)


def plot_f_Pl(temps=[6000, 5000, 4000], N=25, order=4 , scale = 'log', 
              l_min=0.1e-6, l_max=None):
    """here l_min and l_max are in micro meters
    """
    fig, ax = plt.subplots(1, 2)
    
    max_list = []
    
    all_l_min = []
    all_l_max = []

    plot_legend_list1 = []
    plot_label_list1  = []

    plot_legend_list2 = []
    plot_label_list2  = ['binned integral', 'polynomial fit']
    
    c = 0
    for T in temps:
        #estimate maximum of spectrum using Wien's law
        if l_max == None:
            ll = Wien(T)
            l_min = ll * dl_min
            l_max = ll * dl_max
            
        else:
            pass
        
        all_l_max.append(l_max)
        all_l_min.append(l_min)
        
        color = color_list[c]
        c += 1
        
        bins, vals, norm, dl = f_Pl(T=T, N=N, l_min=l_min, l_max=l_max, 
                                    scale=scale)
        
        if scale == 'log':
            x = np.linspace(np.log10(l_min), np.log10(l_max), 50)
            
        elif scale == 'linear':
            x = np.linspace(l_min, l_max, 50)
                
        bins, cdf_list = cdf_Pl(T=T, N=N, l_min=l_min, l = l_max, scale = scale)
        max_list.append(max(vals))
        
        for i in range(N):
            print ('cdf/bin:', round(cdf_list[i], 4),  ', ', bins[i])
        
        #plot PDF
        for i in range(N):
            if scale == 'linear':
                rect = patches.Rectangle(((bins[i]-.5*dl)*1.0e6, 0.), 
                                     dl*1.0e6, 
                                     vals[i]*1.0e-5,
                                     edgecolor= 'w', 
                                     facecolor = color)
            
            elif scale == 'log':
                 rect = patches.Rectangle(((bins[i]-.5*dl)+6, 0.), 
                                          dl, 
                                          vals[i]*1.0e-5,
                                          edgecolor= 'w', 
                                          facecolor = color) 
                 
            ax[0].add_patch(rect)
        
        #perform polynomial fit to cdf
        poly = fit_cdf_Pl(T=T, N=N, l_min=l_min, l_max=l_max, order=order,
                          scale = scale)
        

        if scale == 'linear':
            pl1, = ax[1].plot(bins*1.0e6, cdf_list, color=color, linestyle='None',
              marker='s', markevery=1, markersize=5, label = str(T)+' K')
            
            pl2, = ax[1].plot(x*1.0e6, fitted_cdf(T=T, poly=poly, l=x), 
                      color=color)

        elif scale == 'log':
            pl1, = ax[1].plot(bins+6., cdf_list, color=color, linestyle='None',
              marker='s', markevery=1, markersize=5, label = str(T)+' K')
            
            pl2, = ax[1].plot(x+6, fitted_cdf(T=T, poly=poly, l=x), 
                      color=color)            
            
        if T == temps[0]:
            plot_legend_list2 = [pl1, pl2]
            
        plot_legend_list1.append(pl1)
        plot_label_list1.append(str(T)+' K')

    
    legend1 = ax[0].legend(plot_legend_list1, plot_label_list1, frameon=False)
    legend2 = ax[1].legend(plot_legend_list2, plot_label_list2, frameon=False)
    ax[0].add_artist(legend1)
    ax[1].add_artist(legend2)
    
    ax[0].set_ylim(0., max(max_list)*1.1*1.0e-5)
    
    #ax[1].set_xlim(l_min*1.0e6, l_max*1.0e6)
    ax[1].set_ylim(0., 1.)
    
    if scale == 'linear':
        ax[0].set_xlim(min(all_l_min)*1.0e6, max(all_l_max)*1.0e6)
        ax[0].set_xlabel(r'$\lambda \ [\mu m]$')
        ax[1].set_xlabel(r'$\lambda \ [\mu m]$')


    elif scale == 'log':    
        ax[0].set_xlim(np.log10(min(all_l_min))+6., np.log10(max(all_l_max))+6.)
        ax[0].set_xlabel(r'$log(\lambda) \ [\mu m]$')
        ax[1].set_xlabel(r'$log(\lambda) \ [\mu m]$')
        
        
    ax[0].set_ylabel(r'$PDF_T(\lambda) \ [10^{5}]$')
    
    ax[1].set_ylabel(r'$CDF_T(\lambda)$')

    for axx in ax:
        axx.tick_params(direction='out')

    ax[1].yaxis.set_label_position('right')
    ax[1].yaxis.tick_right()
    
    ax[0].text(1.5, .5, 'N ='+str(N), fontsize=14)
    ax[0].grid()
    ax[1].grid()
    

def sample_Pl(T=None, q=None):
    val = ftool.bisec(a=0., b=1., y=q, f=cdf_Pl)
    
    
def MC2(photons=1000, sigma_a = [.01], sigma_s = [.1], d = 1., nphotons=10,
        n=10, f_bands = [1.0e-6]):
    """d: box height
    """
    Rd = 0 #diffuse reflectance, number of photons leaving the slab through the top
    Tt = 0 #number of photons leaving the slab through the bottom
    A = 0 # number of absorbed photons
    
    
    nb = len(sigma_a)
    
    all_P = [] #final positions of photons that made it out of the top
    
    #define 2d detector geometry
    xmax = 2.
    ymax = 2.
    
    nx = n
    ny = n
    
    dx = 2 * xmax / nx
    dy = 2 * ymax / ny
    
    #initalize detector grid
    x_frame = np.linspace(-xmax, xmax, nx)
    y_frame = np.linspace(-ymax, ymax, ny)
    
    xx, yy = np.meshgrid(x_frame, y_frame)
    
    #set up data frame for each spectral band individually
    dat = np.zeros([nb, nx, ny])
    
    #total extinction coefficient as combination of scattering and absorbtion
    sigma_t = []
    
    for i in range(len(sigma_a)):
        sigma_t.append(sigma_s[i] + sigma_a[i])
    
    for i in range(photons):
        
        P = np.array([0., 0., 0.]) # initial position of photon
        mu = np.array([1.0e-10, 1.0e-10, .99999999]) # initial direction of photon
        
        #set initial light ray weight to unitiy
        w = 1.
        
        #sample Panck function to determine ray wavelength
        q = np.random.uniform(0, 1)
        
        
        #determine wavelength of the photon ray by underlying distribution
        #by default a Planck distribution is assumed
        
        if nb > 1:
            band = int(round(random.uniform(0, nb-1), 0))
        
        else:
            band = 0
        
        while True and w > 0.:     
            q = random.uniform(0, 1)     
            s = -np.log(q)/sigma_t[band]
            
            #update photon position
            P = P + s*mu
            
            #if photon leaves surface through top
            if P[2] > d:
                Rd += w
                all_P.append(P)
                #assign photon to corresponding detector cell
                for i in range(len(x_frame)):
                    for j in range(len(y_frame)):
                        if abs(P[0]-x_frame[i]) < dx/2. \
                        and abs(P[1]-y_frame[j]) < dy/2.:
                            dat[band][i][j] += w
                break
            
            #if photon leaves surface through bottom
            elif P[2] < 0.:
                Tt += w
                break
            
            #if ray is still in the box, compute new direction after scattering
            #and how much of the intensity is lost due to absorption
            else:                
                #absorbtion: reduce weight of light packet according to the
                #on average, a certain fraction of all photons in the ray 
                #will be absorbed
                dw = sigma_a[band] / sigma_t[band]
                
                #on average, at each scattering event, this fraction is
                #removed from the packet
                w -= dw
                
                #the weight cannot be negative
                w = max(0., w)
                
                #scattering: update direction of travel
                #compute azimuthal angle of scattering vector
                phi_s = 2. * np.pi * random.uniform(0, 1)
                
                #compute polar angle of scattering vector
                mu_s = 1. - 2.*random.uniform(0, 1)
                
                #update direction
                matrix = np.array([[mu[0]*mu[2]/np.sqrt(1.-mu[2]**2),
                                    -mu[1]/np.sqrt(1.-mu[2]**2),
                                    mu[0]],
                                   [mu[1]*mu[2]/np.sqrt(1.-mu[2]**2),
                                    mu[0]/np.sqrt(1.-mu[2]**2),
                                    mu[1]],
                                   [-np.sqrt(1.-mu[2]**2),
                                    0.,
                                    mu[2]]])
                
                b = np.array([np.sqrt(1.-mu_s**2)*np.cos(phi_s),
                              np.sqrt(1.-mu_s**2)*np.sin(phi_s),
                              mu_s])

                mu = np.linalg.solve(matrix, b)
                mu = mu/np.sqrt(sum(mu**2))
    
    fig, ax = plt.subplots(1, nb)
    
    
    if not nb == 1:
        for i in range(nb):
            ax[i].imshow(dat[i])
            
    else:
        ax.imshow(dat[0])

    return Rd/photons, Tt/photons
