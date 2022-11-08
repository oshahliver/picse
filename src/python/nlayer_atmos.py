#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 19:29:43 2019

@author: oshah
"""

from PIMPphysicalparams import sigmaSB

import numpy as np
from matplotlib import pyplot as plt
import Atmosphere
import RTE
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PIMPphysicalparams import au, sigmaSB, r_sun, T_sun

def matrix(n = 2, eps_arr = [], abs_arr = [], frac_arr = []):
    """Compute flux matrix for n-layer N-band model. Conventions for parameter
    inputs are:
        
        eps_arr (Nxn array): first index refers to the band and the second
        index refers to the layer
    """
    mat = np.zeros([n + 1, n + 1])

    #iterate over all frequency bands
    for a in range(len(abs_arr)):
        
        #iterate over all layers including the surface
        for i in range(n + 1):
            #print ('-------\ni=', i)
            #set diagonal element for layer i
            mat[i][i] += eps_arr[i]*frac_arr[a][i]
                
            #set elements to the left of diagonal (upward flux compontents)
            #iterate over all layers below layer i
            for j in range(i):
                #print ('\nj1 =', j)
                #print ('eps[j] =', eps_arr[j])
                #print ('frac[a][j] =', frac_arr[a][j])
                k = eps_arr[j]*frac_arr[a][j]
                #print ('k =', k)
                for kk in range(i-j):
                 #   print ('abs =', abs_arr[a][j + 1 + kk])
                    k *= (1. - abs_arr[a][j + 1 + kk])
                                    
                mat[i][j] += k
                #print ('mat[i][j] =', mat[i][j])
        
            #set elements to the right of diagonal (downsward flux components)
            #iterate over all layers above layer i
            for j in range(n - i):
                k = eps_arr[j + i + 1]*frac_arr[a][j]
                #print ('\nj2 =', j)
                #print ('eps[i+j+1] =', eps_arr[i+j+1])
                #print ('frac[a][j] =', frac_arr[a][j])
                #print ('k =', k)
                for kk in range(j - i - 1):
                  #  print ('kk =', kk)
                 #   print ('ind =', j  + kk)
                    k *= (1. - abs_arr[a][i  + kk])
                                    
                mat[i][j+i+1] -=  k
                #print ('mat[i][j] =', mat[i][j])
                
    return mat


def rhs(n = 2, albedo = [], S0_arr = [], abs_arr = []):
    vec = np.zeros([n + 1])
    
    #iterate over all frequency bands
    for nu in range(len(S0_arr)):    
        
        #iterate over all layers including the surface
        for i in range(n + 1):
            k = .25/sigmaSB * S0_arr[nu] * (1. - albedo[nu])
            
            for kk in range(n - i):
                k *= (1. - abs_arr[nu][kk + 1 + i])
                
            vec[i] += k
    
    return vec


def AttenuationShort(abs_arr=[]):
    
    n = len(abs_arr[0]) - 1
    
    opt = np.zeros([n+1])
    
    for i in range(n+1):
        opt[i] = 1.
        for k in range(n-i):
            print (abs_arr[0][k])
            opt[i] *= (1. - abs_arr[0][k + i + 1])
    
    print (opt)
    return -np.log(opt)

def attenuation(eps_arr=[], temp=[], albedo=.0, S0=0.):
    n = len(eps_arr[0]) - 1
    
    mat = matrix(n=n, eps_arr=eps_arr)
    
    vec = rhs(n=n, albedo=[albedo], S0_arr=[S0], eps_arr=eps_arr)
    
    netto = np.zeros([n + 1])
    
    for i in range(n + 1):
        for j in range(n + 1):
            netto[i] += mat[i][j] * temp[j]
    
    print (vec)
    return (4.*sigmaSB)*(netto-vec)/S0
    

def two_band_model(albedo = [0.3, .0], eps_arr = [.95, .78], 
                   abs_arr = [[0., 0.], [.95, .78]], T_star = T_sun, a = 1., 
                   r_star = 1., max_iter = 10, plot = False):
    
    eps0_init, eps1_init = eps_arr
        
    #compute total flux at distance of planet
    flux = sigmaSB*(r_star*r_sun)**2*T_star**4/(a*au)**2
    
    #divide incident flux into short-wave and long-wave band
    frS0s = RTE.cdf_Pl_spec(T=T_star, l=3.0e-6, N=100)
    S0s = flux*frS0s
    S0l = (1. - frS0s)*flux
    
    #first estimate temperature using a 1-band model
    mat = matrix(n = 1, eps_arr = [eps0_init, eps1_init],
                 abs_arr = abs_arr,
                 frac_arr=[[.0, .0], [1., 1.]])
    
    vec = rhs(n=1, albedo=albedo, S0_arr=[flux, 0.], 
              abs_arr=abs_arr)
    
    print ('mat=', mat)
    print ('vec=', vec)
    
    temp = np.linalg.solve(mat, vec)
    temp = temp**(1./4.)
    
    print ('initial temp guess:', temp)
    surface_temp_list = [temp[0]]
    atm_temp_list = [temp[1]]
    
    #start two band model iteration based on initial temperature guess
    iterate = True
    iterationCount = 0
    while iterate:
        
        if iterationCount > max_iter:
            break
        
        iterationCount += 1
        old_surface_temp = temp[0]
        
        print ('------------')
        print ('#', iterationCount)
        print ('------------')        
        
        mat = matrix(n = 1, eps_arr = [eps0_init, eps1_init],
                     abs_arr = abs_arr,
                     frac_arr=[[.0, .0], [1., 1.]])
        
        vec = rhs(n=1, albedo=albedo, S0_arr=[S0s, S0l], 
                  abs_arr=abs_arr)
        
        #print ('mat =', mat)
        #print ('vec =', vec)
        temp = np.linalg.solve(mat, vec)
        temp = temp**(1./4.)
        
        surface_temp_list.append(temp[0])
        atm_temp_list.append(temp[1])
        
        print ('temp =', temp)
        
        dev = abs(old_surface_temp - temp[0])/old_surface_temp
        if dev < .001:
            print ('Convergance reached after', iterationCount,' iterations\n')
            iterate = False
    
    if plot:
        x = np.arange(0, iterationCount + 1)
        
        fig, ax = plt.subplots()
        ax.plot(x, surface_temp_list, marker='s', linestyle = '--', color='g', 
                label='surface')
        ax.plot(x, atm_temp_list, marker='s', linestyle = '--', color='b',
                label ='atmosphere')
        
        ax.set_xlabel(r'$\# \ iterations$')
        ax.set_ylabel(r'$T_{eq} \ [K]$')
        
        major_ticks_x = x
        
        ax.set_xticks(major_ticks_x)
        ax.minorticks_off()
        ax.legend()
        ax.grid(axis='y')
    
    return temp


def specFrac(T=None, frequencies=[], N=200):
    bands = len(frequencies) + 1
    frac = np.zeros([bands])
    frac[0] = RTE.cdf_Pl_spec(T=T, l=frequencies[0]*1.0e-6, N=N)
    
    for f in range(bands - 2):
        frac_left = RTE.cdf_Pl_spec(T=T, l=frequencies[f]*1.0e-6, N=N)
        frac_right = RTE.cdf_Pl_spec(T=T, l=frequencies[f+1]*1.0e-6, N=N)
        
        frac[f+1] = frac_right - frac_left
        
    frac[-1] = 1. - sum(frac)
    return frac
    

def n2_model(albedo = [.3, 0.], eps_arr = [.95, .78], T_star = T_sun,
             r_star = 1, a = 1, abs_arr = [[], []],
             frequencies=[3.], bins=200, plot = False):
    """n-layer 2-band 1d atmpsphere model. The emmitted flux in each layer
    is distributed over the to frequency bands (long and short wave radiation)
    which are by default defined as short: lambda < 3 microns. For each layer,
    the spectral emmissivities have to be estimated from the total emissivity
    and the emissions spectrum at a given layer temperature. Therefore, the 
    layer temperature must be known to find the emissivities which in turn is
    needed to find the temperature. The approach to solve this chicken and egg
    problem is as follows: estimate the temperature in layer i using a i-layer
    model starting at i=0 up to i = n. This means that first, a naked planet
    model is envoked to estimate the surface temperature T0. From this the
    spectral emmissivities are computed using the cdf at T0 from which the
    fluxes at layer 1 can be estimated. Using this emissivities now in a 1 layer
    model gives a first correction to T0 and a first estimate for T1. Using the
    new T0 the emmissivities in layer 1 can be updated. This is iteratively done
    unitl the correction to T0 converges. Then a new layer is added and the 
    proceedure is repeated up to n layers. 
    """
    
    n = len(eps_arr) - 1
    N_freq = len(frequencies) + 1
    
    #compute total flux at distance of planet
    flux = sigmaSB*(r_star*r_sun)**2*T_star**4/(a*au)**2
    print ('incident flux =', flux)
    
    #compute incident flux array for frequency bands
    frS_in = RTE.cdf_Pl_spec(T=T_star, l=frequencies[0]*1.0e-6, N=100)

    S_in_arr = np.zeros([len(frequencies)+1])
    S_in_arr[0] = frS_in*flux
    
    #gather contribution from individual frequencies
    S_in_arr = specFrac(T=T_star, frequencies=frequencies) * flux

    #initially set fraction to 0 for short wave and 1 for long wave
    S_arr = np.zeros([N_freq, n + 1])
    
    for i in range(n + 1):
        S_arr[0][i] = 0.
        S_arr[1][i] = 1.
    
    mat = matrix(n=n, eps_arr = eps_arr,
                 abs_arr = abs_arr,
                 frac_arr=S_arr)
    
    vec = rhs(n=n, albedo=albedo, S0_arr=S_in_arr,
              abs_arr=abs_arr)

    print ('mat =', mat)
    print ('vec =', vec)
    
    temp = np.linalg.solve(mat, vec)
    temp = temp**.25
    '''
    temp = two_band_model(eps_arr = eps_arr, albedo=albedo, T_star=T_star,
                          abs_arr = abs_arr)

    #compute new layer flux array for frequency bands
    for i in range(2):
        frS = RTE.cdf_Pl_spec(T=temp[i], l=frequencies[0]*1.0e-6, N=bins)
    
        S = eps_arr[i]*sigmaSB*temp[i]**4
        
        #the fraction of the total flux is stored, not the flux it self
        S_arr[0][i] = frS
        for f in range(N_freq - 2):
            frS_left = RTE.cdf_Pl_spec(T=temp[i], l=frequencies[f]*1.0e-6,
                                        N=bins)
            frS_right = RTE.cdf_Pl_spec(T=temp[i], l=frequencies[f+1]*1.0e-6,
                                         N=bins)
            frS = frS_right - frS_left
            S_arr[f+1] = S*frS
                         
            S_arr[f+1][i] = frS
            
        S_arr[-1][i] = (1. - frS)
        
    print ('S_arr =', S_arr)
    '''
    if plot:
        fig, ax = plt.subplots()
        
        ax.plot(temp, np.arange(0, n+1), marker='s')
    
    return temp
    

def nN_model(albedo_arr = [], eps_arr = [], T_star=T_sun, albedo=.3,
             r_star=1, a=1, frequencies=[], eps_tot_arr = [], plot=False):
    
    n = len(eps_arr[0]) - 1
        
    for eps in eps_arr:
        if not len(eps) == n + 1:
            print ('WARNING: array missmatch')
    
    #compute total flux at distance of planet
    flux = sigmaSB*(r_star*r_sun)**2*T_star**4/(a*au)**2
    print ('incident flux =', flux)
    
    
    #compute incident flux array for frequency bands
    frS0 = RTE.cdf_Pl_spec(T=T_star, l=frequencies[0]*1.0e-6, N=100)

    S0_arr = np.zeros([len(frequencies)+1])
    S0_arr[0] = frS0*flux
    for f in range(len(frequencies) - 1):
        frS0_left = RTE.cdf_Pl_spec(T=T_star, l=frequencies[f]*1.0e-6,
                                    N=100)
        frS0_right = RTE.cdf_Pl_spec(T=T_star, l=frequencies[f+1]*1.0e-6,
                                     N=100)
        frS0 = frS0_right - frS0_left
        S0_arr[f+1]=flux*frS0
        
    S0_arr[-1]=flux-sum(S0_arr)    
    
    print ('S0_arr =', S0_arr)
    
    '''
    eps_dummy_arr = np.zeros([2, n + 1])
    
    for i in range(n + 1):
        eps_dummy_arr[0][i] = 0.
        eps_dummy_arr[1][i] = eps_tot_arr[i]
    '''
    
    #first estimate temperature using a 1-band model
    vec = rhs(n=n, albedo=albedo_arr, S0_arr=S0_arr,
              eps_arr = eps_arr)
    
    print ('eps_arr =', eps_arr)
    
    mat = matrix(n=n, eps_arr=eps_arr)
    
    print ('vec=', vec)
    #print ('mat=', mat)
    
    temp = np.linalg.solve(mat, vec)
    temp = temp**.25
        
    print ('initial temp guess=', temp)
    
    '''
    #compute new emissivities in each layer
    #iterate over layers
    for i in range(n + 1):
        print ('---------')
        print ('updating layer', i)
        print ('computing T =', temp[i])
        #compute emitted flux for all frequency bands
        frS = RTE.cdf_Pl_spec(T=temp[0], l=frequencies[0]*1.0e-6, N=100)
    
        S_arr = np.zeros([len(frequencies)+1])
        S_arr[0] = eps_arr[0][i]*frS
        
        eps_arr[0][i] = S_arr[0]
        
        #iterate over frequencies
        for f in range(len(frequencies) - 1):
            frS_left = RTE.cdf_Pl_spec(T=temp[i+1], l=frequencies[f]*1.0e-6, 
                                        N=100)
            frS_right = RTE.cdf_Pl_spec(T=temp[i+1], l=frequencies[f+1]*1.0e-6,
                                         N=100)
            frS = frS_right - frS_left
            S_arr[f+1]=eps_tot_arr[i]*frS
            
            eps_arr[f+1][i] = eps_arr[f+1][i]*frS
            
        S_arr[-1]=eps_arr[-1][i]-sum(S_arr)
        
        eps_arr[-1][i] = S_arr[-1]
        
        print ('S_arr=', S_arr)
        print ('sum S_arr=', sum(S_arr))
    
    print ('new eps_arr =', eps_arr)
    
    #compute new temperature profile using new emissivities
    #first estimate temperature using a 1-band model
    vec = rhs(n=n, albedo=albedo_arr, S0_arr=S0_arr, 
              eps_arr = eps_arr)
    
    mat = matrix(n=n, eps_arr=eps_arr)
    
    #print ('mat=', mat)
    print ('vec=', vec)

    temp = np.linalg.solve(mat, vec)
    temp = temp**.25
    '''
    if plot:
        fig, ax = plt.subplots()
        
        ax.plot(temp, np.arange(0, n+1), marker='s')
    
    return temp


def solve_atmos(n = 2, S0 = 1367., albedo = .3, eps_arr = []):
    #initialize flux matrix
    mat = matrix(n = n, eps_arr = eps_arr)
    
    #incoming flux at top of atmosphere
    rhs = .25/sigmaSB * S0 * (1. - albedo)
    
    #initialize rhs vector
    vec = np.zeros([n + 1])
    
    for v in range(n + 1):
        vec[v]  = rhs
        
    x = np.linalg.solve(mat, vec)
    
    return x**.25
    

def iterate_atmos(n = 2, S0 = 1367., albedo = .3):
    
    #initialize and construct atmosphere instance
    atm = Atmosphere.Atmosphere(contents=[1], tempType = 'adiabatic', 
                                eps_r = .5)
    atm.construct()
    atm.plot()
    #extract atmosphere profile
    nlayers = len(atm.layers[0].shells)
    
    temp = np.zeros([nlayers + 1])
    pres = np.zeros([nlayers + 1])
    dens = np.zeros([nlayers + 1])
    height = np.zeros([nlayers + 1])
    
    temp[0] = atm.temp_surface
    pres[0] = atm.pres_surface
    dens[0] = atm.dens_surface
    height[0] = 0.0
    
    lay = atm.layers[0]
    for i in range(nlayers):
        temp[i + 1] = lay.shells[i].temp
        pres[i + 1] = lay.shells[i].pres
        dens[i + 1] = lay.shells[i].dens
        height[i + 1] = lay.shells[i].radius - atm.R_surface
    
    #compute optical depth for each shell
    opt = np.zeros([nlayers + 1])
    opt[0] = 1.
    for i in range(nlayers):
        opt[i + 1] = abs(dens[i] - dens[i + 1])/2.
    
    print ('opt =', opt)
    y = solve_atmos(n = nlayers, eps_arr = opt)
    
    fig, ax = plt.subplots()
    ax.plot(height, y)    
    
    return y