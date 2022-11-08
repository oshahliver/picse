#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 18:54:45 2020

@author: os18o068
"""

import PlanetFactory as plfac
import numpy as np
import pandas as pd
import io
import pickle
import MOI_model
import matplotlib as mpl
import time
import fortfunctions

from PIMPphysicalparams import m_earth, r_earth, MoI_solar, delta_MoI_solar,\
    R_solar, J2_solar, delta_J2_solar, orbit_solar, day, yr, M_solar, T_solar,\
        material_list_fort, material_YMg, material_YSi, material_YO, material_YH, \
        material_YS, mS, mH2O, mFe, mSi, mO, mMg, mH, mEn, mH2O
import Material
import PlanetFort
from PIMPrunparams import color_list
from matplotlib import pyplot as plt
import matplotlib
import functionTools as ftool
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import random
from matplotlib.ticker import MaxNLocator

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath, amssymb}"]

data_dir = '/home/os18o068/Documents/PHD/Projects/Planets/Data/icy_satellites/'+\
    'compare_ice_fractions/'

toolkit = plfac.Toolkit()

ocean_mass_fractions = np.linspace(.0001, .6, 1)
impurity_abundances = np.linspace(0., .5, 1)

hard_ranges = np.array([[.5, 1.],       #Mg# 
                        [1./3., .5],  #Si#
                        [.0, .0],      #M_H2O
                        [0., 1.],       #xi_S
                        [.0, .3],    #xi_Fe
                        [.0, .5],      #X_ice_mean
                        [0., .1]])      #X_ice_0

fnts1 = 8
fnts2 = 10
fnts3 = 12
fnts4 = 16

def convert_X_ice_to_xi_ice(Si_number=None, xi_Fe=None, X_ice=None, 
                            contents = [6,7,1]):
    
    SiMg = Si_number/(1.-Si_number)
    
    #Compute fractions without water
    fractions = fortfunctions.functions.compute_abundance_vector(simg=SiMg, 
                                    femg=1./(1.-xi_Fe)-1.,
                                    n_mats=len(contents),
                                    ymgi=[material_YMg[i-1] for i in contents],
                                    ysii=[material_YSi[i-1] for i in contents],
                                    xih2oi=[0. for i in contents],
                                    xifei=[xi_Fe for i in contents],
                                    xialsii = [0. for i in contents],
                                    xialmgi = [0. for i in contents],
                                    contents = contents,
                                    additional = [.0])
    
    #Compute total normalized mass in the mantle
    m_tilde = sum([fractions[i]*(1.-xi_Fe)*material_YMg[contents[i]-1] 
              for i in range(len(contents))])*mMg+\
            sum([fractions[i]*xi_Fe*material_YMg[contents[i]-1] 
                 for i in range(len(contents))])*mFe +\
            sum([fractions[i]*material_YSi[contents[i]-1] 
                 for i in range(len(contents))])*mSi+\
            sum([fractions[i]*material_YO[contents[i]-1] 
                 for i in range(len(contents))])*mO+\
            sum([fractions[i]*material_YH[contents[i]-1] 
                 for i in range(len(contents))])*mH
    
    #Compute mole fraction of water at given composition
    xi = Material.xi(eta=X_ice, m1=m_tilde, m2=mH2O)
    
    return xi


def compute_max_ocean_frac(Si_number=None, xi_Fe=None, X_ice = None, 
                           contents = [6,7,1], Mg_number = None, 
                           H2O_frac = None, xi_S=None):
    """The total water mass fraction is (M_ice + M_ocean)/M. For a given amount
    of water in the mantle, given by X_ice, the possible value range for the
    ocean mass fraction is therefore limited.
    """
    
    #First convert the ice mass fraction in the mantle to mole fraction for
    #the given composition
    xi_ice = convert_X_ice_to_xi_ice(Si_number = Si_number, xi_Fe = xi_Fe,
                                     X_ice = X_ice, contents = [6,7,1])
    
    #With this, compute the core mass fraction without an ocean. Note that for
    #M_surface = 1 the value of the core mass corresponds to the value of the 
    #core mass fraction
    all_contents = [[2,8], [2,8], [4,5,1], [6,7,1], [1]]
    SiMg = Si_number/(1.-Si_number)
    
    print ('Core mass in compute max ocean')
    core_frac = PlanetFort.ComputeCoreMass(contents = all_contents, 
                                           Mg_number = Mg_number,
                                           Mg_number_mantle = 1.-xi_Fe, 
                                           SiMg = SiMg,
                                           M_surface = 1., 
                                           M_ocean = 0., 
                                           impurity = [xi_ice],
                                           xi_S_core = xi_S)
    
    print ('core_frac =', core_frac)
    #Compute ocean mass fraction
    ocean_frac = (H2O_frac - X_ice*(1.-core_frac))/(1.-X_ice+X_ice*core_frac)
    
    return ocean_frac
    

def compute_max_X_ice(Si_number=None, xi_Fe=None, contents=[6,7,1], Mg_number=None,
                      H2O_frac = None, xi_S=None):
    """The maximum amount of ice in the mantle corresponds to the total amount
    of water in the planet, i.e. to the case where M_ocean = 0.
    """
    #Compute core mass fraction
    all_contents = [[2,8], [2,8], [4,5,1], [6,7,1], [1]]
    SiMg = Si_number/(1.-Si_number)

    core_frac = PlanetFort.ComputeCoreMass(contents = all_contents, 
                                           Mg_number = Mg_number,
                                           Mg_number_mantle = 1.-xi_Fe,
                                           SiMg = SiMg,
                                           M_surface = 1., 
                                           M_ocean = H2O_frac, 
                                           impurity = [0.],
                                           xi_S_core = xi_S)
    
    X_ice_max = H2O_frac/(1.-core_frac)

    return min(X_ice_max, hard_ranges[5][1]), core_frac


def compute_X_ice_slope(M_mantle=None, X_ice=None, X_ice_0=None):
    """
    The slope of the ice fraction in the mantle is given by the mean ice content
    and the ice content at the bottom of the mantle.
    """
    
    #Use linear dependency of the ice mass fraction on the mass
    slope = 2.*(X_ice - X_ice_0)/M_mantle
    
    return slope


def filter_data_convergence(data, accs = [1.0e-4, 1.0e-3, 1.0e-2]):
    """
    Takes a data set as input an filters it using the given accuracies. All
    data points which did not converge according to the criterion will be
    rejected. The parameters which will be checked are:
        
        total mass
        total magnesium number
        ocean mass fraction
        
    """
    for i in range(len(data)):
        M_surface_is = data[i][7]
        M_surface_should = data[i][11]
        
        Mg_number_is = data[i][3]
        Mg_number_should = data[i][8]
        
        ocean_frac_is = data[i][5]
        ocean_frac_should = data[i][10]
        
        reldev_Mg = (Mg_number_is-Mg_number_should)/Mg_number_should
        reldev_mass = (M_surface_is-M_surface_should)/M_surface_should
        reldev_ocean = (ocean_frac_is-ocean_frac_should)/ocean_frac_should

        if ocean_frac_should < .01:
            reldev_ocean = 0.

        if abs(reldev_mass) > accs[0] \
            or abs(reldev_Mg) > accs[1]\
            or abs(reldev_ocean) > accs[2]:
                print ('----')
                print ('not converged')
                print ('reldevs =', reldev_mass, reldev_Mg, reldev_ocean)
                data[i][:] = np.nan
                
    return data


def filter_data_match(data, obj='titan', filter_type=0, acc_R = 1.0e-2,
                      acc_MoI = 2.0e-2):
    """Checks each point of the data set if it matches the properties of the
    given object with the desired precision. The parameters that are probed
    are:
        
        MoI factor
        J2 coefficient
        total Radius
        
    """
    acc = np.sqrt(acc_R**2 + acc_MoI**2)
    
    #Extract the values of the real object
    data_is = [MoI_solar[obj], J2_solar[obj], R_solar[obj]]

    for i in range(len(data)):
        try:
            reldev_R = (data[i][2] - data_is[2])/data_is[2]
        except TypeError:
            reldev_R = 1.0e10
            
        try:
            reldev_MoI = (data[i][0] - data_is[0])/data_is[0]
        except TypeError:
            reldev_MoI = 1.0e10
        
        try:
            reldev_J2 = (data[i][1]*1.0e6 - data_is[1])/data_is[1]
        except TypeError:
            reldev_J2 = 1.0e10
        
        #The total radius and MoI factor are combined to determine a match
        if filter_type == 0:
            reldev = np.sqrt(reldev_R**2 + reldev_MoI**2)
            
            #If the data point does not match the real values, eject it
            if reldev > acc:
                data[i][:] = np.nan
        
        #Total radius and MoI factor are used independently to determine a match
        elif filter_type == 1:
            if abs(reldev_R) < acc_R and abs(reldev_MoI) < acc_MoI:
                pass
            
            else:
                data[i][:] = np.nan

    return data


def run1(Mg_number_should = 0.5, eps_H2O = 0.,
        obj = 'titan'):
    
    M_surface_should = M_solar[obj]
    
    data = np.empty([len(ocean_mass_fractions), len(impurity_abundances), 3])
    
    suff = ftool.param_to_label(params = [Mg_number_should, eps_H2O, M_surface_should],
                                types = ['float', 'float', 'float'],
                                rounds = [1, 2, 4],
                                labels = ['Mg', 'eps_H2O', 'M_surface'])
    
    for i in range(len(ocean_mass_fractions)):
        for j in range(len(impurity_abundances)):
            of = np.log10(ocean_mass_fractions[i])
            xi = impurity_abundances[j]
            
            pl=toolkit.model_hydro(P_center=9.5e10, 
                            match_mass=True, 
                            eps_r=0.25, 
                            predictor_T='none', 
                            ocean=True, 
                            ocean_frac_should=of, 
                            temp_jumps=[0., 0.*1400.*M_surface_should**.75, 
                                        0., 0., 0.], 
                            Si_number_should=1./3., 
                            P_surface_should=1.0e5, 
                            T_surface_should=300., 
                            T_center=4000.,  
                            Fe_number_mantle=0.0, 
                            Mg_number_should=Mg_number_should, 
                            eps_H2O=eps_H2O, 
                            iterationType=1, 
                            M_surface_should=M_surface_should, 
                            predictor_P='linear', 
                            log=False, 
                            sweeps=10, 
                            iteration_limit=25, 
                            subphase_res=32, 
                            xi_Stv=0., 
                            acc_Mg_number=1.0e-3, 
                            acc_T_surface=1.0e-2, 
                            acc_M_surface=1.0e-4, 
                            xi_impurity=xi)
            
            #Compute kf Love number
            kf = MOI_model.k_f(pl.MOI_is)
            
            #Compute J2 coefficient
            omega = 2.*np.pi/(orbit_solar[obj]*day)
            
            J2 = MOI_model.J_2(omega, pl.R_surface_is, pl.M_surface_is, kf)
            
            data[i][j][0] = pl.MOI_is
            data[i][j][1] = J2
            data[i][j][2] = pl.R_surface_is/r_earth


def run2(Mg_number = [.5, .9999], Si_number = [1./3., .5], obj = 'titan', 
         res = 3, impurity_abundances = [0., .1, .2, .3, .4, .5], 
         ocean_mass_fractions = [0.0001]):
    """
    Create sample of models for specified object for different compositions
    and no hydration.
    

    Parameters
    ----------
    Mg_number : TYPE, optional
        DESCRIPTION. Bulk Mg content. The default is [.5, .9999].
    Si_number : TYPE, optional
        DESCRIPTION. Bulk Si content. The default is [1./3., .5].
    obj : TYPE, optional
        DESCRIPTION. Name of the object. The default is 'titan'.
    res : TYPE, optional
        DESCRIPTION. Defines the number of models. The default is 3.
    impurity_abundances : TYPE, optional
        DESCRIPTION. Content of water ice mixed into the mantle.
        The default is [0., .1, .2, .3, .4, .5].
    ocean_mass_fractions : TYPE, optional
        DESCRIPTION. Amount of water in the surface layer.
        The default is [0.0001].

    Returns
    -------
    None.

    """
    N = 2**res
    
    eps_H2O = 0.
    
    M_surface_should = M_solar[obj]
    
    data = np.empty([len(ocean_mass_fractions), len(impurity_abundances), 12, N, N])
    
    suff = '_run2'
    
    print ('N =', N)
    print ('shape =', data.shape)
    
    Mg_numbers = np.linspace(Mg_number[0], Mg_number[1], N)
    Si_numbers = np.linspace(Si_number[0], Si_number[1], N)
    
    for m in range(N):
        for n in range(N):
    
            for i in range(len(ocean_mass_fractions)):
                for j in range(len(impurity_abundances)):
                    of = np.log10(ocean_mass_fractions[i])
                    xi = impurity_abundances[j]
                    
                    pl=toolkit.model_hydro(P_center=9.5e10, 
                                    match_mass=True, 
                                    eps_r=0.25, 
                                    predictor_T='none', 
                                    ocean=True, 
                                    ocean_frac_should=of, 
                                    temp_jumps=[0., 0.*1400.*M_surface_should**.75, 
                                                0., 0., 0.], 
                                    Si_number_should=Si_numbers[m], 
                                    P_surface_should=1.0e5, 
                                    T_surface_should=max(T_solar[obj], 100), 
                                    T_center=4000.,  
                                    Fe_number_mantle=0.0, 
                                    Mg_number_should=Mg_numbers[n], 
                                    eps_H2O=eps_H2O, 
                                    iterationType=0, 
                                    M_surface_should=M_surface_should, 
                                    predictor_P='linear', 
                                    log=False, 
                                    sweeps=10, 
                                    iteration_limit=25, 
                                    subphase_res=32, 
                                    xi_Stv=0., 
                                    acc_Mg_number=1.0e-3, 
                                    acc_T_surface=1.0e-2, 
                                    acc_M_surface=1.0e-4, 
                                    xi_impurity=xi)
                    
                    #Compute kf Love number
                    kf = MOI_model.k_f(pl.MOI_is)
                    
                    #Compute angular velocity
                    omega = 2.*np.pi/(orbit_solar[obj]*day)
                    
                    #Compute J2 gravity coefficient
                    J2 = MOI_model.J_2(omega, pl.R_surface_is, pl.M_surface_is, kf)
                    
                    data[i][j][0][m][n] = pl.MOI_is
                    data[i][j][1][m][n] = J2
                    data[i][j][2][m][n] = pl.R_surface_is/r_earth
                    data[i][j][3][m][n] = pl.Mg_number_is
                    data[i][j][4][m][n] = pl.Si_number_is/(1.-pl.Si_number_is)
                    data[i][j][5][m][n] = pl.ocean_frac_is
                    data[i][j][6][m][n] = xi
                    data[i][j][7][m][n] = pl.M_surface_is/m_earth
                    data[i][j][8][m][n] = pl.Mg_number_should
                    data[i][j][9][m][n] = pl.Si_number_should
                    data[i][j][10][m][n] = pl.ocean_frac_should
                    data[i][j][11][m][n] = pl.M_surface_should/m_earth
    
    
    with io.open(data_dir+obj+suff+'.txt', 'w') as outfile:
        print ('writing to file:', data_dir+obj+suff+'.txt')
        outfile.write('#Array shape: {0}\n'.format(data.shape))
        print ('writing: ', '#Array shape: {0}\n'.format(data.shape))
        for s in data:
            for ss in s:
                for sss in ss:
                    for ssss in sss:

                        np.savetxt(outfile, ssss)
                        
                        outfile.write('# New slice \n')                
                    
                    outfile.write('# New slice \n')    
    
    #output = open(data_dir + 'titan'+suff+'.pkl', 'wb')
    #pickle.dump(data, output)
    #output.close()


def draw(ranges = None):
    vals = []
    
    if ranges == None:
        mids = [ .5*(hr[0]+hr[1]) for hr in hard_ranges]
        ranges = [[m, m] for m in mids]
        print ('ranges =', ranges)
    
    for i in range(len(ranges)):
        
        #Decide whether to draw from above or below
        p = random.random()
        
        #above
        if p > .5:
            interval = [ranges[i][1], hard_ranges[i][1]]
            
        #below
        else:
            interval = [ranges[i][0], hard_ranges[i][0]]            
        
        p = random.random()
        delta = interval[1] - interval[0]
        
        val = interval[0] + p*delta
        vals.append(val)
        
    return vals

        
def run3(obj = 'titan', res = 3, iterations = 2, sampler=0,
         ranges = hard_ranges):
    
    N = 2**res
    eps_H2O = 0.
    M_surface_should = M_solar[obj]
    data = np.empty([N, 15])
    random_vals = np.empty([len(ranges)])    
    
    for it in range(iterations):
        for i in range(N):
    
            #Draw parameters
            if sampler == 0:
                for j in range(len(ranges)):
                    
                    p = random.random()
                    delta = (ranges[j][1]-ranges[j][0])
                    
                    val = ranges[j][0] + p*delta
                    random_vals[j] = val
                    
                #The maximum Mg# is limited by (1-xi_Fe) in the mantle
                p = random.random()
                delta = (ranges[0][1]-random_vals[4] - ranges[0][0])

                val = ranges[0][0] + p*delta
                random_vals[0] = val
                
                #The maximum ice mas4.055489476149217e+22s fraction in the mantle is limited by
                #the total amount of water
        
                X_ice_max, core_frac = compute_max_X_ice(Si_number = random_vals[1],
                                              xi_Fe = random_vals[4],
                                              Mg_number = random_vals[0],
                                              H2O_frac = random_vals[2],
                                              xi_S = random_vals[3])
                
                p = random.random()
                delta = (min(X_ice_max, ranges[5][1]) - ranges[5][0])
                
                val = ranges[5][0] + p*delta
                random_vals[5] = val

                #The maximum ocean mass fraction is limited by the ice content
                #in the mantle
                ocean_frac_max = compute_max_ocean_frac(
                                                        Si_number=random_vals[1],
                                                        xi_Fe = random_vals[4],
                                                        Mg_number = random_vals[0],
                                                        X_ice = random_vals[5],
                                                        H2O_frac = random_vals[2],
                                                        xi_S = random_vals[3])
                
                p = random.random()
                print ('max ocean frac =', ocean_frac_max)
                delta = (min(ocean_frac_max, ranges[2][1]) - ranges[2][0])

                val = ranges[2][0] + p*delta
                random_vals[2] = val                
                
                #Scale down X_ice_0 to be between 0 and X_ice
                random_vals[6] *= random_vals[5]
                
            #Convert wt fraction of water into mole fraction as model input
            xi_ice = convert_X_ice_to_xi_ice(Si_number = random_vals[1],
                                             xi_Fe = random_vals[4],
                                             X_ice = random_vals[5],
                                             contents = [6,7,1])
            
            xi_ice_0 = convert_X_ice_to_xi_ice(Si_number = random_vals[1],
                                             xi_Fe = random_vals[4],
                                             X_ice = random_vals[6],
                                             contents = [6,7,1])
            
            #Now that all the layers are defined the mantle mass can be computed
            M_mantle = (1. - core_frac - random_vals[2])*M_surface_should*m_earth
            
            #From this the slope of the ice mass fraction in the mantle is computed
            X_ice_slope = compute_X_ice_slope(X_ice = random_vals[5],
                                                X_ice_0 = random_vals[6],
                                                M_mantle = M_mantle)
           
            X_impurity_0_layers = [0., 0., 0., random_vals[6], 0.]
            
            pl=toolkit.model_hydro(P_center=9.5e10, 
                            match_mass=True, 
                            eps_r=0.25, 
                            predictor_T='none', 
                            ocean=True, 
                            ocean_frac_should=np.log10(random_vals[2]), 
                            temp_jumps=[0., 0.*1400.*M_surface_should**.75, 
                                        0., 0., 0.], 
                            Si_number_should=random_vals[1], 
                            P_surface_should=1.0e5, 
                            T_surface_should=max(T_solar[obj], 100), 
                            T_center=4000.,  
                            Fe_number_mantle=random_vals[4], 
                            Mg_number_should=random_vals[0], 
                            eps_H2O=eps_H2O, 
                            iterationType=0, 
                            M_surface_should=M_surface_should, 
                            predictor_P='linear', 
                            log=False, 
                            sweeps=10, 
                            iteration_limit=25, 
                            subphase_res=32, 
                            xi_Stv=0., 
                            acc_Mg_number=1.0e-3, 
                            acc_T_surface=1.0e-2, 
                            acc_M_surface=1.0e-4, 
                            X_impurity=random_vals[5],
                            X_impurity_0_layers = X_impurity_0_layers,
                            xi_FeS = random_vals[3])

            print ('core is =', pl.M_core_is*m_earth)
            print ('core =', core_frac*M_surface_should*m_earth)
            print ('mantle is =', pl.layer_properties[3]['indigenous_mass'])
            print ('mantle =', M_mantle)
            print ('X_ice =', random_vals[5])
            print ('X_ice_0 =', random_vals[6])
            print ('xi_ice_0 =', xi_ice_0)
            print ('xi_ice =', xi_ice)
            print ('Mg_number = ', random_vals[0])
            print ('Si_number = ', random_vals[1])
            print ('xi_Fe =', random_vals[4])
            print ('xi_S =', random_vals[3])
            print ('ocean =', random_vals [2])
            print ('X_ice_max =', X_ice_max)
            print ('predict =', random_vals[6] + X_ice_slope*M_mantle)

            #Compute kf Love number
            kf = MOI_model.k_f(pl.MOI_is)
            
            #Compute angular velocity
            omega = 2.*np.pi/(orbit_solar[obj]*day)
            
            #Compute J2 gravity coefficient
            J2 = MOI_model.J_2(omega, pl.R_surface_is, pl.M_surface_is, kf)
            
            #Compute mantle mass
            M_mantle = pl.layer_properties[3]['indigenous_mass']+\
                        pl.layer_properties[2]['indigenous_mass']
                        
            #Compute total mass of ice in the mantle
            M_ice_mantle = M_mantle*random_vals[5]
            
            data[i][0] = pl.MOI_is
            data[i][1] = J2
            data[i][2] = pl.R_surface_is/r_earth
            data[i][3] = pl.Mg_number_is
            data[i][4] = pl.Si_number_is#/(1.-pl.Si_number_is)
            data[i][5] = pl.M_H2O_is
            data[i][6] = random_vals[5] #X_ice
            data[i][7] = pl.M_surface_is/m_earth
            data[i][8] = pl.Mg_number_should
            data[i][9] = pl.Si_number_should
            data[i][10] = random_vals[2] + M_ice_mantle/pl.M_surface_is
            data[i][11] = pl.M_surface_should/m_earth
            data[i][12] = random_vals[4] #xi_Fe
            data[i][13] = pl.xi_S_core
            data[i][14] = (random_vals[5]-random_vals[6])/random_vals[5]
            
            #ranges = extract_parameter_ranges(data, obj=obj, acc=acc)
            
    return data


def extract_parameter_ranges(data, obj='titan', acc_R=1.0e-2,
                             acc_MoI = 2.0e-2, filter_type = 0):
    data = filter_data_convergence(data)
    data = filter_data_match(data, obj=obj, acc_R=acc_R, acc_MoI = acc_MoI,
                             filter_type = filter_type)

    ranges_indices = [3, 4, 5, 13, 12, 6, 14]
    ranges = np.empty([len(hard_ranges), 2])
    
    for i in range(len(hard_ranges)):
        ind = ranges_indices[i]
        ranges[i][0] = np.nanmin(data.T[ind])
        ranges[i][1] = np.nanmax(data.T[ind])
        
    return ranges


def run_sampler(obj = 'titan', 
                res = 3, 
                sampler=0, 
                ranges = [[.5, 1.], [1./3., .5], [.0, .75]],
                acc_R = 1.0e-2,
                acc_MoI = 2.0e-2,
                **kwargs
                ):
    
    N = 2**res
        
    M_surface_should = M_solar[obj]
    
    data = np.empty([N, 14])
        
    random_vals = np.empty([len(ranges)])    
    
    for i in range(N):

        #Draw parameters
        if sampler == 0:
            for j in range(len(ranges)):
                #Decide if min or max is sampled
                p1 = random.random()
                p2 = random.random()
                
                #max
                if p1 > .5:
                    print ('check1')
                    q0 = ranges[j][1]
                    delta = hard_ranges[j][1] - q0
                    
                    val = q0 + p2*delta
                
                #min
                else:
                    print ('check1')
                    q0 = ranges[j][0]
                    delta = hard_ranges[j][0] - q0
                    
                    val = q0 + p2*delta
                                    
                random_vals[j] = val
        
        print ('random vals = ', random_vals)
        
        pl=toolkit.model_hydro(P_center=9.5e10, 
                        match_mass=True, 
                        eps_r=0.25, 
                        predictor_T='none', 
                        ocean=True, 
                        ocean_frac_should=np.log10(random_vals[2]), 
                        temp_jumps=[0., 0.*1400.*M_surface_should**.75, 
                                    0., 0., 0.], 
                        Si_number_should=random_vals[1], 
                        P_surface_should=1.0e5, 
                        T_surface_should=max(T_solar[obj], 100), 
                        T_center=4000.,  
                        Fe_number_mantle=.0, 
                        Mg_number_should=random_vals[0], 
                        eps_H2O=.0, 
                        iterationType=0, 
                        M_surface_should=M_surface_should, 
                        predictor_P='linear', 
                        log=False, 
                        sweeps=10, 
                        iteration_limit=25, 
                        subphase_res=32, 
                        xi_Stv=0., 
                        acc_Mg_number=1.0e-3, 
                        acc_T_surface=1.0e-2, 
                        acc_M_surface=1.0e-4, 
                        xi_impurity = .0,
                        xi_FeS = .0)
        
        #Compute kf Love number
        kf = MOI_model.k_f(pl.MOI_is)
        
        #Compute angular velocity
        omega = 2.*np.pi/(orbit_solar[obj]*day)
        
        #Compute J2 gravity coefficient
        J2 = MOI_model.J_2(omega, pl.R_surface_is, pl.M_surface_is, kf)
        
        data[i][0] = pl.MOI_is
        data[i][1] = J2
        data[i][2] = pl.R_surface_is/r_earth
        data[i][3] = pl.Mg_number_is
        data[i][4] = pl.Si_number_is#/(1.-pl.Si_number_is)
        data[i][5] = pl.ocean_frac_is
        data[i][6] = 0.
        data[i][7] = pl.M_surface_is/m_earth
        data[i][8] = pl.Mg_number_should
        data[i][9] = pl.Si_number_should
        data[i][10] = pl.ocean_frac_should
        data[i][11] = pl.M_surface_should/m_earth
        data[i][12] = 0.
        data[i][13] = 0.
        
    data = filter_data_convergence(data)
    data = filter_data_match(data, obj=obj, acc_R = acc_R, acc_MoI = acc_MoI)
    
    ranges_indices = [3, 4, 5]
    model_ranges = np.empty([3, 2])
    
    #Estimate maximal parameter range based on the data
    for i in range(3):
        ind = ranges_indices[i]
        model_ranges[i][0] = np.nanmin(data.T[ind])
        model_ranges[i][1] = np.nanmax(data.T[ind])

    return data, model_ranges


def check_convergence(obj='titan', N_points = 4, acc_R=1.0e-2, acc_MoI = 2.0e-2,
                      res = 1, write = False, read = True, filter_type = 0):
    all_ranges = []
    ranges = hard_ranges
    
    data = np.empty([0, 15])
    
    t0 = time.time()
    
    for i in range(N_points):
        file = obj+'_output_'+str(i+1)+'.npy'
        if read:
            data = np.load(data_dir + obj.capitalize()+'/'+ file)
        
        else:
            add_data = run3(obj=obj, res=res+i, iterations=1, ranges = hard_ranges)
            
            data = np.append(data, add_data, axis=0)
            
            #write raw data to file before checking for convergence and filtering
            if write:
                np.save(data_dir + obj.capitalize()+'/'+ file, data)
        
        #Extract the maximum parameter range for the current run
        ranges = extract_parameter_ranges(data, obj=obj, acc_R=acc_R,
                                          acc_MoI = acc_MoI, 
                                          filter_type = filter_type)
        all_ranges.append(ranges)
        
    all_ranges = np.array(all_ranges)
    
    fig, ax = plt.subplots(len(hard_ranges), 1, figsize=(6, 8))
    plt.subplots_adjust(hspace=.5)
    
    #Convert from logspace to ocean mass fraction
    all_ranges.T[0][2] = all_ranges.T[0][2]
    all_ranges.T[1][2] = all_ranges.T[1][2]
    
    y_labels = [r'$\rm Mg \#$', r'$\rm Si \#$', r'$M_{\rm H_2 O}/M$', 
                r'$\xi_{\rm S}$', r'$\xi_{\rm Fe}$', r'${\bar{X}_{\rm Ice}}$',
                r'$\frac{{\bar{X}_{\rm Ice}}-X_{\rm Ice, 0}}{\bar{X}_{\rm Ice}}$']
    
    ax[-1].set_xlabel(r'$N \ \rm runs$')
    
    for i in range(len(hard_ranges)):
        ax[i].plot(all_ranges.T[0][i])
        ax[i].plot(all_ranges.T[1][i])
        
        ax[i].set_ylabel(y_labels[i])
        #Force integer ticks for x-axis
        ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))
    
    fig.align_ylabels(ax[:])
        
    bbox_props = dict(edgecolor='k', facecolor='white',
                  linewidth = .5)

    fig.savefig(data_dir+obj.capitalize()+'/'+obj+'_convergence.pdf',
                    format='pdf',
                    bbox_inches='tight')    
    plt.close(fig)
    
    t = time.time()
    ftool.printTime(t-t0, ms=False, where = 'check_convergence')
    
    return all_ranges


def MCMC(obj='titan', 
         N_points = 4, 
         ranges = [[1., .5], [.5, 1./3.], [.75, .00001]],
         res=3,
         **kwargs
         ):
    
    all_ranges = []
    
    for i in range(N_points):
        
        data, ranges = run_sampler(obj=obj, res=res, ranges=ranges, acc=1.0e-1)
        print ('\n ranges =', ranges)
        
        ranges[2][0] = 10**ranges[2][0]
        ranges[2][1] = 10**ranges[2][1]
        
        all_ranges.append(ranges)
    
    all_ranges = np.array(all_ranges)
    
    fig, ax = plt.subplots(3, 1)

    for i in range(3):
        ax[i].plot(all_ranges.T[0][i])
        ax[i].plot(all_ranges.T[1][i])

    bbox_props = dict(edgecolor='k', facecolor='white',
                  linewidth = .5)

    fig.savefig(data_dir+obj.capitalize()+'/'+obj+'_convergence.pdf',
                    format='pdf',
                    bbox_inches='tight')    
    plt.close(fig)
    
    return all_ranges


def plot1(file, Mg = 1.0, obj='titan'):
    
    data, shape = ftool.read_data_file(file)
                      
    data = np.asarray(data).reshape(shape)
    
    data_is = [MoI_solar[obj], J2_solar[obj], R_solar[obj]]

    
    fig, ax = plt.subplots(1, 3, figsize = (10, 3))
    fig.subplots_adjust(hspace=.2, wspace=.4) 
    
    #Compute y-axis limits for each parameter from the data
    y_lims = np.empty([3, 2])
    
    scalings = [1., 1.0e6, 1.]

    for i in range(3):
        low = min(np.min(data.T[i])*scalings[i], data_is[i])
        up = max(np.max(data.T[i])*scalings[i], data_is[i])
        
        delta = .1*(up-low)
        
        low -= delta
        up += delta
        
        y_lims[i][0] = low
        y_lims[i][1] = up
        
        ax[i].set_ylim(y_lims[i][0], y_lims[i][1])
        
    print (y_lims)
    
    
    for j in range(len(impurity_abundances)):
        col = color_list[j]
        
        ax[0].plot(ocean_mass_fractions, data.T[0][j], color=col,
                   label=r'$\xi_{\rm Ice}= \ $'+\
                       str(round(impurity_abundances[j], 2)))
        ax[1].plot(ocean_mass_fractions, data.T[1][j]*scalings[1], color=col)
        ax[2].plot(ocean_mass_fractions, data.T[2][j], color=col)
        
    
    J2 = J2_solar[obj]
    dJ2 = delta_J2_solar[obj]
    
    MoI = MoI_solar[obj]
    dMoI = delta_MoI_solar[obj]
    
    #Plot MoI factor, J2 and radius of Titan
    ax[0].plot([.0, .6], [MoI_solar[obj], MoI_solar[obj]], color='grey')
    ax[1].plot([.0, .6], [J2, J2], color='grey')
    ax[2].plot([.0, .6], [R_solar[obj], R_solar[obj]], color='grey')

    bbox_props = dict(edgecolor='None', facecolor='white',
                  linewidth = .5, pad=.1, alpha=1.)
    
    ax[0].text(.1, MoI, str(MoI), color='grey',
               bbox=bbox_props, va = 'center', fontsize=fnts1)

    ax[1].text(.1, J2, str(J2), color='grey',
               bbox=bbox_props, va = 'center', fontsize=fnts1)

    ax[2].text(.1, R_solar[obj], str(R_solar[obj]), color='grey',
               bbox=bbox_props, va = 'center', fontsize=fnts1)

    for i in range(len(ax)):
        ax[i].tick_params(which='both', right=True, top=True, direction='in')
        ax[i].set_xlim(0., .6)
        
        ax[i].set_xlabel(r'$\rm Ocean \ mass \ fraction$', fontsize=fnts2)
        
    ax[0].set_ylabel(r'$I/(MR^{2})$', fontsize=fnts2)
    ax[1].set_ylabel(r'$J_2 \ \rm [10^{-6}]$', fontsize=fnts2)
    ax[2].set_ylabel(r'$R/R_\oplus$', fontsize=fnts2)
    
    ypos = (y_lims[-1][0]+y_lims[-1][1])/2.
    
    ax[-1].text(ocean_mass_fractions[-1]*1.1, ypos,
                r'$\rm Mg\# ={}$'.format(Mg), 
                rotation=90, 
                fontsize=fnts3,
                va = 'center',
                color='grey',
                weight = 'bold')
    
    bbox_props = dict(edgecolor='k', facecolor='white',
                  linewidth = .5)
    
    ax[2].text(ocean_mass_fractions[-1], np.max(data.T[-1])*1.05, 
               r'$\rm {}$'.format(obj.capitalize()), ha='center',
                  bbox=bbox_props, fontsize=fnts4)
    
    ax[0].legend(fontsize=fnts1, loc = 1)
    
    fig.savefig(file.split('.')[0]+'.pdf', format='pdf', bbox_inches='tight')
    plt.close(fig)
    

def plot2(obj='titan', analyze=True):
    """
    Plot relative difference from observed properties to modeled properties
    for different boundary conditions.

    Parameters
    ----------
    obj : TYPE, optional
        DESCRIPTION. Name of the object under consideration. Possibilities
        are the main solar system bodies. The default is 'titan'.
    analyze : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    """
    edge_colors = np.array([color_list[2], [1,1,1], color_list[0]])
    
    N_col = 101
    
    colors = np.empty([N_col, 3])
    
    colors[0][:] = edge_colors[0][:]
    
    for i in range(N_col-1):
        for j in range(3):
                        
            if i < int(N_col/2):
                ind = [0, 1]
                
            else:
                ind = [1, 2]
            
            grad = np.array(2*(edge_colors[ind[1]][j] - \
                             edge_colors[ind[0]][j])/(N_col-1))
                
            colors[i+1][j] = max(colors[i][j] + grad, .0)
            colors[i+1][j] = min(colors[i+1][j], 1.)
        
    newcmap = mpl.colors.ListedColormap(colors)
    
    file = obj+'_run2.txt'
    
    data, shape = ftool.read_data_file(file)
    print ('shape =', shape)              

    data = np.asarray(data).reshape(shape)        
    
    data_is = [MoI_solar[obj], J2_solar[obj], R_solar[obj]]
    print ("shape =", np.shape(data))
    rows = len(data[0])
    plots = len(data)
    print ("rows =", rows)
    print ("plots =", plots)
    cbar_labels = [r'$\delta I /I$', r'$\delta{J_2}/J_2$', r'$\delta R/R$']
    scalings = [1., 1.0e6, 1.]
    
    vmins = np.empty([3, rows, plots])
    vmaxs = np.empty([3, rows, plots])    
    
    accs = [1.0e-4, 1.0e-3, 1.0e-2]
    
    #Check convergence
    for p in range(plots):
        for r in range(rows):
            for c in range(len(data[p][r])):
                for d in range(4):
                    for i in range(len(data[p][r][0])):
                        for j in range(len(data[p][r][0][i])):
                            M_surface_is = data[p][r][7][i][j]
                            M_surface_should = data[p][r][11][i][j]
                            
                            Mg_number_is = data[p][r][3][i][j]
                            Mg_number_should = data[p][r][8][i][j]
                            
                            ocean_frac_is = data[p][r][5][i][j]
                            ocean_frac_should = data[p][r][10][i][j]

                            reldev_Mg = (Mg_number_is-Mg_number_should)/Mg_number_should
                            reldev_mass = (M_surface_is-M_surface_should)/M_surface_should
                            reldev_ocean = (ocean_frac_is-ocean_frac_should)/ocean_frac_should
                            
                            J2 = data[p][r][1][i][j]*scalings[1]
                            
                            if ocean_frac_should < -2:
                                reldev_ocean = 0.

                            if abs(reldev_mass) > accs[0] \
                                or abs(reldev_Mg) > accs[1]\
                                or abs(reldev_ocean) > accs[2]:
                                    data[p][r][0:][i][j] = None
                        
    
    #Extract value range for each parameter
    for p in range(plots):
        for r in range(rows):
            for i in range(3):
                vals = (data[p][r][i]*scalings[i]-data_is[i])/data_is[i]
                
                vmins[i][r][p] = np.min(vals)
                vmaxs[i][r][p] = np.max(vals)
                
                extrem = max(abs(vmins[i][r][p]), abs(vmaxs[i][r][p]))
                
                vmins[i][r][p] = -min(extrem, .3)
                vmaxs[i][r][p] = min(extrem, .3)
    
    for p in range(plots):
        fig, ax = plt.subplots(rows, 3, figsize=(12, 8))
        fig.subplots_adjust(wspace=.0, hspace=.5)
        
        aspect = 20
        pad_fraction = 10.
        digits = 3
        
        ticks = np.linspace(0, 1, 3)
        xtick_labels = [.5, .75, 1.]#np.linspace(round(np.min(data[0][0][3]), digits), 
                                   #round(np.max(data[0][0][3]), digits), 3)
        ytick_labels = [.5, .75, 1.]#np.linspace(round(np.min(data[0][0][4]), digits), 
                                   #round(np.max(data[0][0][4]), digits), 3)
    
        bbox_props = dict(edgecolor='k', facecolor='white',
                      linewidth = .5)    
        
        ax[0][-1].text(1., 1.2, r'$\rm {}$'.format(obj.capitalize())+\
                       r'$: \ M_{\rm Ocean}/M= \ $'+\
                           r'${}$'.format(round(10**data[p][0][5][0][0],2)), 
                       ha='center', 
                       bbox=bbox_props, fontsize=fnts4)

        for r in range(rows):
            xi_ice = data[p][r][6][0][0]
            ax[r][-1].text(1.75, .5, r'$\xi_{\rm Ice}= \ $'+r'${}$'.format(xi_ice), 
                           color='grey', rotation=90, va='center', fontsize=fnts3)
            
            for i in range(3):
                extrem = max(abs(np.min(vmins[i])), 
                             abs(np.max(vmaxs[i])))
                
                vmin = -extrem
                vmax = extrem
                
                vals = (data[p][r][i]*scalings[i]-data_is[i])/data_is[i]
                
                newcmap = ftool.my_dev_cmap()
                
                im = ax[r][i].imshow(vals, 
                                  aspect='equal', 
                                  origin='lower',
                                  extent = [0, 1, 0, 1],
                                  cmap =newcmap,
                                  vmin = vmin, 
                                  vmax = vmax)
                
                ax[r][i].set_xticks(ticks)
                ax[r][i].set_yticks(ticks)
                ax[r][i].set_xticklabels(xtick_labels)
                ax[r][i].set_yticklabels(ytick_labels)
                
                ax[r][i].set_xlabel(r'$\rm Mg\#$')
                ax[r][i].set_ylabel(r'$\rm Si/Mg$')
                
                ftool.my_colorbar(im, pad=.1, label = cbar_labels[i])
        
        fig.savefig('./'+obj.capitalize()+'/'+file.split('.')[0]+'_'+str(p+1)+'.pdf',
                    format='pdf', 
                    bbox_inches='tight')
        plt.close(fig)        
        
    return data


def plot3(obj = 'titan'):
    
    data = np.load(data_dir+obj + '_run3.npy')
    print (len(data))
    data_is = [MoI_solar[obj], J2_solar[obj], R_solar[obj]]
    
    xmin = [0., 0.5]
    xmax = [.75, 1.]
    
    z = data.T[4]
    
    x_labels = [r'$M_{\rm Ocean}/M$', r'$\rm Mg \#$']
    y_labels = [r'$I/MR^2$', r'$J_2 \ [10^{-6}]$', r'$R/R_\oplus$']
    
    fig, ax = plt.subplots(2, 3, figsize=(10, 4))
    plt.subplots_adjust(wspace=.5, hspace=.5)
    
    cmap = ftool.my_cmap(res=4, start=8)
    
    for i in range(len(ax)):
        for j in range(len(ax[i])):
            ax[i][j].set_xlabel(x_labels[i])
            ax[i][j].set_ylabel(y_labels[j])
            ax[i][j].set_xlim(xmin[i], xmax[i])
            
            ax[i][j].plot([0., 1.], [data_is[j], data_is[j]], color='k')
            
    ax[0][0].scatter(10**data.T[5], data.T[0], c=z, cmap=cmap)
    ax[0][1].scatter(10**data.T[5], data.T[1]*1.0e6, c=z, cmap=cmap)    
    ax[0][2].scatter(10**data.T[5], data.T[2], c=z, cmap=cmap)    

    ax[1][0].scatter(data.T[3], data.T[0], c=z, cmap=cmap)
    ax[1][1].scatter(data.T[3], data.T[1]*1.0e6, c=z, cmap=cmap)    
    sc = ax[1][2].scatter(data.T[3], data.T[2], c=z, cmap=cmap)
    
    bbox_props = dict(edgecolor='k', facecolor='white',
                  linewidth = .5)    
    
    ax[0][-1].text(.5, .35, r'$\rm {}$'.format(obj.capitalize()), 
                   ha='center', 
                   bbox=bbox_props, 
                   fontsize=fnts4)    
    
    colorbar = fig.colorbar(sc, ax=ax, label=r'$\xi_{\rm Ice}$')
    
    fig.savefig('./'+obj.capitalize()+'/MCMC.pdf',
                    format='pdf', 
                    bbox_inches='tight')
    plt.close(fig)
    

def analyze3(obj):
    data = np.load(data_dir+obj + '_run3.npy')
    data_is = [MoI_solar[obj], J2_solar[obj], R_solar[obj]]
    
    indices = []
    accs = [1.0e-4, 1.0e-3, 1.0e-2]

    for i in range(len(data)):
        M_surface_is = data[i][7]
        M_surface_should = data[i][11]
        
        Mg_number_is = data[i][3]
        Mg_number_should = data[i][8]
        
        ocean_frac_is = data[i][5]
        ocean_frac_should = data[i][10]

        reldev_Mg = (Mg_number_is-Mg_number_should)/Mg_number_should
        reldev_mass = (M_surface_is-M_surface_should)/M_surface_should
        reldev_ocean = (ocean_frac_is-ocean_frac_should)/ocean_frac_should
        
        
        if ocean_frac_should < -2:
            reldev_ocean = 0.

        if abs(reldev_mass) > accs[0] \
            or abs(reldev_Mg) > accs[1]\
            or abs(reldev_ocean) > accs[2]:
                print ('not converged')
                data[i][:] = None
        
    
    #Find all models that match constraints with the desired precision
    for i in range(len(data)):
        try:
            reldev_R = (data[i][2] - data_is[2])/data_is[2]
        except TypeError:
            reldev_R = None
            
        try:
            reldev_MoI = (data[i][0] - data_is[0])/data_is[0]
        except TypeError:
            reldev_MoI = None
        
        try:
            reldev_J2 = (data[i][1]*1.0e6 - data_is[1])/data_is[1]
        except TypeError:
            reldev_J2 = None
            
        reldev = np.sqrt(reldev_R**2 + reldev_MoI**2)
                
        if abs(reldev_R) < 1.0e-2 \
            and abs(reldev_MoI) < 2.0e-2 \
            and abs(reldev_J2) < 5.0e-2:
            indices.append(i)
            
    matched_models = np.empty([len(indices), 3])
    
    for i in range(len(indices)):
        dat = data[indices[i]]
        
        matched_models[i][0] = 10**dat[5]
        matched_models[i][1] = dat[3]
        matched_models[i][2] = dat[4]
        
    
    print ('Number of matched models = ', len(indices))
    for i in range(3):
        print ('range = ', np.amin(matched_models.T[i]), np.amax(matched_models.T[i]))
        
        
        
        
            
    
    
 