#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 00:32:12 2021

@author: os18o068
"""

import PlanetFactory as plfac
import numpy as np
from PIMPphysicalparams import m_earth, r_earth
from matplotlib import pyplot as plt
from PIMPrunparams import color_list

toolkit = plfac.Toolkit()
data_dir = '/home/os18o068/Documents/PHD/Projects/Planets/Data/'

#M, Mg#, T_surf
accs = [[1.0e-4, 1.0e-3, 1.0e-4],
        [1.0e-4, 1.0e-3, 1.0e-2],
        [1.0e-5, 1.0e-3, 1.0e-2]]

def delta_T_CMB(M):
    return 1400.*M**.75

def generate_data():
    Mg_numbers = [.2, .7]
    temps = [500, 2000]
    Si_number = .4
    xi_Fe = .25
    
    masses = np.logspace(-1, np.log10(3.), 6)
        
    N_params = 24
    
    data = np.empty([N_params, len(temps), len(Mg_numbers), 
                     len(masses)])
    
    ocean_masses = np.empty([len(temps), len(Mg_numbers), 
                             len(masses)])
    
    for i in range(len(temps)):
        for j in range(len(Mg_numbers)):
            for k in range(len(masses)):
                pl=toolkit.model_hydro(P_center=9.5e10, 
                            match_mass=True, 
                            eps_r=0.25, 
                            predictor_T='none', 
                            ocean=False, 
                            ocean_frac_should=-1.72, 
                            temp_jumps=[0., delta_T_CMB(masses[k]), 
                                        0., 0., 0.], 
                            Si_number_should=.3999, 
                            P_surface_should=1.0e5, 
                            T_surface_should=temps[i], 
                            T_center=4000.,  
                            Fe_number_mantle=.25, 
                            Mg_number_should=Mg_numbers[j], 
                            eps_H2O=1., 
                            iterationType=1, 
                            M_surface_should=masses[k], 
                            predictor_P='linear', 
                            log=False, 
                            sweeps=10, 
                            iteration_limit=20, 
                            subphase_res=32, 
                            xi_Stv=0.,
                            impurity=1,
                            acc_M_surface = accs[0][0],
                            acc_Mg_number = accs[0][1],
                            acc_T_surface = accs[0][2])
    
                dat = pl.M_surface_is/m_earth, pl.R_surface_is/r_earth, pl.Mg_number_is, \
                pl.Si_number_is, pl.M_core_is, pl.T_surface_is, pl.P_surface_is, \
                pl.M_ocean_is, pl.M_H2O_is, pl.MOI_is, pl.M_H2O_core, pl.T_center, \
                pl.P_center, pl.ocean_depth, pl.xi_H_core, pl.xi_H_core_predicted,\
                pl.Mg_number_should, pl.M_surface_should/m_earth, pl.T_surface_should, \
                pl.ocean_frac_should, pl.ocean_frac_is, pl.P_H2_CMB, None, None

                ocean_masses[i][j][k] = dat[8]
                
                for d in range(len(dat)):
                    data[d][i][j][k] = dat[d]
                    
    np.save(data_dir+'data_hyd_2.npy', data)

    for i in range(len(temps)):
        for j in range(len(Mg_numbers)):
            for k in range(len(masses)):
                pl=toolkit.model_hydro(P_center=9.5e10, 
                            match_mass=True, 
                            eps_r=0.25, 
                            predictor_T='none', 
                            ocean=False, 
                            ocean_frac_should=-1.72, 
                            temp_jumps=[0., delta_T_CMB(masses[k]), 
                                        0., 0., 0.], 
                            Si_number_should=.3999, 
                            P_surface_should=1.0e5, 
                            T_surface_should=temps[i], 
                            T_center=4000.,  
                            Fe_number_mantle=.25, 
                            Mg_number_should=Mg_numbers[j], 
                            eps_H2O=0., 
                            iterationType=1, 
                            M_surface_should=masses[k], 
                            predictor_P='linear', 
                            log=False, 
                            sweeps=10, 
                            iteration_limit=20, 
                            subphase_res=32, 
                            xi_Stv=0.,
                            impurity=1,
                            acc_M_surface = accs[1][0],
                            acc_Mg_number = accs[1][1],
                            acc_T_surface = accs[1][2])
    
                dat = pl.M_surface_is/m_earth, pl.R_surface_is/r_earth, pl.Mg_number_is, \
                pl.Si_number_is, pl.M_core_is, pl.T_surface_is, pl.P_surface_is, \
                pl.M_ocean_is, pl.M_H2O_is, pl.MOI_is, pl.M_H2O_core, pl.T_center, \
                pl.P_center, pl.ocean_depth, pl.xi_H_core, pl.xi_H_core_predicted,\
                pl.Mg_number_should, pl.M_surface_should/m_earth, pl.T_surface_should, \
                pl.ocean_frac_should, pl.ocean_frac_is, pl.P_H2_CMB, None, None
                
                for d in range(len(dat)):
                    data[d][i][j][k] = dat[d]
                    
    np.save(data_dir+'data_dry_2.npy', data)

    for i in range(len(temps)):
        for j in range(len(Mg_numbers)):
            for k in range(len(masses)):
                
                of = np.log10(ocean_masses[i][j][k])
                
                pl=toolkit.model_hydro(P_center=9.5e10, 
                            match_mass=True, 
                            eps_r=0.25, 
                            predictor_T='none', 
                            ocean=True, 
                            ocean_frac_should=of, 
                            temp_jumps=[0., delta_T_CMB(masses[k]), 
                                        0., temps[i]-300., 0.], 
                            Si_number_should=.3999, 
                            P_surface_should=1.0e5, 
                            T_surface_should=300., 
                            T_center=4000.,  
                            Fe_number_mantle=.25, 
                            Mg_number_should=Mg_numbers[j], 
                            eps_H2O=0., 
                            iterationType=0, 
                            M_surface_should=masses[k], 
                            predictor_P='linear', 
                            log=False, 
                            sweeps=10, 
                            iteration_limit=20, 
                            subphase_res=32, 
                            xi_Stv=0.,
                            impurity=1,
                            acc_M_surface = accs[2][0],
                            acc_Mg_number = accs[2][1],
                            acc_T_surface = accs[2][2])
    
                dat = pl.M_surface_is/m_earth, pl.R_surface_is/r_earth, pl.Mg_number_is, \
                pl.Si_number_is, pl.M_core_is, pl.T_surface_is, pl.P_surface_is, \
                pl.M_ocean_is, pl.M_H2O_is, pl.MOI_is, pl.M_H2O_core, pl.T_center, \
                pl.P_center, pl.ocean_depth, pl.xi_H_core, pl.xi_H_core_predicted,\
                pl.Mg_number_should, pl.M_surface_should/m_earth, pl.T_surface_should, \
                pl.ocean_frac_should, pl.ocean_frac_is, pl.P_H2_CMB, None, None
                
                for d in range(len(dat)):
                    data[d][i][j][k] = dat[d]
                    
    np.save(data_dir+'data_oce_2.npy', data)
                   
    return data


def plot_data():
    
    #With mole fraction (wrong)
    dry1 = np.load('data_dry_1.npy')
    hyd1 = np.load('data_hyd_1.npy')
    oce1 = np.load('data_oce_1.npy')
    
    #With mass fraction (correct)
    dry2 = np.load('data_dry_2.npy')
    hyd2 = np.load('data_hyd_2.npy')
    oce2 = np.load('data_oce_2.npy')    

    data1 = [dry1, hyd1, oce1]
    data2 = [dry2, hyd2, oce2]
    
    styles = ['-', '--', ':']
    
    fig1, ax1 = plt.subplots(2, 2)
    fig2, ax2 = plt.subplots(2, 2)
    
    ax2[0][0].set_ylabel('rel. error on radius')
    ax2[0][1].set_ylabel('rel. error on water content')

    ax2[1][0].set_ylabel('abs. error on '+r'$\delta R_1/R_1$')
    ax2[1][1].set_ylabel('abs. error on '+r'$\delta R_2/R_2$')
    
    plot_list1 = []
    plot_list2 = []
    
    for i in range(len(dry1[0])):
        for j in range(len(dry1[0][i])):
            for d in range(len(data1)):
                dat1 = data1[d]
                dat2 = data2[d]
                
                dR = (dat1[1][i][j]-dat2[1][i][j])/dat1[1][i][j]
                dH2O = (dat1[8][i][j]-dat2[8][i][j])/dat1[8][i][j]
                
                if dR[-1] > .03:
                    print ('indices =', i, j)
                
                ax1[0][0].plot(dat1[0][i][j], dat1[1][i][j], color=color_list[j],
                           linestyle = styles[d])
                
                ax1[0][1].plot(dat1[0][i][j], dat1[8][i][j], color=color_list[j],
                           linestyle = styles[d])
                
                
                if d == 0:
                    pl1, = ax2[0][0].plot(dat1[0][i][j], dR, color=color_list[j],
                           linestyle = styles[i])
                
                if d == 1:
                    ax2[0][1].plot(dat1[0][i][j], dH2O, color=color_list[j],
                           linestyle = styles[i])

                if i == 0 and d == 0:
                    plot_list1.append(pl1)
                
                if j == 0 and d == 0:
                    plot_list2.append(pl1)
                
            dR11 = (data1[1][1][i][j]-data1[0][1][i][j])/data1[0][1][i][j]
            dR21 = (data1[2][1][i][j]-data1[1][1][i][j])/data1[1][1][i][j]
            dR12 = (data2[1][1][i][j]-data2[0][1][i][j])/data2[0][1][i][j]
            dR22 = (data2[2][1][i][j]-data2[1][1][i][j])/data2[1][1][i][j]
            
            print ('---')
            print (dR21)
            print (dR22)
            
            ddR1 = (dR11-dR12)
            ddR2 = (dR21-dR22)
            
            ax1[1][0].plot(dat1[0][i][j], dR11, color=color_list[j], marker='s')
            ax1[1][1].plot(dat1[0][i][j], dR21, color=color_list[j], marker='s')

            ax2[1][0].plot(dat1[0][i][j], ddR1, color=color_list[j],
                           linestyle = styles[i])
            ax2[1][1].plot(dat1[0][i][j], ddR2, color=color_list[j],
                           linestyle = styles[i])
    
    
    legend1 = ax2[0][0].legend(plot_list1, [r'$\rm Mg\# \ = \ 0.2$', 
                                            r'$\rm Mg\# \ = \ 0.7$'],
                               loc=1)
    ax2[0][0].add_artist(legend1)

    legend2 = ax2[0][0].legend(plot_list2, [r'$\Delta T_{\rm TBL} = 200 \ \rm K$', 
                                            r'$\Delta T_{\rm TBL} = 1700 \ \rm K$'],
                               loc=2)
    ax2[0][0].add_artist(legend2)
    
    for i in range(len(dry2[0])):
        for j in range(len(dry2[0][i])):            
            for d in range(len(data2)):
                dat = data2[d]
                ax1[0][0].plot(dat[0][i][j], dat[1][i][j], color=color_list[j],
                           linestyle = styles[d], marker='o')
                ax1[0][1].plot(dat[0][i][j], dat[8][i][j], color=color_list[j],
                           linestyle = styles[d], marker='o')
                            
            dR1 = (data2[1][1][i][j]-data2[0][1][i][j])/data2[0][1][i][j]
            dR2 = (data2[2][1][i][j]-data2[1][1][i][j])/data2[1][1][i][j]
            
            ax1[1][0].plot(dat[0][i][j], dR1, color=color_list[j], marker='o')
            ax1[1][1].plot(dat[0][i][j], dR2, color=color_list[j], marker='o')