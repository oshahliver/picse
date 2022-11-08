#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 16:45:47 2021

@author: os18o068
"""

from PIMPphysicalparams import R_solar, M_solar, MoI_solar, T_solar

def Calibrate(obj='earth'):
    
    planets = ['mercury', 'venus', 'earth', 'mars']
    
    
    for i in range(len(planets)):
        pl=toolkit.model_hydro(P_center=9.5e10, 
                    match_mass=True, 
                    eps_r=0.25, 
                    predictor_T='none', 
                    ocean=True, 
                    ocean_frac_should=-6, 
                    temp_jumps=[0., 0.*1400.*M_surface_should**.75, 
                                0., 0., 0.], 
                    Si_number_should=Si_number, 
                    P_surface_should=1.0e5, 
                    T_surface_should=300, 
                    T_center=4000.,  
                    Fe_number_mantle=.0, 
                    Mg_number_should=Mg_number, 
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
                    acc_T_surface=1.0e-3, 
                    acc_M_surface=1.0e-4, 
                    xi_impurity = .0,
                    impurity=9,
                    xi_FeS = .18)