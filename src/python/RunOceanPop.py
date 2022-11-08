#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 21:02:11 2020

@author: os18o068
"""

import PlanetFactory as plfac
import numpy as np
from PIMPphysicalparams import m_earth, r_earth
import os
import shutil
import sys
import time
import functionTools as ftool

toolkit = plfac.Toolkit()


def delta_T_CMB(M):
    return 1400.*M**.75


Mg_numbers = np.linspace(.2, .7, 6)
masses = np.logspace(-1, 1, 8)
temps = np.linspace(500., 2000., 4)
ocean_fracs = np.logspace(-4, -1, 4)

planets = []


subsubdir_names = []
overwrite = True

planet_dir = 'ocean_planets'
    #Check if planet output directory exists and create it if not        
try:
    os.mkdir('./'+planet_dir)

except FileExistsError:
    if overwrite:
        shutil.rmtree('./'+planet_dir)
        os.mkdir('./'+planet_dir)            
    
    else:
        print ('WARNING: the output directory already exists!')
        sys.exit()
   
time_i = time.time()
count = 0
for k in range(len(masses)):
    m = masses[k]
    
    for i in range(len(temps)):
        t = temps[i]
        
        for j in range(len(Mg_numbers)):
            mg = Mg_numbers[j]
            
            for l in range(len(ocean_fracs)):
                
                pl = toolkit.model_hydro(P_center=9.5e10, 
                                           match_mass=True, 
                                           eps_r=0.25, 
                                           predictor_T='none', 
                                           ocean=True, 
                                           ocean_frac_should=np.log10(ocean_fracs[l]), 
                                           temp_jumps=[0., delta_T_CMB(m), 
                                                       0., t-300., 0.], 
                                           Si_number_should=.4, 
                                           P_surface_should=1.0e5, 
                                           T_surface_should=300., 
                                           T_center=4000.,  
                                           Fe_number_mantle=.24999, 
                                           Mg_number_should=mg, 
                                           eps_H2O=0., 
                                           iterationType=1, 
                                           M_surface_should=m, 
                                           predictor_P='linear', 
                                           log=False, 
                                           sweeps=10, 
                                           iteration_limit=20, 
                                           subphase_res=32, 
                                           xi_Stv=0.,
                                           acc_M_surface = 1.0e-5,
                                           acc_Mg_number = 1.0e-1,
                                           acc_T_surface = 1.0e-2)
                
                planets.append(pl)
               
                pl.write(loc=planet_dir+'/', out = 'planet_'+str(count))
                
                count += 1
                
time_f=time.time()
ftool.printTime(sec=time_f-time_i, ms=False, where='run hydro curve')