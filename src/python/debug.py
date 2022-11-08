#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 16:31:04 2021

@author: os18o068
"""

import PlanetFactory as plfac;toolkit=plfac.Toolkit()
import numpy as np
from PIMPphysicalparams import R_solar, M_solar, T_solar, axes_solar, \
MoI_solar, mFe, mO, mH, mMg, mS, mSi, m_earth, r_earth

si = .5
obj = 'earth'

of=-10;im=0;ts=1500;m=M_solar[obj];mg=.49;o=False;h2o=1.;fe=.1;dt=3.;it=1;si=si
eps_r=.25;T0=1400.;FeS=.3;Pc=9.5e10

pl3=toolkit.model_hydro(P_center=Pc, match_mass=True, eps_r=eps_r, 
                        predictor_T='none', ocean=o, ocean_frac_should=of, 
                        temp_jumps=[0., T0*m**.75, 300., 0., 0.], 
                        Si_number_should=si, P_surface_should=1.0e5, 
                        T_surface_should=ts, T_center=4000.,  
                        Fe_number_mantle=fe, Mg_number_should=mg, eps_H2O=h2o, 
                        iterationType=it, M_surface_should=m, 
                        predictor_P='linear', log=False, subphase_res=32, 
                        xi_Stv=0., acc_Mg_number=1.0e-3, acc_T_surface=5.0e-4, 
                        acc_M_surface=1.0e-4, xi_FeS=FeS, X_impurity=im, 
                        X_impurity_0_layers=[0., 0., im, im, 0], sweeps=10, 
                        iteration_limit=25, impurity=9, inner_core_frac=.5)

pl3.Plot(spec=obj)
pl3.prt()
