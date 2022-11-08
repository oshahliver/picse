#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 23:37:37 2021

@author: os18o068
"""

from PIMPphysicalparams import m_earth, r_earth
import numpy as np

def update(dat, pl, dir=''):
    dat = pl.M_surface_is/m_earth, pl.R_surface_is/r_earth, pl.Mg_number_is, \
        pl.Si_number_is, pl.M_core_is, pl.T_surface_is, pl.P_surface_is, \
        pl.M_ocean_is, pl.M_H2O_is, pl.MOI_is, pl.M_H2O_core, pl.T_center, \
        pl.P_center, pl.ocean_depth, pl.xi_H_core, pl.xi_H_core_predicted,\
        pl.Mg_number_should, pl.M_surface_should/m_earth, pl.T_surface_should, \
        pl.ocean_frac_should, pl.ocean_frac_is, pl.P_H2_CMB, None, None
        
    np.save(dir, dat)