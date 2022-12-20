#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 08:39:24 2020

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
    return 1400.0 * M ** 0.75


Mg_numbers = np.linspace(0.2, 0.7, 2)
masses = np.logspace(-2, 0, 10)
temps = np.linspace(500.0, 2000.0, 1)

planets = [[], [], []]


subsubdir_names = []
overwrite = True

N = 2

inertias = np.empty([3, len(temps), len(Mg_numbers), N, len(masses)])

for k in range(len(masses)):
    m = masses[k]

    for i in range(len(temps)):
        t = temps[i]

        for j in range(len(Mg_numbers)):
            mg = Mg_numbers[j]

            pl1 = toolkit.model_hydro(
                P_center=9.5e10,
                match_mass=True,
                eps_r=0.25,
                predictor_T="none",
                ocean=False,
                ocean_frac_should=-1.72,
                temp_jumps=[0.0, delta_T_CMB(m), 0.0, 0.0, 0.0],
                Si_number_should=0.4,
                P_surface_should=1.0e5,
                T_surface_should=t,
                T_center=4000.0,
                Fe_number_mantle=0.2599,
                Mg_number_should=mg,
                eps_H2O=0.0,
                iterationType=0,
                M_surface_should=m,
                predictor_P="linear",
                log=False,
                sweeps=10,
                iteration_limit=20,
                subphase_res=32,
                xi_Stv=0.0,
                acc_M_surface=1.0e-4,
                acc_Mg_number=1.0e-1,
                acc_T_surface=1.0e-2,
            )

            pl2 = toolkit.model_hydro(
                P_center=9.5e10,
                match_mass=True,
                eps_r=0.25,
                predictor_T="none",
                ocean=False,
                ocean_frac_should=-1.72,
                temp_jumps=[0.0, delta_T_CMB(m), 0.0, 0.0, 0.0],
                Si_number_should=0.4,
                P_surface_should=1.0e5,
                T_surface_should=t,
                T_center=4000.0,
                Fe_number_mantle=0.2599,
                Mg_number_should=mg,
                eps_H2O=1.0,
                iterationType=1,
                M_surface_should=m,
                predictor_P="linear",
                log=False,
                sweeps=10,
                iteration_limit=20,
                subphase_res=32,
                xi_Stv=0.0,
                acc_M_surface=1.0e-4,
                acc_Mg_number=1.0e-1,
                acc_T_surface=1.0e-4,
            )

            m_oc = pl2.M_H2O_is

            pl3 = toolkit.model_hydro(
                P_center=9.5e10,
                match_mass=True,
                eps_r=0.25,
                predictor_T="none",
                ocean=True,
                ocean_frac_should=np.log10(m_oc),
                temp_jumps=[0.0, delta_T_CMB(m), 0.0, 0.0, 0.0],
                Si_number_should=0.4,
                P_surface_should=1.0e5,
                T_surface_should=300.0,
                T_center=4000.0,
                Fe_number_mantle=0.2599,
                Mg_number_should=mg,
                eps_H2O=0.0,
                iterationType=0,
                M_surface_should=m,
                predictor_P="linear",
                log=False,
                sweeps=10,
                iteration_limit=20,
                subphase_res=32,
                xi_Stv=0.0,
                acc_M_surface=1.0e-5,
                acc_Mg_number=1.0e-1,
                acc_T_surface=1.0e-2,
            )

            planets[0].append(pl1)
            planets[1].append(pl2)
            planets[2].append(pl3)

            for s in range(3):
                inertias[s][i][j][0][k] = planets[s][-1].M_surface_is
                inertias[s][i][j][1][k] = planets[s][-1].MOI_is

np.save("inertias.npy", inertias)
