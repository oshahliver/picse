#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:58:20 2020

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


Fe_number_mantle = 0.24999
CMB_fac = 1.0
TBL_fac = 1.0

suff = "nom_double_check"

planet_dirs = ["planets_hyd_" + suff, "planets_dry_" + suff, "planets_oce_" + suff]

# Create parameter space
Mg_numbers = np.linspace(0.2, 0.7, 2)
xi_Fe = np.linspace(0.0, Fe_number_mantle, 3)
Si_numbers = np.linspace(0.4, 0.5, 1)
xi_Stv = np.linspace(0.0, 0.2, 3)
temps = np.linspace(500, 2000, 2)
masses = np.logspace(-1, np.log10(3.0), 24)


# M, Mg#, T_surf
accs = [[1.0e-4, 1.0e-1, 1.0e-4], [1.0e-4, 1.0e-1, 1.0e-2], [1.0e-5, 1.0e-1, 1.0e-2]]

N_params = 24


data = np.empty([N_params, len(Si_numbers), len(temps), len(Mg_numbers), len(masses)])

ocean_masses = np.empty([len(Si_numbers), len(temps), len(Mg_numbers), len(masses)])

subsubdir_names = []

overwrite = True


for d in range(len(planet_dirs)):
    planet_dir = planet_dirs[d]

    # Check if planet output directory exists and create it if not
    try:
        os.mkdir("./" + planet_dir)

    except FileExistsError:
        if overwrite:
            shutil.rmtree("./" + planet_dir)
            os.mkdir("./" + planet_dir)

        else:
            print("WARNING: the output directory already exists!")
            sys.exit()

    # Now prepare subdirectories for each temperature
    for t in range(len(temps)):
        temp = temps[t]
        subdir_name = "Tsurf_" + str(int(temp))

        os.mkdir("./" + planet_dir + "/" + subdir_name)

        subsubdir_names.append([])

        # Prepare susubdirectories for total Mg#
        for j in range(len(Mg_numbers)):
            Mg_round = round(Mg_numbers[j], 2)
            if round(Mg_round, 1) == Mg_round:
                Mg_string = str(Mg_round) + "0"

            else:
                Mg_string = str(Mg_round)

            subsubdir_name = "Mg_number_" + Mg_string
            subsubdir_name = subsubdir_name.replace(".", "_")
            os.mkdir("./" + planet_dir + "/" + subdir_name + "/" + subsubdir_name)

            subsubdir_names[t].append(subsubdir_name)


t0 = time.time()

# Create hydrated planets
for s in range(len(Si_numbers)):
    for i in range(len(temps)):
        for j in range(len(Mg_numbers)):

            Mg_round = round(Mg_numbers[j], 2)
            if round(Mg_round, 1) == Mg_round:
                Mg_string = str(Mg_round) + "0"

            else:
                Mg_string = str(Mg_round)

            subdir_name = "Tsurf_" + str(int(temps[i]))
            subsubdir_name = "Mg_number_" + Mg_string
            subsubdir_name = subsubdir_name.replace(".", "_")

            dir_name = planet_dirs[0] + "/" + subdir_name + "/" + subsubdir_name + "/"

            for k in range(len(masses)):

                pl = toolkit.model_hydro(
                    P_center=9.5e10,
                    match_mass=True,
                    eps_r=0.25,
                    predictor_T="none",
                    ocean=False,
                    ocean_frac_should=-10,
                    temp_jumps=[0.0, CMB_fac * delta_T_CMB(masses[k]), 0.0, 0.0, 0.0],
                    Si_number_should=Si_numbers[s],
                    P_surface_should=1.0e5,
                    T_surface_should=temps[i],
                    T_center=4000.0,
                    Fe_number_mantle=Fe_number_mantle,
                    Mg_number_should=Mg_numbers[j],
                    eps_H2O=1.0,
                    iterationType=1,
                    M_surface_should=masses[k],
                    predictor_P="linear",
                    log=False,
                    sweeps=10,
                    iteration_limit=20,
                    subphase_res=32,
                    xi_Stv=0.0,
                    acc_M_surface=accs[0][0],
                    acc_Mg_number=accs[0][1],
                    acc_T_surface=accs[0][2],
                )

                dat = (
                    pl.M_surface_is / m_earth,
                    pl.R_surface_is / r_earth,
                    pl.Mg_number_is,
                    pl.Si_number_is,
                    pl.M_core_is,
                    pl.T_surface_is,
                    pl.P_surface_is,
                    pl.M_ocean_is,
                    pl.M_H2O_is,
                    pl.MOI_is,
                    pl.M_H2O_core,
                    pl.T_center,
                    pl.P_center,
                    pl.ocean_depth,
                    pl.xi_H_core,
                    pl.xi_H_core_predicted,
                    pl.Mg_number_should,
                    pl.M_surface_should / m_earth,
                    pl.T_surface_should,
                    pl.ocean_frac_should,
                    pl.ocean_frac_is,
                    pl.P_H2_CMB,
                    None,
                    None,
                )

                ocean_masses[s][i][j][k] = dat[8]

                for d in range(len(dat)):
                    data[d][s][i][j][k] = dat[d]

                if k < 10:
                    num = "00" + str(k)

                elif k >= 10 and k < 100:
                    num = "0" + str(k)

                else:
                    num = str(k)

                # pl.write(loc=dir_name, out = 'planet_'+num)


np.save("planet_output_hyd_" + suff + ".npy", data)

# Create dry planets
for s in range(len(Si_numbers)):
    for i in range(len(temps)):
        for j in range(len(Mg_numbers)):

            Mg_round = round(Mg_numbers[j], 2)
            if round(Mg_round, 1) == Mg_round:
                Mg_string = str(Mg_round) + "0"

            else:
                Mg_string = str(Mg_round)

            subdir_name = "Tsurf_" + str(int(temps[i]))
            subsubdir_name = "Mg_number_" + Mg_string
            subsubdir_name = subsubdir_name.replace(".", "_")

            dir_name = planet_dirs[1] + "/" + subdir_name + "/" + subsubdir_name + "/"

            for k in range(len(masses)):

                pl = toolkit.model_hydro(
                    P_center=9.5e10,
                    match_mass=True,
                    eps_r=0.25,
                    predictor_T="none",
                    ocean=False,
                    ocean_frac_should=-10,
                    temp_jumps=[0.0, CMB_fac * delta_T_CMB(masses[k]), 0.0, 0.0, 0.0],
                    Si_number_should=Si_numbers[s],
                    P_surface_should=1.0e5,
                    T_surface_should=temps[i],
                    T_center=4000.0,
                    Fe_number_mantle=Fe_number_mantle,
                    Mg_number_should=Mg_numbers[j],
                    eps_H2O=0.0,
                    iterationType=0,
                    M_surface_should=masses[k],
                    predictor_P="linear",
                    log=False,
                    sweeps=10,
                    iteration_limit=20,
                    subphase_res=32,
                    xi_Stv=0.0,
                    acc_M_surface=accs[1][0],
                    acc_Mg_number=accs[1][1],
                    acc_T_surface=accs[1][2],
                )

                dat = (
                    pl.M_surface_is / m_earth,
                    pl.R_surface_is / r_earth,
                    pl.Mg_number_is,
                    pl.Si_number_is,
                    pl.M_core_is,
                    pl.T_surface_is,
                    pl.P_surface_is,
                    pl.M_ocean_is,
                    pl.M_H2O_is,
                    pl.MOI_is,
                    pl.M_H2O_core,
                    pl.T_center,
                    pl.P_center,
                    pl.ocean_depth,
                    pl.xi_H_core,
                    pl.xi_H_core_predicted,
                    pl.Mg_number_should,
                    pl.M_surface_should / m_earth,
                    pl.T_surface_should,
                    pl.ocean_frac_should,
                    pl.ocean_frac_is,
                    pl.P_H2_CMB,
                    None,
                    None,
                )

                for d in range(len(dat)):
                    data[d][s][i][j][k] = dat[d]

                if k < 10:
                    num = "00" + str(k)

                elif k >= 10 and k < 100:
                    num = "0" + str(k)

                else:
                    num = str(k)

                for d in range(len(dat)):
                    data[d][s][i][j][k] = dat[d]

                # pl.write(loc=dir_name, out = 'planet_'+num)

np.save("planet_output_dry_" + suff + ".npy", data)

# Create ocean planets
for s in range(len(Si_numbers)):
    for i in range(len(temps)):
        for j in range(len(Mg_numbers)):

            Mg_round = round(Mg_numbers[j], 2)
            if round(Mg_round, 1) == Mg_round:
                Mg_string = str(Mg_round) + "0"

            else:
                Mg_string = str(Mg_round)

            subdir_name = "Tsurf_" + str(int(temps[i]))
            subsubdir_name = "Mg_number_" + Mg_string
            subsubdir_name = subsubdir_name.replace(".", "_")
            dir_name = planet_dirs[2] + "/" + subdir_name + "/" + subsubdir_name + "/"

            for k in range(len(masses)):

                of = np.log10(ocean_masses[s][i][j][k])

                pl = toolkit.model_hydro(
                    P_center=9.5e10,
                    match_mass=True,
                    eps_r=0.25,
                    predictor_T="none",
                    ocean=True,
                    ocean_frac_should=of,
                    temp_jumps=[
                        0.0,
                        CMB_fac * delta_T_CMB(masses[k]),
                        0.0,
                        (temps[i] - 300.0) * TBL_fac,
                        0.0,
                    ],
                    Si_number_should=Si_numbers[s],
                    P_surface_should=1.0e5,
                    T_surface_should=300.0,
                    T_center=4000.0,
                    Fe_number_mantle=Fe_number_mantle,
                    Mg_number_should=Mg_numbers[j],
                    eps_H2O=0.0,
                    iterationType=1,
                    M_surface_should=masses[k],
                    predictor_P="linear",
                    log=False,
                    sweeps=10,
                    iteration_limit=20,
                    subphase_res=32,
                    xi_Stv=0.0,
                    acc_M_surface=accs[2][0],
                    acc_Mg_number=accs[2][1],
                    acc_T_surface=accs[2][2],
                )

                dat = (
                    pl.M_surface_is / m_earth,
                    pl.R_surface_is / r_earth,
                    pl.Mg_number_is,
                    pl.Si_number_is,
                    pl.M_core_is,
                    pl.T_surface_is,
                    pl.P_surface_is,
                    pl.M_ocean_is,
                    pl.M_H2O_is,
                    pl.MOI_is,
                    pl.M_H2O_core,
                    pl.T_center,
                    pl.P_center,
                    pl.ocean_depth,
                    pl.xi_H_core,
                    pl.xi_H_core_predicted,
                    pl.Mg_number_should,
                    pl.M_surface_should / m_earth,
                    pl.T_surface_should,
                    pl.ocean_frac_should,
                    pl.ocean_frac_is,
                    pl.P_H2_CMB,
                    pl.finals["layer_properties"][3]["P_outer"],
                    pl.finals["layer_properties"][3]["T_outer"],
                )

                for d in range(len(dat)):
                    data[d][s][i][j][k] = dat[d]

                if k < 10:
                    num = "00" + str(k)

                elif k >= 10 and k < 100:
                    num = "0" + str(k)

                else:
                    num = str(k)

                for d in range(len(dat)):
                    data[d][s][i][j][k] = dat[d]

                # pl.write(loc=dir_name, out = 'planet_'+num)

np.save("planet_output_oce_" + suff + ".npy", data)

t = time.time()
ftool.printTime(sec=t - t0, ms=False, where="run hydro curve")
