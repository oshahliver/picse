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
    return 1400.0 * M**0.75


Mg_numbers = np.linspace(0.2, 0.7, 6)
masses = np.array([0.1, 0.6, 3.0])
temps = np.linspace(500.0, 2000.0, 4)

planets = [[], [], []]

subsubdir_names = []
overwrite = True

for i in range(len(masses)):
    planet_dir = "M_" + str(round(masses[i], 2))
    planet_dir = planet_dir.replace(".", "_")
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
        for j in range(3):
            type_str = str(j + 1)

            subsubdir_name = "Type_" + type_str
            subsubdir_name = subsubdir_name.replace(".", "_")

            mkdir_name = "./" + planet_dir + "/" + subdir_name + "/" + subsubdir_name

            try:
                os.mkdir(mkdir_name)
                print("creating directory ", mkdir_name)

            except FileExistsError:
                if overwrite:
                    shutil.rmtree(mkdir_name)
                    os.mkdir(mkdir_name)

            subsubdir_names[t].append(subsubdir_name)

time_i = time.time()

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
                Fe_number_mantle=0.24999,
                Mg_number_should=mg,
                eps_H2O=0.0,
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
                Fe_number_mantle=0.24999,
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
                temp_jumps=[0.0, delta_T_CMB(m), 0.0, t - 300.0, 0.0],
                Si_number_should=0.4,
                P_surface_should=1.0e5,
                T_surface_should=300.0,
                T_center=4000.0,
                Fe_number_mantle=0.24999,
                Mg_number_should=mg,
                eps_H2O=0.0,
                iterationType=1,
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

            loc_1 = "./M_" + str(round(m, 2)).replace(".", "_") + "/Tsurf_"
            loc_1 = loc_1 + str(int(t)) + "/Type_1/"

            loc_2 = "./M_" + str(round(m, 2)).replace(".", "_") + "/Tsurf_"
            loc_2 = loc_2 + str(int(t)) + "/Type_2/"

            loc_3 = "./M_" + str(round(m, 2)).replace(".", "_") + "/Tsurf_"
            loc_3 = loc_3 + str(int(t)) + "/Type_3/"

            pl1.write(loc=loc_1, out="planet_" + str(j))
            pl2.write(loc=loc_2, out="planet_" + str(j))
            pl3.write(loc=loc_3, out="planet_" + str(j))

time_f = time.time()
ftool.printTime(sec=time_f - time_i, ms=False, where="run hydro curve")
