#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 18:54:45 2020

@author: os18o068
"""

from pics.interiors import PlanetFactory as plfac
import numpy as np
import pandas as pd
from decimal import Decimal
from pics.interiors import MOI_model
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import sys
import time
from pics.utils import fortfunctions
import matplotlib.patches as patches
import matplotlib.transforms as transforms
from pics.physicalparams import (
    m_earth,
    r_earth,
    MoI_solar,
    delta_MoI_solar,
    R_solar,
    J2_solar,
    delta_J2_solar,
    orbit_solar,
    day,
    yr,
    M_solar,
    T_solar,
    material_list_fort,
    material_YMg,
    material_YSi,
    material_YO,
    material_YH,
    material_YS,
    mS,
    mH2O,
    mFe,
    mSi,
    mO,
    mMg,
    mH,
    mEn,
    mH2O,
    material_YFe,
    Si_earth,
    Mg_earth,
    SFe_earth,
    SFe_solar,
    sigmaSB,
    V_ocean_earth,
)
from pics.materials import material
from pics.utils import readPREM
from pics.interiors import PlanetFort
from pics.runparams import color_list, xi_Fe_mantle_max
from matplotlib import pyplot as plt
from matplotlib import pylab
import matplotlib
from pics.utils import functionTools as ftool
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import random
from matplotlib.ticker import MaxNLocator
from pics.utils import read_exoplanets

matplotlib.rc("text", usetex=True)
matplotlib.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath, amssymb}"]

data_dir = (
    "/home/os18o068/Documents/PHD/Projects/Planets/Data/icy_satellites/"
    + "compare_ice_fractions/"
)

toolkit = plfac.Toolkit()
PREM_data = readPREM.Do()

ocean_mass_fractions = np.linspace(0.0001, 0.6, 1)
impurity_abundances = np.linspace(0.0, 0.5, 1)

hard_ranges = np.array(
    [
        [0.4, 0.7],  # Mg#
        [1.0 / 3.0, 0.5],  # Si#
        [0.0, 1.0e-6],  # M_H2O
        [0.0, 0.5],  # xi_S
        [0.0, 0.3],  # xi_Fe
        [0.0, 1.0e-6],  # X_ice_mean
        [0.0, 0.1],
    ]
)  # X_ice_0

fnts1 = 8
fnts2 = 10
fnts3 = 12
fnts4 = 16

params = [
    "Mg#",
    "R",
    "MoI",
    "Si#",
    "Si# m",
    "xi_Fe",
    "P_CS",
    "T_CS",
    "X_S",
    "X_Si",
    "X_O",
    "xi_FeS",
    "xi_FeSi",
    "xi_FeO",
    "T0_CMB",
    "R_core (km)",
    "R_IC (km)",
    "R_HMI (km)",
    "T_TBL (K)",
    "xi_SiO2 m",
    "xi_FeO m",
    "T_C",
    "P_C",
    "fO2",
    "DeltaT_MTZ",
    "M_H2O",
    "M_H2O_mantle",
    "Fe#_core",
    "core_mass_frac",
    "T_CMB",
    "P_CMB",
    "rho_CMB",
    "T_ICB",
    "P_ICB",
    "rho_ICB",
    "T_MTZ",
    "P_MTZ",
    "rho_MTZ",
    "xi_MgO m",
    "R_MTZ",
    "xi_S_core",
    "xi_S_mantle",
    "S/Fe",
    "Si/Fe",
    "Mg/Fe",
    "S/Si",
    "M",
    "f_IC",
    "ocean_frac",
    "P_surf",
    "w_liq",
    "w_scf",
    "w_sol",
    "h_liq",
    "h_scf",
    "h_sol",
    "w_liq/w_H2O",
    "w_SCF/w_H2O",
    "w_sol/w_H2O",
    "V_liq",
    "V_scf",
    "V_sol",
    "grav",
    "hyd_struc",
    "T_HMI",
    "P_HMI",
    "rho_HMI",
    "bulk_struc",
]

all_param_labels = [
    r"[Mg]/[Mg + Fe]",
    r"$R/R_\oplus$",
    r"$C/MR^2$",
    r"[Si]/[Si + Mg]",
    r"$\rm Si\#_{\rm Mantle}$",
    r"$X_{\rm Fe}^{\rm Mantle}$",
    r"$P_{\rm CS}$",
    r"$T_{\rm CS}$",
    r"$w_{\rm S}$",
    r"$w_{\rm Si}$",
    r"$w_{\rm O}$",
    r"$X_{\rm FeS}^{\rm Core}$",
    r"$X_{\rm FeSi}^{\rm Core}$",
    r"$X_{\rm FeO}^{\rm Core}$",
    r"$T_{0, \rm CMB}$",
    r"$R_{\rm Core}$",
    r"$R_{\rm ICB}$",
    r"$R_{\rm HMI}$",
    r"$T_{\rm TBL}$",
    r"$X_{\rm SiO_2}^{\rm Mantle}$",
    r"$X_{\rm FeO}^{\rm Mantle}$",
    r"$T_{\rm C}$",
    r"$P_{\rm C}$",
    r"$\log(f_{\rm O_2})$",
    r"$\Delta T_{\rm MTZ}$",
    r"$M_{\rm H_2 O}$",
    r"$M_{\rm H_2 O}^{\rm Mantle}$",
    r"$X_{\rm Fe}^{\rm Core}$",
    r"$M_{\rm Core}/M$",
    r"$T_{\rm CMB}$",
    r"$P_{\rm CMB}$",
    r"$\rho_{\rm CMB}$",
    r"$T_{\rm ICB}$",
    r"$P_{\rm ICB}$",
    r"$\rho_{\rm ICB}$",
    r"$T_{\rm MTZ}$",
    r"$P_{\rm MTZ}$",
    r"$\rho_{\rm MTZ}$",
    r"$X_{\rm MgO}^{\rm Mantle}$",
    r"$R_{\rm MTZ}$",
    r"$X_{\rm S}^{\rm Core}$",
    r"$X_{\rm S}^{\rm Mantle}$",
    "[S]/[Fe]",
    "[Si]/[Fe]",
    "[Mg]/[Fe]",
    "[S]/[Si]",
    r"$M/M_\oplus$",
    r"$M_{\rm IC}/M_{\rm Core}$",
    r"$M_{\rm Ocean} / M$",
    r"$P_{\rm Surf}$",
    r"$w_{\rm liq}$",
    r"$w_{\rm SCF}$",
    r"$w_{\rm sol}$" r"$h_{\rm liq}$",
    r"$h_{\rm SCF}$",
    r"$h_{\rm sol}$",
    r"$w_{\rm liq}/w_{\rm H_2O}$",
    r"$w_{\rm SCF}/w_{\rm H_2O}$",
    r"$w_{\rm sol}/w_{\rm H_2O}$",
    r"$V_{\rm liq}$",
    r"$V_{\rm scf}$",
    r"$V_{\rm sol}$",
    r"$g_{\rm Surf}$",
    r"$\rm Hydro \ structure$",
    "T_HMI",
    "P_HMI",
    "rho_HMI",
    "bulk structure",
]


def Mass_n_layer_planet(rho, r):
    mass = 0.0
    for i in range(len(rho)):
        mass += rho[i] * r[i] ** 3

    for i in range(len(rho) - 1):
        mass -= rho[i + 1] * r[i] ** 3

    return 4 / 3 * np.pi * mass


def MoI_n_layer_planet(rho, r):
    moi = 0.0
    for i in range(len(rho)):
        moi += rho[i] * r[i] ** 5

    for i in range(len(rho) - 1):
        moi -= rho[i + 1] * r[i] ** 5

    return 8 / 15 * np.pi * moi


def C(rho, r):
    r = r * r_earth
    return MoI_n_layer_planet(rho, r) / Mass_n_layer_planet(rho, r) / r[-1] ** 2


def plot_N_layer_planet(N=15):
    rho = [1e4, 5e3, 1e3]

    r1 = np.linspace(0.0, 1.0, N)

    r2 = np.linspace(0.0, 1.0, N)
    r3 = 1.0

    data = np.empty([len(r1), len(r2)])

    for i in range(len(r1)):
        print("")
        for j in range(len(r2)):
            r = np.array([r1[i], (1.0 - r1[i]) * r2[j] + r1[i], r3])
            r = r * r_earth
            data[i][j] = MoI_n_layer_planet(rho, r)
            data[i][j] /= Mass_n_layer_planet(rho, r) * r[-1] ** 2
            print(r)
            print(data[i][j])

    fig, ax = plt.subplots(2, 2)

    for i in range(len(r1)):
        ax[0][0].plot(r1, data[i])

    return data


def convert_X_impurity_to_xi_impurity(
    Si_number=None, xi_Fe=None, X=None, contents=[6, 7, 1], m=(2 * mH + mO)
):

    SiMg = Si_number / (1.0 - Si_number)

    # Compute fractions without water
    fractions = fortfunctions.functions.compute_abundance_vector(
        simg=SiMg,
        femg=1.0 / (1.0 - xi_Fe) - 1.0,
        n_mats=len(contents),
        ymgi=[material_YMg[i - 1] for i in contents],
        ysii=[material_YSi[i - 1] for i in contents],
        xih2oi=[0.0 for i in contents],
        xifei=[xi_Fe for i in contents],
        xialsii=[0.0 for i in contents],
        xialmgi=[0.0 for i in contents],
        contents=contents,
        additional=[0.0],
    )

    # Compute total normalized mass in the mantle
    m_tilde = (
        sum(
            [
                fractions[i] * (1.0 - xi_Fe) * material_YMg[contents[i] - 1]
                for i in range(len(contents))
            ]
        )
        * mMg
        + sum(
            [
                fractions[i] * xi_Fe * material_YMg[contents[i] - 1]
                for i in range(len(contents))
            ]
        )
        * mFe
        + sum(
            [fractions[i] * material_YSi[contents[i] - 1] for i in range(len(contents))]
        )
        * mSi
        + sum(
            [fractions[i] * material_YO[contents[i] - 1] for i in range(len(contents))]
        )
        * mO
        + sum(
            [fractions[i] * material_YH[contents[i] - 1] for i in range(len(contents))]
        )
        * mH
        + sum(
            [fractions[i] * material_YS[contents[i] - 1] for i in range(len(contents))]
        )
        * mS
        + sum(
            [fractions[i] * material_YFe[contents[i] - 1] for i in range(len(contents))]
        )
        * mFe
    )

    # Compute mole fraction of impurity at given composition
    xi = material.xi(eta=X, m1=m_tilde, m2=m)

    return xi


def compute_max_ocean_frac(
    Si_number=None,
    xi_Fe=None,
    X_ice=None,
    contents=[6, 7, 1],
    Mg_number=None,
    H2O_frac=None,
    xi_S=None,
):
    """The total water mass fraction is (M_ice + M_ocean)/M. For a given amount
    of water in the mantle, given by X_ice, the possible value range for the
    ocean mass fraction is therefore limited.
    """

    # First convert the ice mass fraction in the mantle to mole fraction for
    # the given composition
    xi_ice = convert_X_impurity_to_xi_impurity(
        Si_number=Si_number, xi_Fe=xi_Fe, X=X_ice, m=mH2O, contents=[6, 7, 1]
    )

    # With this, compute the core mass fraction without an ocean. Note that for
    # M_surface = 1 the value of the core mass corresponds to the value of the
    # core mass fraction
    all_contents = [[2, 8], [2, 8], [4, 5, 1], [6, 7, 1], [1]]
    SiMg = Si_number / (1.0 - Si_number)

    print("Core mass in compute max ocean")
    core_frac = PlanetFort.ComputeCoreMass(
        contents=all_contents,
        Mg_number=Mg_number,
        Mg_number_mantle=1.0 - xi_Fe,
        SiMg=SiMg,
        M_surface=1.0,
        M_ocean=0.0,
        impurity=[xi_ice],
        xi_S_core=xi_S,
    )

    print("core_frac =", core_frac)
    # Compute ocean mass fraction
    ocean_frac = (H2O_frac - X_ice * (1.0 - core_frac)) / (
        1.0 - X_ice + X_ice * core_frac
    )

    return ocean_frac


def compute_max_X_ice(
    Si_number=None,
    xi_Fe=None,
    contents=[6, 7, 1],
    Mg_number=None,
    H2O_frac=None,
    xi_S=None,
):
    """The maximum amount of ice in the mantle corresponds to the total amount
    of water in the planet, i.e. to the case where M_ocean = 0.
    """
    # Compute core mass fraction
    all_contents = [[2, 8], [2, 8], [4, 5, 1], [6, 7, 1], [1]]
    SiMg = Si_number / (1.0 - Si_number)

    core_frac = PlanetFort.ComputeCoreMass(
        contents=all_contents,
        Mg_number=Mg_number,
        Mg_number_mantle=1.0 - xi_Fe,
        SiMg=SiMg,
        M_surface=1.0,
        M_ocean=H2O_frac,
        impurity=[0.0],
        xi_S_core=xi_S,
    )

    X_ice_max = H2O_frac / (1.0 - core_frac)
    return min(X_ice_max, hard_ranges[5][1]), core_frac


def compute_X_impurity_slope(M_mantle=None, X=None, X_0=None):
    """
    The slope of the ice fraction in the mantle is given by the mean ice content
    and the ice content at the bottom of the mantle.
    """

    # Use linear dependency of the ice mass fraction on the mass
    slope = 2.0 * (X - X_0) / M_mantle

    return slope


def filter_data_convergence(data, accs=[1.0e-4, 1.0e-3, 1.0e-2]):
    """
    Takes a data set as input an filters it using the given accuracies. All
    data points which did not converge according to the criterion will be
    rejected. The parameters which will be checked are:

        total mass
        total magnesium number
        ocean mass fraction

    """
    fail_count = 0
    for i in range(len(data)):
        M_surface_is = data[i][7]
        M_surface_should = data[i][11]

        Mg_number_is = data[i][3]
        Mg_number_should = data[i][8]

        ocean_frac_is = data[i][5]
        ocean_frac_should = data[i][10]

        reldev_Mg = (Mg_number_is - Mg_number_should) / Mg_number_should
        reldev_mass = (M_surface_is - M_surface_should) / M_surface_should
        reldev_ocean = (ocean_frac_is - ocean_frac_should) / ocean_frac_should

        if ocean_frac_should < 0.01:
            reldev_ocean = 0.0

        if (
            abs(reldev_mass) > accs[0]
            or abs(reldev_Mg) > accs[1]
            or abs(reldev_ocean) > accs[2]
        ):
            print("----")
            print("not converged")
            print("reldevs =", reldev_mass, reldev_Mg, reldev_ocean)
            data[i][:] = np.nan
            fail_count += 1
    print(fail_count, "planets did not converge.")
    return data, fail_count


def filter_data_match(data, obj="titan", filter_type=1, acc_R=1.0e-2, acc_MoI=2.0e-2):
    """Checks each point of the data set if it matches the properties of the
    given object with the desired precision. The parameters that are probed
    are:

        MoI factor
        J2 coefficient
        total Radius

    """
    acc = np.sqrt(acc_R**2 + acc_MoI**2)

    # Extract the values of the real object
    data_is = [MoI_solar[obj], J2_solar[obj], R_solar[obj]]
    reldev_best_fit_list = np.empty([len(data)])

    for i in range(len(data)):
        try:
            reldev_R = (data[i][2] - data_is[2]) / data_is[2]
        except TypeError:
            reldev_R = 1.0e10

        try:
            reldev_MoI = (data[i][0] - data_is[0]) / data_is[0]
        except TypeError:
            reldev_MoI = 1.0e10

        try:
            reldev_J2 = (data[i][1] * 1.0e6 - data_is[1]) / data_is[1]
        except TypeError:
            reldev_J2 = 1.0e10

        reldev_best_fit = np.sqrt(reldev_R**2 + reldev_MoI**2)
        reldev_best_fit_list[i] = reldev_best_fit

        # The total radius and MoI factor are combined to determine a match
        if filter_type == 0:
            reldev = np.sqrt(reldev_R**2 + reldev_MoI**2)

            # If the data point does not match the real values, eject it
            if reldev > acc:
                data[i][:] = np.nan

        # Total radius and MoI factor are used independently to determine a match
        elif filter_type == 1:
            if abs(reldev_R) < acc_R and abs(reldev_MoI) < acc_MoI:
                pass

            else:
                data[i][:] = np.nan

    best_fit_index = np.nanargmin(reldev_best_fit_list)
    return data, best_fit_index


def run3(obj="titan", res=3, iterations=2, sampler=0, ranges=hard_ranges):

    N = 2**res
    eps_H2O = 1.0
    M_surface_should = M_solar[obj]
    data = np.empty([N, 15])
    random_vals = np.empty([len(ranges)])

    for it in range(iterations):
        for i in range(N):

            # Draw parameters
            if sampler == 0:
                for j in range(len(ranges)):

                    p = random.random()
                    delta = ranges[j][1] - ranges[j][0]

                    val = ranges[j][0] + p * delta
                    random_vals[j] = val

                # The maximum Mg# is limited by (1-xi_Fe) in the mantle
                p = random.random()
                delta = ranges[0][1] - random_vals[4] - ranges[0][0]

                val = ranges[0][0] + p * delta
                random_vals[0] = val

                # The Si# or equivalently Si/Mg is limited by xi_Fe
                Si_number_min = 1.0 / (3.0 - 2 * random_vals[4])
                Si_number_max = 1.0 / (2.0 - random_vals[4])
                p = random.random()
                delta = Si_number_max - Si_number_min
                random_vals[1] = Si_number_min + p * delta

                # The maximum ice mas4.055489476149217e+22s fraction in the mantle is limited by
                # the total amount of water
                X_ice_max, core_frac = compute_max_X_ice(
                    Si_number=random_vals[1],
                    xi_Fe=random_vals[4],
                    Mg_number=random_vals[0],
                    H2O_frac=random_vals[2],
                    xi_S=random_vals[3],
                )

                p = random.random()
                delta = min(X_ice_max, ranges[5][1]) - ranges[5][0]
                print("delta =", delta)
                val = ranges[5][0] + p * delta
                random_vals[5] = val
                print("p/val =", p, val)
                # The maximum ocean mass fraction is limited by the ice content
                # in the mantle
                ocean_frac_max = compute_max_ocean_frac(
                    Si_number=random_vals[1],
                    xi_Fe=random_vals[4],
                    Mg_number=random_vals[0],
                    X_ice=random_vals[5],
                    H2O_frac=random_vals[2],
                    xi_S=random_vals[3],
                )

                p = random.random()

                delta = min(ocean_frac_max, ranges[2][1]) - ranges[2][0]

                val = ranges[2][0] + p * delta
                random_vals[2] = 1e-10  # val

                # Scale down X_ice_0 to be between 0 and X_ice
                random_vals[6] = random_vals[5]

            # Convert wt fraction of water into mole fraction as model input
            xi_ice = convert_X_impurity_to_xi_impurity(
                Si_number=random_vals[1],
                xi_Fe=random_vals[4],
                X=random_vals[5],
                m=mH2O,
                contents=[6, 7, 1],
            )

            xi_ice_0 = convert_X_impurity_to_xi_impurity(
                Si_number=random_vals[1],
                xi_Fe=random_vals[4],
                X=random_vals[6],
                m=mH2O,
                contents=[6, 7, 1],
            )

            # Now that all the layers are defined the mantle mass can be computed
            M_mantle = (1.0 - core_frac - random_vals[2]) * M_surface_should * m_earth

            # From this the slope of the ice mass fraction in the mantle is computed
            X_ice_slope = compute_X_impurity_slope(
                X=random_vals[5], X_0=random_vals[6], M_mantle=M_mantle
            )

            X_impurity_0_layers = [0.0, 0.0, random_vals[6], random_vals[6], 0.0]

            pl = toolkit.model_hydro(
                P_center=9.5e10,
                match_mass=True,
                eps_r=0.25,
                predictor_T="none",
                ocean=False,
                ocean_frac_should=np.log10(random_vals[2]),
                temp_jumps=[
                    0.0,
                    0.0 * 1400.0 * M_surface_should**0.75,
                    0.0,
                    0.0,
                    0.0,
                ],
                Si_number_should=random_vals[1],
                P_surface_should=1.0e5,
                T_surface_should=max(T_solar[obj], 300),
                T_center=4000.0,
                Fe_number_mantle=random_vals[4],
                Mg_number_should=random_vals[0],
                eps_H2O=eps_H2O,
                iterationType=1,
                M_surface_should=M_surface_should,
                predictor_P="linear",
                log=False,
                sweeps=10,
                iteration_limit=25,
                subphase_res=32,
                xi_Stv=0.0,
                acc_Mg_number=1.0e-3,
                acc_T_surface=1.0e-2,
                acc_M_surface=1.0e-4,
                X_impurity=random_vals[5],
                X_impurity_0_layers=X_impurity_0_layers,
                xi_FeS=random_vals[3],
                impurity=1,
            )

            print("core is =", pl.M_core_is * m_earth)
            print("core =", core_frac * M_surface_should * m_earth)
            print("mantle is =", pl.layer_properties[3]["indigenous_mass"])
            print("mantle =", M_mantle)
            print("X_ice =", random_vals[5])
            print("X_ice_0 =", random_vals[6])
            print("xi_ice_0 =", xi_ice_0)
            print("xi_ice =", xi_ice)
            print("Mg_number = ", random_vals[0])
            print("Si_number = ", random_vals[1])
            print("xi_Fe =", random_vals[4])
            print("xi_S =", random_vals[3])
            print("ocean =", random_vals[2])
            print("X_ice_max =", X_ice_max)
            print("predict =", random_vals[6] + X_ice_slope * M_mantle)

            # Compute kf Love number
            kf = MOI_model.k_f(pl.MOI_is)

            # Compute angular velocity
            omega = 2.0 * np.pi / (orbit_solar[obj] * day)

            # Compute J2 gravity coefficient
            J2 = MOI_model.J_2(omega, pl.R_surface_is, pl.M_surface_is, kf)

            # Compute mantle mass
            M_mantle = (
                pl.layer_properties[3]["indigenous_mass"]
                + pl.layer_properties[2]["indigenous_mass"]
            )

            # Compute total mass of ice in the mantle
            M_ice_mantle = M_mantle * random_vals[5]

            data[i][0] = pl.MOI_is
            data[i][1] = J2
            data[i][2] = pl.R_surface_is / r_earth
            data[i][3] = pl.Mg_number_is
            data[i][4] = pl.Si_number_is  # /(1.-pl.Si_number_is)
            data[i][5] = pl.M_H2O_is
            data[i][6] = random_vals[5]  # X_ice
            data[i][7] = pl.M_surface_is / m_earth
            data[i][8] = pl.Mg_number_should
            data[i][9] = pl.Si_number_should
            data[i][10] = random_vals[2] + M_ice_mantle / pl.M_surface_is
            data[i][11] = pl.M_surface_should / m_earth
            data[i][12] = random_vals[4]  # xi_Fe
            data[i][13] = pl.xi_S_core
            data[i][14] = (random_vals[5] - random_vals[6]) / random_vals[5]

    return data


def extract_parameter_ranges(
    data, obj="titan", acc_R=1.0e-2, acc_MoI=2.0e-2, filter_type=0
):
    data, fail_count = filter_data_convergence(data)
    data, best_fit_index = filter_data_match(
        data, obj=obj, acc_R=acc_R, acc_MoI=acc_MoI, filter_type=filter_type
    )

    ranges_indices = [3, 4, 5, 13, 12, 6, 14]
    ranges = np.empty([len(hard_ranges), 2])

    for i in range(len(hard_ranges)):
        ind = ranges_indices[i]
        ranges[i][0] = np.nanmin(data.T[ind])
        ranges[i][1] = np.nanmax(data.T[ind])

    return ranges, fail_count, best_fit_index


def check_convergence(
    obj="titan",
    N_points=4,
    acc_R=1.0e-2,
    acc_MoI=2.0e-2,
    res=1,
    write=False,
    read=True,
    filter_type=0,
):
    all_ranges = []
    ranges = hard_ranges

    data = np.empty([0, 15])

    t0 = time.time()

    fail_count = 0
    for i in range(N_points):
        file = obj + "_hyd_output_" + str(i + 1) + ".npy"
        if read:
            add_data = np.load(data_dir + obj.capitalize() + "/" + file)

        else:
            add_data = run3(obj=obj, res=res + i, iterations=1, ranges=hard_ranges)

            # write raw data to file before checking for convergence and filtering
            if write:
                np.save(data_dir + obj.capitalize() + "/" + file, add_data)

        data = np.append(data, add_data, axis=0)

        # Extract the maximum parameter range for the current run
        ranges, fc, bfi = extract_parameter_ranges(
            data, obj=obj, acc_R=acc_R, acc_MoI=acc_MoI, filter_type=filter_type
        )
        all_ranges.append(ranges)
        fail_count += fc

    print(fail_count, "planets did not converge in total.")
    print("best fit:", data[bfi])

    all_ranges = np.array(all_ranges)

    fig, ax = plt.subplots(len(hard_ranges), 1, figsize=(6, 8))
    plt.subplots_adjust(hspace=0.5)

    # Convert from logspace to ocean mass fraction
    all_ranges.T[0][2] = all_ranges.T[0][2]
    all_ranges.T[1][2] = all_ranges.T[1][2]

    y_labels = [
        r"$\rm Mg \#$",
        r"$\rm Si \#$",
        r"$M_{\rm H_2 O}/M$",
        r"$\xi_{\rm FeS}$",
        r"$\xi_{\rm Fe}$",
        r"${\bar{X}_{\rm Ice}}$",
        r"$\frac{{\bar{X}_{\rm Ice}}-X_{\rm Ice, 0}}{\bar{X}_{\rm Ice}}$",
    ]

    ax[-1].set_xlabel(r"$N \ \rm runs$")

    for i in range(len(hard_ranges)):
        ax[i].plot(all_ranges.T[0][i])
        ax[i].plot(all_ranges.T[1][i])

        ax[i].set_ylabel(y_labels[i])
        # Force integer ticks for x-axis
        ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))

    fig.align_ylabels(ax[:])

    bbox_props = dict(edgecolor="k", facecolor="white", linewidth=0.5)

    fig.savefig(
        data_dir + obj.capitalize() + "/" + obj + "_hyd_convergence.pdf",
        format="pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    t = time.time()
    ftool.printTime(t - t0, ms=False, where="check_convergence")

    return all_ranges


def create_input():
    """Construct all input parameters for the most recent planetary model
    according to some pre-defined conditions. The input parameters are chosen
    from random sampling. Each call of this function generates a random planet
    within the allowed ranges.
    """

    ranges = [
        [0.25, 0.75],  # Mg# tot
        [3e9, 70e9],  # P_CS
        [0.0, 0.5],  # xi_FeS
        [0.0, 0.3],  # xi_FeSi
        [0.0, 0.3],  # xi_FeO
        [1400, 2000],  # T_TBL
        [1200.0, 1600.0],  # T_0 delta_T_core
        [250.0, 350.0],
    ]  # delta_T_MTZ

    # total mass
    p = random.random()
    m = 0.05 + (5.0 - 0.05) * p

    inputs = np.empty([len(ranges) + 1])
    inputs[-1] = m

    while True:
        for j in range(len(ranges)):
            p = random.random()
            slope = ranges[j][1] - ranges[j][0]

            inputs[j] = ranges[j][0] + p * slope

        # Sample core segregation pressure according to mass
        pcsmin = 30e9 * inputs[8] ** 0.9
        pcsmax = 70e9 * inputs[8] ** 0.9

        p = random.random()
        inputs[1] = pcsmin + (pcsmax - pcsmin) * p

        # check if Si# and Fe# are in allowed range and resample if not
        ocmf = [0.0, inputs[2], inputs[3], inputs[4]]
        xi = material.mat2at_core(ocmf, xiH=0.0)
        T_CS = material.T_liquidus_pyrolite(inputs[1])
        fem = material.Fe_number_mantle(inputs[1], T_CS, xi=xi)
        sim = material.Si_number_mantle(inputs[1], T_CS, xi=xi)

        margin = 1e-3

        simmin = material.Si_number_min(1.0 - fem)
        simmax = material.Si_number_max(1.0 - fem)

        femmin = margin
        femmax = xi_Fe_mantle_max * (1.0 + margin)

        # xi_S = Material.xi_S(inputs[i][1], T_CS, xi)
        # xi_S_max = .001

        check = True

        if sum(ocmf) > 1e0:
            check = False

        if fem < femmin or fem > femmax:  # .0005 or fem > xi_Fe_mantle_max*.999:
            check = False

        if sim < simmin * (1.0 + margin) or sim > simmax * (
            1.0 - margin
        ):  # simmin*1.001 or sim > simmax*.999:
            check = False

        if check:
            # print (fem)
            break

    return inputs


def planet_sampler(N=1):
    planets = []
    all_inputs = np.empty([N, 9])
    for i in range(N):
        all_inputs[i] = create_input()

    fail_count = 0
    for i in range(N):
        mg = all_inputs[i][0]
        ocmf = [1.0, all_inputs[i][2], all_inputs[i][3], all_inputs[i][4]]
        icmf = [1.0]
        pcs = all_inputs[i][1]
        ts = all_inputs[i][5]
        T0 = all_inputs[i][6]
        dtmtz = all_inputs[i][7]

        of = -8
        im = 0
        m = all_inputs[i][8]
        o = False
        it = 1
        eps_r = 0.25
        FeS = ocmf[2]
        Pc = 9.5e10
        pl = toolkit.model_hydro(
            P_center=Pc,
            match_mass=True,
            eps_r=eps_r,
            predictor_T="none",
            ocean=o,
            ocean_frac_should=of,
            temp_jumps=[0.0, T0 * m**0.75, dtmtz, 0.0, 0.0],
            Si_number_should=0.5,
            P_surface_should=1e5,
            T_surface_should=ts,
            T_center=4000.0,
            Fe_number_mantle=0.0,
            Mg_number_should=mg,
            eps_H2O=0.0,
            iterationType=it,
            M_surface_should=m,
            predictor_P="linear",
            log=False,
            subphase_res=16,
            xi_Stv=0.0,
            acc_Mg_number=1.0e-3,
            acc_T_surface=5.0e-4,
            acc_M_surface=1.0e-4,
            xi_FeS=FeS,
            X_impurity=im,
            X_impurity_0_layers=[0.0, 0.0, im, im, 0],
            sweeps=10,
            iteration_limit=25,
            impurity=9,
            inner_core_frac=0.5,
            P_CS=pcs,
            outer_core_mat_fracs=ocmf,
            inner_core_mat_fracs=icmf,
        )
        pl.check_convergence()
        if pl.converged:
            planets.append(pl)

        else:
            fail_count += 1

    print("Converged:", int((1.0 - fail_count / N) * 100), "%")

    return planets


def plot_literature_ranges():
    ranges = np.array(
        [
            [0.329, 0.341],
            [0.338, 0.338],
            [0.327, 0.342],
            [0.317, 0.33],
            [0.317, 0.338],
            [0.32, 0.351],
        ]
    )

    labels = [
        "Zhang and Zhang 1995",
        "Aitta 2012",
        "Dumoulin et al. 2017",
        "This study (S-free)",
        "This study (Nominal)",
        "This study (S-rich)",
    ]

    add_planets = ["mercury", "venus", "earth", "mars", "moon"]
    add_cols = [
        (0.5, 0.5, 0.5),
        (0.75, 0.4, 0.2),
        (0.0, 0.0, 0.75),
        (0.75, 0.25, 0.0),
        (0.75, 0.75, 0.75),
    ]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_ylabel(r"$C/MR^2$")
    ax.tick_params(right=True)
    x1 = np.arange(0, len(labels))
    ax.set_xticks(x1)
    ax.set_ylim(0.3, 0.4)
    ax.set_xticklabels(labels, minor=False, rotation=45, ha="right", fontsize=10)

    trafo = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    rect = patches.Rectangle(
        (0, 0.313),
        1,
        0.048,
        linewidth=2,
        edgecolor="none",
        facecolor=np.asarray(add_cols[1]) + np.array([0.2, 0.3, 0.3]),
        transform=trafo,
    )

    # Add the patch to the Axes
    ax.add_patch(rect)

    for i in range(len(add_planets)):
        pl = add_planets[i]
        ax.plot(
            [0, 1], [MoI_solar[pl], MoI_solar[pl]], transform=trafo, color=add_cols[i]
        )
        ax.text(0.1, MoI_solar[pl], pl.capitalize(), va="bottom", color=add_cols[i])
    # ax.arrow(0, .35, 0, .01, head_length = 0.002, head_width = .1, fc = 'k', ec = 'k')
    # ax.text (.1, .35, 'more oxidized', rotation = 90)
    ax.text(0.1, 0.315, "Possible range (Margot et al. 2021)", color=add_cols[1])

    for i in range(len(ranges)):
        ax.plot([i, i], ranges[i], linewidth=4, color="k")

    ax.scatter([1], [0.338], color="k", zorder=10, marker="o")
    fig.savefig(
        "/home/os18o068/Documents/PHD/Abbildungen/MoI_ranges.pdf",
        format="pdf",
        bbox_inches="tight",
    )
    plt.close(fig)


def run_planet(
    res=0,
    solar_SFe=False,
    mass=1.0,
    iterationType=0,
    deltaType=0,
    initial_predictor=0,
    deltas=None,
    core_segregation_model=False,
    ocean=False,
    ocean_mass_fraction=0.1,
    Si_number_mantle=0.5,
    eps_r=0.25,
    inner_core_segregation_model=False,
    T_surface=300.0,
    P_surface=1e5,
    E_tot=1.0,
    L_int=1.0,
    Mg_number=0.5,
):

    pl = toolkit.model_hydro(
        match_mass=True,
        eps_r=eps_r,
        predictor_T="none",
        ocean=ocean,
        ocean_frac_should=np.log10(ocean_mass_fraction),
        temp_jumps=[0.0, 0.0, 0.0, 0.0, 0.0],
        Si_number_mantle=Si_number_mantle,
        P_surface_should=P_surface,
        T_surface_should=T_surface,
        T_center=None,
        P_center=None,
        E_tot_should=E_tot,
        L_int_should=L_int,
        Fe_number_mantle=0.0,
        Mg_number_should=Mg_number,
        eps_H2O=0.0,
        iterationType=iterationType,
        deltaType=deltaType,
        M_surface_should=mass,
        predictor_P="linear",
        log=False,
        subphase_res=16,
        xi_Stv=0.0,
        acc_Mg_number=1.0e-3,
        acc_T_surface=1.0e-3,
        acc_M_surface=1.0e-4,
        X_impurity=0.0,
        X_impurity_0_layers=[0.0, 0.0, 0.0, 0.0, 0],
        sweeps=10,
        iteration_limit=10,
        outer_core_mat_fracs=[1.0, 0.0, 0.0, 0.0],
        inner_core_mat_fracs=[1.0, 0.0],
        initial_predictor=0,
        deltas=deltas,
        core_segregation_model=False,
        inner_core_segregation_model=False,
    )

    return pl


def T_surf(S0, eps=0.0, alpha=0.0):
    return (S0 * (1.0 - alpha) / (4.0 * sigmaSB * (1.0 - eps / 2.0))) ** (1.0 / 4.0)


def S0_ev(t, S0):
    fac = 365 * 24 * 3600 * 1e6
    return S0 * (1.0 + 0.0 * np.sin(2 * np.pi * t / fac / 50))


def run_track(T_range=[600, 200.0], N=5, S0=1360.0, eps=0.0, alpha=0.0, dT=10.0):
    planets1 = []
    times1 = [0.0]
    planets2 = []
    times2 = [0.0]
    dTdt_list = []
    dLdT_list = []
    dTdL_list = []

    fac = 365 * 24 * 3600 * 1e6
    temps = np.linspace(T_range[0], T_range[1], N)

    # Run initial planet
    pl = run_planet(T_surface=T_range[0], P_surface=1e7, ocean_mass_frac=0.1)
    pl.compute_luminosity(S0, eps=eps, alpha=alpha)
    planets1.append(pl)
    planets2.append(pl)
    i = 0
    """
    while True:
        T = planets1[i].T_surface_is
        pl = run_planet(T_surface = T + dT, P_surface = 1e7, ocean_mass_frac = 0.1)
        planets1.append(pl)
        
        dTdL = .25 / (sigmaSB * (1. - eps / 2.) * 4. * np.pi * pl.R_surface_is**2)
        dTdL *= (1. / (sigmaSB * (1. - eps / 2.)) * (pl.luminosity / (4 * np.pi *pl.R_surface_is**2) + .25 * S0 * (1. - alpha)))**(-3./4.)
        dLdT = 16. * np.pi * pl.R_surface_is**2
        dLdT *= sigmaSB * pl.T_surface_is**3 * (1. - eps / 2.)
        dT = -pl.luminosity / dLdT
        
        dE_grav = - (planets1[i+1].E_grav - planets1[i].E_grav)
        dE_int = planets1[i+1].E_int - planets1[i].E_int
        dE = dE_grav + dE_int
        dt = - 2. * dE / (planets1[i+1].luminosity + planets1[i].luminosity)
        times1.append(times1[i] + dt / fac)
        dTdt = (planets1[i+1].T_surface_is - planets1[i].T_surface_is) / dt
        dT = dt * dTdt
        
        dTdt_list.append(dT / dt)
        dLdT_list.append(dLdT)
        dTdL_list.append(dTdL)
        
        i += 1
        
        if i > N:
            break
    """
    i = 0
    while True:
        S0 = S0_ev(times2[i], S0)
        T = planets2[i].T_surface_is - dT * np.sign(planets2[i].luminosity)
        pl = run_planet(T_surface=T, P_surface=1e7, ocean_mass_frac=0.1)
        pl.compute_luminosity(S0, eps=eps, alpha=alpha)
        planets2.append(pl)
        dLdT = 16.0 * np.pi * pl.R_surface_is**2
        dLdT *= sigmaSB * pl.T_surface_is**3 * (1.0 - eps / 2.0)

        dE_grav = -(planets2[i + 1].E_grav - planets2[i].E_grav)
        dE_int = planets2[i + 1].E_int - planets2[i].E_int
        dE = dE_grav + dE_int
        dt = 2.0 * abs(dE / (planets2[i + 1].luminosity + planets2[i].luminosity))
        dt = max(dt, 1.0)
        times2.append(times2[i] + dt)

        # Compute new dT
        dT = (
            min(abs(pl.luminosity / dLdT), abs((times2[i + 1] - times2[i]) / (dt / dT)))
            * 2.0
        )
        # dT *= 0.5
        i += 1
        if abs(pl.luminosity) < 1e14:
            break

        if i > N:
            break

    fig, ax = plt.subplots(1, 3, figsize=(10, 4))
    plt.subplots_adjust(wspace=0.5)
    ax[0].plot(
        np.asarray(times2) / fac,
        [pl.luminosity * 1e-15 for pl in planets2],
        marker="s",
        markerfacecolor="None",
        label="effective",
    )
    ax[0].plot(
        np.asarray(times2) / fac,
        [pl.intrinsic_luminosity * 1e-15 for pl in planets2],
        marker="s",
        markerfacecolor="None",
        label="intrinsic",
    )
    # ax[0].semilogy(times1, [pl.luminosity * 1e-15 for pl in planets1],
    #           marker = 's', markerfacecolor = 'None')
    ax[0].set_xlabel(r"Time [Myr]")
    ax[0].set_ylabel(r"Luminosity [PW]")

    ax[1].plot(
        np.asarray(times2) / fac,
        [pl.R_surface_is / r_earth for pl in planets2],
        marker="s",
        markerfacecolor="None",
    )
    # ax[1].plot(times1, [pl.R_surface_is / r_earth for pl in planets1],
    #           marker = 's', markerfacecolor = 'None')
    ax[1].set_xlabel(r"Time [Myr]")
    ax[1].set_ylabel(r"Radius [$R_\oplus$]")

    ax[2].plot(
        np.asarray(times2) / fac,
        [
            (pl.luminosity / (4.0 * np.pi * sigmaSB * pl.R_surface_is**2))
            ** (1.0 / 4.0)
            for pl in planets2
        ],
        marker="s",
        markerfacecolor="None",
        label="effective",
    )
    ax[2].plot(
        np.asarray(times2) / fac,
        [pl.T_surface_is for pl in planets2],
        marker="s",
        markerfacecolor="None",
        label="intrinsic",
    )
    # ax[2].plot(times1, [pl.T_surface_is for pl in planets1],
    #           marker = 's', markerfacecolor = 'None')
    ax[2].set_xlabel(r"Time [Myr]")
    ax[2].set_ylabel(r"Temperature [K]")
    ax[0].legend()
    ax[2].legend()
    plt.savefig("./time_evolution_test.pdf", format="pdf", bbox_inches="tight")
    print("dTdt =", np.asarray([dTdt_list]) * fac)
    print("dLdT =", np.asarray([dLdT_list]))
    print("dTdL =", np.asarray([dTdL_list]))

    return planets1, times1


def run_sample(
    res=0,
    obj="earth",
    add=True,
    hyd=False,
    earth_like=False,
    solar_SFe=False,
    sulfur_poor=False,
    mass=None,
    ocean_scale="lin",
    sample_mass=False,
    mass_scale="log",
    iterationType=0,
    deltaType=0,
    initial_predictor=0,
    deltas=None,
    core_segregation_model=True,
    ocean=False,
    mass_range=[0.1, 5.0],
    Fe_number_mantle_range=[0.0, 1.0],
    ocean_mass_frac_range=[-3.0, np.log10(0.5)],
    PS_range=[1e5, 1e5],
    T_TBL_range=None,
    Mg_range=None,
    Si_number_mantle=None,
    T_0_delta_T_core_range=[1200.0, 1600.0],
    delta_T_MTZ_range=[0.0, 0.0],
    save_grid=False,
    inner_core_segregation_model=False,
    P_CS_range=[30e9, 70e9],
    scale_P_CS=False,
    log_core_FeO=False,
    log_core_FeSi=False,
    external_temp_index=0,
):
    N = 2**res
    # Fix total mass if it is not to be sampled
    if not sample_mass:
        if mass == None:
            m = M_solar[obj]
        else:
            m = mass

    if hyd:
        h2o = 1.0

    else:
        h2o = 0.0

    if obj == "earth":
        PS = 1e5
        T_TBL_range = [1400.0, 1800.0]
        Mg_range = Mg_earth

    elif obj == "venus":
        PS = 9e6
        Mg_range = [0.25, 0.75]
        T_TBL_range = [1400.0, 1800.0]

    elif obj == "mars":
        PS = 1e4
        Mg_range = [0.47, 0.53]
        T_TBL_range = [1600.0, 1900.0]

    else:
        if T_TBL_range == None:
            T_TBL_range = [1400.0, 2000.0]
        if Mg_range == None:
            Mg_range = [0.2, 0.8]

    if obj == "earth":
        ranges = [
            Mg_range,  # Mg# tot
            P_CS_range,  # P_CS
            [0.0, 0.2],  # xi_FeS
            [0.0, 0.05],  # xi_FeSi
            [0.0, 0.05],  # xi_FeO
            T_TBL_range,  # T_TBL
            T_0_delta_T_core_range,  # T_0 delta_T_core
            delta_T_MTZ_range,  # delta_T_MTZ
            [0.0, 0.0],
        ]

    elif obj == "mars":
        ranges = [
            Mg_range,  # Mg# tot
            P_CS_range,  # P_CS
            [0.3, 0.5],  # xi_FeS
            [5e-6, 2e-4],  # xi_FeSi
            [0.01, 0.1],  # xi_FeO
            T_TBL_range,  # T_TBL
            T_0_delta_T_core_range,  # T_0 delta_T_core
            delta_T_MTZ_range,  # delta_T_MTZ
            [0.0, 0.0],
        ]

    elif obj == "venus":
        ranges = [
            Mg_range,  # Mg# tot
            P_CS_range,  # P_CS
            [0.0, 0.5],  # xi_FeS
            [0.0, 0.3],  # xi_FeSi
            [0.0, 0.2],  # xi_FeO
            T_TBL_range,  # T_TBL
            T_0_delta_T_core_range,  # T_0 delta_T_core
            delta_T_MTZ_range,  # delta_T_MTZ
            [0.0, 0.0],
        ]
    else:
        ranges = [
            Mg_range,  # Mg# tot
            P_CS_range,  # P_CS
            [0.0, 0.0],  # [.0, 1.], #xi_FeS
            [0.0, 0.0],  # [.0, .5], #xi_FeSi
            [0.0, 0.0],  # [.0, .5], #xi_FeO
            T_TBL_range,  # T_TBL
            T_0_delta_T_core_range,  # T_0 delta_T_core
            delta_T_MTZ_range,  # delta_T_MTZ
            ocean_mass_frac_range,
        ]

        if log_core_FeO:
            ranges[4] = [-10, -2]

        if log_core_FeSi:
            ranges[3] = [-12, -4]

    if earth_like:
        ranges[0] = Mg_earth

        if not solar_SFe and not sulfur_poor:
            ranges[2] = [0.08, 0.15]

    if solar_SFe:
        ranges[2] = [0.2, 0.5]

    if sulfur_poor:
        ranges[2] = [0.0, 1e-6]

    if solar_SFe or earth_like:
        ranges[4] = [0.0, 0.2]
        if log_core_FeO:
            ranges[4] = [-10, -1]

    if add:
        try:
            planets = ftool.load_objects(obj + "_cali.pkl")
        except FileNotFoundError:
            print("File does not exist")
            planets = []
    else:
        planets = []

    inputs = np.empty([N, len(ranges)])
    masses = np.empty([N])
    Fe_numbers_mantle = np.empty([N])
    Si_numbers_mantle = np.empty([N])
    put_ocean = np.array([ocean] * N)

    # Prepare all input parameters for N planets
    for i in range(N):
        print("\n")
        # Sample total mass within given mass range
        if sample_mass:
            p = random.random()
            if mass_scale == "lin":
                slope = mass_range[1] - mass_range[0]
                m = mass_range[0] + p * slope

            elif mass_scale == "log":
                slope = np.log10(mass_range[1] / mass_range[0])
                m = np.log10(mass_range[0]) + p * slope
                m = 10**m

            masses[i] = m

        if put_ocean[i]:
            # Probe ocean mass fraction in given range
            p = random.random()
            slope = ocean_mass_frac_range[1] - ocean_mass_frac_range[0]
            of = ocean_mass_frac_range[0] + slope * p
            # of = np.log10(of)

        else:
            of = -10.0
        print("of = ", of)
        # If selected ocean mass fraction is below 0.1% don't model it
        if of < -3:
            of = -10
            put_ocean[i] = False
            print("turn off ocean")

        inputs[i][8] = of

        check_count = 0
        while True:
            for j in range(len(ranges)):
                p = random.random()
                slope = ranges[j][1] - ranges[j][0]
                inputs[i][j] = ranges[j][0] + p * slope

                # For outer core composition previous values restrict current ones
                if j == 3 or j == 4:
                    # print ("inputs before =", inputs[i][2:5])
                    inputs[i][j] = min(
                        inputs[i][j], max(1.0 - sum(inputs[i][2:j]), 0.0)
                    )
                    # print ("inputs after =", inputs[i][2:5])
            # Rescale FeO or FeSi content if logarithmically probed
            if log_core_FeO:
                inputs[i][4] = 10 ** inputs[i][4]

            if log_core_FeSi:
                inputs[i][3] = 10 ** inputs[i][3]

            # The core segregation pressure must be scaled to the total
            # Scaling law P / P_Earth = (M / M_Earth)^eta
            if scale_P_CS:
                eta = 1.0
                inputs[i][1] *= m * (1.0 - 10 ** inputs[i][8]) ** eta

            # check if Si# and Fe# are in allowed range and resample if
            # Note that mat2at_core computese first element of the array (Fe
            # content) as a function of the other elements such that sum = 1.
            ocmf = [0.0, inputs[i][2], inputs[i][3], inputs[i][4]]
            # print ("\nocmf =", ocmf)
            # print ("pcs =", inputs[i][1]*1e-9)
            icmf = [1.0]
            xi = material.mat2at_core(ocmf, xiH=0.0)
            # print ("xi =", xi)
            T_CS = material.T_liquidus_pyrolite(inputs[i][1])
            fem = material.Fe_number_mantle(inputs[i][1], T_CS, xi=xi)
            sim = material.Si_number_mantle(inputs[i][1], T_CS, xi=xi)
            # print ('P_CS, fem, sim =', inputs[i][1]*1e-9, fem, sim)
            Fe_numbers_mantle[i] = fem
            Si_numbers_mantle[i] = sim
            margin = 1e-4
            # oxide_fracs = Material.compute_oxide_fractions(sim, fem)

            simmin = material.Si_number_min(fem)
            simmax = material.Si_number_max(fem)

            femmin = margin
            femmax = xi_Fe_mantle_max * (1.0 + margin)

            # Impose additional constraint for the Earths mantle to improve the
            # statistics by not probing ranges that will not match the PREM
            # data anyways.
            if obj == "earth":
                simmin = max(0.45, simmin)
                simmax = min(0.55, simmax)

                femmin = max(0.01, femmin)
                femmax = min(0.2, femmax)

            # xi_S = Material.xi_S(inputs[i][1], T_CS, xi)
            # xi_S_max = .001

            # print ("P_CS = ", inputs[i][1] * 1e-9)
            # print ('femmin, femmax, fem =', femmin, femmax, fem)
            # print ('simmin, simmax, sim =', simmin, simmax, sim)
            check = True

            if sum(ocmf) > 1e0:
                check = False

            if fem < femmin or fem > femmax:  # .0005 or fem > xi_Fe_mantle_max*.999:
                check = False

            if fem < Fe_number_mantle_range[0] or fem > Fe_number_mantle_range[1]:
                check = False

            if sim < simmin * (1.0 + margin) or sim > simmax * (
                1.0 - margin
            ):  # simmin*1.001 or sim > simmax*.99
                check = False

            if check:
                check_count += 1
            # if sum(ocmf) < .03:
            #   check = False
            # If CS model is inactive any combination within the specified ranges
            # for the input parameters will be accepted by the model. Otherwise
            # the CS model must set check = True.
            if not core_segregation_model or check:
                # print ('----------------')
                break

    print("\n--------------\n")
    # sys.exit()
    # Run N models with the given input parameters
    for i in range(N):
        if sample_mass:
            m = masses[i]
        else:
            pass

        PS = np.log10(PS_range[1]) - np.log10(PS_range[0])
        PS *= random.random()
        PS += np.log10(PS_range[0])
        PS = 10**PS
        mg = inputs[i][0]
        ocmf = [1.0, inputs[i][2], inputs[i][3], inputs[i][4]]
        # print ("ocmf =", ocmf)
        icmf = [1.0]
        pcs = inputs[i][1]  # Material.P_CS(m, P0 = inputs[i][1])
        ts = inputs[i][5]
        # ammas = [inputs[i][6], inputs[i][6], 1.96, 1.26, 1.0]
        T0 = inputs[i][6]
        dtmtz = inputs[i][7]
        of = inputs[i][8]

        # The Fe and Si content in the mantle defined here have no impact on
        # the model if core_segregation_model is switched on because they
        # will be computed self-consistently in the fortran routine during
        # integration. If core_segregation_model is off the values defined here
        # will be taken as the mantle composition in the fortran routine. But
        # the Fe and Si content in the mantle are used to predict P_C, T_C and
        # M_core prior to the iteration. Therefore, in order to optimize the
        # proceedure the correct values should be passed here even if CS model
        # is activated.
        if not core_segregation_model:
            pcs = 0e0  # Fortran can't handle None input so set it to zero
            p = random.random()
            slope = Fe_number_mantle_range[1] - Fe_number_mantle_range[0]
            fem = Fe_number_mantle_range[0] + slope * p
            simmin = material.Si_number_min(fem)
            simmax = material.Si_number_max(fem)

            if Si_number_mantle == None:
                # Sample Si# in mantle within allowed range
                # NOTE: fortplanet updates the Si# number in the layer self-consistently
                # in the minor_bisection subroutine after core-mantle transition is reached.
                # Therefore, setting it here does not actually do anything if CS is on. The true
                # value for the Si# in the mantle is an output from the planet integration
                # and can be accessed afterwards as pl.Si_number_mantle
                p = random.random()
                sim = simmin + (simmax - simmin) * p
                # fem = 0.
                # sim = .5
            else:
                sim = Si_number_mantle

        else:
            sim = Si_numbers_mantle[i]
            fem = Fe_numbers_mantle[i]

        print("Fe mantle in run_sample =", fem)
        im = 0
        eps_r = 0.25
        Pc = None
        Tc = None
        pl = toolkit.model_hydro(
            P_center=Pc,
            match_mass=True,
            eps_r=eps_r,
            predictor_T="none",
            ocean=put_ocean[i],
            ocean_frac_should=of,
            temp_jumps=[0.0, T0 * (m * (1.0 - 10**of)) ** 0.75, dtmtz, 0.0, 0.0],
            Si_number_mantle=sim,
            P_surface_should=PS,
            T_surface_should=ts,
            T_center=Tc,
            Fe_number_mantle=fem,
            Mg_number_should=mg,
            eps_H2O=h2o,
            iterationType=iterationType,
            deltaType=deltaType,
            M_surface_should=m,
            predictor_P="linear",
            log=False,
            subphase_res=16,
            xi_Stv=0.0,
            acc_Mg_number=1.0e-3,
            acc_T_surface=1.0e-3,
            acc_M_surface=1.0e-4,
            X_impurity=im,
            X_impurity_0_layers=[0.0, 0.0, im, im, 0],
            sweeps=5,
            iteration_limit=10,
            impurity=9,
            P_CS=pcs,
            outer_core_mat_fracs=ocmf,
            inner_core_mat_fracs=icmf,
            initial_predictor=initial_predictor,
            deltas=deltas,
            external_temp_index=external_temp_index,
            core_segregation_model=core_segregation_model,
            inner_core_segregation_model=inner_core_segregation_model,
        )
        planets.append(pl)
        pl.prt()

    if not save_grid:
        ftool.save_objects(planets, obj + "_cali.pkl")

    return planets


def run_curve(
    obj="pure_iron",
    N=3,
    eps_r=0.5,
    ocean=False,
    of=-10,
    sim=0.5,
    fem=0.0,
    mg=1e-5,
    ts=300.0,
    PS=1e5,
    iterationType=0,
    deltaType=1,
    massRange=[0.5, 1.5],
):
    planets = []
    masses = np.logspace(np.log10(massRange[0]), np.log10(massRange[1]), N)
    for i in range(N):
        m = masses[i]
        pl = toolkit.model_hydro(
            P_center=5e11,
            match_mass=True,
            eps_r=eps_r,
            predictor_T="none",
            ocean=False,
            contents=[[9]],
            fractions=[[1.0]],
            ocean_frac_should=of,
            temp_jumps=[0.0, 0.0, 0.0, 0.0],
            Si_number_mantle=sim,
            P_surface_should=PS,
            T_surface_should=ts,
            T_center=5000.0,
            Fe_number_mantle=fem,
            Mg_number_should=mg,
            iterationType=iterationType,
            deltaType=deltaType,
            M_surface_should=m,
            predictor_P="linear",
            log=False,
            subphase_res=16,
            xi_Stv=0.0,
            acc_Mg_number=1.0e-3,
            acc_T_surface=1.0e-3,
            acc_M_surface=1.0e-4,
            sweeps=10,
            iteration_limit=20,
            initial_predictor=0,
            core_segregation_model=False,
            inner_core_segregation_model=False,
        )
        planets.append(pl)
        pl.prt()

    ftool.save_objects(planets, obj + "_cali.pkl")

    return planets


def extract_profiles(planets):
    """Takes Planet objects and extracts interior profiles into numpy
    array. Convenction is [Radius, Temperature, Pressure, Density]
    """
    profiles = []

    for pl in planets:
        pl.trim_profiles()
        profiles.append(list(pl.profiles[0:4]))

    print(len(profiles))
    return np.asarray(profiles)


def extract_params(planets):
    """Takes Planet object and extracts all parameters to an array"""
    N = len(planets)
    all_all_out = np.empty([N, len(params)])

    i = -1
    while i < N - 1:
        i += 1
        pl = planets[i]

        pl.update_oxide_fractions()
        pl.extract_bulk_structure()
        pl.complete_layer_properties()

        # Compute mole fractions of oxides in the mantle
        xiFeO = pl.xi_FeO_mantle
        xiSiO2 = pl.xi_SiO2_mantle
        xiMgO = pl.xi_MgO_mantle

        if not abs(xiMgO + xiFeO + xiSiO2 - 1.0) < 0.0001:
            print("WARNING: Inconsistent oxide fractions in the mantle!")

        all_all_out[i][0] = pl.Mg_number_is
        all_all_out[i][1] = pl.R_surface_is / r_earth
        all_all_out[i][2] = pl.MOI_is
        all_all_out[i][3] = pl.Si_number_is
        all_all_out[i][4] = pl.Si_number_mantle
        all_all_out[i][5] = pl.xi_Fe_mantle
        all_all_out[i][6] = pl.P_CS * 1e-9
        all_all_out[i][7] = pl.T_CS
        all_all_out[i][8] = pl.X_all_core[2]  # S
        all_all_out[i][9] = pl.X_all_core[3]  # Si
        all_all_out[i][10] = pl.X_all_core[4]  # O
        all_all_out[i][11] = pl.fractions[1][1]  # FeS
        all_all_out[i][12] = pl.fractions[1][2]  # FeSi
        all_all_out[i][13] = pl.fractions[1][3]  # FeO
        all_all_out[i][14] = pl.temp_jumps[1]  # T0_CMB
        all_all_out[i][15] = pl.layer_properties[1]["R_outer"] / 1000  # R CMB
        all_all_out[i][16] = pl.layer_properties[0]["R_outer"] / 1000  # R ICB
        try:
            all_all_out[i][17] = pl.layer_properties[3]["R_outer"] / 1000  # R HMI
        except IndexError:
            all_all_out[i][17] = None
        all_all_out[i][18] = pl.T_surface_is
        all_all_out[i][19] = xiSiO2  # xi_SiO2 mantle
        all_all_out[i][20] = xiFeO  # xi_FeO mantle
        all_all_out[i][21] = pl.T_center
        all_all_out[i][22] = pl.P_center * 1e-9
        all_all_out[i][23] = pl.logfO2
        all_all_out[i][24] = pl.temp_jumps[2]  # DeltaT_MTZ
        all_all_out[i][25] = pl.M_H2O_is
        all_all_out[i][26] = pl.M_H2O_mantle
        all_all_out[i][27] = pl.xi_all_core[0] / (sum(pl.xi_all_core))
        all_all_out[i][28] = pl.M_core_is / pl.M_surface_is * m_earth
        all_all_out[i][29] = pl.layer_properties[1]["T_outer"]
        all_all_out[i][30] = pl.layer_properties[1]["P_outer"]
        all_all_out[i][31] = pl.layer_properties[1]["rho_outer"]
        all_all_out[i][32] = pl.layer_properties[0]["T_outer"]
        all_all_out[i][33] = pl.layer_properties[0]["P_outer"]
        all_all_out[i][34] = pl.layer_properties[0]["rho_outer"]
        all_all_out[i][35] = pl.layer_properties[2]["T_outer"]
        all_all_out[i][36] = pl.layer_properties[2]["P_outer"]
        all_all_out[i][37] = pl.layer_properties[2]["rho_outer"]
        all_all_out[i][38] = 1.0 - xiSiO2 - xiFeO  # xi_MgO mantle
        all_all_out[i][39] = pl.layer_properties[2]["R_outer"] / 1000
        all_all_out[i][40] = pl.fractions[1][1] / (
            1.0 - pl.fractions[1][0] + sum(pl.fractions[1])
        )
        all_all_out[i][41] = pl.xi_S_mantle
        all_all_out[i][42] = pl.S_count / pl.Fe_count
        all_all_out[i][43] = pl.Si_count / pl.Fe_count
        all_all_out[i][44] = pl.Mg_count / pl.Fe_count
        try:
            all_all_out[i][45] = pl.S_count / pl.Si_count
        except ZeroDivisionError:
            all_all_out[i][45] = np.nan
        all_all_out[i][46] = pl.M_surface_is / m_earth
        all_all_out[i][47] = (
            pl.layer_properties[0]["indigenous_mass"] / m_earth / pl.M_core_is
        )  # inner core mass fraction
        all_all_out[i][48] = 10**pl.ocean_frac_is
        all_all_out[i][49] = pl.P_surface_is

        try:
            all_all_out[i][50] = pl.masses[0]  # liq. water mass fraction wrt. tot. mass
            all_all_out[i][51] = pl.masses[1]
            all_all_out[i][52] = pl.masses[2]
            all_all_out[i][53] = pl.heights[0]
            all_all_out[i][54] = pl.heights[1]
            all_all_out[i][55] = pl.heights[2]
            all_all_out[i][56] = pl.masses[0] / sum(
                pl.masses
            )  # liq. water mass frac. wrt. water
            all_all_out[i][57] = pl.masses[1] / sum(pl.masses)
            all_all_out[i][58] = pl.masses[2] / sum(pl.masses)
            all_all_out[i][59] = pl.vols[0] / V_ocean_earth
            all_all_out[i][60] = pl.vols[1] / V_ocean_earth
            all_all_out[i][61] = pl.vols[2] / V_ocean_earth
        except ZeroDivisionError:
            all_all_out[i][50] = 0.0  # liq. water mass fraction wrt. tot. mass
            all_all_out[i][51] = 0.0
            all_all_out[i][52] = 0.0
            all_all_out[i][53] = 0.0
            all_all_out[i][54] = 0.0
            all_all_out[i][55] = 0.0
            all_all_out[i][56] = 0.0  # liq. water mass frac. wrt. water
            all_all_out[i][57] = 0.0
            all_all_out[i][58] = 0.0
            all_all_out[i][59] = 0.0
            all_all_out[i][60] = 0.0
            all_all_out[i][61] = 0.0
        except AttributeError:
            all_all_out[i][50] = 0.0  # liq. water mass fraction wrt. tot. mass
            all_all_out[i][51] = 0.0
            all_all_out[i][52] = 0.0
            all_all_out[i][53] = 0.0
            all_all_out[i][54] = 0.0
            all_all_out[i][55] = 0.0
            all_all_out[i][56] = 0.0  # liq. water mass frac. wrt. water
            all_all_out[i][57] = 0.0
            all_all_out[i][58] = 0.0
            all_all_out[i][59] = 0.0
            all_all_out[i][60] = 0.0
            all_all_out[i][61] = 0.0
        all_all_out[i][62] = (
            pl.M_surface_is / m_earth / (pl.R_surface_is / r_earth) ** 2
        )  # Surface gravity
        try:
            all_all_out[i][63] = "".join([str(int(hs)) for hs in pl.hydro_structure])
        except AttributeError:
            all_all_out[i][63] = None

        all_all_out[i][64] = pl.layer_properties[3]["T_outer"]
        all_all_out[i][65] = pl.layer_properties[3]["P_outer"]
        all_all_out[i][66] = pl.layer_properties[3]["rho_outer"]
        try:
            all_all_out[i][67] = "".join([str(int(hs)) for hs in pl.bulk_structure])
        except AttributeError:
            all_all_out[i][67] = None

    return all_all_out


def calibrate(
    all_planets,
    obj="earth",
    plot=False,
    reldev_R=0.01,
    reldev_C=0.01,
    MoI_real=None,
    add_matches=False,
    only_outer_core=False,
    R_real=None,
    only_earth_comp=False,
    scatter_only_match=False,
    only_solar_SFe=False,
    plot_label="All models",
):

    try:
        all_matches = ftool.load_objects(obj + "_cali_match.pkl")

    except FileNotFoundError:
        all_matches = []

    planets = []
    N_tot = len(all_planets)
    N = N_tot
    print(N, "planets loaded.")
    not_converged_count = 0

    # Only consider planets which converged
    for i in range(N):
        all_planets[i].check_convergence()

        if not all_planets[i].converged:
            not_converged_count += 1

        else:
            planets.append(all_planets[i])

    N = len(planets)
    print("# converged planets =", N)
    if MoI_real == None:
        MoI_real = MoI_solar[obj]

    if R_real == None:
        R_real = R_solar[obj]

    print("R / MoI = ", R_real, MoI_real)
    deltas_r = np.empty([N, 5])
    deltas_d = np.empty([N, 6])
    match_count = 0
    indices = []

    all_all_out = extract_params(planets)
    # all_all_out = np.empty([N, len(params)])

    i = -1
    while i < N - 1:
        i += 1
        pl = planets[i]

        pl.update_oxide_fractions()
        # pl.compute_sulfur_in_mantle()
        pl.trim_profiles()
        props = pl.layer_properties

        # Compute rel. deviation from MoI and R of the real Planet
        dev_R = (pl.R_surface_is / r_earth - R_real) / R_real
        dev_MoI = (pl.MOI_is - MoI_real) / MoI_real

        # Check if inner core exists
        core_check = True

        # Check if bulk composition is Earth-like and/or S/Fe is solar
        comp_check = True

        go = False

        if obj == "earth":
            transitions = np.array([p["R_outer"] * 1e-3 for p in props])

            densities = np.empty([6])

            # loop over layers
            for l in range(len(props) - 1):
                densities[2 * l] = props[l]["rho_inner"]
                densities[2 * l + 1] = props[l]["rho_outer"]

            # Compute difference from PREM in km
            dr = transitions - readPREM.transition_radii
            dd = (densities - readPREM.layer_densities) / readPREM.layer_densities

            # Check if all layer transitions match with PREM
            count_r = 0
            for j in range(len(dr) - 2):
                deltas_r[i][j] = dr[j]
                if j == 0:
                    err = 200.0
                else:
                    err = 200.0
                if abs(dr[j]) < err:
                    count_r += 1

            # Check if all densities match with PREM
            count_d = 0
            for j in range(len(dd)):
                deltas_d[i][j] = dd[j]
                if j < 2:
                    err = 0.05
                else:
                    err = 0.05
                if abs(dd[j]) < err:
                    count_d += 1

            if (
                count_r == len(dr) - 2
                and abs(dev_R) < reldev_R
                and abs(dev_MoI) < reldev_C
                and count_d == len(dd)
                and core_check
            ):
                go = True

        else:
            if (
                core_check
                and comp_check
                and abs(dev_R) < reldev_R
                and abs(dev_MoI) < reldev_C
            ):
                go = True

        if go:

            # planets[i].Plot(spec='match_test', prem=True)
            match_count += 1
            indices.append(i)

    N_match = len(indices)

    # Gather all relevant parameters for the subsequent analysis
    i = -1
    while i < N_match - 1:
        i += 1
        pl = planets[indices[i]]

        if add_matches:
            all_matches.append(pl)

        # Compute mole fractions of oxides in the mantle
        xiFeO = 1.0 - pl.Si_number_mantle
        xiFeO /= (1.0 - pl.xi_Fe_mantle) / pl.xi_Fe_mantle + 1.0 - pl.Si_number_mantle
        xiSiO2 = (1.0 - xiFeO) * pl.Si_number_mantle
        xiMgO = 1.0 - xiFeO - xiSiO2

        if not abs(xiMgO + xiFeO + xiSiO2 - 1.0) < 0.0001:
            print("WARNING: Inconsistent oxide fractions in the mantle!")

    if add_matches:
        ftool.save_objects(all_matches, obj + "_cali_match.pkl")

    print("-----")
    print(N, "planets converged.")
    print(N_match, "matches found.")
    print(
        not_converged_count,
        "not converged (",
        int(100 * not_converged_count / (N + not_converged_count)),
        "%)",
    )
    print("-----")
    return all_all_out, planets, indices


def filter_data(data, param, val_should, delta, relative=True):
    """Filters the input data using the parameter param. All data points for
    which the value of param is within val_should +- delta are considered. delta is
    either the absolut deviation or relative deviation. This can be set by
    relative = True/False.
    """
    i = -1

    indices = []

    while i < len(data) - 1:
        i += 1

        # Check if condition is met for current data point
        val_is = data[i][param]

        if relative:
            dev = (val_is - val_should) / val_should

        else:
            dev = val_is - val_should

        if abs(dev) <= delta:
            indices.append(i)

    data_dummy = np.empty([len(indices), len(data[0])])

    for i in range(len(indices)):
        data_dummy[i] = data[indices[i]]

    return data_dummy


def plot_correlations(all_all_out, indices, obj="earth", spec="all"):
    N_match = len(indices)
    xi_FeO = np.linspace(0.01, 0.15, 10)
    reldevs = [0.1, 0.2]

    core_devs = np.empty([len(reldevs), len(xi_FeO)])

    # gather plot data
    all_out = np.empty([N_match, len(params)])
    i = -1
    while i < N_match - 1:
        i += 1
        all_out[i][:] = all_all_out[indices[i]][:]
    """
    for i in range(len(reldevs)):
        for j in range(len(xi_FeO)):
            dat = filter_data(all_all_out, 19, xi_FeO[j], reldevs[i])
            d = np.nanmax(dat.T[15]) - np.nanmin(dat.T[15])
            core_devs[i][j] = d
    """

    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    # fig1, ax1 = plt.subplots()

    plt.subplots_adjust(wspace=0.25, hspace=0.25)

    param_indices = [[22, 15], [19, 15]]
    x_labels = [r"$\log{f_{\rm O_2}}$", r"$X_{\rm FeO}^{\rm Mantle} \ \rm [mol\%]$"]
    y_labels = [r"$R_{\rm Core} \ \rm [km]$", r"$R_{\rm Core} \ \rm [km]$"]
    scales = [1.0, 100.0]
    axes_lims = [[[-6, 0], [2800.0, 4400.0]], [[0.0, 15.0], [2800.0, 4400.0]]]
    earth_logfo2 = [-2.8, -1.8]
    earth_xiFeO = [3.4, 7.2]

    x = np.array([earth_logfo2, earth_xiFeO])

    fnts1 = 12
    fnts2 = 14
    """
    for i in range(len(reldevs)):
        ax1.plot(xi_FeO, core_devs[i])
    """
    for i in range(2):
        ax[i].set_xlim(axes_lims[i][0])
        ax[i].set_ylim(axes_lims[i][1])
        ax[i].tick_params(
            which="both", direction="out", right=True, top=True, labelsize=fnts2, pad=10
        )
        ax[i].set_xlabel(x_labels[i], fontsize=fnts2)
        ax[i].set_ylabel(y_labels[i], fontsize=fnts2)
        sc = ax[i].scatter(
            all_out.T[param_indices[i][0]] * scales[i],
            all_out.T[param_indices[i][1]],
            c=all_out.T[2],
            cmap="copper",
        )

        # Create a Rectangle patch
        dx = x[i][1] - x[i][0]

        trans = transforms.blended_transform_factory(ax[i].transData, ax[i].transAxes)

        rect = patches.Rectangle(
            (x[i][0], 0),
            dx,
            1.0,
            linewidth=2,
            edgecolor="none",
            facecolor=(0.75, 0.75, 0.75),
            transform=trans,
            zorder=0,
        )

        # Add the patch to the Axes
        ax[i].add_patch(rect)

        ax[i].text(
            sum(x[i]) / 2.0,
            0.8,
            r"Earth reference",
            transform=trans,
            ha="center",
            size=fnts2,
        )

    inset_ax = inset_axes(
        ax[-1],
        height="40%",
        width="5%",
        loc=1,
        bbox_to_anchor=(-0.1, 0, 0.9, 0.9),
        bbox_transform=ax[-1].transAxes,
        borderpad=0,
    )

    cbar = plt.colorbar(sc, cax=inset_ax)
    cbar.ax.tick_params(labelsize=fnts1)
    cbar.set_label(r"$C/MR^2$", fontsize=fnts1)

    fig.savefig(
        spec + "/" + obj + "_corrs_" + spec + ".png",
        format="png",
        bbox_inches="tight",
        dpi=240,
    )

    # fig1.savefig(spec+'/'+obj+'_core_devs_'+spec+'.png', format='png',
    #        bbox_inches='tight', dpi = 240)


def plot_models(
    all_all_out,
    planets,
    indices,
    y_real=None,
    obj="earth",
    reldev_x=0.01,
    reldev_y=0.01,
    plot_label="All models",
    spec="all",
    x_real=None,
    x_range=0.1,
    y_range=0.1,
    plot_rect=True,
    plot_cross=True,
    mass=None,
    xy=[1, 2],
    xy_scale=["log", "linear"],
    x_lims=[0.1, 5.0],
    y_lims=[0.5, 1.75],
    plot_PT=False,
    profile_x_lims=[0.0, 1.0],
    plot_PREM=True,
    plot_CMB=True,
    profile_y_lims=[[1000.0, 6000.0], [0.0, 400.0], [1.0, 14.0]],
    profile_cbar_loc=3,
    plot_real=True,
    colormap="viridis",
    colormap_arch="jet",
    real_filter=False,
    max_mass_err=0.1,
):

    if plot_PT:
        import phase_transitions_water_Wagner2002 as phase

        fig2, ax2 = phase.plot_the_shit()

    if mass == None:
        mass = M_solar[obj]

    N = len(planets)
    N_match = len(indices)

    # gather plot data
    all_out = np.empty([N_match, len(params)])
    i = -1
    while i < N_match - 1:
        i += 1
        all_out[i][:] = all_all_out[indices[i]][:]

    # Define x and y values of real object
    if y_real == None:
        y_real = MoI_solar[obj]

    if x_real == None:
        x_real = R_solar[obj]

    # Allowed range
    dx = reldev_x * x_real
    dy = reldev_y * y_real

    # Define plot axes ranges
    mins = np.array(
        [
            [0.1, 0.3, 0.0],
            [100.0, 4, 0.0],
            [0, -6, -3],
            [0.0, 0.3, 1400.0],
            [2000.0, 30.0, 0.2],
        ]
    )
    maxs = np.array(
        [
            [0.9, 0.6, 1.0],
            [2000, 9, 1.0],
            [0.3, -4, -1],
            [0.5, 0.6, 1800],
            [8000.0, 70, 0.4],
        ]
    )
    overview_params = [[0, 3, 27], [17, 48, 47], [8, 9, 10], [5, 18, 17], [7, 6, 2]]
    param_scales = [
        ["lin", "lin", "lin"],
        ["lin", "log", "lin"],
        [
            "lin",
            "log",
            "log",
        ],
        ["lin", "lin", "lin"],
        ["lin", "lin", "lin"],
    ]

    cbar_labels = [
        [
            all_param_labels[overview_params[i][j]]
            for j in range(len(overview_params[i]))
        ]
        for i in range(len(overview_params))
    ]

    fnts = 15

    fig, ax = plt.subplots(
        len(overview_params), len(overview_params[0]), figsize=(16, 20)
    )

    plt.subplots_adjust(wspace=0.5, hspace=0.25)
    x, y = xy

    if plot_real:
        arch = read_exoplanets.Archive()

    print("reals = ", x_real, y_real)
    # Plot some parameters in R-MoI plane
    for i in range(len(overview_params)):
        for j in range(len(overview_params[0])):
            ax[i][j].tick_params(labelsize=fnts)
            ax[i][j].set_xlabel(all_param_labels[x], fontsize=fnts)
            ax[i][j].set_ylabel(all_param_labels[y], fontsize=fnts)
            # ax[i][j].set_aspect('equal')
            ax[i][j].set_xscale(xy_scale[0])
            ax[i][j].set_yscale(xy_scale[1])
            ax[i][j].set_xlim(x_lims[0], x_lims[1])
            ax[i][j].set_ylim(y_lims[0], y_lims[1])

            # Overplot real exoplanets with error bars
            if plot_real:
                if real_filter:
                    arch.filter(dM_max=max_mass_err)

                arch.plot(ax=ax[i][j], cmap=colormap_arch)
            # Plot line at real x and y
            if plot_cross:
                print("Plotting cross...")
                trafo = transforms.blended_transform_factory(
                    ax[i][j].transData, ax[i][j].transAxes
                )
                ax[i][j].plot(
                    [x_real, x_real], [0.0, 1.0], color="grey", transform=trafo
                )

                trafo = transforms.blended_transform_factory(
                    ax[i][j].transAxes, ax[i][j].transData
                )
                ax[i][j].plot(
                    [0.0, 1.0], [y_real, y_real], color="grey", transform=trafo
                )

            if param_scales[i][j] == "lin":
                c = all_all_out.T[overview_params[i][j]]
            elif param_scales[i][j] == "log":
                c = np.log10(all_all_out.T[overview_params[i][j]])
            sc = ax[i][j].scatter(
                all_all_out.T[x],
                all_all_out.T[y],
                c=c,
                vmin=mins[i][j],
                vmax=maxs[i][j],
                cmap=colormap,
                zorder=0,
            )

            if plot_rect:
                # Create a Rectangle patch
                x_rect = x_real - dx
                y_rect = y_real - dy
                rect = patches.Rectangle(
                    (x_rect, y_rect),
                    2 * dx,
                    2 * dy,
                    linewidth=2,
                    edgecolor="r",
                    facecolor="none",
                )

                # Add the patch to the Axes
                ax[i][j].add_patch(rect)

            # Add colorbar to the subplot
            ftool.my_colorbar(sc, pad=0.25, label=cbar_labels[i][j])

    # Value range for the colors of the PT profiles
    if plot_PT:
        PT_value_range_cmap = [all_out.T[47].min(), all_out.T[47].max()]

    # plot all exoplanets and selected exoplanets
    fnts_arch = 18
    fig3, ax3 = plt.subplots(1, 2, figsize=(16, 6))
    c = all_all_out.T[48]
    ax3[0].scatter(
        all_all_out.T[x],
        all_all_out.T[y],
        c=c,
        vmin=0.0,
        vmax=1.0,
        cmap=colormap,
        zorder=0,
    )
    sc = ax3[1].scatter(
        all_all_out.T[x],
        all_all_out.T[y],
        c=c,
        vmin=0.0,
        vmax=1.0,
        cmap=colormap,
        zorder=0,
    )

    arch.plot(ax=ax3[0], cmap="jet", fnts1=fnts_arch)
    # cbar_ax = fig3.add_axes([1, 0, 0.01, 1])
    inset_ax_cbar = ax3[1].inset_axes((1.1, 0.0, 0.05, 1))
    cbar = plt.colorbar(sc, cax=inset_ax_cbar)
    cbar.ax.set_ylabel(
        r"${\rm Water \ mass \ fraction} \ w_{\rm H_2O}$", fontsize=fnts_arch
    )
    cbar.ax.tick_params(labelsize=fnts_arch)

    # Plot T-profiles, P-profiles and rho-profiles
    n = 48  # parameter for which profiles should be colored
    c1 = np.min(all_out.T[n])
    c2 = np.max(all_out.T[n])

    # prepare colors
    c = (all_out.T[n] - c1) / (c2 - c1)
    cols = np.array([[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.5, 0.5, 1.0]])
    cols = np.array([[0.0, 0.0, 0.0], [0.5, 0.25, 0.0], [1.0, 0.5, 0.0]])
    cols = np.array([[0.0, 0.0, 0.0], [0.5, 0.125, 0.0], [1.0, 0.25, 0.0]])
    cols = np.array([[0, 0, 0], [0, 0.25, 0.5], [0, 0.5, 1]])
    cmap = ftool.my_dev_cmap(N_col=32, border_colors=cols)
    cmap = plt.cm.get_cmap(colormap)

    sm = plt.cm.ScalarMappable(cmap=colormap)

    y_labels = [
        r"$\rm Temperature \ [K]$",
        r"$\rm Pressure \ [GPa]$",
        r"$\rm Density \ [g/cm^3]$",
    ]

    profile_scales = [1.0, 1e-9, 1e-3]

    fig1, ax1 = plt.subplots(1, 3, figsize=(16, 4))

    # inset_ax = inset_axes(ax1[-1], height="40%", width="5%", loc=profile_cbar_loc,
    #                     borderpad=2, bbox_to_anchor = [.5, 0, 0.5, 1.],
    #                    bbox_transform = ax1[-1].transAxes)
    inset_ax = ax1[-1].inset_axes((0.6, 0.5, 0.05, 0.4))
    plt.subplots_adjust(wspace=0.3)
    cbar = plt.colorbar(sm, cax=inset_ax, ticks=[0.0, 0.5, 1.0])

    ticklabels = np.array(
        [
            all_out.T[n].min(),
            (all_out.T[n].min() + all_out.T[n].max()) / 2.0,
            all_out.T[n].max(),
        ]
    ).round(3)

    cbar.ax.set_yticklabels(ticklabels)
    cbar.ax.tick_params(labelsize=fnts)
    cbar.set_label(r"$w_{\rm H_2O}$", fontsize=fnts)

    try:
        for i in range(3):
            ax1[i].tick_params(
                which="both",
                direction="out",
                top="on",
                right="on",
                labelsize=fnts + 3,
                pad=10,
            )
            ax1[i].set_xlim(profile_x_lims[0], profile_x_lims[1])
            ax1[i].set_ylim(profile_y_lims[i][0], profile_y_lims[i][1])
            ax1[i].set_ylabel(y_labels[i], fontsize=fnts + 3)
            ax1[i].set_xlabel(r"$R/R_\oplus$", fontsize=fnts + 3)

    except ValueError:
        pass

    for i in range(N_match):
        pl = planets[indices[i]]
        pl.Update()
        pl.trim_profiles()

        # Plot Temperature, Pressure and density profile
        for j in range(3):
            c = (all_out.T[n][i] - c1) / (c2 - c1)
            ax1[j].plot(
                pl.profiles[0] / r_earth,
                pl.profiles[j + 1] * profile_scales[j],
                color=cmap(c),  # (c[i], c[i]*.125, 0.), #color_list[j*3],
                linewidth=0.5,
            )

        # plot PT-profiles
        if plot_PT:
            col = ftool.extract_col_from_cmap(
                cmap, PT_value_range_cmap, 10**pl.ocean_frac_is
            )
            ax2.semilogy(pl.profiles[1], pl.profiles[2], color=col)
            ax2.scatter(
                pl.layer_properties[3]["T_outer"],
                pl.layer_properties[3]["P_outer"],
                color="r",
            )

            inset_ax = inset_axes(ax2, height="40%", width="5%", loc=4, borderpad=10)
            plt.subplots_adjust(wspace=0.3)
            cbar = plt.colorbar(sm, cax=inset_ax, ticks=[0.0, 0.5, 1.0])

            n = 47

            ticklabels = np.array(
                [
                    all_out.T[n].min(),
                    (all_out.T[n].min() + all_out.T[n].max()) / 2.0,
                    all_out.T[n].max(),
                ]
            ).round(3)

            cbar.ax.set_yticklabels(ticklabels)
            cbar.ax.tick_params(labelsize=fnts)
            cbar.set_label(r"$\rm Ocean \ mass \ fraction$", fontsize=fnts)

    if plot_CMB:
        # Indicated ranges for the core-mantle-boundary
        box_col = (0.75, 0.75, 0.75)
        for j in range(3):
            minsx = all_out.T[15].min() / r_earth * 1000.0
            maxsx = all_out.T[15].max() / r_earth * 1000.0

            minsy = all_out.T[28 + j].min() * profile_scales[j]
            maxsy = all_out.T[28 + j].max() * profile_scales[j]

            trans1 = transforms.blended_transform_factory(
                ax1[j].transData, ax1[j].transAxes
            )

            trans2 = transforms.blended_transform_factory(
                ax1[j].transAxes, ax1[j].transData
            )

            rect1 = patches.Rectangle(
                (minsx, 0),
                width=maxsx - minsx,
                height=1,
                transform=trans1,
                color=box_col,
                alpha=1.0,
                linestyle=None,
            )

            rect2 = patches.Rectangle(
                (0, minsy),
                width=1,
                height=maxsy - minsy,
                transform=trans2,
                color=box_col,
                alpha=1.0,
                linestyle=None,
            )

            ax1[j].add_patch(rect1)
            ax1[j].add_patch(rect2)

    # Plot PREM density profile
    if plot_PREM:
        ax1[2].plot(
            PREM_data[0] * 1000 / r_earth,
            PREM_data[1] / 1000,
            color=(0.5, 0.5, 0.5),
            linestyle="--",
            linewidth=3,
        )

        ax1[2].text(0.6, 10, "PREM", color=(0.5, 0.5, 0.5), fontsize=fnts + 4)

    bbox_props = dict(edgecolor="k", facecolor="white", linewidth=0.5)

    # ax1[2].text(1.1, 1., plot_label,
    #           fontsize=24,
    #          transform=ax1[-1].transAxes,
    #         bbox=bbox_props,
    #        horizontalalignment = 'right')

    # fig1.savefig(spec+'/'+obj+'_'+spec+'.pdf', format='pdf',
    #            bbox_inches='tight')

    fig1.savefig(spec + "/" + obj + ".png", format="png", bbox_inches="tight", dpi=240)

    if plot_PT:
        x, y = 0.6, 0.1
        trafo = transforms.blended_transform_factory(ax2.transAxes, ax2.transAxes)
        ax2.scatter([x], [y], color="r", transform=trafo, s=75)
        ax2.text(
            x + 0.02,
            y,
            r"$\rm Ocean-Manlte \ interphase$",
            va="center",
            transform=trafo,
            fontsize=fnts,
        )
        fig2.savefig(
            spec + "/" + obj + "_" + spec + "_PT.png",
            format="png",
            bbox_inches="tight",
            dpi=240,
        )
    fig.savefig(
        spec + "/" + obj + "_calibration.png",
        format="png",
        bbox_inches="tight",
        dpi=240,
    )

    fig3.savefig(
        spec + "/" + obj + "_data_selection.png",
        format="png",
        bbox_inches="tight",
        dpi=240,
    )


def check_statistics(
    data, N=10, red=2, obj="earth", write=False, spec="all", N_all=None
):

    N_sub = int(len(data) / red)

    all_ranges = np.empty([N, len(data[0]), 2])
    RMSD = np.empty([len(data[0]), 2])

    # Loop over number of random samples of given size
    for i in range(N):
        sub_data = np.empty([N_sub, len(data[0])])
        indices = np.arange(0, len(data))
        ind = create_subset(indices, red=red)

        for j in range(N_sub):
            sub_data[j] = data[ind[j]]

        # Loop over individual parameters to be extracted
        for j in range(len(data[0])):
            all_ranges[i][j][0] = np.nanmin(sub_data.T[j])
            all_ranges[i][j][1] = np.nanmax(sub_data.T[j])

    columns = ["{}_min".format(p) for p in params]
    for i in range(len(params)):
        columns.insert(2 * i + 1, "{}_max".format(params[i]))

    df = pd.DataFrame([a.flatten() for a in all_ranges], columns=columns)

    # Apporach I
    std_low = np.empty([len(data[0])])
    std_high = np.empty([len(data[0])])

    mean_low = np.empty([len(data[0])])
    mean_high = np.empty([len(data[0])])

    for j in range(len(data[0])):
        mean_low[j] = np.mean(all_ranges.T[0][j])
        std_low[j] = ftool.std_corr(all_ranges.T[0][j])
        mean_high[j] = np.mean(all_ranges.T[1][j])
        std_high[j] = ftool.std_corr(all_ranges.T[1][j])

        # print (params[j]+':', round(std_low[j]/mean_low[j]*100, 3),'%',
        #      round(std_high[j]/mean_high[j]*100, 3),'%')

    digits = 4

    if write:
        print("writing to file...")

        # Write ranges to file
        path = spec + "/" + obj + "_ranges.txt"

        print("path =", path)
        df.to_pickle("{a}/{b}_ranges.pkl".format(a=spec, b=obj))

        with open(path, "w") as f:
            f.write(
                str(len(data))
                + " of "
                + str(N_all)
                + " planets ("
                + str(int(len(data) / N_all * 100))
                + "%)\n"
            )
            for i in range(len(params)):
                print("param = ", params[i], np.nanmin(data.T[i]), np.nanmax(data.T[i]))
                # Hydro structure
                if i == len(params) - 1:
                    string = params[i] + ":"
                    string = string + "." * (14 - len(params[i]))
                    string = string + ", ".join(set([str(int(d)) for d in data.T[i]]))

                    # sys.exit()
                else:
                    m1 = ftool.fancyround(abs(np.nanmin(data.T[i])), digits=digits)
                    m2 = ftool.fancyround(abs(np.nanmax(data.T[i])), digits=digits)
                    d1 = abs(std_low[i] / mean_low[i] * np.nanmin(data.T[i]))
                    d1 = ftool.fancyround(abs(d1), 1)

                    d2 = abs(std_high[i] / mean_high[i] * np.nanmax(data.T[i]))
                    d2 = ftool.fancyround(abs(d2), 1)

                    # If std is less than 0.1% artificially set it to zero.
                    if abs(d1 - m1) / m1 < 1e-3:
                        d1 = 0.0

                    if abs(d2 - m2) / m2 < 1e-3:
                        d2 = 0.0

                    s1 = np.nanmin(data.T[i]) / abs(np.nanmin(data.T[i]))
                    s2 = np.nanmax(data.T[i]) / abs(np.nanmax(data.T[i]))
                    if np.isnan(s1):
                        s1 = 1.0
                    if np.isnan(s2):
                        s2 = 1.0

                    string = params[i] + ":"
                    string = string + "." * (14 - len(params[i]))
                    string = string + "{:.3e}".format(Decimal(str(m1 * s1)))
                    string = string + "("
                    string = string + "{:.1e}".format(Decimal(str(d1)))
                    string = string + ") to "
                    string = string + "{:.3e}".format(Decimal(str(m2 * s2)))
                    string = string + "("
                    string = string + "{:.1e}".format(Decimal(str(d2)))
                    string = string + ")\n"
                f.write(string)

    return df  # all_ranges, std_low/mean_low, std_high/mean_high, RMSD


def create_subset(indices, red=2):
    """
    Take array of indices and randomly select a subset of int(N/red) indices

    """
    N = len(indices)

    N_sub = int(N / red)

    ind = []

    if red > 1:
        for i in range(N_sub):

            while True:
                p = random.random()
                ii = int(p * (N - 1))
                if ii in ind:
                    pass

                else:
                    ind.append(ii)
                    break

        return ind

    else:
        return indices


def estimate_maximal_uncertainty(N=10, obj="earth", spec="all"):

    delta0 = 0.05

    pairs = np.empty([N, 2])
    all_ranges = np.empty([N, len(params), 2])

    for i in range(N):
        p1 = random.random()
        p2 = random.random()

        C = MoI_solar[obj] * (1.0 - delta0 + p1 * 2 * delta0)
        deltaC = 0.02 + p2 * 0.03

        pairs[i][0] = C
        pairs[i][1] = deltaC

        all_all_out, planets, indices = calibrate(
            obj=obj, plot=False, reldev_C=deltaC, MoI_real=C
        )

        data = np.empty([len(indices), len(params)])
        for j in range(len(indices)):
            data[:] = all_all_out[j][:]

        print("len data =", len(data))

        for j in range(len(data[0])):
            all_ranges[i][j][0] = np.nanmin(data.T[j])
            all_ranges[i][j][1] = np.nanmax(data.T[j])

    for i in range(N):
        for j in range(len(all_ranges[0])):
            print(params[j], ":", all_ranges[i][j][0], all_ranges[i][j][1])

    N_col = 3
    N_row = 3
    fig, ax = plt.subplots(N_row, N_col, figsize=(10, 8))
    plt.subplots_adjust(wspace=0.75, hspace=0.5)

    all_diffs = all_ranges.T[1] - all_ranges.T[0]
    all_diffs = all_diffs.T

    colorbar_labels = [
        [r"$\delta \rm Mg\#$", r"$\delta \xi_{\rm FeO}$", r"$\delta X_{\rm Si}$"],
        [
            r"$\delta R_{\rm Core} \ [\rm km]$",
            r"$\delta R_{\rm ICB}  \ [\rm km]$",
            r"$\delta X_{\rm O}$",
        ],
        [r"$\delta T_{\rm C}$", r"$\delta P_{\rm C}$", r"$\delta \log(f_{\rm O_2})$"],
    ]

    which_params = [[0, 5, 9], [15, 16, 10], [20, 21, 22]]

    for i in range(N_row):
        for j in range(N_col):
            ax[i][j].set_xlabel(r"$C/MR^2$")
            ax[i][j].set_ylabel(r"$\delta C/C$")

    # Plot the parameter ranges for all [C, delta C] pairs
    for i in range(N):
        for row in range(N_row):
            for col in range(N_col):
                p = which_params[row][col]
                sc = ax[row][col].scatter(
                    pairs[i][0],
                    pairs[i][1],
                    c=all_diffs[i][p],
                    cmap="copper",
                    vmin=np.nanmin(all_diffs.T[p]),
                    vmax=np.nanmax(all_diffs.T[p]),
                )

                ftool.my_colorbar(sc, pad=0.25, label=colorbar_labels[row][col])

    print("pairs = ", pairs)

    fig.savefig(spec + "/" + obj + "_stat.pdf", format="pdf", bbox_inches="tight")
    fig.savefig(
        spec + "/" + obj + "_stat.png", format="png", bbox_inches="tight", dpi=240
    )


def plot_uncertainties(all_all_out, indices, obj="earth", res=6, N=10, spec="all"):
    """
    Computes variability in the inferred parameter ranges for res subsamples
    and plots them as a function of the relative sample size.
    """
    N_match = len(indices)
    max_ranges = np.empty([len(params), 2])
    # gather plot data
    data = np.empty([N_match, len(params)])

    i = -1
    while i < N_match - 1:
        i += 1
        data[i][:] = all_all_out[indices[i]][:]

    # Compute min and max for each parameter of the entire data set
    for i in range(len(params)):
        max_ranges[i][0] = np.nanmin(data.T[i])
        max_ranges[i][1] = np.nanmax(data.T[i])

    if res < 2:
        print("WARNING: res < 2 given in plot_uncertainty.")

    labels = [
        [
            r"$\sigma_j(\rm Mg\#) \ [\%]$",
            r"$\sigma_j(\rm X_{\rm FeO}^{\rm Mantle}) \ [\%]$",
        ],
        [r"$\sigma_j(R_{\rm Core}) \ [\%]$", r"$\sigma_j(R_{\rm ICB}) \ [\%]$"],
        [
            r"$\sigma_j(X_{\rm Si}^{\rm Core}) \ [\%]$",
            r"$\sigma_j(X_{\rm O}^{\rm Core}) \ [\%]$",
        ],
        [r"$\sigma_j(T_{\rm C}) \ [\%]$", r"$\sigma_j(P_{\rm C}) \ [\%]$"],
    ]
    which_params = [[0, 5], [15, 16], [9, 10], [19, 20]]

    all_std = np.empty([res, 2, len(params)])
    all_rmsd = np.empty([res, 2, len(params)])

    # Loop over number of sample sizes
    for i in range(res):

        print("Starting sweep ", i)

        # Compute fraction of data used for sample
        red = 2 ** (res - i)

        ranges, low, high, rmsd = check_statistics(
            data, N=N, red=red, write=True, obj=obj, spec=spec, N_all=len(all_all_out)
        )

        for j in range(len(low)):
            all_std[i][0][j] = low[j]
            all_std[i][1][j] = high[j]
            all_rmsd[i][0][j] = rmsd[j][0]
            all_rmsd[i][1][j] = rmsd[j][1]

    fig, ax = plt.subplots(len(which_params), len(which_params[0]), sharex=True)
    plt.subplots_adjust(wspace=0.4)
    fig.align_ylabels(ax[:])

    ax[-1][0].set_xlabel(r"$j$")
    ax[-1][1].set_xlabel(r"$j$")

    bound_cols = ["k", "g"]
    bound_line_styles = ["-", "--"]
    # x_range = np.arange(res)
    # x_range = x_range[res:0:-1]
    for i in range(len(which_params)):
        for j in range(len(which_params[i])):
            ax[i][j].set_ylabel(labels[i][j])
            ax[i][j].plot(
                np.arange(res),
                100 * all_std.T[which_params[i][j]][0],
                marker="x",
                label="lower bound",
                color="k",
                linestyle=bound_line_styles[0],
            )
            ax[i][j].plot(
                np.arange(res),
                100 * all_std.T[which_params[i][j]][1],
                marker="x",
                label="upper bound",
                color="k",
                linestyle=bound_line_styles[1],
            )

            ax[i][j].xaxis.set_ticks(np.linspace(0, res - 1, res))
            ax[i][j].set_xticklabels(np.arange(res)[::-1])

    ax[0][0].legend()

    fig.savefig(
        spec + "/" + obj + "_conv_" + spec + ".pdf", format="pdf", bbox_inches="tight"
    )

    return max_ranges


def fit_sample(all_all_out):
    """Perform simple MR fit of a sample of planets. By default the fittin
    function is of the form log(R/R_E) = a*log(M/M_E) + b*log(T_surf/300 K)
    """
    N = len(params)

    masses = all_all_out.T[45]
    temps = all_all_out.T[17]
    radii = all_all_out.T[1]
    c1 = 0.0
    c2 = 2000.0

    fig, ax = plt.subplots()
    sc = ax.scatter(
        np.log10(masses), np.log10(radii), c=temps, cmap="jet", vmin=c1, vmax=c2
    )
    cbar = plt.colorbar(sc)
    cmap = plt.cm.get_cmap("jet")

    def fct(theta, M, T):
        a, b, c, d, e = theta
        return a + b * np.log10(M) + c * (T / 1000) + d * np.log10(M) * (T / 1000)

    def diff(theta):
        return fct(theta, masses, temps) - np.log10(radii)

    # Create pandas data frame of data
    df = pd.DataFrame(all_all_out, columns=all_param_labels)

    theta0 = [0.5, 0.0, 0.0, 0.0, 0.0]

    res = least_squares(diff, theta0)

    x = np.linspace(0.1, 10)
    t = np.linspace(1.0, 2000, 10)
    dat = np.empty([len(t), len(x)])
    for i in range(len(t)):
        for j in range(len(x)):
            dat[i][j] = fct(res.x, x[j], t[i])

    for i in range(len(t)):
        c = (t[i] - c1) / (c2 - c1)
        ax.plot(np.log10(x), dat[i], color=cmap(c))

    cbar.set_label(r"$T_{\rm eq}$")
    ax.set_xlabel(r"$\log{M/M_\oplus}$")
    ax.set_ylabel(r"$\log{R/R_\oplus}$")
    fig.savefig("./fit_test.png", format="png", dpi=240)

    return res


def fit_grid(all_all_out, obj="earth", spec="all"):

    import seaborn as sb

    N = len(params)

    df = pd.DataFrame(all_all_out, columns=all_param_labels)
    corrMatrix = df.corr()
    s = corrMatrix.unstack()

    for i in [21, 20, 27]:
        ind = np.argpartition(np.abs(corrMatrix.to_numpy()[i]), -10)[-10:]
        print(params[i] + ": params =", [params[ii] for ii in ind])

    # sb_plot = sb.heatmap(corrMatrix, annot = True)

    l = len(all_param_labels)
    ticks = np.linspace(0.5 / l, 1.0 - 0.5 / l, l)

    fig, ax = plt.subplots(figsize=(12, 12))
    ax.set_xticks(ticks)
    ax.set_xticklabels(all_param_labels, rotation=90)

    ax.set_yticks(ticks)
    ax.set_yticklabels(all_param_labels)

    im = ax.imshow(
        df.corr(), extent=[0.0, 1.0, 0.0, 1.0], aspect="equal", origin="lower"
    )
    ftool.my_colorbar(im, pad=0.25, label="correlation")

    fig.savefig(obj + "_corr_" + spec + ".png", format="png", bbox_inches="tight")

    thetas = [[2.0, 0.0, 0.0, 0.0], [3000.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0]]

    logs = [[True, True], [True, False], [True, False]]

    xyz = [[-1, 21, 0], [-1, 20, 0], [-1, 27, 0]]

    fig, ax = plt.subplots(1, len(xyz), figsize=(12, 3))
    plt.subplots_adjust(wspace=0.5)

    for i in range(len(xyz)):
        x_ = np.linspace(-1, 0.5)
        z_ = np.linspace(0.25, 0.75, 6)

        x_data = all_all_out.T[xyz[i][0]]
        y_data = all_all_out.T[xyz[i][1]]
        z_data = all_all_out.T[xyz[i][2]]

        if logs[i][0]:
            x_data = np.log10(x_data)

        if logs[i][1]:
            y_data = np.log10(y_data)

        def yy(theta, x, z):
            return theta[0] + theta[1] * x + theta[2] * z + theta[3] * x * z

        def fct(theta):
            return yy(theta, x_data, z_data) - y_data

        def driver(x, x_data, y_data):
            ynew = fct(x_data, *x)
            yerr = np.sum((ynew - y_data) ** 2)
            return yerr

        theta0 = thetas[i]
        res = least_squares(fct, theta0)
        print(res.x)

        x, y, z = xyz[i]

        ax[i].set_xlabel(all_param_labels[x])
        ax[i].set_ylabel(all_param_labels[y])

        sc = ax[i].scatter(x_data, y_data, c=z_data, cmap="jet")

        ftool.my_colorbar(sc, pad=0.25, label=all_param_labels[z])

        n = len(z_)
        colors = pylab.cm.jet(np.linspace(0, 1, n))

        for j in range(len(z_)):
            ax[i].plot(x_, yy(res.x, x_, z_[j]), color=colors[j])

    fig.savefig(obj + "_fit_" + spec + ".png", format="png", bbox_inches="tight")
