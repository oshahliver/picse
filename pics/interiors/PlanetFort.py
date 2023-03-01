# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 14:05:05 2018

@author: os18o068
"""

from pics.utils.calibration_tools import regression as regr
import matplotlib.ticker
from matplotlib import pyplot as plt
from matplotlib.ticker import (
    MultipleLocator,
    FormatStrFormatter,
    AutoMinorLocator,
    LogLocator,
    FixedLocator,
)

from pics.materials import phase_transitions_water_Wagner2002 as phase
from pics.physicalparams import (
    T0_list,
    K0_list,
    K0prime_list,
    rho0_list,
    aP_list,
    molar_mass_list,
    Rgas,
    mH,
    mO,
    kB,
    G,
    EOS_type,
    material_list,
    EOSstring_list,
    NA,
    P_earth,
    m_earth,
    r_earth,
    material_plot_list,
    solar_system,
    abbrevations_solar,
    P_zero,
    T_zero,
    mFe,
    mMg,
    mSi,
    mH2O,
    Mg_number_solar,
    temperature_jumps,
    mS,
    material_list_fort,
    material_YMg,
    material_YSi,
    material_YO,
    material_YH,
    material_YS,
    material_YFe,
    mS,
    sigmaSB,
)

from pics.runparams import (
    eps_Psurf,
    eps_Mtot,
    eps_layer,
    param_colors,
    plot_units,
    suffix,
    eps_Mg_number,
    plot_params,
    eps_Tsurf,
    color_list,
    layerCodes,
)

import numpy as np
import time
from pics.utils.function_tools import functionTools as ftool
import sys
from tabulate import tabulate
from decimal import Decimal
import astropy.table
from astropy.io import ascii
import warnings
from pics.utils.plot_tools.plot_tools import plotTools
from pics.utils import fortplanet
from pics.utils import fortfunctions
from pics.materials import material
from pics.utils import readPREM

fortplanet.wrapper.load_eos_tables(table_dir="{}/data/EoS_tables/".format("."))
warnings.filterwarnings("ignore")


def refine_hydro_structure(lower, upper, res=1):
    R, T, P, d, M = [[l, u] for l, u in zip(lower, upper)]
    dR = R[1] - R[0]
    dT = T[1] - T[0]
    dP = P[1] - P[0]
    dM = M[1] - M[0]

    slope_T = dT / dR
    slope_P = dP / dR
    slope_M = dM / dR

    N = 2**res + 1
    temps = np.linspace(T[0], T[1], N)
    pres = np.linspace(P[0], P[1], N)
    rads = np.linspace(R[0], R[1], N)
    phases = [phase.phase(P=P[0], T=T[0]), phase.phase(P=P[1], T=T[1])]

    for t, p, r in zip(temps, pres, rads):
        ph = phase.phase(P=p, T=t)
        if ph == phases[0]:
            pass
        else:
            m = M[0] + slope_M * (r - R[0])
            return r, t, p, m


def convert_X_ice_to_xi_ice(Si_number=None, xi_Fe=None, X_ice=None, contents=[6, 7, 1]):

    SiMg = Si_number / (1.0 - Si_number)

    # Compute fractions without water
    fractions = fortfunctions.functionspy.compute_abundance_vector(
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

    print(fractions)

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

    # Compute mole fraction of water at given composition
    xi = material.xi(eta=X_ice, m1=m_tilde, m2=mH2O)

    return xi


def convert_X_impurity_to_xi_impurity(
    Si_number=None, xi_Fe=None, X=None, contents=[6, 7, 1], m=(2 * mH + mO)
):

    SiMg = Si_number / (1.0 - Si_number)

    # Compute fractions without water
    fractions = fortfunctions.functionspy.compute_abundance_vector(
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


def compute_core_mass(
    contents=None,
    Mg_number=None,
    M_surface=1.0,
    Mg_number_mantle=None,
    SiMg=None,
    M_ocean=0.0,
    xi_H_core=0.0,
    impurity=[],
    xi_S_core=0.0,
    n=3,
    xi_all_core=[],
    M_IC=0.0,
):
    """Computes the core mass of a planet at given total mass, composition and
    value for Mg#
    """
    """
    print ('\nCompute Core Mass:')
    print ('contents = ', contents)
    print ('Mg_number =', Mg_number)
    print ('M_surface =', M_surface)
    print ('Mg_number_mantle = ', Mg_number_mantle)
    print ('SiMg =', SiMg)
    print ('M_ocean =', M_ocean)
    print ('xi_H_core =', xi_H_core)
    print ('impurity =', impurity)
    print ('xi_S_core =', xi_S_core)
    """
    # print (contents, Mg_number, M_surface, Mg_number_mantle, SiMg, M_ocean, xi_H_core)
    Mg_number_mantle = min(Mg_number_mantle, 0.9999999999)
    FeMg_mantle = (1.0 - Mg_number_mantle) / Mg_number_mantle
    FeMg = (1.0 - Mg_number) / Mg_number

    # Compute the fractions in the mantle
    fractions = fortfunctions.functionspy.compute_abundance_vector(
        simg=SiMg,
        femg=FeMg_mantle,
        n_mats=len(contents[n]),
        ymgi=[material_YMg[i - 1] for i in contents[n]],
        ysii=[material_YSi[i - 1] for i in contents[n]],
        xih2oi=[0.0 for i in contents[n]],
        xifei=[1.0 - Mg_number_mantle for i in contents[n]],
        xialsii=[0.0 for i in contents[n]],
        xialmgi=[0.0 for i in contents[n]],
        contents=contents[n],
        additional=impurity,
    )
    # Count
    # Note that the indices are shifted by one because the original definition
    # of the arrays comes from the fortran code.
    Q1 = sum(
        [
            fractions[i] * Mg_number_mantle * material_YMg[contents[n][i] - 1]
            for i in range(len(contents[n]))
        ]
    )

    # Compute total normalized mass in the mantle
    Q2 = (
        sum(
            [
                fractions[i] * Mg_number_mantle * material_YMg[contents[n][i] - 1]
                for i in range(len(contents[n]))
            ]
        )
        * mMg
        + sum(
            [
                fractions[i]
                * (1.0 - Mg_number_mantle)
                * material_YMg[contents[n][i] - 1]
                for i in range(len(contents[n]))
            ]
        )
        * mFe
        + sum(
            [
                fractions[i] * material_YSi[contents[n][i] - 1]
                for i in range(len(contents[n]))
            ]
        )
        * mSi
        + sum(
            [
                fractions[i] * material_YO[contents[n][i] - 1]
                for i in range(len(contents[n]))
            ]
        )
        * mO
        + sum(
            [
                fractions[i] * material_YH[contents[n][i] - 1]
                for i in range(len(contents[n]))
            ]
        )
        * mH
    )

    Q3 = sum(
        [
            fractions[i] * (1.0 - Mg_number_mantle) * material_YMg[contents[n][i] - 1]
            for i in range(len(contents[n]))
        ]
    )

    # Q4 = 1. - xi_H_core
    # Q5 = (mFe + xi_S_core*mS)*(1.-xi_H_core) + xi_H_core*mH
    # core_frac = (1.-M_ocean/M_surface)*(Q1/Q2-Mg_number*(Q1/Q2+Q3/Q2))/\
    #             (Mg_number*(Q4/Q5-Q1/Q2-Q3/Q2)+Q1/Q2)

    Q4 = xi_all_core[0]
    m_core = [mFe, mH, mS, mSi, mO]
    Q5 = sum([m_core[i] * xi_all_core[i] for i in range(len(xi_all_core))])
    # print ('Q =', Q1, Q2, Q3, Q4, Q5)
    core_frac = 1.0 - M_ocean / M_surface
    core_frac *= Q1 / Q2 - Mg_number * (Q1 / Q2 + Q3 / Q2)
    core_frac /= Mg_number * (Q4 / Q5 - Q1 / Q2 - Q3 / Q2) + Q1 / Q2
    # print ('core mass old =', core_frac * M_surface)
    core_frac = 1.0 - M_ocean / M_surface
    core_frac *= Q3 / Q2 - Q1 / Q2 * FeMg
    core_frac += M_IC / M_surface * (1.0 / mFe - Q4 / Q5)
    core_frac /= Q3 / Q2 - Q4 / Q5 - FeMg * Q1 / Q2
    # print ('')
    return core_frac * M_surface


class Planet:
    def __init__(self, run_params, planetary_params, **kwargs):
        pass

    def restore_params(self):
        pass


class Planet:
    def __init__(
        self,
        P_center=1.0e12,
        T_center=3000,
        R_seed=0.123,
        contents=[[2]],
        tempType=1,
        layerType=1,
        fractions=[[1.0]],
        P_surface_should=P_earth,
        majorConstraint="P_surface",
        layerradii=[],
        layerConstraint=[1],
        add_atmosphere=False,
        eps_r=0.432,
        layermasses=[],
        layerpres=[],
        M_surface_should=1.0,
        T_surface_should=300,
        R_surface_should=1.0,
        Si_number_layers=[],
        Fe_number_layers=[],
        rhoType=1,
        echo=False,
        vapor_stop=False,
        Mg_number_should=Mg_number_solar,
        Si_number_should=0.0,
        Mg_number_layer=False,
        subphase_res=32,
        modelType=0,
        temp_jumps=[0.0, 0.0, 0.0, 0.0, 0.0],
        gammas_layer=[1.36, 1.36, 1.96, 1.26, 1.0],  # 2.43
        layertemps=[0, 0, 0, 0, 0],
        T_zero=[-1000, -1000, -1000.0, -1000.0, -1000.0],
        M_core_should=None,
        ocean_frac_should=-10.0,
        adiabatType=1,
        q=[0.91, 0.91, 2.5, 2.9, 1.0],  # 0.489
        silence=False,
        measureTime=False,
        eps_T_zero=0.0,
        omega=0.0,
        # xi_Stv=0.,
        xi_S_core=0.0,
        P_CS=50e9,
        core_segregation_type=0,
        inner_core_segregation_type=0,
        E_tot_should=1.0,
        L_int_should=1.0,
        L_eff_should=1.0,
        Si_number_mantle=0.4,
        Mg_number_mantle=1.0,
        #  run_params,
        #  planetary_params,
        **kwargs,
    ):

        self.measureTime = measureTime
        self.core_segregation_type = core_segregation_type
        self.inner_core_segregation_type = inner_core_segregation_type
        # print ('initializing planet ...')
        self.status = "Shadow of the Future"

        # if set to True construction is stopped if vapor phase is reached
        self.vapor_stop = vapor_stop
        self.M_tot_should = None
        self.R_tot_should = None
        self.M_tot_is = None
        self.R_tot_is = None
        self.R_seed = R_seed
        self.rho_center = None
        self.P_center = P_center
        self.T_center = T_center
        self.P_surface_should = P_surface_should
        self.M_surface_should = M_surface_should * m_earth
        self.P_space_should = P_surface_should
        self.P_space_is = None
        self.T_surface_should = T_surface_should
        self.R_surface_should = R_surface_should * r_earth
        self.contents = contents
        self.temp_jumps = temp_jumps
        # self.materials = [[material_list[i] for i in con] for con in contents]
        self.fractions = fractions
        self.v_esc_surface = 0.0
        self.gravity = 0.0
        self.Mg_number_is = 0.0
        self.Mg_number_should = Mg_number_should
        self.Si_number_is = 0.0
        self.Si_number_should = Si_number_should
        self.Si_number_layers = Si_number_layers
        self.Fe_number_layers = Fe_number_layers
        self.echo = echo
        self.M_H2O_should = 0.0
        self.M_H2O_is = 0.0
        self.M_H2O_hidden_is = 0.0
        self.delta_M_H2O = 0.0
        self.delta_M_H2O_hidden = 0.0
        self.H2O_count = 0.0
        self.Fe_count = 0.0
        self.Si_count = 0.0
        self.Al_count = 0.0
        self.H_count = 0.0
        self.O_count = 0.0
        self.Mg_count = 0.0
        self.S_count = 0.0
        self.density_inversion = False
        self.match = True
        self.converged = False
        self.vapor_reached = False
        self.N_shells = 0
        self.N_shells_layers = np.empty([len(contents)])
        self.test_variable = 0.0
        self.delta_Mg_count = 0.0
        self.Mg_number_layer = Mg_number_layer
        # self.subphase_refined = False
        self.subphase_res = subphase_res
        # self.subphase_refine_inhibit = False
        self.gammas_layer = gammas_layer
        self.adiabatType = adiabatType
        self.eps_T_zero = eps_T_zero
        # self.xi_Stv = xi_Stv
        self.simple_densities = []
        self.simple_MoI = None
        self.M_outer_mantle = 0.0
        self.M_ocean_is = 0.0
        self.M_DWR_is = 0.0
        self.M_ocean_should = M_surface_should * 10**ocean_frac_should
        self.M_core_should = M_core_should
        self.M_core_is = 0.0
        self.ocean_frac_should = ocean_frac_should
        self.ocean_frac_is = -10
        self.T_zero = T_zero
        self.omega = omega
        self.ocean_depth = 0.0
        self.rho_mean = 0.0
        self.T_CS = 0.0
        self.xi_Fe_mantle = None
        self.Si_number_mantle = Si_number_mantle
        self.Mg_number_mantle = Mg_number_mantle
        self.P_CS = P_CS
        self.logfO2 = None
        self.xi_S_mantle = 0e0
        self.E_grav = 0e0
        self.E_int = 0e0
        self.E_tot_is = 0e0
        self.E_tot_should = E_tot_should * G * m_earth**2 / r_earth * 3.0 / 5.0
        self.L_int_should = (
            L_int_should * 4.0 * np.pi * r_earth**2 * sigmaSB * 300**4
        )
        self.L_int_is = 0e0
        self.L_eff_is = 0e0
        self.L_eff_should = (
            L_eff_should * 4.0 * np.pi * r_earth**2 * sigmaSB * 300**4
        )

        # defines the physical model that is invoked for the planet
        # e.g Brucite dissociation model without silicates: modelType = 0
        # " " " with hydrated Olivine: modelType = 1
        self.modelType = modelType
        self.mantle_exists = False
        self.inner_core_exists = False
        self.outer_core_exists = False

        self.layermasses = layermasses
        self.layerradii = layerradii
        self.layertemps = layertemps
        self.layerpres = layerpres
        self.d0 = []
        self.q = q
        self.xi_H_core = 0.0
        self.xi_FeO_mantle = 0.0
        self.xi_S_core = xi_S_core
        self.P_H2_CMB = 0.0
        self.M_H2O_core = 0.0
        self.M_H2O_mantle = 0.0
        self.M_ocean = 0.0
        self.M_mantle_is = 0.0

        try:
            self.xi_all_core = material.mat2at_core(
                xi=self.fractions[1], xiH=self.xi_H_core
            )

            self.X_all_core = material.at2wt_core(self.xi_all_core)

        except IndexError:
            self.xi_all_core = None
            self.X_all_core = None

        # It is important that the inner core fraction is accounted for properly
        # If a BC are to met via toolkit.iterate(), the inner core fraction
        # must be updated each time the core mass is updated otherwise the
        # iteration will in many cases crash and no solution can be found
        # print ('inner core frac in Planet.__init__():', inner_core_frac)
        """
        if not inner_core_frac == None:
            print ('check1')
            self.inner_core_frac = inner_core_frac
            self.layermasses[0] = self.layermasses[1] * inner_core_frac
            #sys.exit()
            #self.layermasses[1] = self.M_core_should*(1.-self.inner_core_frac)
        
        else:
            print ('check2')
            try:
                M_core =  self.layermasses[0] + self.layermasses[1]
                self.inner_core_frac = self.layermasses[0] / M_core * 0.
                M_core = 0.
            except IndexError:
                self.inner_core_frac = inner_core_frac * 0.
       
        """
        M_core = 0.0
        self.inner_core_frac = 0.0
        self.eps_r = eps_r

        if len(self.contents) > 0:
            self.differentiated = True

        else:
            self.differentiated = False

        self.shellIterationCount = 0
        self.layerIterationCount = 0
        self.changebisec = True
        self.shellIteration = True
        self.layerIteration = True

        self.tempType = tempType
        self.layerType = layerType
        self.rhoType = rhoType
        self.majorConstraint = majorConstraint
        self.layerConstraint = layerConstraint
        self.layers = []
        self.lay = 0
        self.test = False
        self.atmosphere = None
        self.add_atmosphere = add_atmosphere
        self.have_atmosphere = False
        self.majorComplete = False
        self.layerComplete = False
        self.layer_properties = [
            {
                "P_outer": None,
                "T_outer": None,
                "rho_outer": None,
                "rho_inner": None,
                "indigenous_mass": 0.0,
                "R_outer": None,
                "n_shells": 0,
            }
            for i in range(len(self.contents))
        ]

        self.silence = silence
        # if no fractions have been specified, distribute the components evenly
        # over the mixture
        if len(fractions) == 0:
            if not self.silence:
                print("setting up uniform fractions")
            for c in range(len(self.contents)):
                frac = 1.0 / len(self.contents[c])
                self.fractions.append([frac for i in self.contents[c]])

        """
        if len(self.layermasses) != len(self.contents):
            if not self.silence:
                print ('\nWARNING: number of layers and layermasses does not match')
                print ('given layermasses:', len(self.layermasses))
                print ('given layers: ', len(self.contents))
        """

        self.status = "seed"

        self.M_surface_is = 0.0
        self.T_surface_is = T_center
        self.P_surface_is = P_center
        self.R_surface_is = R_seed
        self.exit_code = 0
        self.direction = None
        self.constraintValue_is = None
        self.constraintValue_should = None

        # gather all initial properties in a dictionary
        self.initials = {
            "P_center": self.P_center,
            "T_center": self.T_center,
            "R_seed": self.R_seed,
            "contents": self.contents,
            "tempType": self.tempType,
            "layerType": self.layerType,
            "fractions": self.fractions,
            "P_surface_should": self.P_surface_should,
            "majorConstraint": self.majorConstraint,
            "layerConstraint": self.layerConstraint,
            "add_atmosphere": self.add_atmosphere,
            "eps_r": self.eps_r,
            "layermasses": self.layermasses,
            "layerradii": self.layerradii,
            "differentiated": self.differentiated,
            "M_surface_should": self.M_surface_should / m_earth,
            "T_surface_should": self.T_surface_should,
            "R_surface_should": self.R_surface_should / r_earth,
            "rhoType": self.rhoType,
            "echo": self.echo,
            "Mg_number_should": self.Mg_number_should,
            "Si_number_should": self.Si_number_should,
            "Si_number_layers": self.Si_number_layers,
            "Fe_number_layers": self.Fe_number_layers,
            "vapor_stop": self.vapor_stop,
            "Mg_number_layer": self.Mg_number_layer,
            "modelType": self.modelType,
            "temp_jumps": self.temp_jumps,
            "gammas_layer": self.gammas_layer,
            "M_ocean_should": self.M_ocean_should,
            "ocean_frac_should": self.ocean_frac_should,
            "layertemps": self.layertemps,
            "T_zero": self.T_zero,
            "M_core_should": self.M_core_should,
            "subphase_res": self.subphase_res,
            "layerpres": self.layerpres,
            #'xi_Stv':self.xi_Stv,
            "omega": self.omega,
            "xi_S_core": self.xi_S_core,
            "xi_all_core": self.xi_all_core,
            "X_all_core": self.X_all_core,
            "P_CS": self.P_CS,
            "core_segregation_type": self.core_segregation_type,
            "inner_core_segregation_type": self.inner_core_segregation_type,
            "E_tot_should": self.E_tot_should / (G * m_earth**2 / r_earth * 3 / 5),
            "L_int_should": self.L_int_should
            / (4 * np.pi * r_earth**2 * sigmaSB * 300**4),
            "L_eff_should": self.L_eff_should
            / (4 * np.pi * r_earth**2 * sigmaSB * 300**4),
        }

        self.finals = {
            "P_surface_is": self.P_surface_is,
            "M_surface_is": self.M_surface_is,
            "T_surface_is": self.T_surface_is,
            "R_surface_is": self.R_surface_is,
            "Mg_number_is": self.Mg_number_is,
            "Si_number_is": self.Si_number_is,
            "M_H2O_is": self.M_H2O_is,
            "M_H2O_should": self.M_H2O_should,
            "M_H2O_hidden_is": self.M_H2O_hidden_is,
            "density_inversion": self.density_inversion,
            "match": self.match,
            "vapor_reached": self.vapor_reached,
            "layer_properties": self.layer_properties,
            "N_shells": self.N_shells,
            "M_ocean_is": self.M_ocean_is,
            "ocean_frac_is": self.ocean_frac_is,
            "M_core_is": self.M_core_is,
            "M_DWR_is": self.M_DWR_is,
            "xi_H_core": self.xi_H_core,
            "xi_FeO_mantle": self.xi_FeO_mantle,
            "M_H2O_core": self.M_H2O_core,
            "ocean_depth": self.ocean_depth,
            "N_shells_layers": self.N_shells_layers,
            "M_H2O_mantle": self.M_H2O_mantle,
            "M_ocean": self.M_ocean,
            "rho_mean": self.rho_mean,
            "M_mantle_is": self.M_mantle_is,
            "P_CS": self.P_CS,
            "T_CS": self.T_CS,
            "xi_Fe_mantle": self.xi_Fe_mantle,
            "Si_number_mantle": self.Si_number_mantle,
            "mantle_exists": self.mantle_exists,
            "inner_core_exists": self.inner_core_exists,
            "outer_core_exists": self.outer_core_exists,
            "inner_core_frac": self.inner_core_frac,
            "E_tot_is": self.E_tot_is,
            "L_int_should": self.L_int_is,
            "L_eff_should": self.L_eff_is,
        }

    def ComputeCoreMass(self, n=3):
        # contents = None, Mg_number=None, M_surface = 1.,
        #                 Mg_number_mantle = None, SiMg = None, M_ocean=0.,
        #                 xi_H_core=0., impurity=[], xi_S_core=0., n=3,
        #                 xi_all_core = [], M_IC = 0.):
        """Computes the core mass of a planet at given total mass, composition and
        value for Mg#
        """
        """
        print ('\nCompute Core Mass:')
        print ('contents = ', contents)
        print ('Mg_number =', Mg_number)
        print ('M_surface =', M_surface)
        print ('Mg_number_mantle = ', Mg_number_mantle)
        print ('SiMg =', SiMg)
        print ('M_ocean =', M_ocean)
        print ('xi_H_core =', xi_H_core)
        print ('impurity =', impurity)
        print ('xi_S_core =', xi_S_core)
        """
        # print (contents, Mg_number, M_surface, Mg_number_mantle, SiMg, M_ocean, xi_H_core)
        self.Mg_number_mantle = min(self.Mg_number_mantle, 0.9999999999)
        FeMg_mantle = (1.0 - self.Mg_number_mantle) / self.Mg_number_mantle
        FeMg = (1.0 - self.Mg_number_should) / self.Mg_number_should
        SiMg = self.Si_number_mantle / (1.0 - self.Si_number_mantle)

        # Compute the fractions in the mantle
        fractions = fortfunctions.functionspy.compute_abundance_vector(
            simg=SiMg,
            femg=FeMg_mantle,
            n_mats=len(self.contents[n]),
            ymgi=[material_YMg[i - 1] for i in self.contents[n]],
            ysii=[material_YSi[i - 1] for i in self.contents[n]],
            xih2oi=[0.0 for i in self.contents[n]],
            xifei=[1.0 - self.Mg_number_mantle for i in self.contents[n]],
            xialsii=[0.0 for i in self.contents[n]],
            xialmgi=[0.0 for i in self.contents[n]],
            contents=self.contents[n],
        )
        # Count
        # Note that the indices are shifted by one because the original definition
        # of the arrays comes from the fortran code.
        Q1 = sum(
            [
                fractions[i]
                * self.Mg_number_mantle
                * material_YMg[self.contents[n][i] - 1]
                for i in range(len(self.contents[n]))
            ]
        )

        # Compute total normalized mass in the mantle
        Q2 = (
            sum(
                [
                    fractions[i]
                    * self.Mg_number_mantle
                    * material_YMg[self.contents[n][i] - 1]
                    for i in range(len(self.contents[n]))
                ]
            )
            * mMg
            + sum(
                [
                    fractions[i]
                    * (1.0 - self.Mg_number_mantle)
                    * material_YMg[self.contents[n][i] - 1]
                    for i in range(len(self.contents[n]))
                ]
            )
            * mFe
            + sum(
                [
                    fractions[i] * material_YSi[self.contents[n][i] - 1]
                    for i in range(len(self.contents[n]))
                ]
            )
            * mSi
            + sum(
                [
                    fractions[i] * material_YO[self.contents[n][i] - 1]
                    for i in range(len(self.contents[n]))
                ]
            )
            * mO
            + sum(
                [
                    fractions[i] * material_YH[self.contents[n][i] - 1]
                    for i in range(len(self.contents[n]))
                ]
            )
            * mH
        )

        Q3 = sum(
            [
                fractions[i]
                * (1.0 - self.Mg_number_mantle)
                * material_YMg[self.contents[n][i] - 1]
                for i in range(len(self.contents[n]))
            ]
        )

        # Q4 = 1. #- xi_H_core
        # Q5 = (mFe + xi_S_core*mS)#*(1.-xi_H_core) + xi_H_core*mH
        # core_frac = (1.-M_ocean/self.M_surface_should)*(Q1/Q2-self.Mg_number_should*(Q1/Q2+Q3/Q2))/\
        #             (self.Mg_number_should*(Q4/Q5-Q1/Q2-Q3/Q2)+Q1/Q2)

        Q4 = self.xi_all_core[0]
        m_core = [mFe, mH, mS, mSi, mO]
        Q5 = sum(
            [m_core[i] * self.xi_all_core[i] for i in range(len(self.xi_all_core))]
        )
        # print ('Q =', Q1, Q2, Q3, Q4, Q5)
        core_frac = 1.0 - M_ocean / self.M_surface_should
        core_frac *= Q1 / Q2 - self.Mg_number_should * (Q1 / Q2 + Q3 / Q2)
        core_frac /= self.Mg_number_should * (Q4 / Q5 - Q1 / Q2 - Q3 / Q2) + Q1 / Q2
        # print ('core mass old =', core_frac * M_surface)
        core_frac = 1.0 - M_ocean / self.M_surface_should
        core_frac *= Q3 / Q2 - Q1 / Q2 * FeMg
        core_frac += M_IC / self.M_surface_should * (1.0 / mFe - Q4 / Q5)
        core_frac /= Q3 / Q2 - Q4 / Q5 - FeMg * Q1 / Q2
        # print ('')

    def prt(self, digits=4, **kwargs):
        print("=======================================")
        print("Planet properties:")
        print("\n--------------------------------")
        print("Layer details:")
        print("--------------------------------")
        for i in range(len(self.layers)):
            print("\nlayer: ", i)
            first_shell = self.layers[i].shells[0]
            for j in range(len(self.contents[i])):
                print(
                    round(first_shell.fractions[j] * 100, digits),
                    "%",
                    material_list[self.contents[i][j]],
                    " ph =",
                    self.layers[i].shells[-1].mix.mix[j].phase,
                )

            print(
                "layer mass [M_earth]: ",
                round(self.layers[i].indigenous_mass / m_earth, digits),
            )
            print(
                "outer radius [R_earth]:",
                round(self.layers[i].radius / r_earth, digits),
            )
            print("outer P [GPa]: ", round(self.layers[i].pres * 1.0e-9, digits))
            print("outer T [K]: ", round(self.layers[i].temp, digits))

        try:
            print("\n--------------------------------")
            print("Major parameters:")
            print("--------------------------------")
            print(
                "R_surface_is [R_earth]:",
                round(self.R_surface_is / r_earth, digits),
                "\nM_surface_is [M_earth]:",
                round(self.M_surface_is / m_earth, digits),
                "\nrho mean [gcc]:",
                round(self.rho_mean / 1000, digits),
                "\nT_surface_is [K]:",
                round(self.T_surface_is, digits),
                "\nT_surface_should [K]:",
                round(self.T_surface_should, digits),
                "\nT_center [K]:",
                round(self.T_center, digits),
                "\nP_surface_is [bar]:",
                round(self.P_surface_is * 1.0e-5, digits),
                "\nP_center [GPa]:",
                round(self.P_center * 1.0e-9, digits),
                "\nMOI factor:",
                round(self.MOI_is, digits),
                "\n\nMg_number_should:",
                round(self.Mg_number_should, digits),
                "\nMg_number_is:",
                round(self.Mg_number_is, digits),
                "\nSi_number_should:",
                round(self.Si_number_should, digits),
                "\nSi_number_is:",
                round(self.Si_number_is, digits),
                #'\nxi_H_core:', round(self.xi_H_core, digits),
                #'\nxi_H_core_predicted:', round(self.xi_H_core_predicted, digits),
                #'\nP_H2_CMB [MPa]:', round(self.P_H2_CMB*1.0e-6, digits),
                "\nL_int_should [W]:",
                round(
                    self.L_int_should / (4 * np.pi * r_earth**2 * sigmaSB * 300**4),
                    digits,
                ),
                "\nL_int_is [W]:",
                round(
                    self.L_int_is / (4 * np.pi * r_earth**2 * sigmaSB * 300**4),
                    digits,
                ),
                "\nE_tot_should [G M_E²/R_E]:",
                round(self.E_tot_should / (G * m_earth**2 / r_earth * 3 / 5), digits),
                "\nE_tot_is [G M_E²/R_E]:",
                round(self.E_tot_is / (G * m_earth**2 / r_earth * 3 / 5), digits),
            )

            try:
                print(
                    "M_H2O [wt%]:",
                    round(
                        (self.H2O_count + self.H_count / 2.0)
                        * mH2O
                        / self.M_surface_is
                        * 100,
                        digits,
                    ),
                )
                print(
                    "M_H2O core [wt%]:",
                    round(self.H_count / 2.0 * mH2O / self.M_surface_is * 100, digits),
                )
                print(
                    "M_H2O mantle [wt%]:",
                    round(self.H2O_count * mH2O / self.M_surface_is * 100, digits),
                )
                print("Ocean frac is:", round(10**self.ocean_frac_is, digits))
                print("Ocean frac should:", round(10**self.ocean_frac_should, digits))
                print(
                    "Core mass frac:",
                    round(self.M_core_is / self.M_surface_is * m_earth, digits),
                )
                print(
                    "Core radius [km]:",
                    round(self.layer_properties[1]["R_outer"] / 1000, digits),
                )
            except ZeroDivisionError:
                print("M_H2O [wt%]: NaN")

        except TypeError:
            print("WARNING: Type Error in Planet.prt()")

        print("\n--------------------------------")
        print("Layer overview:")
        print("--------------------------------")

        dat = []
        for i in range(len(self.layer_properties)):
            lay = self.layer_properties[i]
            material_str = ""
            for c in range(len(self.contents[i])):
                frac = str(round(self.fractions[i][c] * 100, 1)) + "% "
                material_str += frac + material_list_fort[self.contents[i][c] - 1]
                if c < len(self.contents[i]) - 1:
                    material_str += ", "

            dat.append(
                [
                    i,
                    material_str,
                    ftool.scinot(lay["R_outer"] / r_earth, digits=digits),
                    ftool.scinot(lay["indigenous_mass"] / m_earth, digits=digits),
                    ftool.scinot(lay["P_outer"] * 1.0e-9, digits=digits),
                    ftool.scinot(lay["T_outer"], digits=digits),
                    ftool.scinot(lay["rho_outer"], digits=digits),
                ]
            )

        tabl = tabulate(
            dat,
            headers=[
                "Layer",
                "Contents",
                "R [R_e]",
                "m [M_e]",
                "P [GPa]",
                "T [K]",
                "rho [kg m-3]",
            ],
        )

        print()
        print(f"{tabl}")
        print()

    def complete_layer_properties(self):
        for i in range(len(self.layer_properties)):
            self.layer_properties[i].update({"x_SiO2": self.xi_SiO2_mantle})
            self.layer_properties[i].update({"x_MgO": self.xi_MgO_mantle})
            self.layer_properties[i].update({"x_FeO": self.xi_FeO_mantle})
            self.layer_properties[i].update(
                {
                    "mass_fraction": self.layer_properties[i]["indigenous_mass"]
                    / self.M_surface_is
                }
            )

            if i > 0:
                self.layer_properties[i]["P_inner"] = self.layer_properties[i - 1][
                    "P_outer"
                ]
                self.layer_properties[i][
                    "T_inner"
                ] = None  # self.layer_properties[i - 1]['T_outer']
                self.layer_properties[i]["R_inner"] = self.layer_properties[i - 1][
                    "R_outer"
                ]

            else:
                self.layer_properties[i]["P_inner"] = self.P_center
                self.layer_properties[i]["T_inner"] = self.T_center
                self.layer_properties[i]["R_inner"] = self.R_seed

        for i in range(len(self.cleaned_layer_properties)):
            self.cleaned_layer_properties[i].update({"x_SiO2": self.xi_SiO2_mantle})
            self.cleaned_layer_properties[i].update({"x_MgO": self.xi_MgO_mantle})
            self.cleaned_layer_properties[i].update({"x_FeO": self.xi_FeO_mantle})
            self.cleaned_layer_properties[i].update(
                {
                    "mass_fraction": self.cleaned_layer_properties[i]["indigenous_mass"]
                    / self.M_surface_is
                }
            )

            if i > 0:
                self.cleaned_layer_properties[i][
                    "P_inner"
                ] = self.cleaned_layer_properties[i - 1]["P_outer"]
                self.cleaned_layer_properties[i][
                    "T_inner"
                ] = None  # self.layer_properties[i - 1]['T_outer']
                self.cleaned_layer_properties[i][
                    "R_inner"
                ] = self.cleaned_layer_properties[i - 1]["R_outer"]

            else:
                self.cleaned_layer_properties[i]["P_inner"] = self.P_center
                self.cleaned_layer_properties[i]["T_inner"] = self.T_center
                self.cleaned_layer_properties[i]["R_inner"] = self.R_seed

    def extract_bulk_structure(self, trim=False):
        # Extract bulk layer strucutre
        # Offset by 3 because 0,1,2 = liquid, supercrit., solid water
        # ommit last shell because hydrosphere is treated seperately below
        # layer only exists if it contains more than the one base shell
        self.bulk_structure = []
        self.cleaned_layer_properties = []
        for i in range(len(self.N_shells_layers) - 1):
            k = i + 3
            if self.N_shells_layers[i] > 1:
                self.bulk_structure.append(k)
                self.cleaned_layer_properties.append(self.layer_properties[i])

        # Extract sublayers of hydrosphere
        self.extract_hydro_structure(trim=trim)

        # Add hydro structure to bulk structure
        self.bulk_structure += self.hydro_structure
        self.complete_layer_properties()

    def extract_hydro_structure(self, trim=False):
        phases = ["liquid", "supercritical", "solid", "vapor", "gas"]
        combos = [["liquid"], ["supercritical"], ["solid"], ["vapor"], ["gas"]]

        if trim:
            self.trim_profiles()

        # Extract hydrosphere structure types
        phase_change = False
        # Start at bottom of the hydrosphere
        T, P = self.layer_properties[3]["T_outer"], self.layer_properties[3]["P_outer"]
        ph = phase.phase(T=T, P=P)  # H2O phase
        self.phases = [ph]

        # enclosed masses at each phase transition
        self.dms = [self.M_surface_is * (1.0 - 10**self.ocean_frac_is)]
        # radius at each phase transition
        self.locs = [self.layer_properties[3]["R_outer"]]
        # pressure and temperature
        self.temps = [self.layer_properties[3]["T_outer"]]
        self.pres = [self.layer_properties[3]["P_outer"]]
        self.dens_out = [self.layer_properties[3]["rho_outer"]]
        self.dens_in = [self.layer_properties[3]["rho_inner"]]
        # volume
        self.dvols = [self.locs[0] ** 3 * 4 / 3 * np.pi]
        M0 = 0.0
        # Loop over all shells in the hydrosphere
        for i in range(len(self.profiles[0])):
            R, T, P, d, M = self.profiles[0:5, i]
            if R >= self.layer_properties[3]["R_outer"]:
                ph = phase.phase(T=T, P=P)

                # If no phase change do nothing
                if ph == self.phases[-1]:
                    pass
                # New phase region reached
                else:
                    lower = self.profiles[0:5, i - 1]
                    upper = self.profiles[0:5, i]
                    R, T, P, M = refine_hydro_structure(lower, upper, res=4)
                    self.phases.append(ph)
                    self.locs.append(R)
                    self.temps.append(T)
                    self.pres.append(P)
                    self.dens_in.append(self.dens_out[-1])
                    self.dens_out.append(d)
                    self.dvols.append(4 / 3 * np.pi * R**3)
                    self.dms.append(M)

                    # Update cleaned layer properties
                    self.cleaned_layer_properties.append(
                        {
                            "P_outer": self.pres[-1],
                            "T_outer": self.temps[-1],
                            "rho_outer": self.dens_out[-1],
                            "rho_inner": self.dens_in[-1],
                            "indigenous_mass": self.dms[-1] - self.dms[-2],
                            "mass_fraction": (self.dms[-1] - self.dms[-2])
                            / self.M_surface_is,
                            "R_outer": self.locs[-1],
                            "n_shells": None,
                        }
                    )

        self.locs.append(self.R_surface_is)
        self.temps.append(self.T_surface_is)
        self.pres.append(self.P_surface_is)
        self.dms.append(self.M_surface_is)
        self.dvols.append(self.R_surface_is**3 * 4 / 3 * np.pi)

        # Update cleaned layer properties
        self.cleaned_layer_properties.append(
            {
                "P_outer": self.pres[-1],
                "T_outer": self.temps[-1],
                "rho_outer": self.dens_out[-1],
                "rho_inner": self.dens_in[-1],
                "indigenous_mass": self.dms[-1] - self.dms[-2],
                "mass_fraction": (self.dms[-1] - self.dms[-2]) / self.M_surface_is,
                "R_outer": self.locs[-1],
                "n_shells": None,
            }
        )

        # mass and radius increments for each phase region
        hdiffs = [(self.locs[i + 1] - self.locs[i]) for i in range(len(self.locs) - 1)]
        mdiffs = [(self.dms[i + 1] - self.dms[i]) for i in range(len(self.dms) - 1)]
        voldiffs = [
            (self.dvols[i + 1] - self.dvols[i]) for i in range(len(self.dvols) - 1)
        ]

        # print ('\nphases =', phases)
        # print ('heights =', diffs)

        # self.structure_types = [list(t) for t in set(tuple(element) for element in all_phases)]
        # print ('\nstructure types =')
        # for s in structure_types: print (s)
        self.heights = []
        self.masses = []
        self.vols = []
        self.hydro_structure = [phases.index(self.phases[0])]

        for j in range(len(self.phases) - 1):
            if self.phases[j] != self.phases[j + 1]:
                self.hydro_structure.append(phases.index(self.phases[j + 1]))

        # Compute cumulative height, mass and volume of different phase regions
        for o in range(3):
            height = 0.0
            mass = 0.0
            vol = 0.0
            for j in range(len(self.phases)):
                # Phase change¨
                if self.phases[j] == phases[o]:
                    height += hdiffs[j]
                    mass += mdiffs[j]
                    vol += voldiffs[j]
            self.heights.append(height)
            self.masses.append(mass / self.M_surface_is)
            self.vols.append(vol)

    def update_oxide_fractions(self):
        # Compute mole fractions of oxides in the mantle
        if self.xi_Fe_mantle > 0.0:
            self.xi_FeO_mantle = 1.0 - self.Si_number_mantle
            self.xi_FeO_mantle /= (
                (1.0 - self.xi_Fe_mantle) / self.xi_Fe_mantle
                + 1.0
                - self.Si_number_mantle
            )

        else:
            self.xi_FeO_mantle = 0.0

        self.xi_SiO2_mantle = (1.0 - self.xi_FeO_mantle) * self.Si_number_mantle
        self.xi_MgO_mantle = 1.0 - self.xi_FeO_mantle - self.xi_SiO2_mantle

        # Update oxygen fugacity
        self.logfO2 = material.logfO2(
            self.P_CS, self.xi_all_core[0], self.xi_FeO_mantle
        )

    def compute_sulfur_in_mantle(self, **kwargs):
        # Compute mole fraction of sulfur in mantle from S partitioning
        self.xi_S_mantle = material.xi_S(self.T_CS, self.P_CS, self.xi_all_core)

        # Update total sulfur content
        N_mantle = self.xi_FeO_mantle * (mO + mFe)
        N_mantle += self.xi_MgO_mantle * (mO + mMg)
        N_mantle += self.xi_SiO2_mantle * (2 * mO + mSi)
        N_mantle = self.M_mantle_is * m_earth / N_mantle

        self.S_count += N_mantle * self.xi_S_mantle

    def compute_luminosity(self, S0, eps=0.0, alpha=0.0):
        self.L_eff_is = sigmaSB * self.T_surface_is**4 * (1.0 - eps / 2.0)
        self.L_eff_is -= 0.25 * S0 * (1.0 - alpha)
        self.L_eff_is *= 4 * np.pi * self.R_surface_is**2
        self.L_int_is = (
            sigmaSB * self.T_surface_is**4 * 4.0 * np.pi * self.R_surface_is**2
        )

    def update_core_mass(self, **kwargs):
        """Estimates the core mass fraction for a given composition
        and boundary conditions after integration with the now known value
        for xi_H_core that would be required to match the Mg#
        """
        cm = compute_core_mass(
            Mg_number=self.Mg_number_should,
            M_surface=self.M_surface_should / m_earth,
            Mg_number_mantle=1.0 - self.Fe_number_layers[2],
            M_ocean=10 ** (self.ocean_frac_should) * self.M_surface_should / m_earth,
            SiMg=self.Si_number_should / (1.0 - self.Si_number_should),
            contents=self.contents,
            xi_H_core=self.xi_H_core_predicted,
            impurity=[self.xi_impurity],
        )
        # xi_S_core = self.xi_S_core)

        # print ('new core mass frac =', cm/self.M_surface_should*m_earth)

        self.initials["layermasses"][0] = cm * self.initials["inner_core_frac"]
        self.initials["layermasses"][1] = cm - self.initials["layermasses"][0]

    def Update(self):
        # Check for inner core
        self.inner_core_exists = True
        if self.layer_properties[0]["R_outer"] < 1e3:
            self.inner_core_exists = False

        # Check for outer core
        self.outer_core_exists = True
        if (
            self.layer_properties[1]["R_outer"] - self.layer_properties[0]["R_outer"]
            < 1e3
        ):
            self.outer_core_exists = False

    def Update_finals(self, **kwargs):
        # gather final properties in a dictionary
        self.finals = {
            "P_surface_is": self.P_surface_is,
            "M_surface_is": self.M_surface_is / m_earth,
            "T_surface_is": self.T_surface_is,
            "R_surface_is": self.R_surface_is / r_earth,
            "Mg_number_is": self.Mg_number_is,
            "Si_number_is": self.Si_number_is,
            "M_H2O_is": self.M_H2O_is,
            "M_H2O_should": self.M_H2O_should,
            "M_H2O_hidden_is": self.M_H2O_hidden_is,
            "density_inversion": self.density_inversion,
            "match": self.match,
            "vapor_reached": self.vapor_reached,
            "layer_properties": self.layer_properties,
            "N_shells": self.N_shells,
            "M_ocean_is": self.M_ocean_is,
            "ocean_frac_is": self.ocean_frac_is,
            "M_core_is": self.M_core_is,
            "M_DWR_is": self.M_DWR_is,
            "MOI_is": self.MOI_is,
            "xi_H_core": self.xi_H_core,
            "xi_FeO_mantle": self.xi_FeO_mantle,
            "M_H2O_core": self.M_H2O_core,
            "ocean_depth": self.ocean_depth,
            "N_shells_layers": self.N_shells_layers,
            "M_H2O_mantle": self.M_H2O_mantle,
            "M_ocean": self.M_ocean,
            "rho_mean": self.rho_mean,
            "M_mantle_is": self.M_mantle_is,
            "P_CS": self.P_CS,
            "T_CS": self.T_CS,
            "xi_Fe_mantle": self.xi_Fe_mantle,
            "Si_number_mantle": self.Si_number_mantle,
            "mantle_exists": self.mantle_exists,
            "inner_core_exists": self.inner_core_exists,
            "outer_core_exists": self.outer_core_exists,
            "inner_core_frac": self.inner_core_frac,
            "E_tot_is": self.E_tot_is,
            "L_int_is": self.L_int_is,
            "L_eff_is": self.L_eff_is,
        }

    def Update_initials(self, **kwargs):
        """ """
        if self.M_surface_should == float("inf"):
            print("HERE infinity in Planet.update_initials()")

        # gather all initial properties in a dictionary
        self.initials = {
            "P_center": self.P_center,
            "T_center": self.T_center,
            "R_seed": self.R_seed,
            "contents": self.contents,
            "tempType": self.tempType,
            "layerType": self.layerType,
            "fractions": self.fractions,
            "P_surface_should": self.P_surface_should,
            "majorConstraint": self.majorConstraint,
            "layerConstraint": self.layerConstraint,
            "add_atmosphere": self.add_atmosphere,
            "eps_r": self.eps_r,
            "layermasses": self.layermasses,
            "layerradii": self.layerradii,
            "differentiated": self.differentiated,
            "M_surface_should": self.M_surface_should / m_earth,
            "T_surface_should": self.T_surface_should,
            "R_surface_should": self.R_surface_should / r_earth,
            "rhoType": self.rhoType,
            "echo": self.echo,
            "Mg_number_should": self.Mg_number_should,
            "Si_number_should": self.Si_number_should,
            "Si_number_layers": self.Si_number_layers,
            "Fe_number_layers": self.Fe_number_layers,
            "vapor_stop": self.vapor_stop,
            "Mg_number_layer": self.Mg_number_layer,
            "modelType": self.modelType,
            "temp_jumps": self.temp_jumps,
            "gammas_layer": self.gammas_layer,
            "M_ocean_should": self.M_ocean_should,
            "ocean_frac_should": self.ocean_frac_should,
            "layertemps": self.layertemps,
            "T_zero": self.T_zero,
            "M_core_should": self.M_core_should,
            "subphase_res": self.subphase_res,
            "layerpres": self.layerpres,
            #'xi_Stv':self.xi_Stv,
            "omega": self.omega,
            #'xi_S_core':self.xi_S_core,
            "xi_all_core": self.xi_all_core,
            "X_all_core": self.X_all_core,
            "P_CS": self.P_CS,
            "core_segregation_type": self.core_segregation_type,
            "inner_core_segregation_type": self.inner_core_segregation_type,
            "E_tot_should": self.E_tot_should / (G * m_earth**2 / r_earth * 3 / 5),
            "L_int_should": self.L_int_should
            / (4 * np.pi * r_earth**2 * sigmaSB * 300**4),
            "L_eff_should": self.L_eff_should
            / (4 * np.pi * r_earth**2 * sigmaSB * 300**4),
        }

    def Reset(self, **kwargs):
        """Resets all planetary properties to the initial state. This allows
        to use the exact same specifications for re-constructing the planet.
        """

        # delete old properties to release memory
        for lay in self.layers:
            del lay.shells
            del lay

        # print ('layermasses in Planet.Reset():', self.initials['layermasses'])
        self.__init__(**self.initials)
        # print ('layermasses in Planet.Reset():', self.initials['layermasses'])
        if self.initials["M_surface_should"] > 1.0e30:
            print("HERE large after Planet.reset() in Planet")
            print("large value=", self.initials["M_surface_should"])

        if self.M_surface_should == float("inf"):
            print("HERE infinity in Planet.reset() after")

    def Construct(self, print_time=False, fortran=True, echo=False, **kwargs):
        """N sets the number of shells in the case that no surface constrain
        type has been specified. Otherwise N will be ignored.
        """
        # Gather layer dims for fortplanet routine
        layer_dims = []
        conts = []
        fracs = []

        # print ('xi_all_core before construct =', self.xi_all_core)
        for i in range(len(self.contents)):
            layer_dims.append(len(self.contents[i]))
            for j in range(len(self.contents[i])):
                conts.append(self.contents[i][j])
                fracs.append(self.fractions[i][j])

        print("fractions before construct =", self.fractions)
        print("contents before construct =", self.contents)
        print("xi_all_core before construct =", self.xi_all_core)
        print("x_all_core before construct =", self.X_all_core)

        # print ('fracs, conts new =', fracs, conts)
        (
            self.M_surface_is,
            self.R_surface_is,
            self.P_surface_is,
            self.T_surface_is,
            self.Mg_number_is,
            self.Si_number_is,
            self.Fe_count,
            self.Si_count,
            self.Mg_count,
            self.O_count,
            self.H2O_count,
            self.H_count,
            self.S_count,
            self.ocean_frac_is,
            self.MOI_is,
            layer_props_dummy,
            self.profiles,
            self.N_shells,
            self.N_shells_layers,
            fractions_dummy,
            self.xi_Fe_mantle,
            self.Si_number_mantle,
            self.mantle_exists,
            self.inner_core_exists,
            self.outer_core_exists,
            self.E_grav,
            self.E_int,
        ) = fortplanet.wrapper.create_planet(
            t_center=self.T_center,
            p_center=self.P_center,
            contents=conts,
            fractions=fracs,
            layer_dims=layer_dims,
            layer_masses=self.layermasses,
            layer_pres=self.layerpres,
            layer_radii=self.layerradii,
            layer_constraints=self.layerConstraint,
            temp_jumps=self.temp_jumps,
            q_layers=self.q,
            gammag0_layers=self.gammas_layer,
            si_number_layers=self.Si_number_layers,
            fe_number_layers=self.Fe_number_layers,
            r_seed=self.R_seed,
            adiabattype=self.adiabatType,
            rhotype=self.rhoType,
            layertype=self.layerType,
            temptype=self.tempType,
            eps_r=self.eps_r,
            # eps_h2o=self.eps_H2O,
            # eps_al=self.eps_Al,
            p_surface_should=self.P_surface_should,
            t_surface_should=self.T_surface_should,
            spin=self.omega,
            # xi_h_core_predicted=self.xi_H_core_predicted,
            subphase_res=self.subphase_res,
            # xi_stv=self.xi_Stv,
            # x_impurity_0_layers=self.X_impurity_0_layers,
            # x_impurity_slope_layers=self.X_impurity_slope_layers,
            xi_all_core=self.xi_all_core,
            x_all_core=self.X_all_core,
            pres_core_segregation=self.P_CS,
            core_segregation_type=self.core_segregation_type,
            m_ocean_should=self.M_ocean_should,
            m_surface_should=self.M_surface_should,
            mg_number_should=self.Mg_number_should,
            inner_core_segregation_type=self.inner_core_segregation_type,
            r_surface_should=1.0,
            m_core_should=1.0,
        )

        self.status = "constructed"
        for i in range(len(layer_props_dummy)):

            self.layer_properties[i]["P_outer"] = layer_props_dummy[i][0]
            self.layer_properties[i]["T_outer"] = layer_props_dummy[i][1]
            self.layer_properties[i]["rho_outer"] = layer_props_dummy[i][2]
            self.layer_properties[i]["rho_inner"] = layer_props_dummy[i][3]
            self.layer_properties[i]["indigenous_mass"] = layer_props_dummy[i][4]
            self.layer_properties[i]["mass_fraction"] = (
                layer_props_dummy[i][4] / self.M_surface_is
            )
            self.layer_properties[i]["R_outer"] = layer_props_dummy[i][5]
            self.layer_properties[i]["n_shells"] = self.N_shells_layers[i]

            if i > 0:
                self.layer_properties[i]["P_inner"] = self.layer_properties[i - 1][
                    "P_outer"
                ]
                self.layer_properties[i][
                    "T_inner"
                ] = None  # self.layer_properties[i - 1]['T_outer']
                self.layer_properties[i]["R_inner"] = self.layer_properties[i - 1][
                    "R_outer"
                ]

            else:
                self.layer_properties[i]["P_inner"] = self.P_center
                self.layer_properties[i]["T_inner"] = self.T_center
                self.layer_properties[i]["R_inner"] = self.R_seed

        # Check for inner core
        self.inner_core_exists = True
        if self.layer_properties[0]["R_outer"] < 1e3:
            self.inner_core_exists = False

        # Check for outer core
        self.outer_core_exists = True
        if (
            self.layer_properties[1]["R_outer"] - self.layer_properties[0]["R_outer"]
            < 1e3
        ):
            self.outer_core_exists = False

        # If no mantle exists don't update the fractions as they would be set
        # to zero in the mantle for the next iteration in all subsequent layers
        # which leads to errors in the computation of the metal-silicate partitioning
        # and sucks in general. I don't know if this really fixes the current
        # sweep or if it just avoids a segfault but worst case scenario is that
        # this one planet went to waste as it didn't properly converge.
        if self.mantle_exists:
            for i in range(len(self.fractions)):
                for j in range(len(self.fractions[i])):
                    self.fractions[i][j] = fractions_dummy[i][j]

        print("fractions after construct =", self.fractions)
        print("contents after construct =", self.contents)

        print("M_surface is = ", self.M_surface_is)
        print("R_surface is = ", self.R_surface_is)
        print("T_surface is = ", self.T_surface_is)
        print("P_surface is = ", self.P_surface_is)
        print("Mg_number is = ", self.Mg_number_is)
        print("Si_number is = ", self.Si_number_is)

        # print ('fractions after construct =', self.fractions)
        self.MOI_is = self.MOI_is / self.M_surface_is / self.R_surface_is**2
        self.M_core_is = (
            self.layer_properties[0]["indigenous_mass"]
            + self.layer_properties[1]["indigenous_mass"]
        ) / m_earth
        self.M_mantle_is = (
            self.layer_properties[2]["indigenous_mass"]
            + self.layer_properties[3]["indigenous_mass"]
        ) / m_earth

        self.M_H2O_core = self.H_count / 2.0 * mH2O / self.M_surface_is
        self.M_H2O_is = (self.H2O_count + self.H_count / 2.0) * mH2O / m_earth
        self.M_ocean = 10**self.ocean_frac_is * self.M_surface_is / m_earth
        self.M_H2O_mantle = self.M_H2O_is - self.M_H2O_core - self.M_ocean
        self.inner_core_frac = self.layer_properties[0]["indigenous_mass"] / m_earth
        self.inner_core_frac /= self.M_core_is
        self.L_int_is = (
            4 * np.pi * sigmaSB * self.T_surface_is**4 * self.R_surface_is**2
        )
        self.E_tot_is = self.E_grav + self.E_int

        """
        #Update inner core mass fraction
        try:
            self.inner_core_frac = self.layer_properties[0]['indigenous_mass']/m_earth /\
                self.M_core_is
        except ZeroDivisionError:
            self.inner_core_frac = 0.
        """
        try:
            self.ocean_depth = (
                self.layer_properties[4]["R_outer"]
                - self.layer_properties[3]["R_outer"]
            ) / r_earth
        except IndexError:
            self.ocean_depth = 0.0
        # Update mean density of the planet
        self.rho_mean = self.M_surface_is / (4.0 / 3.0 * np.pi * self.R_surface_is**3)
        # Compute mol and wt fractions of the core composition
        self.xi_all_core = material.mat2at_core(
            xi=self.fractions[1], xiH=self.xi_H_core
        )
        self.X_all_core = material.at2wt_core(self.xi_all_core)
        # Compute mol fractions of different oxides in the mantle and the
        # oxygen fugacity of the mantle
        self.update_oxide_fractions()
        # self.complete_layer_properties()
        self.Update_finals()
        self.Update_initials()

    def arange_profiles(self):
        """Aranges the profile data for each individual layer in a separate
        array for processing.
        """

        profs = []

        for k in range(len(self.profiles)):
            c = 0
            profs.append([])
            for i in range(len(self.contents)):
                profs[k].append([])
                for j in range(self.N_shells_layers[i]):
                    profs[k][i].append(self.profiles[k][c])

                    c += 1

        self.layer_profiles = np.array(profs)

    def generate_simple_model(self, res=0, equi="mass"):
        """Takes the averages of each layer and constructs a simple N-layer model
        of the planet. This allows to compare the more rigorous structure model
        with simple models that are often used in the literature to model e.g.
        the properties of the icy satellites (s.a. radius, MoI factor, core mass)
        """
        N = len(self.layer_properties)

        # Take each layer and compute the mass averaged density
        for i in range(N):
            if i > 0:
                R_in = self.layer_properties[i - 1]["R_outer"]

            else:
                R_in = 0.0

            R_out = self.layer_properties[i]["R_outer"]

            M = self.layer_properties[i]["indigenous_mass"]

            rho = 3.0 / 4.0 * M / np.pi / (R_out**3 - R_in**3)

            if np.isnan(rho):
                rho = 0.0

            self.simple_densities.append(rho)

            print("---")
            print("R_in/Rout = ", R_in, R_out)
            print("rho mean =", rho)
            print("mass =", M)

        # Compute new MoI factor using the averaged densities
        self.simple_MoI = 0.0

        for i in range(N):
            if i > 0:
                R_in = self.layer_properties[i - 1]["R_outer"]

            else:
                R_in = 0.0

            R_out = self.layer_properties[i]["R_outer"]

            dMoI = (
                8.0 / 15.0 * np.pi * self.simple_densities[i] * (R_out**5 - R_in**5)
            )
            print("dMoI =", dMoI)
            self.simple_MoI += dMoI

        self.simple_MoI = self.simple_MoI / (self.M_surface_is * self.R_surface_is**2)

    def average_profiles(self, N, equiv="mass"):
        """


        Parameters
        ----------
        N : TYPE
            Number of sublayers each layer is devided in.
        equiv : TYPE, optional
            Parameter that is used to subdivide the layers. The default is 'mass'.

        Returns
        -------
        Average density for each sublayer.

        """

    def trim_profiles(self):
        """The profiles output from the Fortran routine contains a sequence
        of zeros at the end. These elements are being removed here. The
        convention of the profiles is:

            Radius (m), Temperature (K), Pressure (Pa), Density (kg/m3),
            Mass (kg), Gravity (m/s2), Gravitational energy (J)
        """
        off = len(self.contents) - 1
        prof = np.empty([len(self.profiles), self.N_shells + off])

        i = 0
        while True:
            for p in range(len(self.profiles)):
                prof[p][i] = self.profiles[p][i]

            i += 1

            if i == self.N_shells + off:
                break

        self.profiles = prof

    def plot(
        self,
        scatter=False,
        layerColoring=False,
        axis=[],
        x_axis="radius",
        save=True,
        file_name="planet_structure",
        spec="",
        file_dir="./",
        prem=False,
        **kwargs,
    ):

        if axis == []:
            fig, axis = plt.subplots(2, 3, sharex=True, figsize=(10, 4))
            fig.subplots_adjust(hspace=0.1, wspace=0.3)

        color = "k"
        lwdth = 0.5
        self.trim_profiles()
        X_Si = self.X_all_core[3]
        X_S = self.X_all_core[2]
        X_O = self.X_all_core[4]
        T_melt_Fe = material.T_melt_Fe(self.profiles[2], X_Si, X_O, X_S)
        T_melt_MgSiO3 = material.T_liquidus_pyrolite(self.profiles[2])
        T_melt_CaSiO3 = material.T_melt_MgSiO3(self.profiles[2])

        axis[0][0].plot(
            self.profiles[0] * 1e-3, self.profiles[1], linewidth=lwdth, color=color
        )
        axis[0][1].plot(
            self.profiles[0] * 1e-3,
            self.profiles[2] * 1.0e-9,
            linewidth=lwdth,
            color=color,
        )
        axis[0][2].plot(
            self.profiles[0] * 1e-3,
            self.profiles[3] * 1.0e-3,
            linewidth=lwdth,
            color=color,
        )
        axis[1][0].plot(
            self.profiles[0] * 1e-3,
            self.profiles[4] / m_earth,
            linewidth=lwdth,
            color=color,
        )
        axis[1][1].plot(
            self.profiles[0] * 1e-3, self.profiles[5], linewidth=lwdth, color=color
        )
        axis[1][2].plot(
            self.profiles[0] * 1e-3,
            self.profiles[7] / (G * m_earth**2) * r_earth,
            linewidth=lwdth,
            color=color,
        )

        # Plot PREM
        if prem:
            x, y = readPREM.Do()

            axis[0][2].plot(x * 1000 / r_earth, y / 1000, color="b")
            axis[0][2].text(0.6, 9.0, "PREM", color="b")
        """  
        axis[0][0].plot(self.profiles[0]/r_earth, T_melt_Fe, linewidth=lwdth,
                        color = 'red',
                        label=r'$T_{\rm melt}^{\rm Fe-S-Si-O}$')
        axis[0][0].plot(self.profiles[0]/r_earth, T_melt_MgSiO3, linewidth=lwdth,
                        color = 'green',
                        label=r'$T_{\rm liquidus}^{\rm Pyrolite} $')
        axis[0][0].plot(self.profiles[0]/r_earth, T_melt_CaSiO3, linewidth=lwdth,
                        color = 'grey',
                        label=r'$T_{\rm melt}^{\rm MgSiO_3}$')
        
        axis[0][0].legend(loc=3, fontsize=6)
        """
        M_H2O_core = round(self.H_count / 2.0 * mH2O / self.M_surface_is * 100, 3)
        M_H2O_mantle = round(self.H2O_count * mH2O / self.M_surface_is * 100, 3)

        ax = 1
        d = 0.1
        pos0 = -1.75
        axis[ax][0].text(
            0.05,
            -1.75,
            r"$\rm Mg \# = \ $" + str(round(self.Mg_number_is, 3)),
            transform=axis[0][0].transAxes,
        )
        axis[ax][0].text(
            0.05,
            pos0 - d,
            r"$\rm Si \# = \ $" + str(round(self.Si_number_is, 3)),
            transform=axis[0][0].transAxes,
        )
        axis[ax][0].text(
            0.05,
            pos0 - 2 * d,
            r"$\xi_{\rm Fe} = \ $" + str(round(self.xi_Fe_mantle, 4)),
            transform=axis[0][0].transAxes,
        )
        # axis[ax][0].text(.05, pos0-3*d,
        #                 r'$\xi_{\rm S} = \ $'+str(round(self.xi_S_core,3)),
        #               transform = axis[0][0].transAxes)
        axis[ax][0].text(
            0.05,
            pos0 - 4 * d,
            r"$T_{\rm S} = \ $" + str(round(self.T_surface_is, 1)),
            transform=axis[0][0].transAxes,
        )
        axis[ax][1].text(
            1.55,
            pos0,
            r"$M_{\rm H_2 O}/M = \ $"
            + str(round(self.M_H2O_is * 100, 3))
            + r"$\ \rm wt\%$",
            transform=axis[0][0].transAxes,
        )
        axis[ax][1].text(
            1.55,
            pos0 - d,
            r"$M_{\rm H_2 O, Core}/M = \ $" + str(M_H2O_core) + r"$\ \rm wt\%$",
            transform=axis[0][0].transAxes,
        )
        axis[ax][1].text(
            1.55,
            pos0 - 2 * d,
            r"$M_{\rm H_2 O, Mantle}/M = \ $" + str(M_H2O_mantle) + r"$\ \rm wt\%$",
            transform=axis[0][0].transAxes,
        )
        axis[ax][1].text(
            1.55,
            pos0 - 3 * d,
            r"$X_{\rm Ice}= \ $" + str(0.0) + r"$\ \rm wt\%$",
            transform=axis[0][0].transAxes,
        )
        axis[ax][1].text(
            1.55,
            pos0 - 4 * d,
            r"$X^{\prime}_{\rm Ice}= \ $" + str(0.0) + r"$\ {\rm wt\%} M_\oplus^{-1}$",
            transform=axis[0][0].transAxes,
        )

        axis[ax][2].text(
            3.05,
            pos0 - 0 * d,
            r"$M/M_\oplus= \ $" + str(round(self.M_surface_is / m_earth, 4)),
            transform=axis[0][0].transAxes,
        )

        axis[ax][2].text(
            3.05,
            pos0 - 1 * d,
            r"$R/R_\oplus= \ $" + str(round(self.R_surface_is / r_earth, 4)),
            transform=axis[0][0].transAxes,
        )
        axis[ax][2].text(
            3.05,
            pos0 - 2 * d,
            r"$I/MR^2= \ $" + str(round(self.MOI_is, 4)),
            transform=axis[0][0].transAxes,
        )

        y_labels = [
            [
                r"$\rm Temperature \ [K]$",
                r"$\rm Pressure \ [GPa]$",
                r"$\rm Density \ [g \ cm^{-3}]$",
            ],
            [
                r"$M/M_\oplus$",
                r"$\rm Gravity \ [kg \ m \ s^{-2}]$",
                r"$E_{\rm grav} / \ \left(GM_\oplus^2 / R_\oplus\right)$",
            ],
        ]

        for i in range(len(axis)):
            for j in range(len(axis[i])):
                axis[i][j].set_ylabel(y_labels[i][j])
                axis[i][j].tick_params(direction="in", top=True, right=True)
                if i == 1:
                    axis[i][j].set_xlabel(r"$R / \rm (km)$")

        fig.savefig(
            file_dir + file_name + "_" + spec, bbox_inches="tight", format="pdf"
        )

        """
        title_list = [r'$\rmPressure \ [GPa]$', r'$\rmMass \ [M_\oplus]$', 
                      r'$\rmDensity \ [10^3 \ kg/m^3]$', r'$\rmTemperature  \ [K]$',
                      r'$\rmGravity \ [m/s^2]$', r'$\rmv_{esc} \ [km/s]$']
        
        radius_list_dummy = np.array([[shell.radius/r_earth 
                                       for shell in layer.shells] 
                                        for layer in self.layers])
        
        
        #collect data for all layers
        plot_list_dummy = [[[shell.pres*1.0e-9 for shell in layer.shells] 
                                            for layer in self.layers ],
                    [[shell.mass/m_earth for shell in layer.shells] 
                                            for layer in self.layers ],
                    [[shell.dens*1.0e-3 for shell in layer.shells] 
                                            for layer in self.layers ],
                    [[shell.temp for shell in layer.shells] 
                                            for layer in self.layers ],
                    [[shell.gravity for shell in layer.shells]
                                            for layer in self.layers],
                    [[shell.v_esc*1.0e-3 for shell in layer.shells]
                                            for layer in self.layers]
                    ]
                    
        #collect thermal velocity data for water vapor to overplot it on
        #the graph for the escape velocity for comparison
        
        v_th_H2O_dummy = [[shell.v_th_H2O*1.0e-3 for shell in layer.shells]
                                            for layer in self.layers]
        
        #extract data of individual layers to one array for each parameter
        plot_list = []
        radius_list = []
        v_th_H2O_list = []
        
        for p in range(len(plot_list_dummy)):
            plot_list.append([])
            
        for p in range(len(plot_list_dummy)):
            for l in range(len(self.layers)):
                for i in range(len(plot_list_dummy[p][l])):
                    plot_list[p].append(plot_list_dummy[p][l][i])
            
        for l in range(len(self.layers)):
            for i in range(len(radius_list_dummy[l])):
                    radius_list.append(radius_list_dummy[l][i])
                    v_th_H2O_list.append(v_th_H2O_dummy[l][i])
                    
        ax = [axis[0][0], axis[0][1], axis[1][0], axis[1][1], axis[0][2], 
              axis[1][2]]
        
        ax[2].set_xlabel(r'$\rmRadius  \ [R_\oplus]$')
        ax[3].set_xlabel(r'$\rmRadius \ [R_\oplus]$')
        ax[5].set_xlabel(r'$\rmRadius \ [R_\oplus]$')
        #for axx in ax:
         #   axx.grid(which='both', axis='both')
        
        for i in range(len(ax)):
            if x_axis == 'radius':
                ax[i].plot(radius_list, plot_list[i], color=param_colors[i], zorder=3)

            elif x_axis == 'pres':
                ax[i].semilogx(plot_list[0], plot_list[i], color=param_colors[i], zorder=3)

            ax[i].tick_params(right=True, top=True, direction='in', which='both')
            ax[i].tick_params(which='major', axis='both', length=8)
            ax[i].grid(True, which='both', zorder=0, 
              color=plot_params['gridcolor'],
              alpha=plot_params['gridalpha'])
            
            ax[i].set_facecolor(plot_params['backgroundcol'])
            ax[i].grid(which='major', axis='both', linewidth=2, zorder=0)
            ax[i].xaxis.set_minor_locator(AutoMinorLocator())
            ax[i].yaxis.set_minor_locator(AutoMinorLocator())
            
            for j in range(len(radius_list)):
                if scatter:
                    ax[i].scatter(radius_list[j], plot_list[i][j])
                else:
                    ax[i].plot(radius_list[j], plot_list[i][j])
                  
            ax[i].set_ylabel(title_list[i])
        
        ax[-1].plot(radius_list, v_th_H2O_list, color=param_colors[-1],
          linestyle = '--')
        
        #save figure as image file if turned on
        if save:
            fig.savefig(path+filename+'.'+suffix)
            
        fig2, ax2 = plt.subplots()
        ax2.loglog(np.asarray(plot_list[3]), np.asarray(plot_list[0])*1.0e9)
        """

    def write(self, out="planet", loc="./"):
        """Creates data table of the planetary structure and writes it to ascii
        file. By default, the data is organized as follows:

        R           M           T       P       rho         v_esc       g
        (r_earth)   (m_earth)   (K)     (GPa)   (kg/m3)     (km/s)      (m/s2)
        -----------------------------------------------------------------------
        R_center    M_center    ...
        ...         ...         ...
        R_tot       M_tot       ...

        The meta data contains the properties that are stored in Planet.initials
        and Planet.finas by default
        """
        self.trim_profiles()

        # gather all planetary parameters as meta data
        meta = {"initials": self.initials, "finals": self.finals}

        names = [
            "R (r_earth)",
            "T (K)",
            "P (GPa)",
            "rho (kg/m3)",
            "M (m_earth)",
            "g (m/s2)",
        ]

        data_table = astropy.table.Table(self.profiles.T, meta=meta)
        ascii.write(
            data_table,
            loc + out + suffix["planet"],
            overwrite=True,
            format="ecsv",
            names=names,
        )

    def check_convergence(
        self,
        accs={
            "M_surface": 1.0e-4,
            "Mg_number": 1.0e-3,
            "T_surface": 1.0e-2,
            "ocean_frac": 1.0e-2,
            "x_H_core": 1.0e-2,
        },
        check_params=["M_surface", "Mg_number", "T_surface", "ocean_frac"],
    ):
        """Uses total mass, Mg#, ocean mass fraction, xi_H, and surface temperature
        to check if the planet has properly converged to the desired solution.
        """
        """
        self.reldev_M = (self.M_surface_is-self.M_surface_should)/self.M_surface_should
        self.reldev_Mg = (self.Mg_number_is-self.Mg_number_should)/self.Mg_number_should      
        self.reldev_T = (self.T_surface_is-self.T_surface_should)/self.T_surface_should
        self.reldev_ocean = (self.ocean_frac_is-self.ocean_frac_should)/self.ocean_frac_should
        
        try:
            self.reldev_H = (self.xi_H_core - self.xi_H_core_predicted)/self.xi_H_core_predicted
        except ZeroDivisionError:
            self.reldev_H = 0.
        
        #for very small oceans neglect uncertainties
        if self.ocean_frac_should < -3 and self.ocean_frac_is < -3:
            self.reldev_ocean = 0.
        
        #for very small Mg# neglect uncertainties
        if self.Mg_number_should < 1e-3 and self.Mg_number_is < 1e-3:
            self.reldev_Mg = 0e0
        
        reldevs = [self.reldev_M, 
                   self.reldev_Mg, 
                   self.reldev_T, 
                   self.reldev_ocean, 
                   self.reldev_H]
        """
        # print ('reldevs =', reldevs)
        self.converged = True
        for param in check_params:
            val_is = self.finals["{}_is".format(param)]
            val_should = self.initials["{}_should".format(param)]
            reldev = abs(val_should - val_is) / val_should
            # print ("{}: reldev".format(param), reldev)
            # For very small ocean mass fractions neglect uncertainty
            if param == "ocean_frac":
                if val_is < -2 and val_should < -2:
                    reldev = 0.0

            # For very small magnesium numbers neglect uncertainties
            if param == "Mg_number":
                if val_is < 1e-3 and val_should < 1e-3:
                    reldev = 0.0

            if abs(reldev) > accs[param]:
                self.converged = False

    def correct_density_profiles(self):
        """
        If a layer is skipped it is possible that the ghost shell contains
        spurious values for the density that are lower than the density of the
        subsequent shells. In this case, ommit the ghost shell by overwriting it
        with the parameter values of the neighbourhing shell.
        """

        i = 0

        df = self.data_table.to_pandas()

        N_shells = len(df["rho (kg/m3)"])

        while True:
            if i == N_shells - 2:
                break

            else:
                d = df["rho (kg/m3)"][i]
                dd = df["rho (kg/m3)"][i + 1]

                if d < dd:
                    self.data_table[i][3] = dd

            i += 1

    def load(self, loc="./", file_name="planet"):
        """ """
        # read ascii file for eos table
        self.data_table = ascii.read(loc + file_name + suffix["planet"], format="ecsv")

        self.initials = self.data_table.meta["initials"]
        self.finals = self.data_table.meta["finals"]
