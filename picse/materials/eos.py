#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 14:26:56 2019

@author: oshah
"""
import time
import sys
from picse.utils.function_tools import PIMPtools
from picse.materials import eosTablesUse as eosTab
import numpy as np
from matplotlib import pyplot as plt

# from picse.materials import Mazevet2018EOS as mazevet
import matplotlib.transforms as transforms
from picse.materials.equations_of_state import Feistel2006EOS as feistel

# from picse.materials import eoswater
from picse.materials import phaseCheck

# from picse.materials import hydration
# from picse.materials import brucite_phase as bruce
from picse.utils.function_tools import functionTools as ftool
from picse.physicalparams import (
    T0_list,
    K0_list,
    K0prime_list,
    rho0_list,
    aP_list,
    molar_mass_list,
    Rgas,
    mH,
    kB,
    EOS_type,
    material_list,
    EOSstring_list,
    NA,
    gas_molar_mass_list,
    gas_EOS_type,
    a_VDW_list,
    b_VDW_list,
    K0prime_list_Mg2SiO4,
    phase_list_Mg2SiO4,
    material_plot_list,
    K0_list_Brucite,
    K0prime_list_Brucite,
    rho0_list_Brucite,
    phase_list_Brucite,
    q_list,
    n_list,
    thetaD0_list,
    gamma0_list,
    dq_pv,
    dK0_pv,
    dK0prime_pv,
    dgamma0_pv,
    dthetaD0_pv,
    dd0_pv,
    dq_ppv,
    dK0_ppv,
    dK0prime_ppv,
    dgamma0_ppv,
    dthetaD0_ppv,
    dd0_ppv,
    aT_list,
    bT_list,
    cT_list,
    aP_list_Mg2SiO4,
    aT_list_Mg2SiO4,
    bT_list_Mg2SiO4,
    cT_list_Mg2SiO4,
    daTdH2O_list_Mg2SiO4,
    dbTdH2O_list_Mg2SiO4,
    EOSstring_list,
)

from picse.runparams import color_list, grid_color, pressure_switch_H2O
from matplotlib.ticker import MultipleLocator
from matplotlib import cm
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import ImageGrid, AxesGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.mpl_axes import Axes

# Mixing rule for material mixtures
def mix1(densities, fractions, **kwargs):
    return 1.0 / sum([fractions[i] / densities[i] for i in range(len(fractions))])


# Birch-Murnaghan (BM)
def P_BM(eta=0, K=0, K0prime=0, **kwargs):
    return (
        3.0
        / 2.0
        * K
        * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
        * (1.0 + 3.0 / 4.0 * (K0prime - 4.0) * (eta ** (2.0 / 3.0) - 1.0))
    )


def rho_BM(
    ll=0,
    d0=0,
    T=0,
    P=0,
    dmin=None,
    dmax=None,
    T0=300,
    aT=None,
    bT=None,
    cT=None,
    **kwargs,
):
    d0 = rho0_list[ll]
    K0 = K0_list[ll]
    d0T = PIMPtools.rho0T(T, T0, d0, ll, aT, bT, cT)
    # linear slope of K0(T-T0)
    aP = aP_list[ll]

    # ambient temperature
    T0 = T0_list[ll]
    K = K0 + aP * (T - T0)
    K0prime = K0prime_list[ll]
    d = d0T * ftool.bisec(
        f=P_BM,
        whicharg="eta",
        K=K,
        K0prime=K0prime,
        a=dmin / d0T,
        b=dmax / d0T,
        y=P,
        identity="rho_BM()",
        **kwargs,
    )

    return d


def dPdrho_BM(ll=0, eta=0, K=0, K0prime=0, d0T=0, **kwargs):
    return (
        3.0
        / 2.0
        * K
        / d0T
        * (
            (7.0 / 3.0 * eta ** (4.0 / 3.0) - 5.0 / 3.0 * eta ** (2.0 / 3.0))
            * (1.0 - 3.0 / 4.0 * (4.0 - K0prime) * (eta ** (2.0 / 3.0) - 1.0))
            - 1.0
            / 2.0
            * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
            * (4.0 - K0prime)
            * eta ** (-1.0 / 3.0)
        )
    )


# Mie-Grüneisen-Debye (MGD)
def dP_MGD(
    eta=0,
    ll=5,
    T=300,
    gamma0=None,
    thetaD0=None,
    rho0=None,
    q=None,
    molmass=None,
    nparticles=None,
    **kwargs,
):
    d = eta * rho0
    tempType = ftool.checkKey(name="tempType", **kwargs)

    if tempType == "isothermal":
        deltaP = 0.0

    else:
        E1 = PIMPtools.Eth(d, T, thetaD0, rho0, gamma0, q, molmass, nparticles, ll)
        E2 = PIMPtools.Eth(d, 300.0, thetaD0, rho0, gamma0, q, molmass, nparticles, ll)
        deltaP = (E1 - E2) * PIMPtools.gamma(q, gamma0, rho0, d)

    return deltaP


def P_MGD(
    eta=0.0,
    K0=0.0,
    K0prime=0.0,
    ll=5,
    T=300,
    gamma0=None,
    thetaD0=None,
    rho0=None,
    q=None,
    molmass=None,
    nparticles=None,
    **kwargs,
):
    d = eta * rho0
    tempType = ftool.checkKey(name="tempType", **kwargs)

    if tempType == "isothermal":
        deltaP = 0.0

    else:
        # call routine to compute thermal energy contribution
        E1 = PIMPtools.Eth(d, T, thetaD0, rho0, gamma0, q, molmass, nparticles, ll)
        E2 = PIMPtools.Eth(d, 300.0, thetaD0, rho0, gamma0, q, molmass, nparticles, ll)
        deltaP = (E1 - E2) * PIMPtools.gamma(q, gamma0, rho0, d)

    P = (
        3.0
        / 2.0
        * K0
        * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
        * (1.0 + 3.0 / 4.0 * (K0prime - 4.0) * (eta ** (2.0 / 3.0) - 1.0))
    )
    # print("P =", P * 1e-5)
    # print("delta P =", deltaP * 1e-5)
    return P + deltaP


def rho_MGD(
    eta=0.0,
    K0=0,
    K0prime=0.0,
    ll=5,
    T=300,
    gamma0=None,
    thetaD0=None,
    rho0=None,
    q=None,
    molmass=None,
    P=0.0,
    nparticles=None,
    dmin=None,
    dmax=None,
    **kwargs,
):

    rho = rho0 * ftool.bisec(
        f=P_MGD,
        whicharg="eta",
        K0=K0,
        K0prime=K0prime,
        ll=ll,
        T=T,
        y=P,
        a=dmin / rho0,
        b=dmax / rho0,
        q=q,
        molmass=molmass,
        nparticles=nparticles,
        gamma0=gamma0,
        rho0=rho0,
        thetaD0=thetaD0,
        **kwargs,
    )
    return rho


def dPdrho_MGD(
    eta=0,
    K0=0,
    K0prime=0,
    ll=0,
    T=300,
    gamma0=None,
    thetaD0=None,
    rho0=None,
    q=None,
    molmass=None,
    nparticles=None,
    **kwargs,
):
    d0 = rho0

    ddeltaPth = ftool.deriv(
        f=dP_MGD,
        x0=eta,
        whicharg="eta",
        T=T,
        gamma0=gamma0,
        thetaD0=thetaD0,
        rho0=rho0,
        q=q,
        molmass=molmass,
        nparticles=nparticles,
        **kwargs,
    )
    return (
        3.0
        / 2.0
        * K0
        / d0
        * (
            (7.0 / 3.0 * eta ** (4.0 / 3.0) - 5.0 / 3.0 * eta ** (2.0 / 3.0))
            * (1.0 - 3.0 / 4.0 * (4.0 - K0prime) * (eta ** (2.0 / 3.0) - 1.0))
            - 1.0
            / 2.0
            * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
            * (4.0 - K0prime)
            * eta ** (-1.0 / 3.0)
        )
        + ddeltaPth
    )


# Vinet (Vinet)
def P_Vinet(eta=0, K0=0, K0prime=0):
    return (
        3.0
        * K0
        * eta ** (2.0 / 3.0)
        * (1.0 - eta ** (-1.0 / 3.0))
        * np.exp(3.0 / 2.0 * (K0prime - 1.0) * (1.0 - eta ** (-1.0 / 3.0)))
    )


def dPdrho_Vinet(eta=0, K0=0, K0prime=0, ll=0):
    d0 = rho0_list[ll]
    return (
        K0
        / d0
        * np.exp(3.0 / 2.0 * (K0prime - 1.0) * (1.0 - eta ** (-1.0 / 3.0)))
        * (
            2.0 * eta ** (-1.0 / 3.0) * (1.0 - eta ** (-1.0 / 3.0))
            + eta ** (-2.0 / 3.0)
            * 3.0
            / 2.0
            * eta ** (-2.0 / 3.0)
            * (1.0 - eta ** (-1.0 / 3.0))
            * (K0prime - 1.0)
        )
    )


# Belonoshko (Bel)
def P_Bel(
    eta=0,
    K0=0,
    K0prime=0,
    ll=0,
    T=300,
    q=None,
    gamma0=None,
    rho0=None,
    molmass=None,
    **kwargs,
):
    d = eta * rho0
    tempType = ftool.checkKey(name="tempType", **kwargs)

    if tempType == "isothermal":
        deltaP = 0.0
    else:
        deltaP = (
            3.0
            * Rgas
            * PIMPtools.gamma(q, gamma0, rho0, d)
            * (T - T0_list[ll])
            * (d / molmass)
        )

    return (
        3.0
        / 2.0
        * K0
        * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
        * (1.0 + 3.0 / 4.0 * (K0prime - 4.0) * (eta ** (2.0 / 3.0) - 1.0))
        + deltaP
    )


def rho_Bel(
    rho0=None,
    T=None,
    P=None,
    K0=None,
    K0prime=None,
    ll=9,
    q=None,
    molmass=None,
    gamma0=None,
    dmin=None,
    dmax=None,
    **kwargs,
):
    rho = rho0 * ftool.bisec(
        f=P_Bel,
        whicharg="eta",
        K0=K0,
        K0prime=K0prime,
        a=dmin / rho0,
        b=dmax / rho0,
        y=P,
        ll=ll,
        T=T,
        q=q,
        molmass=molmass,
        gamma0=gamma0,
        rho0=rho0,
        **kwargs,
    )

    return rho


def dPdrho_Bel(
    eta=0,
    K0=0,
    K0prime=0,
    ll=0,
    T=300,
    q=None,
    gamma0=None,
    rho0=None,
    molmass=None,
    **kwargs,
):
    d0 = rho0
    d = eta * d0
    tempType = ftool.checkKey(name="tempType", **kwargs)
    if tempType == "isothermal":
        delta = 0.0

    else:
        delta = (
            3.0
            * Rgas
            * PIMPtools.gamma(q, gamma0, rho0, d)
            * (T - T0_list[ll])
            / molmass
        )

    return (
        3.0
        / 2.0
        * K0
        / d0
        * (
            (7.0 / 3.0 * eta ** (4.0 / 3.0) - 5.0 / 3.0 * eta ** (2.0 / 3.0))
            * (1.0 - 3.0 / 4.0 * (4.0 - K0prime) * (eta ** (2.0 / 3.0) - 1.0))
            - 1.0
            / 2.0
            * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
            * 0.5
            * (4.0 - K0prime)
            * eta ** (-1.0 / 3.0)
        )
        + delta
    )


# Vinet-Rydberg (VR)
def P_VR(
    eta=0.0,
    K0=0.0,
    K0prime=0.0,
    ll=5,
    T=300,
    gamma0=None,
    theta0=None,
    rho0=None,
    q=None,
    molmass=None,
    **kwargs,
):

    T = T - T0_list[ll]

    d = eta * rho0
    # print ('d =', d)
    P0 = (
        3.0
        * K0
        * (eta ** (2.0 / 3.0) - eta ** (1.0 / 3.0))
        * np.exp((3.0 / 2.0 * K0prime - 3.0 / 2.0) * (1.0 - eta ** (-1.3)))
    )

    theta = theta0 * np.exp(gamma0 / q * (1.0 - eta ** (-q)))

    print("theta =", theta)
    print("P0 =", P0)
    Pth = 3.0 * Rgas * (theta / (np.exp(theta / T) - 1.0))
    # print ('gamma =', PIMPtools.gamma(q, gamma0, rho0, d))
    Pth *= PIMPtools.gamma(q, gamma0, rho0, d) * n_list[ll] * d / molmass
    print("Pth =", Pth)

    g = -0.884
    e0 = 198e-6
    e = e0 * eta ** (-g)

    Pe = 3.0 / 2.0 * n_list[ll] * Rgas * e * g * T ** 2 * d / molmass
    print("Pe =", Pe)
    return P0 + Pth + Pe


def rho_VR(
    rho0=None,
    T=None,
    P=None,
    K0=None,
    K0prime=None,
    ll=9,
    q=None,
    molmass=None,
    gamma0=None,
    dmin=None,
    dmax=None,
    **kwargs,
):
    rho = rho0 * ftool.bisec(
        f=P_VR,
        whicharg="eta",
        K0=K0,
        K0prime=K0prime,
        a=dmin / rho0,
        b=dmax / rho0,
        y=P,
        ll=ll,
        T=T,
        q=q,
        molmass=molmass,
        gamma0=gamma0,
        rho0=rho0,
        **kwargs,
    )

    return rho


def dPdrho_VR(
    eta=0.0,
    K0=0.0,
    K0prime=0.0,
    ll=5,
    T=300,
    gamma0=None,
    theta0=None,
    rho0=None,
    q=None,
    molmass=None,
    **kwargs,
):

    e = np.exp((3.0 / 2.0 * K0prime - 3.0 / 2.0) * (1.0 - eta ** (-1.0 / 3.0)))
    dP0drho = (
        3.0 * K0 * (2.0 / 3.0 * eta ** (-1.0 / 3.0) - 1.0 / 3.0 * eta ** (-2.0 / 3.0))
    )
    dP0drho *= e
    dP0drho += (
        3.0
        * K0
        * (eta ** (2 / 3.0) - eta ** (-1.0 / 3.0))
        * (-1.0 / 3.0 * eta ** (-4.0 / 3.0))
        * (3 / 2.0 * K0prime - 3.0 / 2.0)
        * e
    )
    dP0drho = dP0drho / rho0

    dPthdrho = (
        3.0 * n_list[ll] * Rgas * molmass * gamma0 * theta0 / (np.exp(theta0 / T) - 1.0)
    )
    dPthdrho *= eta ** (-q) / rho0 * (1.0 - q)
    dPthdrho = dPthdrho / rho0

    return dP0drho + dPthdrho


# Ichikawa
def P_Ichi(
    eta=0.0,
    K0=0.0,
    K0prime=0.0,
    ll=5,
    T=300,
    gamma0=None,
    thetaD0=None,
    rho0=None,
    q=None,
    molmass=None,
    **kwargs,
):

    d = eta * rho0

    # print ('eta =', eta)
    P0 = (
        3.0
        * K0
        * (eta ** (2.0 / 3.0) - eta ** (1.0 / 3.0))
        * np.exp((3.0 / 2.0 * K0prime - 3.0 / 2.0) * (1.0 - eta ** (-1.0 / 3.0)))
    )

    # print (np.exp((3./2.*K0prime-3./2.)*(1.-eta**(-1./3.))))
    # print ('K0 =', K0)

    a = 1.0
    b = 0.35
    e0 = 0.314e-4
    g = -0.4

    gamma = gamma0 * (1.0 + a * (eta ** (-b) - 1.0))
    Eth = 3.0 * n_list[ll] * Rgas * (T + e0 * eta ** (-g) * T ** 2)
    Eth0 = 3.0 * n_list[ll] * Rgas * (T0_list[ll] + e0 * eta ** (-g) * T0_list[ll] ** 2)

    V = molmass / d / eta

    Pth = gamma / V * (Eth - Eth0)
    # print ('P0 =', P0*1e-9)
    # print ('Pth =', Pth*1e-9)
    return P0 + Pth


def rho_Ichi(
    rho0=None,
    T=None,
    P=None,
    K0=None,
    K0prime=None,
    ll=9,
    q=None,
    molmass=None,
    gamma0=None,
    dmin=None,
    dmax=None,
    **kwargs,
):

    rho = rho0 * ftool.bisec(
        f=P_Ichi,
        whicharg="eta",
        K0=K0,
        K0prime=K0prime,
        a=dmin / rho0,
        b=dmax / rho0,
        y=P,
        ll=ll,
        T=T,
        q=q,
        molmass=molmass,
        gamma0=gamma0,
        rho0=rho0,
        **kwargs,
    )

    return rho


def dPdrho_Ichi(
    eta=0.0,
    K0=0.0,
    K0prime=0.0,
    ll=5,
    T=300,
    gamma0=None,
    theta0=None,
    rho0=None,
    q=None,
    molmass=None,
    **kwargs,
):

    e = np.exp((3.0 / 2.0 * K0prime - 3.0 / 2.0) * (1.0 - eta ** (-1.0 / 3.0)))
    dP0drho = (
        3.0 * K0 * (2.0 / 3.0 * eta ** (-1.0 / 3.0) - 1.0 / 3.0 * eta ** (-2.0 / 3.0))
    )
    dP0drho *= e
    dP0drho += (
        3.0
        * K0
        * (eta ** (2 / 3.0) - eta ** (-1.0 / 3.0))
        * (-1.0 / 3.0 * eta ** (-4.0 / 3.0))
        * (3 / 2.0 * K0prime - 3.0 / 2.0)
        * e
    )
    dP0drho = dP0drho / rho0

    dPthdrho = (
        3.0 * n_list[ll] * Rgas * molmass * gamma0 * theta0 / (np.exp(theta0 / T) - 1.0)
    )
    dPthdrho *= eta ** (-q) / rho0 * (1.0 - q)
    dPthdrho = dPthdrho / rho0

    return dP0drho + dPthdrho


# van der Waal equation for fluids
def P_VDW(d=0.0, T=0.0, ll=None, **kwargs):
    mm = gas_molar_mass_list[ll]
    a = a_VDW_list[ll]
    b = b_VDW_list[ll]
    c1 = d * Rgas * T / (mm - d * b)
    c2 = d ** 2 * a / mm ** 2
    print("c1 =", c1)
    print("c2 =", c2)
    return c1 - c2


def rho_VDW(P=0.0, T=0.0, ll=None, **kwargs):
    return ftool.bisec(
        f=P_VDW, whicharg="d", T=T, a=1.0e-8, b=10.0, ll=ll, y=P, **kwargs
    )


def dPdrho_VDW(d=0.0, T=0.0, ll=None, **kwargs):
    mm = gas_molar_mass_list[ll]
    a = a_VDW_list[ll]
    b = b_VDW_list[ll]
    return (
        Rgas * T / (mm - b * d)
        - 2.0 * d * a / mm ** 2
        + b * d * Rgas * T / (mm - b * d) ** 2
    )


def dTdP_VDW(d=0.0, T=0.0, ll=None, **kwargs):
    mm = gas_molar_mass_list[ll]
    b = b_VDW_list[ll]
    return (mm - b * d) / (d * Rgas)


# ideal gas law (IGL)
def P_IGL(d=0, T=300, ll=None):
    mm = gas_molar_mass_list[ll]
    return d * kB * T * NA / mm


def rho_IGL(P=0, T=300, ll=None):
    mm = gas_molar_mass_list[ll]
    return P * mm / NA / kB / T


def dPdrho_IGL(T=300, ll=0):
    mm = gas_molar_mass_list[ll]
    return kB * T * NA / mm


def dTdP_IGL(d=0.0, ll=0):
    mm = gas_molar_mass_list[ll]
    return mm / (d * kB * NA)


# IAPWS + Mazevet for pure water
def rho_water_eos(P=None, T=None, **kwargs):
    raise NotImplementedError("This feature is not available yet.")

    dens = eoswater.density(P=P, T=T)
    return dens


def dPdrho_water(P=None, T=None, d=None, **kwargs):
    """This computes dPdrho_S for water"""
    raise NotImplementedError("This feature is not available yet.")

    dPdrho = eoswater.dPdrho_T(P=P, d=d, T=T)
    return dPdrho


def rho0_enstatite(Fe_number=None):
    """Fe_number is in mole fraction (NOT mol%)
    Fe_number = 1 -> pure Fe2Si2O6
    Fe_number = 0 -> pure Mg2Si2O6
    """
    delta_rho0_enstatite = rho0_list[4] - rho0_list[2]

    return rho0_list[2] + Fe_number * delta_rho0_enstatite


def rho0_perovskite(Fe_number=None):
    """Fe_number is in mole fraction (NOT mol%)
    Fe_number = 1 -> pure FeSiO3
    Fe_number = 0 -> pure MgSiO3
    """
    delta_rho0_perovskite = rho0_list[8] - rho0_list[6]

    return rho0_list[6] + Fe_number * delta_rho0_perovskite


def rho0_wuestite(Fe_number=None):
    """Fe_number is in mole fraction (NOT mol%)
    Fe_number = 1 -> pure FeSiO3
    Fe_number = 0 -> pure MgSiO3
    """
    delta_rho0_wuestite = rho0_list[7] - rho0_list[5]

    return rho0_list[5] + Fe_number * delta_rho0_wuestite


# def rho0_silicates(
#     SiMg=None, Fe_number=None, X_H2O=0.0, T=None, P=None, phase=None, lay=1
# ):
#     """SiMg and Fe_number are in mole fraction (NOT mol%) but X_H2O is in wt%"""
#     # start form rho0 for hydrated Olivine
#     rho0_silicates = hydration.rho0(X_H2O, 100 * Fe_number, phase)

#     rho0_En = rho0_enstatite(Fe_number)
#     rho0_Ol = hydration.rho0(0.0, 100 * Fe_number, phase)

#     delta = rho0_En - rho0_Ol

#     return rho0_silicates + 2 * delta * (SiMg - 0.5)


def K0_wuestite(Fe_number=None, phase=None):
    """Fe_number is in mole fraction (NOT mol%)
    Fe_number = 1 -> pure FeSiO3
    Fe_number = 0 -> pure MgSiO3
    """
    delta_K0_wuestite = K0_list_Brucite[phase] - K0_list[7]

    return K0_list_Brucite[phase] - Fe_number * delta_K0_wuestite


def K0_silicates(
    SiMg=None, Fe_number=None, X_H2O=0.0, T=None, P=None, phase=None, lay=1
):
    """SiMg and Fe_number are in mole fraction (NOT mol%) but X_H2O is in wt%"""
    if lay == 1:
        # start from rho0 for hydrated Wuestite
        K0_silicates = K0_list_Brucite[phase]

        # for perovskite no data on how Fe changes K0 or K0' are used
        K0_Perov = K0_list[6]

        # for wuestite K0 and K0' change if iron is added
        K0_Ws = K0_wuestite(Fe_number=Fe_number, phase=phase)

        delta = K0_Perov - K0_Ws

        return K0_silicates + delta * SiMg

    elif lay == 2:
        # start form rho0 for hydrated Olivine
        K0_silicates = hydration.K0(X_H2O, 100 * Fe_number, phase)

        # for enstatite no data on how Fe changes K0 or K0' are used
        K0_En = K0_list[2]
        K0_Ol = hydration.K0(0.0, 100 * Fe_number, phase)

        delta = K0_En - K0_Ol

        return K0_silicates + 2 * delta * (SiMg - 1 / 2)


def K0prime_silicates(
    SiMg=None, Fe_number=None, X_H2O=0.0, T=None, P=None, phase=None, lay=1
):
    """SiMg and Fe_number are in mole fraction (NOT mol%) but X_H2O is in wt%"""
    if lay == 1:
        # start from rho0 for hydrated Periclase = Brucite
        K0_prime_silicates = K0prime_list_Brucite[phase]

        # for perovskite no data on how Fe changes K0 and K0' are used
        K0_prime_Perov = K0prime_list[6]

        ##for wuestite no data on how Fe changes K0 and K0' are used
        K0_prime_Ws = K0prime_list_Brucite[phase]

        delta = K0_prime_Perov - K0_prime_Ws

        return K0_prime_silicates + delta * SiMg

    elif lay == 2:
        # start from rho0 for hydrated Olivine
        K0_prime_silicates = K0prime_list[1]

        K0_prime_En = K0prime_list[2]
        K0_prime_Ol = K0prime_list[1]

        delta = K0_prime_En - K0_prime_Ol

        return K0_prime_silicates + 2 * delta * (SiMg - 0.5)


def compute(
    ll=0,
    what="all",
    P=None,
    T=None,
    d=None,
    dmin=0.001,
    dmax=3e4,
    Tmin=10.0,
    Tmax=1.0e5,
    whichEOS=None,
    Fe_number=0.0,
    X_H2O=0.0,
    saturation=False,
    phase=None,
    table=False,
    **kwargs,
):
    """
    Employs the EOS for the given material and computes the specified
    quantity at a given conditions.

    Parameters:

    ll (int, optional): The material signature. Defaults to 0.
    P (float, optional): The pressure in Pa. Defaults to None.
    T (float, optional): The temperature in K. Defaults to None.
    d (float, optional): The density in kg/m³. Defaults to None.

    Outputs:
        Tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float], Optional[int], Optional[float], Optional[float], Optional[float]]:
        - Element 1: Density. Defaults to None.
        - Element 2: Temperature. Defaults to None.
        - Element 3: Pressure. Defaults to None.
        - Element 4: Pressure derivative of temperature at constant entropy. Defaults to None.
        - Element 5: Density derivative of pressure at constant temperature. Defaults to None.
        - Element 6: Phase tag. Defaults to None.
        - Element 7: Water equivalent content for hydrated or hydrogenated compounds. Defaults to None.
        - Element 8: Thermal expansion. Defaults to None.
        - Element 9: Al content for silicates. Defaults to None.
    """

    Fe_number = min(Fe_number, 100.0 - 1.0e-10)
    FeMg = Fe_number / (100.0 - Fe_number)

    if ll == 9 and phase == 1:
        ll = 19

    xi_Al = 0.0

    # Initiate pressure derivative with respect to density
    dPdrho = None

    # Initiate shear modulus
    G_shear = None

    # Initiate thermal expansion coefficient
    alpha = None

    # initial adiabatic gradient
    dTdP_S = None

    if not table:
        stopper = False

        try:

            K0 = K0_list[ll]
            K0prime = K0prime_list[ll]
            d0 = rho0_list[ll]
            aT = aT_list[ll]
            bT = bT_list[ll]
            cT = cT_list[ll]
            aP = aP_list[ll]
            T0 = T0_list[ll]
            gamma0 = gamma0_list[ll]
            thetaD0 = thetaD0_list[ll]
            q = q_list[ll]
            nparticles = n_list[ll]
            molmass = molar_mass_list[ll]
            # linear slope of K0(T-T0)
            aP = aP_list[ll]

            # ambient temperature
            T0 = T0_list[ll]

            # for Mg2SiO4 (Olivine) different phases are taken into account
            # hydration effects are taken into account
            # effect of iron content are taken into account
            if ll == 1 or ll == 12 or ll == 13 or ll == 14:
                # determine phase but only if no phase has been specified
                # and a pressure is given. If a phase is passed as argument
                # the Compute() methode will NOT check if the phase is correct
                # at the given TP, but just use the EOS for the given phase and
                # therefore possibly extrapolate into a different phase
                if not P == None and phase == None:
                    phase = hydration.phaseCheck(T=T, P=P)

                else:

                    # different Olivine phases
                    if ll == 12:
                        phase = 0

                    elif ll == 13:
                        phase = 1

                    elif ll == 14:
                        phase = 2

                # compute saturated water content saturation=True
                # this overwrites any value of X_H2O that is passed to the function
                # Note that the value is computed in wt%
                if saturation:

                    # The EoS of Ol missbehaves at low/high T and high P in the
                    # hydrated case and is hence truncated here.
                    T = max(T, 300.0)
                    T = min(T, 5000.0)

                    # P = min(2.5e10, P)

                    X_H2O = hydration.X_H2O_sat(
                        P=P, T=T, Fe_number=Fe_number, phase=phase
                    )
                    # print("X_H2O =", X_H2O)

                # extract parameters for the BM eos for the phase at hand
                # compute K0 using linear dependence on water content in wt% and
                # iron number in mol%
                K0 = hydration.K0(X_H2O, Fe_number, phase)

                # pressure derivatives of isothermal bulk moduli are assumed to
                # be constant due to lack of reliable experimental constrains
                K0prime = K0prime_list_Mg2SiO4[phase]
                # ambient density scales linearly with iron and water content
                d0 = hydration.rho0(X_H2O, Fe_number, phase)

                aP = aP_list_Mg2SiO4[phase]
                aT = max(
                    aT_list_Mg2SiO4[phase] + X_H2O * daTdH2O_list_Mg2SiO4[phase], 1.0e-7
                )
                bT = max(
                    bT_list_Mg2SiO4[phase] + X_H2O * dbTdH2O_list_Mg2SiO4[phase],
                    1.0e-10,
                )
                cT = cT_list_Mg2SiO4[phase]

            # for Mg(OH)2 (Brucite) different phases are taken into account
            elif ll == 11:
                raise NotImplementedError(
                    "Mg(OH)2 EoS is currently not available here."
                )
                stopper = False

                if phase == None:
                    phase = bruce.phaseCheck(T=T, P=P)

                if phase > 0:
                    ll = 5
                    stopper = False

                else:
                    X_H2O = 30.1

                K0 = K0_list_Brucite[phase]
                K0prime = K0prime_list_Brucite[phase]
                d0 = rho0_list_Brucite[phase]
                # return d, phase, X_H2O

            elif ll == 5:
                stopper = False

                if phase == None:
                    raise NotImplementedError(
                        "Mg(OH)2 EoS is currently not available here. Specify the phase argument to 0 to avoid this error."
                    )
                    phase = bruce.phaseCheck(T=T, P=P)

                # phase 0 is Mg(OH) which does not exist in the anhydrous case
                if not saturation:
                    if phase == 0:
                        phase = 2

                K0 = K0_list_Brucite[phase]
                K0prime = K0prime_list_Brucite[phase]
                d0 = rho0_list_Brucite[phase]
                """
                #check if dissociated phase is present. If Mg(OH)2 -> MgO + H2O,
                #take the pure MgO eos as water is assumed to be transported up
                if not P == None and phase == None:
                    phase = bruce.phaseCheck(T=T, P=P)       
                    if phase == 2:
                        ll = 5
                    
                K0 = K0_list[ll]
                K0prime = K0prime_list[ll]
                d0 = rho0_list[ll]
                """

                """
                elif ll == 9:
                    stopper = True
                    if what == 'dens' or what != 'dens':
                        d, phase = tabIron.interpolate(order=1, x=T, y=P, param=2)
                        return d, phase, X_H2O
                """

            elif ll == 19 or ll == 9:
                Fe_number = min(Fe_number, 99.9999999)

                # Fits to Thompson et al 2018 for hydrogen content
                d0 += -1.62e3 * (100.0 - Fe_number) / Fe_number
                K0 += -40.2e9 * (100.0 - Fe_number) / Fe_number

            # Adjust EOS parameters for perovskite according to Fe content
            # The adopted values are taken for 13 mol% Fe and are here
            # scaled to the specified iron content (Sun 2018)
            elif ll == 6:
                if phase == None:
                    phase = phaseCheck.getPhase(ll=ll, P=P, T=T)

                if phase == 0:
                    q += Fe_number * dq_pv
                    q = max(q, 0.0)
                    thetaD0 += Fe_number * dthetaD0_pv
                    gamma0 += Fe_number * dgamma0_pv

                    d0 += Fe_number * dd0_pv
                    K0 += Fe_number * dK0_pv
                    K0prime += Fe_number * dK0prime_pv

                    # Adopt fixed saturation value of 0.1 wt% H2O in Pv and
                    # neglect effects on EoS parameters as they are very small
                    # (e.g. Jacobsen 2010, Ohtani 2015)
                    if saturation:
                        X_H2O = 0.1

                # Adjust EOS parameters for post-perovskite according to Fe content
                # and water content
                # Hydration is taken from Townsend 2015
                # Fe and Fe free is taken from Sun 2018
                if phase == 1:
                    """
                    ll=15
                    K0 = K0_list[ll]
                    K0prime = K0prime_list[ll]
                    d0 = rho0_list[ll]
                    #aT = aT_list[ll]
                    #bT = bT_list[ll]
                    #cT = cT_list[ll]
                    #aP = aP_list[ll]
                    T0 = T0_list[ll]
                    gamma0 = gamma0_list[ll]
                    thetaD0 = thetaD0_list[ll]
                    q = q_list[ll]
                    T0 = T0_list[ll]
                    nparticles = n_list[ll]
                    molmass = molar_mass_list[ll]
                    """
                    # linear slope of K0(T-T0)
                    # aP = aP_list[ll]
                    K0 = K0_list[15]
                    d0 = rho0_list[15]
                    q = q_list[15]
                    thetaD0 = thetaD0_list[15]
                    gamma0 = gamma0_list[15]

                    q += Fe_number * dq_pv
                    q = max(q, 0.0)
                    thetaD0 += Fe_number * dthetaD0_pv
                    gamma0 += Fe_number * dgamma0_pv

                    d0 += Fe_number * dd0_pv
                    K0 += Fe_number * dK0_pv
                    K0prime += Fe_number * dK0prime_pv

                    """
                    q += Fe_number*dq_ppv
                    q = max(q, 0.0)
                    thetaD0 += Fe_number*dthetaD0_ppv
                    gamma0 += Fe_number*dgamma0_ppv
                    d0 += dd0_ppv
                    K0 += Fe_number*dK0_ppv
                    """

                    # Effect of temperature on derivative of Bulk modulus
                    # Fit to Lin 2014 Fig 4
                    # K0prime +=  2.27e-4*(T-300)

                    # acconut for hydration
                    if saturation:
                        X_H2O = 3.0  # in wt % taken from Townsend 2015
                        K0 -= (9.6 - 0.0183 * P * 1.0e-9) * X_H2O
                        d0 -= 0.0037 * d0 * X_H2O

            elif ll == 16 or ll == 17:
                if saturation:
                    X_H2O = hydration.XH2O_Stv(T=T, P=P, AlSi=0.0) * 100.0
                    d0 -= d0 * 0.0051 * X_H2O

            # K0 = K0_list[ll]
            # K0prime = K0prime_list[ll]
            # d0 = rho0_list[ll]
            # temperature dependence of zero pressure density
            d0T = PIMPtools.rho0T(T, T0, d0, ll, aT, bT, cT)

        # for ideal gases no such parameters are specified
        except IndexError:
            pass

        # if no EOS type is specified use the one predefined in physicalparams.py
        # for the given material
        if whichEOS == None:
            whichEOS = EOS_type[ll]

        # Birch-Murnaghan (BM3)
        if whichEOS == 0 and not stopper:
            try:
                K = K0 + aP * (T - T0)

            except TypeError:
                K = K0

            if what == "pres":
                eta = d / d0T
                P = P_BM(eta=eta, K=K, K0prime=K0prime, **kwargs)

            elif what == "dens":
                d = rho_BM(ll=ll, P=P, dmin=dmin, dmax=dmax, T=T)

            elif what == "dPdrho":
                eta = d / d0T
                dPdrho = dPdrho_BM(
                    eta=eta, K=K, K0prime=K0prime, ll=ll, d0T=d0T, **kwargs
                )

            elif what == "all":
                # Compute density
                d = d0T * ftool.bisec(
                    f=P_BM,
                    whicharg="eta",
                    K=K,
                    K0prime=K0prime,
                    a=dmin / d0T,
                    b=dmax / d0T,
                    y=P,
                    identity="eos.compute()",
                    acc=1.0e-5,
                    **kwargs,
                )

                # Compute relative density change if hydrated and adjust
                # unhydrous density accordingly
                d = d * hydration.delta_rho(P=P, T=T, ll=ll, X_H2O=X_H2O)

                eta = d / d0T

                # Compute dPdrho at given density and temperature
                dPdrho = dPdrho_BM(
                    eta=eta, K=K, K0prime=K0prime, ll=ll, d0T=d0T, **kwargs
                )

                # Compute thermal expansion coefficient
                # Note, for some reason acc < 1.0e-3 leads to numerical
                # issues in the derivative. But alpha does not need to be
                # calculated very accurately anyways as it has negligible effect
                # on the adiabatic gradient and the relative error induced by the
                # physical error on P(T,rho) is probably larger than 1.0e-3
                drhodT = ftool.deriv(
                    f=rho_BM,
                    whicharg="T",
                    P=P,
                    dmax=dmax,
                    dmin=dmin,
                    ll=ll,
                    x0=T,
                    T0=T0,
                    aT=aT,
                    bT=bT,
                    cT=cT,
                    acc=1.0e-3,
                    **kwargs,
                )

                alpha = -1.0 / d * drhodT

        # Mie-Grüneisen-Debye (MGD)
        elif whichEOS == 1 and not stopper:

            if what == "pres":
                eta = d / d0
                P = P_MGD(
                    eta=eta,
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    T=T,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    **kwargs,
                )

            elif what == "dens":
                d = d0 * ftool.bisec(
                    f=P_MGD,
                    whicharg="eta",
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    T=T,
                    y=P,
                    a=dmin / d0,
                    b=dmax / d0,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    **kwargs,
                )

                d = d * hydration.delta_rho(P=P, T=T, ll=ll, X_H2O=X_H2O)

            elif what == "temp":
                eta = d / d0
                T = ftool.bisec(
                    f=P_MGD,
                    whicharg="T",
                    K0=K0,
                    K0prime=K0prime,
                    a=Tmin,
                    b=Tmax,
                    y=P,
                    ll=ll,
                    eta=eta,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    **kwargs,
                )

            elif what == "all":
                # Compute density
                d = rho_MGD(
                    T=T,
                    P=P,
                    dmin=dmin,
                    dmax=dmax,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    ll=ll,
                    K0=K0,
                    K0prime=K0prime,
                )
                d = d * hydration.delta_rho(P=P, T=T, ll=ll, X_H2O=X_H2O)

                # Compute dPdrho at given density and temperature
                eta = d / d0
                dPdrho = dPdrho_MGD(
                    eta=eta,
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    **kwargs,
                )

                drhodT = ftool.deriv(
                    f=rho_MGD,
                    whicharg="T",
                    P=P,
                    dmax=dmax,
                    dmin=dmin,
                    ll=ll,
                    x0=T,
                    acc=1.0e-3,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    K0=K0,
                    K0prime=K0prime,
                )

                alpha = -1.0 / d * drhodT

        # Vinet (Vinet)
        elif whichEOS == 2 and not stopper:

            if what == "pres":
                eta = d / d0
                P = P_Vinet(eta=eta, K0=K0, K0prime=K0prime, **kwargs)

            elif what == "dens":
                d = d0 * ftool.bisec(
                    f=P_Vinet,
                    whicharg="eta",
                    K0=K0,
                    K0prime=K0prime,
                    a=dmin / d0,
                    b=dmax / d0,
                    y=P,
                    **kwargs,
                )

            elif what == "temp":
                eta = d / d0
                T = ftool.bisec(
                    f=P_Vinet,
                    whicharg="T",
                    K0=K0,
                    K0prime=K0prime,
                    a=Tmin,
                    b=Tmax,
                    y=P,
                    ll=ll,
                    eta=eta,
                    **kwargs,
                )

            elif what == "all":
                # Compute density
                d = d0 * ftool.bisec(
                    f=P_Vinet,
                    whicharg="eta",
                    K0=K0,
                    K0prime=K0prime,
                    a=dmin / d0,
                    b=dmax / d0,
                    y=P,
                    **kwargs,
                )

                # Compute dPdrho at given density and temperature
                eta = d / d0
                dPdrho = dPdrho_Vinet(eta=eta, K0=K0, K0prime=K0prime, ll=ll)

        # Belonoshko (Bel)
        elif whichEOS == 3 and not stopper:
            if what == "pres":
                eta = d / d0
                print(" Bel:", eta, d0)
                print(" T, rho =", T, d)
                P = P_Bel(
                    eta=eta,
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    T=T,
                    q=q,
                    molmass=molmass,
                    gamma0=gamma0,
                    rho=d0,
                    **kwargs,
                )

            elif what == "dens":
                d = d0 * ftool.bisec(
                    f=P_Bel,
                    whicharg="eta",
                    K0=K0,
                    K0prime=K0prime,
                    a=dmin / d0,
                    b=dmax / d0,
                    y=P,
                    ll=ll,
                    T=T,
                    q=q,
                    molmass=molmass,
                    gamma0=gamma0,
                    rho0=d0,
                    **kwargs,
                )

            elif what == "temp":
                eta = d / d0
                T = ftool.bisec(
                    f=P_Bel,
                    whicharg="T",
                    K0=K0,
                    K0prime=K0prime,
                    a=Tmin,
                    b=Tmax,
                    y=P,
                    ll=ll,
                    eta=eta,
                    q=q,
                    molmass=molmass,
                    gamma0=gamma0,
                    rho0=d0,
                    **kwargs,
                )

            elif what == "all":
                # Compute density
                d = rho_Bel(
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    q=q,
                    dmin=dmin,
                    dmax=dmax,
                    gamma0=gamma0,
                    rho0=d0,
                    molmass=molmass,
                    T=T,
                    P=P,
                )

                # Compute dPdrho at given density and temperature
                eta = d / d0
                dPdrho = dPdrho_Bel(
                    eta=eta,
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    T=T,
                    q=q,
                    molmass=molmass,
                    gamma0=gamma0,
                    rho0=d0,
                )

                alpha = (
                    -1.0
                    / d
                    * ftool.deriv(
                        f=rho_Bel,
                        whicharg="T",
                        P=P,
                        dmax=dmax,
                        dmin=dmin,
                        ll=ll,
                        x0=T,
                        acc=1.0e-3,
                        q=q,
                        molmass=molmass,
                        nparticles=nparticles,
                        gamma0=gamma0,
                        rho0=d0,
                        K0=K0,
                        K0prime=K0prime,
                    )
                )

                phase = 0

        # Water EOS
        elif whichEOS == 4 and not stopper:
            raise NotImplementedError("This feature is not available yet.")

            if what == "pres":
                P = eoswater.pressure(d=d, t=T, p0=P, **kwargs)
                phase = eoswater.phase(t=T, p=P, **kwargs)

            elif what == "dens":
                d = rho_water_eos(P=P, T=T, **kwargs)
                phase = eoswater.phase(T=T, P=P, dens=d, **kwargs)

            elif what == "temp":
                T = eoswater.temperature(d=d, p=P, **kwargs)
                phase = eoswater.phase(t=T, p=P, **kwargs)

            elif what == "all":
                # Compute density, phase and dPdrho
                (
                    d,
                    dTdP_S,
                    dPdrho,
                    alpha,
                    cp_spec,
                    s_spec,
                    u_spec,
                    phase,
                ) = eoswater.together(T=T, P=P, ph=phase)

        # Ideal gas law (IGL)
        elif whichEOS == 5:  # for H and He gas in the atmosphere
            if what == "pres":
                P = P_IGL(d=d, T=T, ll=ll)

            elif what == "dens":
                d = rho_IGL(P=P, T=T, ll=ll)

        # incompressible substances
        elif whichEOS == 6:
            if what == "dens":
                print(
                    "WARNING: Incompressible substance in eos.Compute() not available"
                )

        # Vinet_Rydberg
        elif whichEOS == 7:
            if what == "pres":
                eta = d / d0
                P = P_VR(
                    eta=eta,
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    T=T,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    theta0=thetaD0,
                    **kwargs,
                )

            elif what == "dens":
                d = d0 * ftool.bisec(
                    f=P_VR,
                    whicharg="eta",
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    T=T,
                    y=P,
                    a=dmin / d0,
                    b=dmax / d0,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    theta0=thetaD0,
                    **kwargs,
                )

                d = d * hydration.delta_rho(P=P, T=T, ll=ll, X_H2O=X_H2O)

            elif what == "temp":
                eta = d / d0
                T = ftool.bisec(
                    f=P_VR,
                    whicharg="T",
                    K0=K0,
                    K0prime=K0prime,
                    a=Tmin,
                    b=Tmax,
                    y=P,
                    ll=ll,
                    eta=eta,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    theta0=thetaD0,
                    **kwargs,
                )

            elif what == "all":
                # Compute density
                d = rho_VR(
                    T=T,
                    P=P,
                    dmin=dmin,
                    dmax=dmax,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    theta0=thetaD0,
                    ll=ll,
                    K0=K0,
                    K0prime=K0prime,
                )

                d = d * hydration.delta_rho(P=P, T=T, ll=ll, X_H2O=X_H2O)

                # Compute dPdrho at given density and temperature
                eta = d / d0
                dPdrho = dPdrho_VR(
                    eta=eta,
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    theta0=thetaD0,
                    **kwargs,
                )

                drhodT = ftool.deriv(
                    f=rho_VR,
                    whicharg="T",
                    P=P,
                    dmax=dmax,
                    dmin=dmin,
                    ll=ll,
                    x0=T,
                    acc=1.0e-3,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    theta0=thetaD0,
                    K0=K0,
                    K0prime=K0prime,
                    eta=eta,
                )

                alpha = -1.0 / d * drhodT

        # Ichi dude
        elif whichEOS == 8:
            if what == "pres":
                eta = d / d0
                print("d, d0 =", d, d0)
                P = P_Ichi(
                    eta=eta,
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    T=T,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    **kwargs,
                )

            elif what == "dens":
                d = d0 * ftool.bisec(
                    f=P_Ichi,
                    whicharg="eta",
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    T=T,
                    y=P,
                    a=dmin / d0,
                    b=dmax / d0,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    **kwargs,
                )
                d = d * hydration.delta_rho(P=P, T=T, ll=ll, X_H2O=X_H2O)

            elif what == "temp":
                eta = d / d0
                T = ftool.bisec(
                    f=P_Ichi,
                    whicharg="T",
                    K0=K0,
                    K0prime=K0prime,
                    a=Tmin,
                    b=Tmax,
                    y=P,
                    ll=ll,
                    eta=eta,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    **kwargs,
                )

            elif what == "all":
                # Compute density
                d = rho_Ichi(
                    T=T,
                    P=P,
                    dmin=dmin,
                    dmax=dmax,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    ll=ll,
                    K0=K0,
                    K0prime=K0prime,
                )

                d = d * hydration.delta_rho(P=P, T=T, ll=ll, X_H2O=X_H2O)

                # Compute dPdrho at given density and temperature
                eta = d / d0
                dPdrho = dPdrho_Ichi(
                    eta=eta,
                    K0=K0,
                    K0prime=K0prime,
                    ll=ll,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    **kwargs,
                )

                drhodT = ftool.deriv(
                    f=rho_Ichi,
                    whicharg="T",
                    P=P,
                    dmax=dmax,
                    dmin=dmin,
                    ll=ll,
                    x0=T,
                    acc=1.0e-3,
                    q=q,
                    molmass=molmass,
                    nparticles=nparticles,
                    gamma0=gamma0,
                    rho0=d0,
                    thetaD0=thetaD0,
                    K0=K0,
                    K0prime=K0prime,
                    eta=eta,
                )

                alpha = -1.0 / d * drhodT

    else:
        pass
        """
        if what == 'dens':
            d = tabs[ll].interpolate(param=2, x=T, y=P)
        """
    gamma = 1.2
    # KS = K*(1.+gamma *alpha*T)
    # print ('KS =', KS*1e-9)
    # print ('dTdP =', gamma/KS)
    return d, T, P, dTdP_S, dPdrho, phase, X_H2O / 100.0, alpha, xi_Al


def dTdP_ad(
    ll=None, T=None, phase=None, P=None, Fe_number=0.0, saturation=False, gamma=1,
):
    stuff = compute(
        ll=ll, T=T, P=P, Fe_number=Fe_number, saturation=saturation, phase=phase
    )

    KT = -stuff[0] * stuff[4]

    return gamma * T / (1.0 + gamma * stuff[7] * T) / KT

    # def PlotCmap(
    #     ll=1,
    #     Tmin=300,
    #     Tmax=2400,
    #     Pmin=0,
    #     Pmax=30,
    #     res=4,
    #     params=[0],
    #     Fe_numbers=[0.0, 25.0],
    # ):
    #     """Plots parameter as function of pressure and temperature in a colormap.
    #     Pressure and temperature ranges are defined via Pmin, Pmax and Tmin, Tmax
    #     in units of GPa and K. The PT-grid is given by a 2**res x 2**res array.
    #     """

    #     # Compute grid resolution
    #     N = 2**res

    #     ticklabel_fontsize = 8
    #     label_fontsize = 10

    #     nxticks = 4
    #     nyticks = 6

    #     param_labels = [
    #         r"$\rho \ [\rm kg /m^3]$",
    #         "T",
    #         "P",
    #         r"$d T/dP_S$",
    #         r"$d P/d \rho \ [\rm GPa \ m³ /kg]$",
    #         "phase",
    #         r"$\rm log(X_{\rm H_2 O}) \ \rm [wt \%]$",
    #         r"$\alpha_T \ [10^{-5} K^{-1}]$",
    #         r"$\xi_{\rm Al}$",
    #     ]

    #     cmaps = ["rainbow", "viridis", "seismic", "PuOr"]
    #     params_scalings = [1.0, 1.0, 1.0, 1.0e-9, 1.0, 1.0, 1.0, 1.0e5, 1.0]
    #     data_lims = [[3000.0, 4200.0], [-2, 1], [0.0, 10.0], [0, 10]]

    #     border_colors = [
    #         [0.1, 0, 0],
    #         [0.2, 0, 0],
    #         [0.25, 0, 0],
    #         [0.5, 0, 0],
    #         # [.6, 0, 0],
    #         # [.75, 0, 0],
    #         # [.9, 0, 0],
    #         [1, 0, 0],
    #         [1, 0.25, 0],
    #         [1, 0.5, 0],
    #         [1, 1, 0],
    #         [0.75, 0.5, 0],
    #         # [.75, .25, 0],
    #         [0.7, 0.3, 0],
    #         [0.6, 0.4, 0],
    #         [0.5, 0.5, 0],
    #         [0, 1, 0],
    #         [0, 0, 1],
    #         [1, 1, 1],
    #     ]

    #     border_colors = list(
    #         reversed(
    #             [
    #                 (0.7, 0.2, 0.6),
    #                 (0.6, 0.4, 1.0),
    #                 (0.2, 0.4, 1.0),
    #                 (0.2, 0.6, 0.9),
    #                 (0.2, 0.8, 0.8),
    #                 (0.2, 0.8, 0.4),
    #                 (0.6, 0.8, 0.4),
    #                 (0.6, 0.6, 0.2),
    #                 (0.8, 0.4, 0.2),
    #                 (1.0, 0.2, 0.2),
    #                 (1.0, 0.5, 0.5),
    #                 (0.7, 0.7, 0.7),
    #                 (0.5, 0.5, 0.5),
    #                 (0.2, 0.2, 0.2),
    #                 (0.0, 0.0, 0.0),
    #             ]
    #         )
    #     )

    #     bin_refine_exponent = 4
    #     n_additional_bins = 2**bin_refine_exponent
    #     colors = []

    #     xticks = np.linspace(0, N - 1, nxticks)
    #     yticks = np.linspace(0, N - 1, nyticks)

    #     xtick_labels = np.linspace(Tmin, Tmax, nxticks, dtype=int)
    #     ytick_labels = np.linspace(Pmin, Pmax, nyticks, dtype=int)

    #     for i in range(len(border_colors) - 1):
    #         colors.append(border_colors[i])
    #         for r in range(n_additional_bins):
    #             colors.append(
    #                 np.asarray(border_colors[i])
    #                 + (r + 1)
    #                 / (n_additional_bins + 1)
    #                 * (np.asarray(border_colors[i + 1]) - np.asarray(border_colors[i]))
    #             )

    #     colors.append(border_colors[-1])

    #     newcmap = mpl.colors.ListedColormap(colors)

    #     cmaps[2] = newcmap
    #     cmaps = [newcmap, newcmap, newcmap, newcmap]

    #     x = np.linspace(Tmin, Tmax, N)
    #     y = np.linspace(Pmin, Pmax, N)

    #     xx, yy = np.meshgrid(x, y)

    #     zz = np.zeros([len(Fe_numbers), len(params), N, N])

    #     # Compute data
    #     for f in range(len(Fe_numbers)):
    #         for i in range(N):
    #             for j in range(N):

    #                 eos = Compute(
    #                     ll=ll,
    #                     T=x[j],
    #                     P=y[i] * 1.0e9,
    #                     saturation=True,
    #                     Fe_number=Fe_numbers[f],
    #                 )

    #                 dat = [eos[p] for p in params]

    #                 for k in range(len(params)):
    #                     zz[f][k][i][j] = dat[k] * params_scalings[params[k]]

    #     nrows = len(params)
    #     ncols = len(Fe_numbers)

    #     fig = plt.figure(figsize=(10, 8))
    #     axrows = []

    #     i = 0
    #     for r in range(nrows):

    #         axcols = AxesGrid(
    #             fig,
    #             (nrows, 1, r + 1),
    #             nrows_ncols=(1, ncols),
    #             axes_pad=0.25,
    #             share_all=True,
    #             label_mode="L",
    #             cbar_mode="edge",
    #             cbar_location="right",
    #             cbar_size="7%",
    #             cbar_pad="10%",
    #         )

    #         axrows.append(axcols)

    #         for c in range(ncols):
    #             ax = axcols[c]
    #             ax.tick_params(right=True, top=True)
    #             im = ax.imshow(
    #                 zz[c][r],
    #                 origin="lower",
    #                 cmap=cmaps[r],
    #                 vmin=data_lims[r][0],
    #                 vmax=data_lims[r][1],
    #             )

    #             if r < nrows - 1:
    #                 ax.tick_params(labelbottom=False)

    #             else:
    #                 ax.set_xlabel(r"$T \ \rm [K]$", fontsize=label_fontsize)

    #             if c == 0:
    #                 ax.set_ylabel(r"$P \ \rm [GPa]$", fontsize=label_fontsize)

    #             if r == 0:
    #                 ax.set_title(
    #                     r"$\rm {a} \ mol \% \ Fe$".format(a=str(Fe_numbers[c])),
    #                     fontsize=label_fontsize,
    #                 )

    #         for rax in axcols.cbar_axes:
    #             rax.colorbar(im)
    #             rax.tick_params(
    #                 labelright=True,
    #                 labelleft=False,
    #                 left=False,
    #                 right=True,
    #                 labelsize=ticklabel_fontsize,
    #             )

    #             rax.set_ylabel(param_labels[params[r]], fontsize=label_fontsize)

    #     for r in range(nrows):
    #         for c in range(ncols):
    #             ax = axrows[r][c]

    #             ax.set_xticks(xticks)
    #             ax.set_xticklabels(xtick_labels, fontsize=ticklabel_fontsize)

    #             ax.set_yticks(yticks)
    #             ax.set_yticklabels(ytick_labels, fontsize=ticklabel_fontsize)

    # plt.savefig('/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/XH2O_Ol.pdf',
    #           format='pdf', bbox_inches='tight')
    # plt.close(fig)
    """
    #Plot data
    fig, ax = plt.subplots(4, len(Fe_numbers)+1, sharex=True, sharey = True,
                 )
    cbar_ax = []
    
    fig = plt.figure()
    gs = fig.add_gridspec(4, len(Fe_numbers))

    if len(Fe_numbers) == 1:
        ax = [ax]
    
    
    ax = []
    
    for i in range(4):
        ax.append([])
        for f in range(len(Fe_numbers)):
           
            #ax[i].append(plt.subplot2grid((8,3), (i, f)))
            ax[i].append(fig.add_subplot(gs[i, f]))
    
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    
    for f in range(len(Fe_numbers)):
        for i in range(4):
            
            if i == 0:
                lt = True
                
            else:
                lt = False
            
            if i == 3:
                ax[i][f].set_xlabel(r'$T \ [\rm K]$')

            if f == 0:
                ax[i][f].set_ylabel(r'$P \ [\rm GPa]$')
            
            ax[i][f].tick_params(top=True, right=True, labeltop=lt)
            
            im = ax[i][f].imshow(zz[f][i], origin='lower', 
                   extent=[Tmin, Tmax, Pmin, Pmax],
                   aspect='auto', cmap = cmap[i])
    """


# def plotFancy(
#     ll=1,
#     N=10,
#     P=[10, 12],
#     isotherms=np.linspace(1e3, 5e3, 3),
#     mark="",
#     Fe_number=100,
#     phase=0,
#     presScale=1e-5,
#     **kwargs,
# ):

#     try:
#         ax = kwargs["ax"]
#     except KeyError:
#         fig, ax = plt.subplots()
#     cols = ["r", "g", "b"]
#     pres = np.logspace(P[0], P[1], N)
#     dens = np.array(
#         [
#             [
#                 Compute(
#                     what="all",
#                     ll=ll,
#                     P=p,
#                     T=iso,
#                     Fe_number=Fe_number,
#                     phase=phase,
#                     **kwargs,
#                 )[0]
#                 for p in pres
#             ]
#             for iso in isotherms
#         ]
#     )

#     for i in range(len(dens)):
#         ax.semilogx(pres * presScale, dens[i] * 1e-3, color=cols[i])


# def plotOverview(mats=[1, 9], fnts=14, N=20):
#     fig, ax = plt.subplots(1, len(mats), figsize=(9, 3))
#     plt.subplots_adjust(wspace=0.25)
#     temps = np.linspace(1e3, 5e3, 3)
#     pres = [8, 12]
#     presScale = 11
#     for i in range(len(mats)):
#         plotFancy(
#             ll=mats[i],
#             ax=ax[i],
#             fnts=fnts,
#             isotherms=temps,
#             P=pres,
#             N=N,
#             presScale=10 ** (-presScale),
#         )
#         trafo = transforms.blended_transform_factory(ax[i].transAxes, ax[i].transAxes)

#         a = material_plot_list[mats[i]]
#         b = EOSstring_list[EOS_type[mats[i]]]
#         label = f"${a} \ ({b})$"
#         ax[i].text(0.1, 0.9, label, transform=trafo, fontsize=fnts)
#         ax[i].set_xlabel(r"$\rm Pressure \ [Mbar]$", fontsize=fnts)
#         ax[i].tick_params(labelsize=fnts)
#         ax[i].set_xticks(
#             np.logspace(pres[0] - presScale, pres[1] - presScale, pres[1] - pres[0] + 1)
#         )
#         ax[i].set_xlim(10 ** (pres[0] - presScale), 10 ** (pres[1] - presScale))

#     ax[0].set_ylabel(r"$\rm Density \ [g / cm^3]$", fontsize=fnts)
#     plots = ax[0].lines
#     u = r"\rm K"
#     ax[0].legend(plots, [f"${T:.0f} \ {u}$" for T in temps], fontsize=fnts, loc=3)
#     fig.savefig(
#         "/home/os18o068/Documents/PHD/Abbildungen/eos.pdf",
#         format="pdf",
#         bbox_inches="tight",
#     )
#     plt.close(fig)


# def plotTable(
#     N=2,
#     off=0,
#     params=[2],
#     ll=1,
#     cmap="jet",
#     scalings=[1e-3, 1e-8, 1e5],
#     diff=False,
#     starts=[500, 1e5],
#     ends=[3000, 1e11],
#     fnts=14,
# ):
#     import eosTablesTest

#     indexTrafo = {2: 0, 1: 3, 3: 3, 4: 4, 5: 7}
#     mpl.rc("text", usetex=True)
#     mpl.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath, amssymb}"]

#     fig, ax = plt.subplots(len(params), N, figsize=(9, 9))

#     plt.subplots_adjust(wspace=0.25, hspace=0.5)
#     labels = [
#         r"$\rho \ \rm [g \ cm^{-3}]$",
#         r"$ dP/d\rho \ \rm [10^8 \ Pa \ g^{-1} \ cm^{3}]$",
#         r"$\alpha_{th} \ \rm [10^{-5} \ K^{-1}]$",
#     ]
#     paramLabels = [r"\rho", r"(dP/d\rho)", r"\alpha_{\rm th}"]
#     extent = [starts[0], ends[0], starts[1], ends[1]]
#     aspect = (extent[1] - extent[0]) / (extent[3] - extent[2])

#     # Prepare tables
#     tables = []
#     for i in range(N):
#         tab = eosTablesTest.Table(
#             starts=starts, ends=ends, scalings=["lin", "lin"], alphas=[i + off, i + off]
#         )
#         tab.construct_single_axes()
#         tab.construct_grid()
#         tab.generate(ll=ll, Fe_number=0.0)
#         tables.append(tab)

#     NN = 2 ** (N + off + 1)
#     # Compute relative error on table interpolation
#     if diff:
#         datDiff = []

#         for n in range(N):
#             tab = tables[n]
#             x = np.linspace(starts[0] * 1.001, ends[0] * 0.999, NN)
#             y = np.linspace(starts[1] * 1.001, ends[1] * 0.999, NN)

#             datAcc = np.empty([len(x), len(y), len(params)])
#             datInt = np.empty([len(x), len(y), len(params)])

#             for i in range(len(x)):
#                 for j in range(len(y)):
#                     # Compute accurate data from EoS
#                     da = Compute(
#                         ll=ll, what="all", P=y[j], T=x[i], Fe_number=0.0, phase=0
#                     )

#                     # Compute data from table interpolation
#                     di = tab.interpolate(vals=[x[i], y[j]], params=params)
#                     print("di =", di)

#                     datInt[i][j][:] = di
#                     for ind in range(len(params)):
#                         datAcc[i][j][ind] = da[indexTrafo[params[ind]]]

#             # Compute realtive difference of accurate and interpolated data
#             datDiff.append(np.log10(abs(datAcc - datInt) / datAcc))

#     for r in range(len(params)):
#         for i in range(N):
#             ax[r][i].tick_params(labelsize=fnts)
#             tab = tables[i]
#             if diff:

#                 im = ax[r][i].imshow(
#                     datDiff[i].T[r],
#                     vmin=-5,
#                     vmax=-1,
#                     cmap=cmap,
#                     origin="lower",
#                     extent=extent,
#                     aspect=aspect,
#                 )

#                 xTicks = np.linspace(starts[0], ends[0], 3)
#                 yTicks = np.linspace(starts[1], ends[1], 5)
#                 if scalings[0] == "log":
#                     xString = f"$log({tab.param_strings[0]})$"
#                 else:
#                     xString = tab.param_strings[0]

#                 if scalings[1] == "log":
#                     fore = r"\rm log"
#                     yString = f"${fore}({tab.param_strings[1]})$"
#                 else:
#                     yString = tab.param_strings[1]

#                 ax[r][i].set_xticklabels([f"${t:.0f}$" for t in xTicks], fontsize=fnts)
#                 ax[r][i].set_yticklabels(
#                     [f"${t:.0f}$" for t in yTicks * 1e-9], fontsize=fnts
#                 )

#                 ax[r][i].set_xlabel(xString, fontsize=fnts)
#                 ax[r][i].set_ylabel(yString, fontsize=fnts)
#                 ax[r][i].set_xticks(xTicks)
#                 ax[r][i].set_yticks(yTicks)
#                 if i == N - 1:
#                     axin = ax[r][i].inset_axes([1.1, 0.0, 0.05, 1.0])
#                     cbar = plt.colorbar(im, cax=axin)
#                     cbar.set_label(
#                         f"$\log(|\Delta {paramLabels[r]}|/{paramLabels[r]})$",
#                         fontsize=fnts,
#                     )
#                     cbar.ax.tick_params(labelsize=fnts)

#             else:
#                 if i == N - 1:
#                     addCbar = True
#                 else:
#                     addCbar = False
#                 tab.plot_slice(
#                     params=[params[r]],
#                     scalings=["lin", "lin"],
#                     ax=ax[r][i],
#                     fig=fig,
#                     label=labels[r],
#                     addCbar=addCbar,
#                     cmap=cmap,
#                     scale=scalings[r],
#                     fnts=fnts,
#                 )

#             # Add label for resolution
#             trafo = transforms.blended_transform_factory(
#                 ax[r][i].transAxes, ax[r][i].transAxes
#             )
#             res = 2 ** (i + off + 1)
#             ax[r][i].text(
#                 0.1,
#                 0.9,
#                 r"$\rm Res = {a}\times{b}$".format(a=res, b=res),
#                 transform=trafo,
#                 color="white",
#                 fontsize=fnts,
#             )

#     if diff:
#         fileName = "eos_grid_diff"
#     else:
#         fileName = "eos_grid"

#     fig.savefig(
#         r"/home/os18o068/Documents/PHD/Abbildungen/{fn}.pdf".format(fn=fileName),
#         format="pdf",
#         bbox_inches="tight",
#     )
#     plt.close(fig)


# def Plot(
#     ll=1, N=5, P_min=1.0e5, P_max=1.0e12, mark="", phase=0, Fe_number=0.0, **kwargs
# ):
#     """Generates simple phase diagram over a certain temperature and pressure
#     range accounting for all the implemented phases for the given material
#     """

#     # define temperature range
#     pres = np.logspace(np.log10(P_min), np.log10(P_max), N)
#     temp = np.logspace(2, 4, N)

#     isotherms = np.linspace(1000, 3000, 2)

#     dens = np.array(
#         [
#             [
#                 Compute(
#                     what="all",
#                     ll=ll,
#                     P=p,
#                     T=iso,
#                     Fe_number=Fe_number,
#                     phase=phase,
#                     **kwargs,
#                 )[0]
#                 for p in pres
#             ]
#             for iso in isotherms
#         ]
#     )

#     if ll == 1:
#         isochores = np.linspace(rho0_list[ll], rho0_list[ll] * 1.5, 5)

#     else:
#         isochores = np.linspace(rho0_list[ll], rho0_list[ll] * 1.5, 5)

#     """
#     #to compute the pressure the phase region must be known as it can not
#     #be determined from the phase diagram which requires a PT point as input
#     pres_plot = np.array([[Compute(what='pres', ll=ll, T=t, d=iso,
#                                    phase = phase, P=1.0e9,
#                                    Fe_number=Fe_number, **kwargs)[0]
#                         for t in temp] for iso in isochores])
#     """
#     # compute phase transition pressures for temperature range
#     phase_pres_alpha = hydration.P_alpha_beta_fit(temp) * 1.0e9
#     phase_pres_beta = hydration.P_beta_gamma_fit(temp) * 1.0e9

#     plot_list1 = []
#     plot_list2 = []

#     fnts = 14
#     labelsize = 12
#     lwdt = 2

#     # create two subfigures to plot isotherms and isochores seperately
#     fig, ax = plt.subplots(1, 2, sharey=True)

#     for i in range(len(dens)):
#         (plot,) = ax[0].semilogy(
#             dens[i], pres, color=color_list[i], linewidth=lwdt, marker=mark
#         )
#         plot_list1.append(plot)
#     """
#     for i in range(len(pres_plot)):
#         plot, = ax[1].loglog(temp, pres_plot[i], color=color_list[i],
#                   linewidth=lwdt, marker=mark)
#         plot_list2.append(plot)
#     """
#     """
#     #plot phase curves for Mg2SiO4
#     if ll == 1:
#         ax[1].plot(temp, phase_pres_alpha, color='grey', linestyle='--')
#         ax[1].plot(temp, phase_pres_beta, color='grey', linestyle='--')
#     """

#     legend1 = ax[0].legend(
#         plot_list1, [str(iso) + " K" for iso in isotherms], loc=4, framealpha=0.5
#     )

#     legend2 = ax[1].legend(
#         plot_list1,
#         [str(int(iso)) + r"$ \ kg \ m^{-3}$" for iso in isochores],
#         loc=3,
#         framealpha=0.5,
#     )

#     ax[0].add_artist(legend1)
#     ax[1].add_artist(legend2)
#     ax[1].tick_params(labelright=True)
#     ax[0].set_ylabel(r"$P \ [Pa]$", fontsize=labelsize)
#     ax[0].set_xlabel(r"$\rho \ [kg \ m^{-3}]$", fontsize=labelsize)
#     ax[1].set_xlabel(r"$T \ [K]$", fontsize=labelsize)

#     ml = MultipleLocator(100)
#     ax[0].xaxis.set_minor_locator(ml)
#     fig.suptitle(r"$ phase \ diagram \  $" + material_plot_list[ll])

#     if ll == 1:
#         ax[1].text(3100, 1.1e9, r"$\alpha$", color="grey", fontsize=fnts)
#         ax[1].text(3100, 1.9e10, r"$\beta$", color="grey", fontsize=fnts)
#         ax[1].text(3100, 5.0e10, r"$\gamma$", color="grey", fontsize=fnts)

#     for axx in ax:
#         axx.tick_params(
#             which="both",
#             direction="in",
#             right=True,
#             top=True,
#             pad=10,
#             labelsize=labelsize,
#         )
#         axx.set_ylim(P_min, P_max)
#         # axx.grid(color=grid_color, linewidth=2)
#         # axx.grid(which='minor', color=grid_color)

#     # ax[0].set_xlim(2500, 5000)
#     # ax[1].set_xlim(100, 10000)
#     fig.savefig(
#         "/home/os18o068/Documents/PHD/Abbildungen/eos.pdf",
#         format="pdf",
#         bbox_inches="tight",
#     )
#     plt.close(fig)


"""
fig, ax = plt.subplots()
ax.set_xlabel('Density [kg m/3]')
ax.set_ylabel('Pressure [Pa]')
ax.grid(zorder=0)

temp = 300
mat = 9
dens = rho0_list[mat]
print ('dPdrho =',g(x=dens))
pres = Compute(ll=mat, d=dens,T=temp, what='pres', P=2.0e5)

test = ftool.integrate(dydx=g, start=dens, y=[pres], h=1., T=temp, noisy=False,
                       whicharg='x', plot=True, end=dens+100., N=100, axis=ax,
                       ll=mat, temperature_type = 'ambient isothermal')
print ('pres:',pres)
print ('dens',Compute(what='dens', ll=mat, P=pres, T=temp, acc=1.0e-8))

xx=np.linspace(dens, dens+100, 20)
yy = []
for x in xx:
    yy.append(Compute(ll=mat, d=x, T=temp, what='pres',
                      temperature_type = 'ambient isothermal'))

ax.scatter(xx,yy, zorder=2)

def f(x=0):
    return Compute(ll=mat, what='pres',d=x, P=1.0e12, T=temp)

test = ftool.deriv(f=f, x0=dens, whicharg='x', acc = 1.0e-6, noisy = False, 
                   plot = False, eps=1.0e-3, brash = False)

pres = Compute(ll=mat, what='pres',d=dens, P=1.0e12, T=temp)
print ('test=',test)
print ('pres =',pres)
"""
