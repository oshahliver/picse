#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 14:26:56 2019

@author: oshah
"""

import time
import sys
from pics.utils import PIMPtools
import numpy as np
from matplotlib import pyplot as plt
from pics.utils import functionTools as ftool
from pics.physicalparams import (
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
)

# Mixing rule for material mixtures
def mix1(densities, fractions):
    return 1.0 / sum([fractions[i] / densities[i] for i in range(len(fractions))])


# Birch-Murnaghan (BM)
def P_BM(eta=0, K=0, K0prime=0):
    return (
        3.0
        / 2.0
        * K
        * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
        * (1.0 + 3.0 / 4.0 * (K0prime - 4.0) * (eta ** (2.0 / 3.0) - 1.0))
    )


def dPdrho_BM(ll=0, eta=0, K=0, K0prime=0, d0T=0):
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


def dP_MGD(eta=0, ll=5, T=300, **kwargs):
    d = eta * rho0_list[ll]
    tempType = ftool.checkKey(name="tempType", **kwargs)

    if tempType == "isothermal":
        deltaP = 0.0
    else:
        E1 = PIMPtools.Eth(d, T, ll, **kwargs)
        E2 = PIMPtools.Eth(d, 300.0, ll, **kwargs)
        deltaP = (E1 - E2) * PIMPtools.gamma(d, ll)

    return deltaP


def P_MGD(eta=0.0, K0=0.0, K0prime=0.0, ll=5, T=300, **kwargs):
    d = eta * rho0_list[ll]
    tempType = ftool.checkKey(name="tempType", **kwargs)

    if tempType == "isothermal":
        deltaP = 0.0
    else:
        # call routine to compute thermal energy contribution
        E1 = PIMPtools.Eth(d, T, ll, **kwargs)
        E2 = PIMPtools.Eth(d, 300.0, ll, **kwargs)
        deltaP = (E1 - E2) * PIMPtools.gamma(d, ll)

    return (
        3.0
        / 2.0
        * K0
        * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
        * (1.0 + 3.0 / 4.0 * (K0prime - 4.0) * (eta ** (2.0 / 3.0) - 1.0))
        + deltaP
    )


def dPdrho_MGD(eta=0, K0=0, K0prime=0, ll=0, T=300):
    d0 = rho0_list[ll]

    ddeltaPth = ftool.deriv(f=dP_MGD, x0=eta, whicharg="eta", T=T)
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
def P_Bel(eta=0, K0=0, K0prime=0, ll=0, T=300, **kwargs):
    d = eta * rho0_list[ll]
    tempType = ftool.checkKey(name="tempType", **kwargs)
    if tempType == "isothermal":
        deltaP = 0.0
    else:
        deltaP = (
            3.0
            * Rgas
            * PIMPtools.gamma(d, ll)
            * (T - T0_list[ll])
            * (d / molar_mass_list[ll])
        )
    return (
        3.0
        / 2.0
        * K0
        * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
        * (1.0 + 3.0 / 4.0 * (K0prime - 4.0) * (eta ** (2.0 / 3.0) - 1.0))
        + deltaP
    )


def dPdrho_Bel(eta=0, K0=0, K0prime=0, ll=0, T=300, **kwargs):
    d0 = rho0_list[ll]
    d = eta * d0
    tempType = ftool.checkKey(name="tempType", **kwargs)
    if tempType == "isothermal":
        delta = 0.0
    else:
        delta = (
            3.0
            * Rgas
            * PIMPtools.gamma(d, ll)
            * (T - T0_list[ll])
            / molar_mass_list[ll]
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


# van der Waal EOS for fluids (VDW)
def P_VDW(d=0.0, T=0.0, ll=None, **kwargs):
    mm = gas_molar_mass_list[ll]
    a = a_VDW_list[ll]
    b = b_VDW_list[ll]
    return d * Rgas * T / (mm - d * b) - d ** 2 * a / mm ** 2


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
    if P < 1.0e10:
        return P * mm / NA / kB / T
    else:
        return 86.0


def dPdrho_IGL(T=300, ll=0):
    mm = gas_molar_mass_list[ll]
    return kB * T * NA / mm


def dTdP_IGL(d=0.0, ll=0):
    mm = gas_molar_mass_list[ll]
    return mm / (d * kB * NA)
