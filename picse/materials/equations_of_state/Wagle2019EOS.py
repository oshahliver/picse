# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 11:01:20 2023

@author: oshahliver
"""

import numpy as np


# Define all relevant constants
# NOTE. for details and units see supplementary in Wagle & Steinle-Neumann 2019
V0 = 7.12
rhoe = 7040
pe = 11.59769255256658
Ke = -0.13091598231783
X = 2.1293455089530795
xi0 = 4.891e-3
xi1 = 2.899e-5
xi = 0.849
s = 0.79
n = 2
Vref = 7.12
Tref = 2000
m = 55.845
kB = 1.380649e-23
NA = 6.02214e23
hbar = 1.05457e-34
alist = [
    [-7.023e2, 2.838e1, -4.243e1],
    [-4.989e1, 4.872e2, -1.1013e2],
    [8.605e3, 9.465e1, 6.302e2],
    [2.836e4, -7.78677e3],
]


def F(V=None, T=None):
    """
    Compute the Helmholtz energy of liquid iron according to eq. 1 of
    Wagle & Steinle-Neumann 2019.

    Parameters
    ----------
    V: FLOAT
        Specific volume in cm^3 per mol.
    T: FLOAT
        Temperature in K.

    Returns
    -------
    FLOAT
        Helmholtz energy in kJ per mol.
    """
    return F_el(V, T) + F_XS(V, T) + F_ig(V, T)  # + F_corr(V)


def F_ig(V, T):
    """
    Ideal gas contribution to the Helmholtz energy of liquid iron according to
    eq. 2 of Walge & Steinle-Neumann 2019.

    Parameters
    ----------
    V : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    res : TYPE
        DESCRIPTION.

    """
    res = (
        -1e-3
        * kB
        * NA
        * T
        * (
            np.log((m * kB * T / (2 * np.pi * hbar ** 2)) ** (3 / 2) * 1e-6 * V / NA)
            + 1
        )
    )
    return res


def F_el(V, T):
    """
    Electronic contribution to Helholtz energy  according to eq. 3 of
    Walge & Steinle-Neumann 2019.

    Parameters
    ----------
    V : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    x = V / V0
    return -1e-3 * 0.5 * xi0 * x ** xi * (T ** 2 - xi1 / 3 * T ** 3)


def F_XS(V, T):
    """
    Compute the excess term to the Helmholtz energy  according to eq. 6 of
    Walge & Steinle-Neumann 2019.

    Parameters
    ----------
    V : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    fac : TYPE
        DESCRIPTION.

    """
    fac = 0.0

    f = 1 / n * ((Vref / V) ** (n / 3.0) - 1.0)
    theta = (T / Tref) ** s - 1.0
    for i in range(4):
        for j in range(3):
            if i + j < 5:
                fac += (
                    alist[i][j]
                    / (np.math.factorial(i) * np.math.factorial(j))
                    * f ** i
                    * theta ** j
                )
            else:
                pass

    return fac


def F_corr(V):
    """
    Compute the low-pressure correction term to the Helmholtz energy  according to
    eq. 7 of Walge & Steinle-Neumann 2019.

    Parameters
    ----------
    V : TYPE
        DESCRIPTION.

    Returns
    -------
    FLOAT
        The correction term to the Helmholtz energy in kJ per mol.

    """
    rho = 55845 / V
    rho0 = rhoe * (1 - Ke / (pe * (X + 1))) ** (-1 / X)
    p0 = pe * (rhoe / rho0) ** (-X - 1) * np.exp((X + 1) / X * ((rhoe / rho0) ** X - 1))
    return p0 / (rho0 * (X + 1)) * (1 - np.exp((X + 1) / X * (1 - (rho / rho0) ** X)))


def P_corr(V):
    """
    Compute the low-pressure correction term to the pressure according to
    eq. 10 of Walge & Steinle-Neumann 2019.

    Parameters
    ----------
    V : TYPE
        DESCRIPTION.

    Returns
    -------
    FLOAT
        The correction term to the Helmholtz energy in kJ per mol.

    """
    rho = 55845 / V
    rho0 = rhoe * (1 - Ke / (pe * (X + 1))) ** (-1 / X)
    p0 = pe * (rhoe / rho0) ** (-X - 1) * np.exp((X + 1) / X * ((rhoe / rho0) ** X - 1))
    return p0 * (rho / rho0) ** (1 + X) * np.exp((X + 1) / X * (1 - (rho / rho0) ** X))
