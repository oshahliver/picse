#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 18:25:13 2019

@author: oshah

Implementation of French et al. 2015 for the EoS for ice VII and X of
pure water between 300 K and 1000 K at elevated pressures P > 0.1 GPa.
It is not applicable to the liquid phase that occures in this temperautre
range.

"""

import numpy as np
import PIMPtools
import functionTools as ftool
from scipy.integrate import quad
from PIMPphysicalparams import kBcgs, kB, mH, mO, NA, Rgas

# Fit parameters from French 2015

TABLE1_PBE = [78.0913, -172.712, 18.095, -0.717051, 130.541, 67.2915]

TABLE2_alpha = [6.869192e-3, -3.919234e-3, 6.790631e-5]

TABLE2_gamma = [
    [-2.604037e-2, 2.104060e-2, -1.515928e-3],
    [6.721305e-2, -1.158057e-1, 5.234025e-2],
    [-1.178112e-1, 1.754266e-1, -6.213482e-2],
    [1.206828e-1, -1.910820e-1, 7.712887e-2],
]

TABLE3_b = [2.08, 0.272, 0.096, 0.788, 2.54e-5]


a0, a1, a2, a3, a4, a5 = TABLE1_PBE
b0, b1, b2, b3, b4 = TABLE3_b
k0 = 3.0
# ------------------------------------------------------------------------------
def D(z=None, n=3, a=0.0):
    """Computes the Debye function in the interval [a,b]"""

    def integrand(x):
        return x ** n / (np.exp(x) - 1.0)

    return n / z ** n * quad(integrand, a, z)[0]


def ue(d=None):
    """Electronic ground state contribution to free energy"""
    return (
        a0 + a1 * d + a2 * d ** 2 + a3 * d ** 3 + a4 * np.log(d) + a5 * np.log(d) ** 2
    )


def un(d=None):
    """Nuclear ground state contribution to free energy"""
    return b0 + b1 * d + b2 * d ** 2 + b3 * np.exp(-b4 * d ** 10)


def ft(d=None, T=None):
    """Thermal contribution to free energy"""
    K1 = 0.0
    K2 = 0.0
    for i in range(1):
        Ti = 700.0 * 2.0 ** i
        for k in range(3):
            kk = k - 4
            # print (i, kk, TABLE2_alpha[k])
            K1 += (
                TABLE2_alpha[k]
                * T
                * (3.0 * np.log(1.0 - np.exp(-Ti / T)) - D(z=Ti / T, n=3))
                * d ** (kk / k0)
            )

    for j in range(4):
        Tj = 1000.0 * 2.0 ** j
        for k in range(3):
            kk = k - 4
            # print (j, kk, TABLE2_gamma[j][k])
            K2 += (
                TABLE2_gamma[j][k] * T * np.log(1.0 - np.exp(-Tj / T)) * d ** (kk / k0)
            )

    return K1 + K2


def f(d=None, T=None):
    """Free energy using equation (6)"""
    return ue(d=d) + un(d=d) + ft(d=d, T=T)


def dfedd(d=None, T=None):
    """Electronic ground state contribution obtained by differentiation of
    equation (9)
    """
    return a1 + 2 * a2 * d + 3 * a3 * d ** 2 + a4 / d + 2 * a5 * np.log(d) / d


def dfndd(d=None, T=None):
    """Nuclear ground state contribution obtained by differentiation of
    equation (15)
    """
    return b1 + 2 * b2 * d - 10 * b4 * d ** 9 * b3 * np.exp(-b4 * d ** 10)


def dftdd(d=None, T=None):
    """First density derivative of thermal contribution"""
    K1 = 0.0
    K2 = 0.0
    for i in range(1):
        Ti = 700.0 * 2.0 ** i
        for k in range(3):
            kk = k - 4
            # print (i, kk, TABLE2_alpha[k])
            K1 += (
                TABLE2_alpha[k]
                * T
                * (3.0 * np.log(1.0 - np.exp(-Ti / T)) - D(z=Ti / T, n=3))
                * d ** (kk / k0 - 1.0)
                * kk
                / k0
            )

    for j in range(4):
        Tj = 1000.0 * 2.0 ** j
        for k in range(3):
            kk = k - 4
            # print (j, kk, TABLE2_gamma[j][k])
            K2 += (
                TABLE2_gamma[j][k]
                * T
                * np.log(1.0 - np.exp(-Tj / T))
                * d ** (kk / k0 - 1.0)
                * kk
                / k0
            )

    return K1 + K2


def dfeddd(d=None, T=None):
    """Second density derivative of electronic contributiion"""
    return 2 * a2 + 6 * a3 * d - a4 / d ** 2 - 2 * a5 / d ** 2 * (np.log(d) - 1.0)


def dfnddd(d=None, T=None):
    """Second density derivative of nuclear contribution"""
    return 2 * b2 - 10.0 * b3 * b4 * d ** 8 * np.exp(-b4 * d ** 10) * (
        9.0 - 10.0 * b4 * d ** 10
    )


def dftddd(d=None, T=None):
    """Second density derivative of thermal contribution"""
    K1 = 0.0
    K2 = 0.0

    for i in range(1):
        Ti = 700.0 * 2.0 ** i
        for k in range(3):
            kk = k - 4
            factor = kk / k0 * (kk / k0 - 1.0) * d ** (kk / k0 - 2.0)
            # print (i, kk, TABLE2_alpha[k])
            K1 += (
                TABLE2_alpha[k]
                * T
                * (3.0 * np.log(1.0 - np.exp(-Ti / T)) - D(z=Ti / T, n=3))
                * factor
            )

    for j in range(4):
        Tj = 1000.0 * 2.0 ** j
        for k in range(3):
            kk = k - 4
            factor = kk / k0 * (kk / k0 - 1.0) * d ** (kk / k0 - 2.0)
            # print (j, kk, TABLE2_gamma[j][k])
            K2 += TABLE2_gamma[j][k] * T * np.log(1.0 - np.exp(-Tj / T)) * factor

    return K1 + K2


def dftdT(d=None, T=None):
    """specific entropy obtained from equation (13)"""
    K1 = 0.0
    K2 = 0.0
    for i in range(1):
        Ti = 700.0 * 2.0 ** i
        for k in range(3):
            kk = k - 4
            # print (i, kk, TABLE2_alpha[k])
            K1 += (
                TABLE2_alpha[k]
                * (4.0 * D(Ti / T) - 3.0 * np.log(1.0 - np.exp(-Ti / T)))
                * d ** (kk / k0)
            )

    for j in range(4):
        Tj = 1000.0 * 2.0 ** j
        for k in range(3):
            kk = k - 4
            # print (j, kk, TABLE2_gamma[j][k])
            K2 += (
                TABLE2_gamma[j][k]
                * (
                    Tj / T * 1.0 / (np.exp(Tj / T) - 1.0)
                    - np.log(1.0 - np.exp(-Tj / T))
                )
                * d ** (kk / k0)
            )

    # Return specific entropy in kJ/gK
    return -(K1 + K2)


def dftdT2(d=None, T=None):
    """T-derivative of specific entropy obtained from equation (13)"""
    K1 = 0.0
    K2 = 0.0
    for i in range(1):
        Ti = 700.0 * 2.0 ** i
        for k in range(3):
            kk = k - 4
            # print (i, kk, TABLE2_alpha[k])
            K1 += (
                TABLE2_alpha[k]
                * (
                    -12.0 / T * (Ti / T / (np.exp(Ti / T) - 1.0) - D(Ti / T))
                    + 3.0 * np.exp(-Ti / T) * Ti / T ** 2
                )
                * d ** (kk / 3.0)
            )

    for j in range(4):
        Tj = 1000.0 * 2.0 ** j
        for k in range(3):
            kk = k - 4
            # print (j, kk, TABLE2_gamma[j][k])
            K2 += (
                TABLE2_gamma[j][k]
                * (
                    -Tj / T ** 2 * 1.0 / (np.exp(Tj / T) - 1.0)
                    + Tj ** 2 / T ** 3 * (np.exp(Tj / T - 1.0)) ** (-2) * np.exp(Tj / T)
                    + np.exp(-Tj / T) * Tj / T / (1.0 - np.exp(-Tj / T))
                )
                * d ** (kk / 3.0)
            )

    # Return specific entropy in kJ/gK
    return K1 + K2


def dftdT2_test(d=None, T=None):
    res = ftool.deriv(f=dftdT, whicharg="T", x0=T, d=d, acc=1.0e-5)
    return res


def dftdT_test(d=None, T=None):
    res = ftool.deriv(f=ft, whicharg="T", x0=T, d=d, acc=1.0e-5)
    return res


def dfdd(d=None, T=None):
    return dfedd(d=d, T=T) + dfndd(d=d, T=T) + dftdd(d=d, T=T)


def dfddd(d=None, T=None):
    return dfeddd(d=d, T=T) + dfnddd(d=d, T=T) + dftddd(d=d, T=T)


def s_spec(P=None, T=None, d=None):
    """Compute specific entropy in J/kg K as function
    Input units are SI
    Output unit is J/kg K
    """
    if d == None and not P == None:
        d = Density(P=P, T=T)

    # Convert density to gcc
    d *= 1.0e-3

    return -dftdT(d=d, T=T) * 1.0e6


def cv_spec(d=None, T=None, P=None):
    """Compute specific isochoric heat capacity"""
    if d == None and not P == None:
        d = Density(T=T, P=P) * 1.0e-3

    # Factor of 1.0e6 to convert from kJ/gK to J/kg K
    return -T * dftdT2_test(d=d, T=T) * 1.0e6


def cp_spec(d=None, T=None, P=None):
    """Compute specific isobaric heat capacity"""
    if d == None and not P == None:
        d = Density(T=T, P=P)

    dsdT_P = ftool.deriv(f=s_spec, whicharg="T", x0=T, P=P, d=d)

    return T * dsdT_P


def u_spec(d=None, T=None, P=None):
    """Compute specific internal energy"""
    if d == None and not P == None:
        d = Density(T=T, P=P)

    else:
        d *= 1.0e-3

    # Return in J/kg, f is in kJ/g
    return (f(d=d, T=T) - T * dftdT(d=d, T=T)) * 1.0e6


def h(d=None, T=None):
    """specific Enthalpy"""

    return (f(d=d, T=T) - T * dftdT(d=d, T=T) + d * dftdd(d=d, T=T)) * 1.0e6


def alpha_th(P=None, T=None, d=None):
    """Compute thermal expansion coefficient in 1/K"""
    if d == None and not P == None:
        d = Density(P=P, T=T)

    drhodT = ftool.deriv(f=Density, whicharg="T", x0=T, P=P, acc=1.0e-4)

    return -1.0 / d * drhodT


def dTdP_S(d=None, T=None, P=None):
    """Compute adiabatic gradient"""
    if d == None and not P == None:
        d = Density(T=T, P=P)

    alpha = alpha_th(T=T, P=P, d=d)
    c = cp_spec(d=d, T=T, P=P)

    return alpha * T / (d * c)


def dPdrho_T(d=None, T=None, P=None):
    """Compute isothermal density derivative of pressure in J/kg"""
    if d == None and not P == None:
        d = Density(T=T, P=P)

    # Convert density into gcc
    d *= 1.0e-3

    return (2 * d * dfdd(d=d, T=T) + d ** 2 * dfddd(d=d, T=T)) * 1.0e6


def Pressure(d=None, T=None):
    """Total pressure obtained by differentiation of equation (6)"""
    # Convert density input in gcc
    d *= 1.0e-3

    # Compute pressure in GPa
    ptot = d ** 2 * (dfndd(d=d) + dfedd(d=d) + dftdd(d=d, T=T))

    # Return pressure in Pa
    return ptot * 1.0e9


def Density(P=None, T=None):
    """Compute density in kg/m3"""
    # Increase efficiency and reliability by adjusting bisection range to
    # pressure region at hand
    if P > 1.0e10:
        a = 1500

    else:
        a = 950.0

    d = ftool.bisec(
        f=Pressure,
        a=a,
        b=8000.0,
        y=P,
        whicharg="d",
        T=T,
        identity="French dens",
        acc=1.0e-10,
        limit=100,
    )

    return d


def together(T=None, P=None, d=None):
    """Compute all relevant parameters for planetary modelling in one go"""
    if d == None and not P == None:
        d = Density(T=T, P=P)

    dPdrho = dPdrho_T(d=d, T=T, P=P)
    alpha = alpha_th(d=d, P=P, T=T)
    c = cp_spec(d=d, T=T, P=P)
    s = s_spec(d=d, T=T, P=P)
    u = u_spec(d=d, T=T, P=P)
    dTdP_S = max(alpha * T / (d * c), 0.0)

    return d, dTdP_S, dPdrho, alpha, c, s, u
