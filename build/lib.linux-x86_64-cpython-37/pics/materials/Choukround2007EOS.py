#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:58:33 2020

@author: oshah
"""

import numpy as np
from matplotlib import pyplot as plt
import time
import sys
import functionTools as ftool

c_list = [
    [4215.3, -2.2584, 3.2069e-2],  # liquid
    [2100, 0, 0],  # ice Ih
    [2500, 0, 0],  # ice III
    [2500, 0, 0],  # ice V
    [2600, 0, 0],
]  # ice VI


V0_list = [9.921e-4, 1.083e-3, 8.686e-4, 8.02e-4, 7.38e-4]
Pref_list = [0.1, 6.11657e-4, 355.0, 618.4, 2216.0]
Tref_list = [180.0, 273.16, 256.43, 272.73, 356.15]

aT_list = [
    [1.0, 0.0, 2.963e-9, 3.17],
    [1.0, -1.0e-4, -5.0e-7, 2.0],
    [1.0, -9.95e-5, -1.02e-6, 2.0],
    [1.0, -8.33e-5, -4.95e-6, 2.0],
    [1.0, -1.8e-4, 3.82e-8, 2.0],
]

aP_list = [
    [0.7518, -3.07e-5, -3.6e-19, 5.2],
    [1.0, -9.84e-5, -2.146e-7, 2.0],
    [1.0, -2.43e-4, -5.26e-7, 2.0],
    [1.0, -2.07e-5, 2.42e-8, 2.0],
    [1.0, 0.0, 0.0, 2.0],
]

P0_list = [2500.0, 6.11657e-4, 355.0, 618.4, 2216.0]


def xiT(T=None, ph=0):
    Tref = Tref_list[ph]
    aT1, aT2, aT3, aT4 = aT_list[ph]
    return aT1 + aT2 * (T - Tref) + aT3 * (T - Tref) ** aT4


def xiP(P=None, ph=0):
    Pref = Pref_list[ph]
    aP1, aP2, aP3, aP4 = aP_list[ph]
    return aP1 + aP2 * (P - Pref) + aP3 * (P - Pref) ** aP4


def dxiTdT(T=None, ph=0):
    Tref = Tref_list[ph]
    aT1, aT2, aT3, aT4 = aT_list[ph]
    return aT2 + aT3 * aT4 * (T - Tref) ** (aT4 - 1)


def dxiTdT_test(T=None, P=None, ph=None):

    return ftool.deriv(f=xiT, x0=T, whicharg="T")


def dxiPdP(P=None, ph=0):
    Pref = Pref_list[ph]
    aP1, aP2, aP3, aP4 = aP_list[ph]
    return aP2 + aP3 * aP4 * aP3 * (P - Pref) ** (aP4 - 1)


def ddxiTddT(T=None, ph=0):
    Tref = Tref_list[ph]
    aT1, aT2, aT3, aT4 = aT_list[ph]
    return aT3 * aT4 * (aT4 - 1) * (T - Tref) ** (aT4 - 2)


def ddxiPddP(P=None, ph=0):
    Pref = Pref_list[ph]
    aP1, aP2, aP3, aP4 = aP_list[ph]
    return aP3 * aP4 * (aP4 - 1) * (P - Pref) ** (aP4 - 2)


def dVdT(P=None, T=None, ph=0):
    print("xiP =", xiP(P=P, ph=ph))
    print("dxiTdT=", dxiTdT(T=T, ph=ph))
    print("dxiTdT test =", dxiTdT_test(T=T, ph=ph))
    return V0_list[ph] * xiP(P=P, ph=ph) * dxiTdT(T=T, ph=ph)


def V(T=None, P=None, ph=0):
    V0 = V0_list[ph]
    a = xiT(T=T, ph=ph)
    b = xiP(P=P, ph=ph)
    return V0 * a * b


def cp0(T=None, P=None, ph=0):
    c0, c1, c2 = c_list[ph]
    return c0 + c1 * T + c2 / T**2


def alpha(T=None, P=None, ph=0):
    # Convert pressure into MPa for input
    P *= 1.0e-6
    a = V(T=T, P=P, ph=ph)
    b = dVdT(P=P, T=T, ph=ph)
    print("aa, bb=", a, b)
    return 1.0 / a * b


def cp(T=None, P=None, ph=0):
    # Convert pressure into MPa for input
    P *= 1.0e-6
    Pref = Pref_list[ph]
    aP1, aP2, aP3, aP4 = aP_list[ph]
    V0 = V0_list[ph]

    integral = (
        aP1 * (P - Pref)
        + 0.5 * aP2 * (P - Pref) ** 2
        + 1.0 / (aP4 + 1.0) * (P - Pref) ** (aP4 + 1.0)
    )

    return cp0(T=T, ph=ph) - T * V0 * ddxiTddT(T=T, ph=ph) * integral


def Density(T=None, P=None, ph=0):
    # Convert pressure into MPa for input
    P *= 1.0e-6
    return 1.0 / V(T=T, P=P, ph=ph)


def Plot(res=3, ph=0):

    N = 2**res
    pres = np.linspace(1.0e7, 1.0e8, N)
    temp = np.linspace(100, 400, N)

    xx, yy = np.meshgrid(temp, pres)

    dens = np.zeros([N, N])

    for i in range(N):
        for j in range(N):
            dens[i][j] = Density(P=pres[j], T=temp[i], ph=ph)

    fig, ax = plt.subplots()

    xmin = 100
    xmax = 300
    ymin = 1.0e7
    ymax = 1.0e8

    aspect = (xmax - xmin) / (ymax - ymin)

    im = ax.imshow(dens.T, origin="lower", cmap="Reds", extent=[xmin, xmax, ymin, ymax])
    ax.set_aspect(aspect)
