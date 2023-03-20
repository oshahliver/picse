#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:55:11 2022

@author: os18o068
"""

import numpy as np
import pics.utils.function_tools.functionTools as ftool
import matplotlib.pyplot as plt
from pics.physicalparams import R_solar, M_solar, m_earth, r_earth

# Parameters for conductive part of ice sheet and silicate mantle
T0_list = [110.0, 273]
rho_list = [930.0, 3300]
k_list = [2.6, 4]
qb_list = [40e-3, 40e-3]
Qi_list = [2e12, 1.3e11]


def temp(
    r=None,
    rb=10e5,
    M=M_solar["europa"] * m_earth,
    aa=R_solar["europa"] * r_earth,
    layer=0,
):

    """Eq. 2 from Cammarano et al. 2006"""

    T0, rho, k, qb, Qi = (
        T0_list[layer],
        rho_list[layer],
        k_list[layer],
        qb_list[layer],
        Qi_list[layer],
    )
    H = Qi / M
    res = T0 + rho * H / (6 * k) * (aa ** 2 - r ** 2)
    res += (rho * H * rb ** 3 / (3 * k) - qb * rb ** 2 / k) * (1 / aa - 1 / r)
    return res


obj = "europa"
f = 1.0025

# Radius at bottom of radiative ice shell (free parameter in the model)
rb = R_solar[obj] * r_earth - 20e3
test = ftool.bisec(
    f=temp,
    whicharg="r",
    a=1e5,
    b=1e7,
    layer=0,
    y=273 * 0.95,
    noisy=False,
    M=M_solar["europa"] * m_earth,
    aa=R_solar["europa"] * r_earth * f,
    rb=rb,
    acc=1e-5,
)

print("test (km) =", (R_solar[obj] * r_earth - test) * 1e-3)
print(
    "temp =",
    temp(
        r=test,
        M=M_solar["europa"] * m_earth,
        aa=R_solar["europa"] * r_earth * f,
        rb=rb,
        layer=0,
    )
    / 0.95,
)

x = np.linspace(rb, R_solar[obj] * r_earth)
y = temp(x, rb=rb, layer=0)

plt.plot(x * 1e-3, y)

rb = R_solar[obj] * r_earth - 730e3

x = np.linspace(rb, R_solar[obj] * r_earth)
y = temp(x, rb=rb, layer=1)

plt.plot(x * 1e-3, y)
