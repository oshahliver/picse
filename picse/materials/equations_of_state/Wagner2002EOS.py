# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 15:19:16 2018

@author: os18o068

This module incorporates the IAPWS equation of state (EOS) for pure water in
vapor and liquid phase (Wagner 2002, IAPWS 2009)
Range over which the EOS is constrained by experimental data:
    teperature: 0-273.16 K
    pressure: 0-210 MPa
"""

from matplotlib import pyplot as plt
import matplotlib.ticker
from matplotlib import rc
from picse.utils.function_tools import functionTools as ftool
from picse.materials import phase_transitions_water_Wagner2002 as phase
from picse.physicalparams import Rgas
import numpy as np
import iapws

# import phase_transitions_water_Wagner2002 as phase
Rgasspec = 0.46151805e3
hbar = 1.0545718e-34
EHartree = 4.359745e-18  # Hartree energy in joule
kB = 1.38065e-23
NA = 6.0221409e23
mFe = 55.845e-3  # molar mass in kg
mO = 15.999e-3
mH = 1.00794e-3
mMg = 24.305e-3
mS = 32.065e-3
mSi = 28.0855e-3
mH2O = 2 * mH + mO

degrees_of_freedom = 3 + 3 * 3 - 6.0
heat_ratio = 1.0 + 2.0 / degrees_of_freedom

T_triple = 273.16
P_triple = 611.655

T_critical = 647.096
P_critical = 22.064e6
rho_critical = 322.0

# ideal gas part (eq. 5 in IAPWS 2009, eq. 6.5 in Wagner 2002)
n0_list = [
    -8.32044648201,
    6.6832105268,
    3.00632,
    0.012436,
    0.97315,
    1.27950,
    0.96956,
    0.24873,
]

gamma0_list = [
    None,
    None,
    None,
    1.28728967,
    3.53734222,
    7.74073708,
    9.244378,
    27.5075105,
]

# residual part (eq. 6 in IAPWS 2009, eq. 6.6 in Wagner 2002)
nr_list = [
    0.12533547935523e-1,
    0.78957634722828e1,
    -0.87803203303561e1,
    0.31802509345418e0,
    -0.26145533859358e0,
    -0.78199751687981e-2,
    0.88089493102134e-2,
    -0.66856572307965e0,
    0.20433810950965e0,
    -0.66212605039687e-4,
    -0.19232721156002e0,
    -0.25709043003438e0,
    0.16074868486251e0,
    -0.40092828925907e-1,
    0.39343422603254e-6,
    -0.75941377088144e-5,
    0.56250979351888e-3,
    -0.15608652257135e-4,
    0.11537996422951e-8,
    0.36582165144204e-6,
    -0.13251180074668e-11,
    -0.62639586912454e-9,
    -0.10793600908932e0,
    0.17611491008752e-1,
    0.22132295167546e0,
    -0.40247669763528e0,
    0.58083399985759e0,
    0.49969146990806e-2,
    -0.31358700712549e-1,
    -0.74315929710341e0,
    0.47807329915480e0,
    0.20527940895948e-1,
    -0.13636435110343e0,
    0.14180634400617e-1,
    0.83326504880713e-2,
    -0.29052336009585e-1,
    0.38615085574206e-1,
    -0.20393486513704e-1,
    -0.16554050063734e-2,
    0.19955571979541e-2,
    0.15870308324157e-3,
    -0.16388568342530e-4,
    0.43613615723811e-1,
    0.34994005463765e-1,
    -0.76788197944621e-1,
    0.22446277332006e-1,
    -0.62689710414685e-4,
    -0.55711118565645e-9,
    -0.19905718354408e0,
    0.31777497330738e0,
    -0.11841182425981e0,
    -0.31306260323435e2,
    0.31546140237781e2,
    -0.25213154341695e4,
    -0.14874640856724e0,
    0.31806110878444e0,
]

ar_list = [3.5, 3.5]

br_list = [0.85, 0.95]

aar_list = [0.32, 0.32]

bbr_list = [0.2, 0.2]

ccr_list = [28, 32]

ddr_list = [700, 800]

cr_list = [
    None,
    None,
    None,
    None,
    None,
    None,
    None,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    3,
    4,
    6,
    6,
    6,
    6,
    None,
    None,
    None,
]

dr_list = [
    1,
    1,
    1,
    2,
    2,
    3,
    4,
    1,
    1,
    1,
    2,
    2,
    3,
    4,
    4,
    5,
    7,
    9,
    10,
    11,
    13,
    15,
    1,
    2,
    2,
    2,
    3,
    4,
    4,
    4,
    5,
    6,
    6,
    7,
    9,
    9,
    9,
    9,
    9,
    10,
    10,
    12,
    3,
    4,
    4,
    5,
    14,
    3,
    6,
    6,
    6,
    3,
    3,
    3,
]

tr_list = [
    -0.5,
    0.875,
    1.0,
    0.5,
    0.75,
    0.375,
    1,
    4,
    6,
    12,
    1,
    5,
    4,
    2,
    13,
    9,
    3,
    4,
    11,
    4,
    13,
    1,
    7,
    1,
    9,
    10,
    10,
    3,
    7,
    10,
    10,
    6,
    10,
    10,
    1,
    2,
    3,
    4,
    8,
    6,
    9,
    8,
    16,
    22,
    23,
    23,
    10,
    50,
    44,
    46,
    50,
    0,
    1,
    4,
]

alphar_list = [20, 20, 20]

betar_list = [150, 150, 250, 0.3, 0.3]

gammar_list = [1.21, 1.21, 1.25]

epsilonr_list = [1, 1, 1]


def Delta(the, B, delta, a):
    return the ** 2 + B * ((delta - 1.0) ** 2) ** a


def theta(tau, A, delta, beta):
    return (1.0 - tau) + A * ((delta - 1.0) ** 2) ** (1.0 / (2.0 * beta))


def psi(tau, C, D, delta):
    return np.exp(-C * (delta - 1.0) ** 2 - D * (tau - 1.0) ** 2)


def dDeltadtau(the, B, delta, a):
    return -2 * the


def dDeltaddelta(the, A, B, beta, delta, a):
    return 2 * the * A / beta * (delta - 1) ** (1.0 / beta - 1.0) + 2 * a * B * (
        delta - 1.0
    ) ** (2.0 * a - 1)


def dDeltaddelta_2(the, A, B, beta, delta, a):
    return (
        2.0 * A ** 2.0 / beta ** 2.0 * (delta - 1.0) ** (2.0 / beta - 2.0)
        + 2.0
        * the
        * A
        / beta
        * (1.0 / beta - 1.0)
        * (delta - 1.0) ** (1.0 / beta - 2.0)
        + 2.0 * a * B * (2.0 * a - 1.0) * (delta - 1.0) ** (2.0 * a - 2.0)
    )


def dpsidtau(tau, C, D, delta):
    return -2.0 * D * (tau - 1.0) * psi(tau, C, D, delta)


def dpsiddeltadtau(tau, C, D, delta):
    return 4 * C * D * (delta - 1) * (tau - 1) * psi(tau, C, D, delta)


def dDeltaddeltadtau(the, A, B, beta, delta):
    return -2.0 * A / beta * (delta - 1.0) ** (1.0 / beta - 1.0)


def dpsiddelta(tau, C, D, delta):
    return -2.0 * C * (delta - 1.0) * psi(tau, C, D, delta)


def dpsiddelta_2(tau, C, D, delta):
    return -2 * C * psi(tau, C, D, delta) * (1 - (delta - 1) ** 2 * 2 * C)


def dDeltadtau_2():
    return 2.0


def dpsidtau_2(tau, C, D, delta):
    return 2 * D * psi(tau, C, D, delta) * ((tau - 1) ** 2 * 2 * D - 1)


def phi0(d=None, t=None, **kwargs):
    delta = d / rho_critical
    tau = T_critical / t
    n1 = n0_list[0]
    n2 = n0_list[1]
    n3 = n0_list[2]
    k = 0

    for i in range(5):
        k += n0_list[i + 3] * np.log(1.0 - np.exp(-gamma0_list[i + 3] * tau))

    return np.log(delta) + n1 + n2 * tau + n3 * np.log(tau) + k


def dphi0dtau(d=None, T=None):
    delta = d / rho_critical
    tau = T_critical / T
    n1 = n0_list[0]
    n2 = n0_list[1]
    n3 = n0_list[2]
    k = 0
    for i in range(5):
        k += (
            n0_list[i + 3]
            * gamma0_list[i + 3]
            * (1.0 / (1.0 - np.exp(-gamma0_list[i + 3] * tau)) - 1.0)
        )

    return n2 + n3 / tau + k


def dphi0ddelta(d=None, T=None):
    return rho_critical / d


def dphi0dtau_2(d=None, T=None):
    delta = d / rho_critical
    tau = T_critical / T
    n1 = n0_list[0]
    n2 = n0_list[1]
    n3 = n0_list[2]
    k = 0
    for i in range(5):
        k += (
            n0_list[i + 3]
            * gamma0_list[i + 3] ** 2
            * np.exp(-gamma0_list[i + 3] * tau)
            * (1 - np.exp(-gamma0_list[i + 3] * tau)) ** (-2)
        )

    return -n3 / tau ** 2 - k


def phir(d=None, t=None, **kwargs):
    delta = d / rho_critical
    tau = T_critical / t
    K1 = 0.0
    for i in range(7):
        K1 += nr_list[i] * delta ** (dr_list[i]) * tau ** (tr_list[i])

    K2 = 0.0
    for i in range(51 - 8 + 1):
        ii = i + 7
        K2 += (
            nr_list[ii]
            * delta ** (dr_list[ii])
            * tau ** (tr_list[ii])
            * np.exp(-(delta ** (cr_list[ii])))
        )

    K3 = 0.0
    for i in range(54 - 52 + 1):
        ii = i + 51
        K3 += (
            nr_list[ii]
            * delta ** (dr_list[ii])
            * tau ** (tr_list[ii])
            * np.exp(
                -alphar_list[i] * (delta - epsilonr_list[i]) ** 2
                - betar_list[i] * (tau - gammar_list[i]) ** 2
            )
        )

    theta1 = theta(tau, 0.32, delta, 0.3)
    theta2 = theta(tau, 0.32, delta, 0.3)
    K4 = nr_list[54] * Delta(theta1, 0.2, delta, 3.5) ** (br_list[0]) * delta * psi(
        tau, 28, 700, delta
    ) + nr_list[55] * Delta(theta2, 0.2, delta, 3.5) ** (br_list[1]) * delta * psi(
        tau, 32, 800, delta
    )

    return K1 + K2 + K3 + K4


def dphirdtau(d=None, T=None):
    delta = d / rho_critical
    tau = T_critical / T
    K1 = 0.0
    for i in range(7):
        K1 += (
            nr_list[i] * tr_list[i] * delta ** (dr_list[i]) * tau ** (tr_list[i] - 1.0)
        )

    K2 = 0.0
    for i in range(51 - 8 + 1):
        ii = i + 7
        K2 += (
            nr_list[ii]
            * tr_list[ii]
            * delta ** (dr_list[ii])
            * tau ** (tr_list[ii] - 1.0)
            * np.exp(-(delta ** (cr_list[ii])))
        )

    K3 = 0.0
    for i in range(54 - 52 + 1):
        ii = i + 51
        K3 += (
            nr_list[ii]
            * delta ** (dr_list[ii])
            * tau ** (tr_list[ii])
            * np.exp(
                -alphar_list[i] * (delta - epsilonr_list[i]) ** 2
                - betar_list[i] * (tau - gammar_list[i]) ** 2
            )
            * (tr_list[ii] / tau - 2.0 * betar_list[i] * (tau - gammar_list[i]))
        )

    theta1 = theta(tau, 0.32, delta, 0.3)
    theta2 = theta(tau, 0.32, delta, 0.3)
    Delta1 = Delta(theta1, 0.2, delta, 3.5)
    Delta2 = Delta(theta2, 0.2, delta, 3.5)
    psi1 = psi(tau, 28, 700, delta)
    psi2 = psi(tau, 32, 800, delta)

    K4 = nr_list[54] * delta * (
        Delta1 ** (br_list[0] - 1) * psi1 * dDeltadtau(theta1, 0.2, delta, 3.5)
        + Delta1 ** (br_list[0]) * dpsidtau(tau, 28, 700, delta)
    ) + nr_list[55] * delta * (
        Delta2 ** (br_list[1] - 1) * psi2 * dDeltadtau(theta2, 0.2, delta, 3.5)
        + Delta2 ** (br_list[1]) * dpsidtau(tau, 32, 800, delta)
    )

    return np.real(K1 + K2 + K3 + K4)


def dphirddelta(d=None, T=None):
    delta = d / rho_critical
    tau = T_critical / T
    K1 = 0.0
    for i in range(7):
        K1 += nr_list[i] * dr_list[i] * delta ** (dr_list[i] - 1) * tau ** (tr_list[i])

    K2 = 0.0
    for i in range(51 - 8 + 1):
        ii = i + 7
        K2 += (
            nr_list[ii]
            * np.exp(-(delta ** cr_list[ii]))
            * (
                delta ** (dr_list[ii] - 1)
                * tau ** tr_list[ii]
                * (dr_list[ii] - cr_list[ii] * delta ** cr_list[ii])
            )
        )

    K3 = 0.0
    for i in range(54 - 52 + 1):
        ii = i + 51
        K3 += (
            nr_list[ii]
            * delta ** (dr_list[ii])
            * tau ** (tr_list[ii])
            * np.exp(
                -alphar_list[i] * (delta - epsilonr_list[i]) ** 2
                - betar_list[i] * (tau - gammar_list[i]) ** 2
            )
            * (dr_list[ii] / delta - 2.0 * alphar_list[i] * (delta - epsilonr_list[i]))
        )

    theta1 = theta(tau, 0.32, delta, 0.3)
    theta2 = theta(tau, 0.32, delta, 0.3)
    Delta1 = Delta(theta1, 0.2, delta, 3.5)
    Delta2 = Delta(theta2, 0.2, delta, 3.5)
    psi1 = psi(tau, 28, 700, delta)
    psi2 = psi(tau, 32, 800, delta)

    dpsiddelta1 = dpsiddelta(tau, 28, 700, delta)
    dpsiddelta2 = dpsiddelta(tau, 32, 800, delta)

    dDeltaddelta1 = dDeltaddelta(theta1, 0.32, 0.2, 0.3, delta, 3.5)
    dDeltaddelta2 = dDeltaddelta(theta2, 0.32, 0.2, 0.3, delta, 3.5)

    K4 = nr_list[54] * (
        Delta1 ** br_list[0] * (psi1 + delta * dpsiddelta1)
        + Delta1 ** (br_list[0] - 1) * dDeltaddelta1 * delta * psi1
    ) + nr_list[55] * (
        Delta2 ** br_list[1] * (psi2 + delta * dpsiddelta2)
        + Delta2 ** (br_list[1] - 1) * dDeltaddelta2 * delta * psi2
    )

    return np.real(K1 + K2 + K3 + K4)


def dphirddelta_2(d=None, T=None):
    delta = d / rho_critical
    tau = T_critical / T
    K1 = 0.0
    for i in range(7):
        K1 += (
            nr_list[i]
            * dr_list[i]
            * (dr_list[i] - 1)
            * delta ** (dr_list[i] - 2)
            * tau ** (tr_list[i])
        )

    K2 = 0.0
    for i in range(51 - 8 + 1):
        ii = i + 7
        K2 += (
            nr_list[ii]
            * np.exp(-(delta ** cr_list[ii]))
            * (
                delta ** (dr_list[ii] - 2)
                * tau ** tr_list[ii]
                * (
                    (dr_list[ii] - cr_list[ii] * delta ** cr_list[ii])
                    * (dr_list[ii] - 1 - cr_list[ii] * delta ** cr_list[ii])
                    - cr_list[ii] ** 2 * delta ** cr_list[ii]
                )
            )
        )

    K3 = 0.0
    for i in range(54 - 52 + 1):
        ii = i + 51
        K3 += (
            nr_list[ii]
            * tau ** (tr_list[ii])
            * np.exp(
                -alphar_list[i] * (delta - epsilonr_list[i]) ** 2
                - betar_list[i] * (tau - gammar_list[i]) ** 2
            )
            * (
                -2 * alphar_list[i] * delta ** dr_list[ii]
                + 4
                * alphar_list[i] ** 2
                * delta ** dr_list[ii]
                * (delta - epsilonr_list[i]) ** 2
                - 4
                * dr_list[ii]
                * alphar_list[i]
                * delta ** (dr_list[ii] - 1)
                * (delta - epsilonr_list[i])
                + dr_list[ii] * (dr_list[ii] - 1) * delta ** (dr_list[ii] - 2)
            )
        )

    theta1 = theta(tau, 0.32, delta, 0.3)
    theta2 = theta(tau, 0.32, delta, 0.3)
    Delta1 = Delta(theta1, 0.2, delta, 3.5)
    Delta2 = Delta(theta2, 0.2, delta, 3.5)
    psi1 = psi(tau, 28, 700, delta)
    psi2 = psi(tau, 32, 800, delta)
    dpsiddelta1 = dpsiddelta(tau, 28, 700, delta)
    dpsiddelta2 = dpsiddelta(tau, 32, 800, delta)
    dDeltaddelta1 = dDeltaddelta(theta1, 0.32, 0.2, 0.3, delta, 3.5)
    dDeltaddelta2 = dDeltaddelta(theta2, 0.32, 0.2, 0.3, delta, 3.5)

    dpsiddelta11 = dpsiddelta_2(tau, 28, 700, delta)
    dpsiddelta22 = dpsiddelta_2(tau, 32, 800, delta)

    dDeltaddelta11 = dDeltaddelta_2(theta1, 0.32, 0.2, 0.3, delta, 3.5)
    dDeltaddelta22 = dDeltaddelta_2(theta2, 0.32, 0.2, 0.3, delta, 3.5)

    K4 = nr_list[54] * (
        Delta1 ** br_list[0] * (2 * dpsiddelta1 + delta * dpsiddelta11)
        + 2
        * dDeltaddelta1
        * br_list[0]
        * br_list[0]
        * Delta1 ** (br_list[0] - 1)
        * (psi1 + delta * dpsiddelta1)
        + delta
        * psi1
        * (
            br_list[0]
            * (br_list[0] - 1)
            * Delta1 ** (br_list[0] - 2)
            * dDeltaddelta1 ** 2
            + br_list[0] * Delta1 ** (br_list[0] - 1) * dDeltaddelta11
        )
    ) + nr_list[55] * (
        Delta2 ** br_list[1] * (2 * dpsiddelta2 + delta * dpsiddelta22)
        + 2
        * dDeltaddelta2
        * br_list[1]
        * br_list[1]
        * Delta2 ** (br_list[1] - 1)
        * (psi2 + delta * dpsiddelta2)
        + delta
        * psi2
        * (
            br_list[1]
            * (br_list[1] - 1)
            * Delta2 ** (br_list[1] - 2)
            * dDeltaddelta2 ** 2
            + br_list[1] * Delta2 ** (br_list[1] - 1) * dDeltaddelta22
        )
    )

    return np.real(K1 + K2 + K3 + K4)


def dphirdtau_2(d=None, T=None):
    delta = d / rho_critical
    tau = T_critical / T
    K1 = 0.0
    for i in range(7):
        K1 += (
            nr_list[i]
            * tr_list[i]
            * (tr_list[i] - 1)
            * delta ** dr_list[i]
            * tau ** (tr_list[i] - 2)
        )

    K2 = 0.0
    for i in range(51 - 8 + 1):
        ii = i + 7
        K2 += (
            nr_list[ii]
            * tr_list[ii]
            * (tr_list[ii] - 1)
            * delta ** dr_list[ii]
            * tau ** (tr_list[ii] - 2)
            * np.exp(-(delta ** cr_list[ii]))
        )

    K3 = 0.0
    for i in range(54 - 52 + 1):
        ii = i + 51
        K3 += (
            nr_list[ii]
            * delta ** (dr_list[ii])
            * tau ** (tr_list[ii])
            * np.exp(
                -alphar_list[i] * (delta - epsilonr_list[i]) ** 2
                - betar_list[i] * (tau - gammar_list[i]) ** 2
            )
            * (
                (tr_list[ii] / tau - 2 * betar_list[i] * (tau - gammar_list[i])) ** 2
                - tr_list[ii] / tau ** 2
                - 2 * betar_list[i]
            )
        )

    theta1 = theta(tau, 0.32, delta, 0.3)
    theta2 = theta(tau, 0.32, delta, 0.3)
    Delta1 = Delta(theta1, 0.2, delta, 3.5)
    Delta2 = Delta(theta2, 0.2, delta, 3.5)
    psi1 = psi(tau, 28, 700, delta)
    psi2 = psi(tau, 32, 800, delta)

    K4 = nr_list[54] * delta * (
        (
            2 * br_list[0] * Delta1 ** (br_list[0] - 1)
            + 4
            * theta1 ** 2
            * br_list[0]
            * (br_list[0] - 1)
            * Delta1 ** (br_list[0] - 2)
        )
        * psi1
        + 4
        * theta1
        * br_list[0]
        * Delta1 ** (br_list[0] - 1)
        * 2
        * 700
        * (tau - 1)
        * psi1
        + Delta1 ** br_list[0] * (2 * 700 * (tau - 1) ** 2 - 1) * 2 * 700 * psi1
    ) + nr_list[55] * delta * (
        (
            2 * br_list[1] * Delta2 ** (br_list[1] - 1)
            + 4
            * theta2 ** 2
            * br_list[1]
            * (br_list[1] - 1)
            * Delta2 ** (br_list[1] - 2)
        )
        * psi2
        + 4
        * theta2
        * br_list[1]
        * Delta2 ** (br_list[1] - 1)
        * 2
        * 800
        * (tau - 1)
        * psi2
        + Delta2 ** br_list[1] * (2 * 800 * (tau - 1) ** 2 - 1) * 2 * 800 * psi2
    )

    return np.real(K1 + K2 + K3 + K4)


def dphirddeltadtau(d=None, T=None):
    delta = d / rho_critical
    tau = T_critical / T
    K1 = 0.0
    for i in range(7):
        K1 += (
            nr_list[i]
            * dr_list[i]
            * tr_list[i]
            * delta ** (dr_list[i] - 1)
            * tau ** (tr_list[i] - 1)
        )

    K2 = 0.0
    for i in range(51 - 8 + 1):
        ii = i + 7
        K2 += (
            nr_list[ii]
            * tr_list[ii]
            * delta ** (dr_list[ii] - 1)
            * tau ** (tr_list[ii] - 1)
            * (dr_list[ii] - cr_list[ii] * delta ** cr_list[ii])
            * np.exp(-(delta ** cr_list[ii]))
        )

    K3 = 0.0
    for i in range(54 - 52 + 1):
        ii = i + 51
        K3 += (
            nr_list[ii]
            * delta ** (dr_list[ii])
            * tau ** (tr_list[ii])
            * np.exp(
                -alphar_list[i] * (delta - epsilonr_list[i]) ** 2
                - betar_list[i] * (tau - gammar_list[i]) ** 2
            )
            * (
                dr_list[ii] / delta
                - 2
                * alphar_list[i]
                * (delta - epsilonr_list[i])
                * (tr_list[ii] / tau - 2 * betar_list[i] * (tau - gammar_list[i]))
            )
        )

    theta1 = theta(tau, 0.32, delta, 0.3)
    theta2 = theta(tau, 0.32, delta, 0.3)

    Delta1 = Delta(theta1, 0.2, delta, 3.5)
    Delta2 = Delta(theta2, 0.2, delta, 3.5)

    psi1 = psi(tau, 28, 700, delta)
    psi2 = psi(tau, 32, 800, delta)

    dDeltaddelta1 = dDeltaddelta(theta1, 0.32, 0.2, 0.3, delta, 3.5)
    dDeltaddelta2 = dDeltaddelta(theta2, 0.32, 0.2, 0.3, delta, 3.5)

    K4 = nr_list[54] * (
        Delta1 ** br_list[0]
        * (
            -2 * 700 * (tau - 1) * psi1
            + delta * 4 * 28 * 700 * (delta - 1) * (tau - 1) * psi1
        )
        - delta
        * br_list[0]
        * Delta1 ** (br_list[0] - 1)
        * dDeltaddelta1
        * 2
        * 700
        * (tau - 1)
        * psi1
        - 2
        * theta1
        * br_list[0]
        * Delta1 ** (br_list[0] - 1)
        * (psi1 - delta * 2 * 700 * (delta - 1) * psi1)
        + delta
        * psi1
        * 0.32
        * br_list[0]
        * 2
        / betar_list[0]
        * Delta1 ** (br_list[0] - 1)
        * (delta - 1)
        * (
            (delta - 1) ** (1 / betar_list[0] - 2)
            - 2
            * theta1
            * br_list[0]
            * (br_list[0] - 1)
            * Delta1 ** (br_list[0] - 1)
            * dDeltaddelta1
        )
    ) + nr_list[55] * (
        Delta2 ** br_list[1]
        * (
            -2 * 800 * (tau - 1) * psi2
            + delta * 4 * 32 * 800 * (delta - 1) * (tau - 1) * psi2
        )
        - delta
        * br_list[1]
        * Delta2 ** (br_list[1] - 1)
        * dDeltaddelta2
        * 2
        * 800
        * (tau - 1)
        * psi2
        - 2
        * theta2
        * br_list[1]
        * Delta2 ** (br_list[1] - 1)
        * (psi2 - delta * 2 * 800 * (delta - 1) * psi2)
        + delta
        * psi2
        * 0.32
        * br_list[1]
        * 2
        / betar_list[1]
        * Delta2 ** (br_list[1] - 1)
        * (delta - 1)
        * (
            (delta - 1) ** (1 / betar_list[1] - 2)
            - 2
            * theta2
            * br_list[1]
            * (br_list[1] - 1)
            * Delta2 ** (br_list[1] - 1)
            * dDeltaddelta2
        )
    )

    return np.real(K1 + K2 + K3 + K4)


def f(d=None, T=None, P=None):
    """Specific Helmholtz free energy in J/kg"""
    if d == None and not P == None:
        d = density(P=P, T=T)

    return Rgasspec * T * (phi0(d=d, t=T) + phir(d=d, t=T))


def s_spec(T=None, d=None, P=None):
    """Specific entropy"""
    if d == None and not P == None:
        d = density(P=P, T=T)

    # s = ftool.deriv(f=f, whicharg='t', x0=T, d=d, acc=1.0e-4)
    tau = T_critical / T
    s = Rgasspec * (
        tau * (dphi0dtau(d=d, T=T) + dphirdtau(d=d, T=T))
        - phi0(d=d, t=T)
        - phir(d=d, t=T)
    )
    return s


def h_spec(T=None, d=None, P=None):
    """Specific enthalpy"""
    if d == None and not P == None:
        d = density(P=P, T=T)

    delta = d / rho_critical
    tau = T_critical / T
    h = (
        Rgasspec
        * T
        * (
            1
            + tau * (dphi0dtau(d=d, T=T) + dphirdtau(d=d, T=T))
            + delta * dphirddelta(d=d, T=T)
        )
    )
    return h


def cV_spec(T=None, d=None, P=None):
    """Specific isochoric heat capacity"""
    if d == None and not P == None:
        d = density(P=P, T=T)

    tau = T_critical / T
    return -Rgasspec * tau ** 2 * (dphi0dtau_2(T=T, d=d) + dphirdtau_2(d=d, T=T))


def cp_spec(T=None, d=None, P=None):
    """Specific isobaric heat capacity"""
    if d == None and not P == None:
        d = density(P=P, T=T)

    # Ideal gas law as IAPWS gives spurious behaviour at low P
    if P < 1.0e3:
        cp = heat_ratio / (heat_ratio - 1.0) * Rgas / (mO + 2.0 * mH)

    elif P < 1.0e4 and T > 800.0:
        cp = heat_ratio / (heat_ratio - 1.0) * Rgas / (mO + 2.0 * mH)

    else:
        delta = d / rho_critical
        tau = T_critical / T
        drdelta = dphirddelta(d=d, T=T)

        cp = Rgasspec * (
            -(tau ** 2) * (dphi0dtau_2(d=d, T=T) + dphirdtau_2(d=d, T=T))
            + (1 + delta * drdelta - delta * tau * dphirddeltadtau(d=d, T=T)) ** 2
            / (1 + 2 * delta * drdelta + delta ** 2 * dphirddelta_2(d=d, T=T))
        )

    # cp = ftool.deriv(f=h_spec, whicharg='T', x0=T, d=d, acc=1.0e-6, plot=True)
    return cp


def u_spec(T=None, d=None, P=None):
    """Specific internal energy"""
    if d == None and not P == None:
        d = density(P=P, T=T)

    tau = T_critical / T
    u = Rgasspec * T * tau * (dphi0dtau(d=d, T=T) + dphirdtau(d=d, T=T))
    return u


def pressure(d=None, T=None, *args, **kwargs):
    # return d**2 * ftool.deriv(x0 = d, t = t, whicharg = 'd', f = f,
    #                          identity = 'wagner.pres', acc=1.0e-2, **kwargs)
    delta = d / rho_critical
    P = d * Rgasspec * T * (1.0 + delta * dphirddelta(d=d, T=T))
    return P


def dPdrho_T(T=None, P=None, d=None):
    if d == None and not P == None:
        d = density(P=P, T=T)

    delta = d / rho_critical

    K1 = dphirddelta(d=d, T=T)
    K2 = dphirddelta_2(d=d, T=T)

    return Rgasspec * T * (1.0 + delta * K2) + d * Rgasspec * T * (
        1.0 / rho_critical * K1 + delta / rho_critical * K2
    )


def alpha_th_p(T=None, P=None, d=None):
    """Compute isobaric thermal expansion coefficient in 1/K"""
    # Ideal gas law as IAPWS gives spurious behaviour at low P
    if P < 1.0e3:
        alfav = 1.0 / T

    elif P < 1.0e4 and T > 800.0:
        alfav = 1.0 / T

    else:
        alfav = iapws.IAPWS95(T=T, P=P * 1.0e-6).alfav

        if alfav == None:
            alfav = iapws.IAPWS95(T=T, P=P * 1.0e-6 * 0.99).alfav

    # Note: computing the thermal expansion via it's definition and using
    # drho/dT_P from IAPWS gives a factor 1000 smaller results but the numeric
    # values are in very good agreement. There seems to be a mistake in the
    # units in the IAPWS routine for drhodT_P. It seems that the iapws gives
    # it in g/m3K instead of kg/m3K as claimed in the web documentation
    return alfav


def dTdP_S(T=None, P=None, d=None, cp=None, alpha=None):
    if d == None:
        d = density(T=T, P=P)

    if cp == None:
        cp = cp_spec(d=d, T=T, P=P)

    if alpha == None:
        alpha = alpha_th_p(T=T, P=P)

    # print (1/d*(1.3-1)/(1.3*Rgas)*mH2O)

    return T * alpha / (d * cp)


def gamma_G(T=None, P=None, d=None):
    """Grueneisen parameter"""

    if d == None and not P == None:
        d = density(T=T, P=P)

    return alpha_th_p(T=T, P=P) * d * dPdrho_T(T=T, d=d) / cV_spec(d=d, T=T) / d


def density(P=None, T=None, **kwargs):

    # Ideal gas law as IAPWS gives spurious behaviour at low P
    if P < 1.0e3:
        d = P * (2 * mH + mO) / (kB * T * NA)

    elif P < 1.0e4 and T > 800.0:
        d = P * (2 * mH + mO) / (kB * T * NA)

    else:
        try:
            d = iapws.IAPWS95(T=T, P=P * 1.0e-6).rho

        except OverflowError:
            d = None
            print("T/P=", T, "K / ", P * 1.0e-9, "GPa")
            # sys.exit()

    return d


def together(P=None, T=None):
    try:
        d = density(P=P, T=T)

        delta = d / rho_critical
        tau = T_critical / T

        alpha = alpha_th_p(T=T, P=P, d=d)

        K1 = dphirddelta(d=d, T=T)
        K2 = dphirdtau(d=d, T=T)
        K3 = dphirddelta_2(d=d, T=T)
        K4 = dphirddeltadtau(d=d, T=T)
        K5 = dphi0dtau_2(d=d, T=T)
        K6 = dphirdtau_2(d=d, T=T)
        K7 = dphi0dtau(d=d, T=T)

        dPdrho_T = Rgasspec * T * (1.0 + delta * K1) + d * Rgasspec * T * (
            1.0 / rho_critical * K1 + delta / rho_critical * K3
        )

        cp = Rgasspec * (
            -(tau ** 2) * (K5 + K6)
            + (1.0 + delta * K1 - delta * tau * K4) ** 2
            / (1.0 + 2.0 * delta * K1 + delta ** 2 * K3)
        )

        # test = (heat_ratio-1.)/heat_ratio*T/P
        # testest = dTdP_S(T=T, P=P)
        dTdPS = T * alpha / (d * cp)
        # print (T, P, 100*(test-dTdPS)/dTdPS, testest)
        s = Rgasspec * (tau * (K7 + K2) - phi0(d=d, t=T) - phir(d=d, t=T))

        u = Rgasspec * T * tau * (dphi0dtau(d=d, T=T) + dphirdtau(d=d, T=T))

    except TypeError:
        print("Type error:")

        d = None
        dTdPS = None
        dPdrho_T = None
        alpha = None
        cp = None
        s = None
        u = None

    return d, dTdPS, dPdrho_T, alpha, cp, s, u


"""

#print ('dev = ',round((compute('dens', p = 600.0e6, t = 600.) - 1005.92)/1005.92*100,4),'%')
#print ('dev = ',round((compute('dens', p = 400.0e6, t = 1000.) - 676.11)/676.11*100,4),'%')
#print ('dev = ',round((compute('dens', p = 50.0e6, t = 300.) - 1017.85)/1017.85*100,4),'%')
#print (compute('pres', d = 997., t = 277.2)*1.0e-5)

fig, ax = plt.subplots(1, 2)
ax1, ax2 = ax

testtemp = 2000.
testpres = pressure(d = 911., T = testtemp, eps = 1.0e-6, plot = False, axis = ax1)
testdens = density(P = testpres, T = testtemp, eps = 1.0e-6, plot = False, axis = ax1)

print (density(P = testpres, T = testtemp))
print (pressure(d = 911., T = testtemp))

print ('dens =', testdens)
print ('pres =', testpres)

print ('pres dev =', round((testdens - 911.)/911. * 100, 6), '%')

ax1.scatter(testdens, testpres, color = 'k', marker = 'x', zorder = 1000)


temp_list = [300., 400., 500., 600.]

dens_list = np.logspace(np.log10(10.), np.log10(8000.), 50)

for temp in temp_list:
    pres_list = []
    for dens in dens_list:
        pres_list.append(pressure(d = dens, T = temp, acc = 1.0e-8))
    ax1.loglog(dens_list, pres_list, label = 'T = '+ str(int(temp))+' K')

dens_list = [930., 935., 940., 950., 975., 1000.]
temp_list = np.linspace(250., 800., 100)

for dens in dens_list:
    pres_list = []
    for temp in temp_list:
        p = (pressure(d = dens, T = temp))
        check = phase.phaseCheck(p, temp)
        if check == 'solid':
            pres_list.append(None)
        else:
            pres_list.append(p)
    ax2.semilogy(temp_list, pres_list, label = r'$\rho = $'+str(round(dens,2)),
                 linewidth = 2, zorder = 2)
    
tosolid_list = phase.solidLine(log = True)
tovapor_list = phase.vaporLine(log = True)

ax2.plot(tosolid_list[0], tosolid_list[1], color = 'grey', linestyle = '--', zorder = 1)
ax2.scatter(tosolid_list[2], tosolid_list[3], color =  'grey', zorder = 1)
ax2.plot(tovapor_list[0], tovapor_list[1], color = 'grey', linestyle = '--', zorder = 1)
ax2.scatter(phase.T_critical, phase.P_critical, color = 'grey', zorder = 1)

font1 = {'family': 'sans',
        'color':  'grey',
        'weight': 'bold',
        'size': 10,
        }

font2= {'family': 'sans',
        'color':  'k',
        'weight': 'normal',
        'size': 14,
        }

ax2.text(170., 1.0e6, s = 'solid', 
         color = 'grey', fontdict = font1)
ax2.text(300., 1.0e6, s = 'liquid', 
         color = 'grey', fontdict = font1)
ax2.text(550., 1.0e6, s = 'vapor', 
         color = 'grey', fontdict = font1)
ax1.text(975., 1.0e8, s = 'IAPWS', fontdict = font2)

#phi_list = [[phir(d = dens, t = temp)*1.0e6 for dens in dens_list]
 #           for temp in temp_list]

#phi_list = [[f/abs(f) * np.log10(abs(f)) for f in phi_list[i]] 
#            for i in range(len(phi_list))]
#w = 0
#for f in phi_list:
#    t = temp_list[w]
#    ax2.plot(dens_list, f, label = 'T ='+str(int(t))+' K')
#    w += 1
    
#phase.vaporLine(plot = True, axis = ax2)  

#ff = ftool.bisec(a = 10., b = 2000., y = 0., f = phir, whicharg = 'd', t = 500.,
#            noisy = False, plot = True, axis = ax2)

#print ('ff =', ff)

#ax2.plot([0., 2000.], [0, 0], color = 'grey', linestyle = ':')
#ax2.set_xlabel(r'$\rm Density \ [kg \ m^{-3}]$')

ax1.set_title(r'$Water \ ice \ isotherms$')
ax1.set_xlabel(r'$\rm Density \ [kg \ m^{-3}]$')
ax1.set_ylabel(r'$\rm Pressure \ [Pa]$')
ax1.set_ylim(2.0e2, 2.0e9)

ax2.set_title(r'$Water \ ice \ isochors$')
ax2.set_xlabel(r'$\rmTemperature \ [K]$')
ax2.set_ylabel(r'$\rm Pressure \ [Pa]$')
ax2.yaxis.set_label_position('right')
ax2.yaxis.tick_right()
ax2.set_xlim(100, 1000)
ax2.set_ylim(1.0e1, 2.0e9)

ax1.tick_params(which = 'both', axis = 'both', direction = 'in', top = 'on', 
			right = 'on')
ax1.tick_params(which = 'major', axis = 'y', length = 10)
ax1.tick_params(which = 'minor', axis = 'y', length = 5)
ax1.tick_params(which = 'major', axis = 'x', length = 5)


ax2.tick_params(which = 'both', axis = 'both', direction = 'in', top = 'on',
			left = 'on')
ax2.tick_params(which = 'major', axis = 'y', length = 10)
ax2.tick_params(which = 'minor', axis = 'y', length = 5)
ax2.tick_params(which = 'major', axis = 'x', length = 5)

ax1.set_axisbelow(True)
ax2.set_axisbelow(True)
ax1.grid(which = 'both', color = (.8, .8, .8))
ax2.grid(which = 'both', color = (.8, .8, .8))


ax1.legend()
ax2.legend()
plt.show()
"""
