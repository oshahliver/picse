# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 13:01:56 2023

@author: oshahliver
"""

import numpy as np
from scipy.integrate import quad
from picse.utils.function_tools import functionTools as ftool
from matplotlib import pyplot as plt

# Fo, Wds, Rwd, Akm, Prv, PPrv, Au, MgO, Pt
U0_list = [-2200.69, -2173.0, -2163.3, -1472.5, -1524.0, -1447.9, 0.0, -609.571, 0]
V0_list = [43.67, 40.54, 39.5, 24.45, 26.35, 24.2, 10.215, 11.248, 9.091]
K0_list = [127.4, 169.0, 187.4, 252.0, 215.3, 253.7, 167.0, 160.3, 275.0]
Kp_list = [4.3, 4.14, 3.98, 4.38, 4.91, 4.0, 5.9, 4.25, 5.43]
k_list = [5, 5, 5, 5, 5, 5, 2, 5, 2]
theta01_list = [949, 921, 929, 943, 995, 943, 178, 747, 184]
theta02_list = [348, 393, 414, 417, 451, 417, 84, 399, 137]
m1_list = [10.5, 10.5, 10.5, 7.5, 7.5, 7.5, 1.5, 3, 1.5]
m2_list = [10.5, 10.5, 10.5, 7.5, 7.5, 7.5, 1.5, 3, 1.5]
gamma0_list = [1.066, 1.185, 1.21, 1.70, 2, 1.67, 2.918, 1.53, 2.77]
gamma_inf_list = [0, 0, 0, 0, 0.26, 0, 0.66, 0.624, 0.43]
beta_list = [2.225, 2.10, 1.35, 3.0, 1.41, 2.22, 2.406, 2.115, 2.26]
e0_list = [0, 0, 0, 0, 0, 0, 6.1, None, 79]
g_list = [0, 0, 0, 0, 0, 0, 0.66, None, 0.26]
a0_list = [-15.9, -15.9, -15.9, -15.9, -15.9, -15.9, None, -15.9, None]
m_list = [None, None, None, None, None, None, None, 4.48, None]
n_list = [7, None, None, None, None, None, None, 2, None]
R = 8.31446
T0 = 298.15


def F_wds(V, T):
    return get_free_energy(V, T, 0)


def get_free_energy(V, T, ll):

    kwargs = {
        "U0": U0_list[ll],
        "T0": T0,
        "V0": V0_list[ll],
        "kp": Kp_list[ll],
        "k": k_list[ll],
        "K0": K0_list[ll],
        "m1": m1_list[ll],
        "m2": m2_list[ll],
        "n": 1,
        "gamma0": gamma0_list[ll],
        "gamma_inf": gamma_inf_list[ll],
        "beta": beta_list[ll],
        "m": m_list[ll],
        "g": g_list[ll],
        "e0": e0_list[ll],
        "a0": a0_list[ll],
        "theta01": theta01_list[ll],
        "theta02": theta02_list[ll],
    }

    return F(V, T, **kwargs)


def F(
    V,
    T,
    U0=None,
    T0=None,
    V0=None,
    kp=None,
    k=None,
    K0=None,
    m1=None,
    m2=None,
    n=None,
    a0=None,
    m=None,
    beta=None,
    gamma0=None,
    gamma_inf=None,
    theta01=None,
    theta02=None,
    g=None,
    e0=None,
):
    theta1 = theta(
        V, theta0=theta01, gamma_inf=gamma_inf, gamma0=gamma0, beta=beta, V0=V0
    )
    theta2 = theta(
        V, theta0=theta02, gamma_inf=gamma_inf, gamma0=gamma0, beta=beta, V0=V0
    )
    return (
        U0
        + E0(V, V0=V0, kp=kp, k=k, K0=K0)
        + F_th(V, T, m1=m1, m2=m2, theta1=theta1, theta2=theta2)
        - F_th(V, T0, m1=m1, m2=m2, theta1=theta1, theta2=theta2)
        + F_anh(V, T, n=n, a0=a0, m=m, V0=V0)
        - F_anh(V, T0, n=n, a0=a0, m=m, V0=V0)
        + F_e(V, T, V0=V0, g=g, n=n, e0=e0)
        - F_e(V, T0, V0=V0, g=g, n=n, e0=e0)
    )


def E0(V, V0=None, kp=None, k=None, K0=None):
    def integrand(V):
        xx = (V / V0) ** (1 / 3)
        eta = 3 * kp / 2 - k + 0.5
        return 3 * K0 * xx ** -k * (1 - xx) * np.exp(eta * (1 - xx))

    x0 = 1e-7
    x1 = V
    result, error = quad(integrand, x0, x1)
    return result


def F_th(V, T, m1=None, m2=None, ll=None):
    m1 = m1_list[ll]
    m2 = m2_list[ll]
    theta1 = theta(
        V,
        theta0=theta01_list[ll],
        gamma_inf=gamma_inf_list[ll],
        gamma0=gamma0_list[ll],
        V0=V0_list[ll],
        beta=beta_list[ll],
    )
    theta2 = theta(
        V,
        theta0=theta02_list[ll],
        gamma_inf=gamma_inf_list[ll],
        gamma0=gamma0_list[ll],
        V0=V0_list[ll],
        beta=beta_list[ll],
    )

    a1 = m1 * R * T * np.log(1 - np.exp(-theta1 / T))
    a2 = m2 * R * T * np.log(1 - np.exp(-theta2 / T))
    return a1 + a2


def F_anh(V, T, n=None, a0=None, m=None, V0=None):
    x = V / V0
    return -3 / 2 * n * R * a0 * x ** m * T ** 2


def F_e(V, T, V0=None, g=None, n=None, e0=None):
    x = V / V0
    return -3 / 2 * n * R * e0 * x ** g * T ** 2


def gamma(V, V0=None, beta=None, gamma0=None, gamma_inf=None):
    x = V / V0
    return gamma_inf + (gamma0 - gamma_inf) * x ** beta


def theta(V, theta0=None, gamma_inf=None, gamma0=None, beta=None, V0=None):
    x = V / V0
    return (
        theta0
        * x ** (-gamma_inf)
        * np.exp((gamma0 - gamma_inf) / beta * (1 - x ** beta))
    )


def P0(V, ll):
    """
    Compute pressure at room temperature in GPa

    Parameters
    ----------
    V : TYPE
        DESCRIPTION.
    ll : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    X = (V / V0_list[ll]) ** (1 / 3)
    eta = 3 * Kp_list[ll] / 2 - k_list[ll] + 1 / 2
    return 3 * K0_list[ll] * X ** (-k_list[ll]) * (1 - X) * np.exp(eta * (1 - X))


def P_th(V, T, ll):
    """
    Compute thermal pressure in MPa

    Parameters
    ----------
    V : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    ll : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    gam = gamma(
        V,
        V0=V0_list[ll],
        beta=beta_list[ll],
        gamma0=gamma0_list[ll],
        gamma_inf=gamma_inf_list[ll],
    )
    theta1 = theta(
        V,
        theta0=theta01_list[ll],
        gamma_inf=gamma_inf_list[ll],
        gamma0=gamma0_list[ll],
        V0=V0_list[ll],
        beta=beta_list[ll],
    )
    theta2 = theta(
        V,
        theta0=theta02_list[ll],
        gamma_inf=gamma_inf_list[ll],
        gamma0=gamma0_list[ll],
        V0=V0_list[ll],
        beta=beta_list[ll],
    )
    print(theta1, theta2, gam)
    print(np.exp(theta1 / T))
    a1 = m1_list[ll] * R * gam / V * (theta1 / (np.exp(theta1 / T) - 1))
    a2 = m2_list[ll] * R * gam / V * (theta2 / (np.exp(theta2 / T) - 1))

    return a1 + a2
    # return -ftool.der1(F_th, V, which="V", T=T, ll=ll)


def C_V(V, T, ll):
    theta1 = theta(
        V,
        theta0=theta01_list[ll],
        gamma_inf=gamma_inf_list[ll],
        gamma0=gamma0_list[ll],
        V0=V0_list[ll],
        beta=beta_list[ll],
    )
    theta2 = theta(
        V,
        theta0=theta02_list[ll],
        gamma_inf=gamma_inf_list[ll],
        gamma0=gamma0_list[ll],
        V0=V0_list[ll],
        beta=beta_list[ll],
    )
    a1 = (
        m1_list[ll]
        * R
        * ((theta1 / T) ** 2 * np.exp(theta1 / T) / (np.exp(theta1 / T) - 1) ** 2)
    )
    a2 = (
        m2_list[ll]
        * R
        * ((theta2 / T) ** 2 * np.exp(theta2 / T) / (np.exp(theta2 / T) - 1) ** 2)
    )
    return a1 + a2


def K_T0(V, ll):
    X = (V / V0_list[ll]) ** (1 / 3)
    eta = 3 * Kp_list[ll] / 2 - k_list[ll] + 1 / 2
    return (
        K0_list[ll]
        * X ** (-k_list[ll])
        * np.exp(eta * (1 - X))
        * (X + (1 - X) * (eta * X + k_list[ll]))
    )


def K_T(V, T, ll):
    x = V / V0_list[ll]
    gam = gamma(
        V,
        V0=V0_list[ll],
        beta=beta_list[ll],
        gamma0=gamma0_list[ll],
        gamma_inf=gamma_inf_list[ll],
    )
    pth = P_th(V, T, ll)
    q = (
        beta_list[ll]
        * x ** beta_list[ll]
        * (gamma0_list[ll] - gamma_inf_list[ll])
        / gam
    )
    K_th = pth * (1 + gam - q) - gam ** 2 * T * C_V(V, T, ll) / V
    K_Ta = P_a(V, T, 7) * (1 - m_list[7])
    print(K_Ta)
    return K_T0(V, ll) + K_th * 1e-3  # + K_Ta * 1e-3


def P_a(V, T, ll):
    x = V / V0_list[ll]
    a = a0_list[ll] * x ** m_list[ll]
    Ea = (m1_list[ll] + m2_list[ll]) / 2 * R * a * T ** 2
    return m_list[ll] / V * Ea


def pressure(V, T, ll):
    """
    Compute pressure in GPa.

    Parameters
    ----------
    V : TYPE
        DESCRIPTION.
    T : TYPE
        DESCRIPTION.
    ll : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return P0(V, ll) + P_th(V, T, ll) * 1e-3


def plot():
    temps = np.linspace(-100, 100)
    x = 1.0
    pres = np.array([P_th(x, t, 0) for t in temps])

    plt.plot(temps, pres)
