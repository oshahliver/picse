import numpy as np
from pics.utils.function_tools import functionTools
from pics.physicalparams import (mO, mSi, mFe, mMg, mS)

def mantle_comp(P_CS, ocmf=[]):
    xi = mat2at_core(ocmf, xiH=0.0)
    T_CS = T_liquidus_pyrolite(P_CS)
    fem = Fe_number_mantle(P_CS, T_CS, xi=xi)
    sim = Si_number_mantle(P_CS, T_CS, xi=xi)

    return fem, sim

def epsilon_ki(T, k, i):
    m = [mO * 1e-3, mSi * 1e-3]
    return eki[k][i] * m[i] / 0.242 * 1873.0 / T - m[i] / 55.85 + 1.0


def partition_coefficient_i(T, P, i, xi):
    """Compute metal-silicate partition coefficient according to
    Fischer et al. 2015 accounting for non-ideal effects in multi-component
    fluids. Currently only for O and Si. The effect of S on KD is not has not
    been investigated by these authors and no other sources were found.
    """
    P *= 1e-9
    KD = a_KDi[i] + b_KDi[i] / T + c_KDi[i] * P / T
    eps = epsilon_ki(T, i, i)
    KD += eps * np.log(1.0 - xi[i]) / 2.303
    for k in range(2):
        if k == i:
            pass
        else:
            eps = epsilon_ki(T, k, i)
            delta = 1.0 / 2.303 * eps * xi[k]
            delta *= 1.0 + np.log(1.0 - xi[k]) / xi[k] - 1.0 / (1.0 - xi[i])
            KD += delta
    for k in range(2):
        if k == i:
            pass
        else:
            eps = epsilon_ki(T, k, i)
            delta = 1.0 / 2.303 * eps * xi[k] ** 2 * xi[i]
            delta *= (
                1.0 / (1 - xi[i])
                + 1.0 / (1.0 - xi[k])
                + xi[i] / (2.0 * (1.0 - xi[i]) ** 2)
                - 1.0
            )
            KD -= delta
    return 10**KD


def partition_coefficient(T, P, ll):
    """Compute metal-silicate partition coefficient according to
    Fischer et al. 2015. Currently only for O and Si
    """
    P *= 1e-9
    return 10 ** (a_KD[ll] + b_KD[ll] / T + c_KD[ll] * P / T)


def KD_S(T, P, xi, which="suer"):
    """Partition coefficient for S from Boujibar 2014 or Suer et al. 2017"""
    xiFe, xiH, xiS, xiSi, xiO = xi
    xiC = 0.0
    xiNi = 0.0
    P *= 1e-9

    xiFeO = xi_FeO(P, T, xi)
    xiSiO2 = xi_SiO2(P, T, xi)
    xiMgO = 1.0 - xiFeO - xiSiO2
    CS = 10 ** (-5.704 + 3.15 * xiFeO + 0.12 * xiMgO)

    if which == "bijou":
        b = 405.0
        c = 136.0
        d = 32.0
        e = 181.0
        f = 305.0
        g = 30.2
        h = 1.13
        i = 10.7
        j = 31.4
        k = -3.72

        KD = np.log10(xiFeO) - np.log10(CS) + b / T + c * P / T
        KD += d * np.log10(1.0 - xiSi) + e * (np.log10(1.0 - xiSi)) ** 2
        KD += f * (np.log10(1.0 - xiSi)) ** 3 + g * np.log10(1.0 - xiC)
        KD += h * np.log10(1.0 - xiFe) + i * np.log10(1.0 - xiNi)
        KD += j * np.log10(1.0 - xiO) + k

    elif which == "suer":
        KD = (
            -3.3
            + 3000 / T
            + 33 * P / T
            + np.log(xiFeO)
            - np.log(CS)
            + 14 * np.log(1.0 - xiO)
        )

    return 10**KD


def Margules_eq(P, T, A, B, C):
    """
    Implements the of the Margules equation.

    Parameters:
    P (float): Pressure in Pa.
    T (float): Temperature in K.
    A (float): First Margules parameter.
    B (flaot): Second Margules parameter.
    C (float): Third Margules parameter.

    Returns:
    float: Margules coefficient.
    """
    return A - B * T - C * P / 1e5


def W_Fe_FeO(P, T):
    """
    Interaction parameters for the activity coefficient of FeO from
    Frost et al. 2010. Includes smoothing from non-ideal to ideal regime from
    Schaefer et al. 2017.

    Parameters
    ----------
    P : FLOAT
        Pressure in Pa.
    T : FLOAT
        Temperature in K.

    Returns
    -------
    Float.

    """
    A = 83307.0
    B = 8.978
    C = 0.09
    return Margules_eq(P, T, A, B, C)


def W_FeO_Fe(P, T):
    """
    Interaction parameter for the activity coefficient of Fe from
    Frost et al. 2010. Includes smoothing from non-ideal to ideal regime from
    Schaefer et al. 2017.
    """
    A = 135943.0
    B = 31.122
    C = 0.059
    return Margules_eq(P, T, A, B, C)


def W_Fe_FeO_prime(P, T, T0, P0, a, b):
    A = 83307.0
    B = 8.978
    C = 0.09
    return -B * T_liquidus_pyrolite_prime(P) - C


def W_FeO_Fe_prime(P, T, T0, P0, a, b):
    A = 135943.0
    B = 31.122
    C = 0.059
    return -B * T_liquidus_pyrolite_prime(P) - C


def coefs_Fe_FeO(P, xiFe):
    T = T_liquidus_pyrolite(P)
    T_trans = T_liquidus_pyrolite(P_FeO_Fe_trans)

    W1 = W_FeO_Fe(P_FeO_Fe_trans, T_trans)
    W2 = W_Fe_FeO(P_FeO_Fe_trans, T_trans)
    dW1 = W_FeO_Fe_prime(P_FeO_Fe_trans, T_trans, 1940, 0.0, 29e9, 1.0 / 1.9)
    dW2 = W_Fe_FeO_prime(P_FeO_Fe_trans, T_trans, 1940, 0.0, 29e9, 1.0 / 1.9)

    K = np.exp((1 - xiFe) ** 2 * (W2 + 2.0 * (W1 - W2) * xiFe) / (Rgas * T_trans))
    K -= 1.0
    T_prime_trans = T_liquidus_pyrolite_prime(P_FeO_Fe_trans)
    g_prime = (
        1.0 / (Rgas * T_trans) * (1.0 - xiFe) ** 2 * (dW2 + 2.0 * (dW1 - dW2) * xiFe)
    )
    g_prime += (
        -1.0
        / (Rgas * T_trans**2)
        * T_prime_trans
        * (1.0 - xiFe) ** 2
        * (W2 + 2.0 * (W1 - W2) * xiFe)
    )
    g_prime *= np.exp(
        (1.0 - xiFe) ** 2 * (W2 + 2.0 * (W1 - W2) * xiFe) / (Rgas * T_trans)
    )
    g_prime = ftool.deriv(
        f=gamma_Fe, x0=P_FeO_Fe_trans, whicharg="P", xiFe=xiFe, acc=1e-3
    )
    Z = -g_prime / K
    return K, Z


def coefs_FeO_Fe(P, xiFeO):
    T = T_liquidus_pyrolite(P)
    T_trans = T_liquidus_pyrolite(P_FeO_Fe_trans)

    W1 = W_FeO_Fe(P_FeO_Fe_trans, T_trans)
    W2 = W_Fe_FeO(P_FeO_Fe_trans, T_trans)
    dW1 = W_FeO_Fe_prime(P_FeO_Fe_trans, T_trans, 1940, 0.0, 29e9, 1.0 / 1.9)
    dW2 = W_Fe_FeO_prime(P_FeO_Fe_trans, T_trans, 1940, 0.0, 29e9, 1.0 / 1.9)

    K = np.exp((1 - xiFeO) ** 2 * (W1 + 2.0 * (W2 - W1) * xiFeO) / (Rgas * T_trans))
    K -= 1.0
    T_prime_trans = T_liquidus_pyrolite_prime(P_FeO_Fe_trans)
    g_prime = (
        1.0 / (Rgas * T_trans) * (1.0 - xiFeO) ** 2 * (dW1 + 2.0 * (dW2 - dW1) * xiFeO)
    )
    g_prime += (
        -1.0
        / (Rgas * T_trans**2)
        * T_prime_trans
        * (1.0 - xiFeO) ** 2
        * (W1 + 2.0 * (W2 - W1) * xiFeO)
    )
    g_prime *= np.exp(
        (1.0 - xiFeO) ** 2 * (W1 + 2.0 * (W2 - W1) * xiFeO) / (Rgas * T_trans)
    )
    g_prime = ftool.deriv(
        f=gamma_FeO, x0=P_FeO_Fe_trans, whicharg="P", xiFeO=xiFeO, acc=1e-3
    )
    Z = -g_prime / K
    return K, Z


def W_FeO_MgO(P):
    return 11e3 + 0.011 * P / 1e5

def gamma_FeO(P, x_FeO=0.0):
    """Interaction parameter of FeO from Frost et al. 2010.

    Parameters:
    P (float): Pressure in Pa.
    x_FeO (float, optional): Mole fraction of FeO. Defaults to 0.

    Returns:
    float: Interaction parameter.
    """
    T = T_liquidus_pyrolite(P)
    W1 = W_FeO_Fe(P, T)
    W2 = W_Fe_FeO(P, T)
    r = (1.0 - x_FeO) ** 2 * (W1 + 2.0 * (W2 - W1) * x_FeO)
    r /= Rgas * T
    return np.exp(r)


def gamma_Fe(P, x_Fe=None):
    """Interaction parameter of Fe from Frost et al. 2010.

    Parameters:
    P (float): Pressure in Pa.
    x_Fe (float, optional): Mole fraction of Fe. Defaults to 0.

    Returns:
    float: Interaction parameter.
    """
    T = T_liquidus_pyrolite(P)
    W1 = W_FeO_Fe(P, T)
    W2 = W_Fe_FeO(P, T)
    r = (1.0 - x_Fe) ** 2 * (W2 + 2.0 * (W1 - W2) * x_Fe)
    r /= Rgas * T
    return np.exp(r)


def gamma_FeO_smoothed(P, xiFeO):
    xiFeO = min(xiFeO, 0.99999)
    if P < P_FeO_Fe_trans:
        return gamma_FeO(P=P, xiFeO=xiFeO)

    else:
        K, Z = coefs_FeO_Fe(P, xiFeO)
        return max(1.0 + K * np.exp(-Z * (P - P_FeO_Fe_trans)), 1.0)


def gamma_Fe_smoothed(P, xiFe):
    xiFe = min(xiFe, 0.99999)
    if P < P_FeO_Fe_trans:
        return gamma_Fe(P=P, xiFe=xiFe)

    else:
        K, Z = coefs_Fe_FeO(P, xiFe)
        return max(1.0 + K * np.exp(-Z * (P - P_FeO_Fe_trans)), 1.0)



def xi_S(P, T, xi):
    KD = KD_S(T, P, xi)
    return xi[2] / KD


def xi_FeO(P, T, xi):
    """
    Computes mol fraction of iron oxide in the mantle assuming chemical
    equilibrium between mantle and core. P, T at CMB and xiFe and xiO in core.
    """
    xiFe, xiH, xiS, xiSi, xiO = xi
    KD = partition_coefficient_i(T, P, 1, [xi[3], xi[4]])
    # print ('KD_O =', KD)
    return xiFe * xiO / KD


def xi_SiO2(P, T, xi):
    """
    Computes mol fraction of SiO2 in the mantle assuming chemical
    equilibrium between mantle and core. P, T at CMB and xiFe, xiSi and xiO
    in core.
    """
    xi[0] = 1.0 - (sum(xi) - xi[0])
    xiFe, xiH, xiS, xiSi, xiO = xi
    xiFeO = xi_FeO(P, T, xi)
    KD_Si = partition_coefficient_i(T, P, 0, [xi[3], xi[4]])
    # print ('xiFeO =', xiFeO)
    # print ('KD_Si =', KD_Si)

    return xiSi * xiFeO / (xiFe * KD_Si)


def logfO2(P, xiFe, xiFeO):
    """
    Computes oxygen fugacity according to Rose-Weston et al. 2009.
    The activity coefficients for Fe in the metal phase and FeO in the silicate
    phase are computed using a similar smoothing as in Schaefer et al. 2017 and
    along the pyrolite liquidus curve from Andrault 2011.
    """
    gamma_Fe_ = gamma_Fe_smoothed(P, xiFe)
    gamma_FeO_ = gamma_FeO_smoothed(P, xiFeO)
    aFe = gamma_Fe_ * xiFe
    aFeO = gamma_FeO_ * xiFeO
    return 2 * np.log10(aFeO / aFe)


def Fe_number_mantle(P, T, xi):
    """
    Computes Fe# in the mantle assuming chemical equilibrium between mantle
    and core. P, T at CMB and xi = xiFe, xiH, xiS, xiSi, xiO in core.
    """
    xiFe, xiH, xiS, xiSi, xiO = xi
    xiSiO2 = xi_SiO2(P, T, xi)
    xiFeO = xi_FeO(P, T, xi)
    result = xiFeO / (1.0 - xiSiO2)
    # print ('P, T, xi =', P, T, xi)
    # print ('xiSiO2 =', xiSiO2)
    # print ('xiFeO =', xiFeO)
    # if result <= 0. or result > xi_Fe_mantle_max:
    #    result = xi_Fe_mantle_max
    return result


def Si_number_mantle(P, T, xi):
    """
    Computes Si# in the mantle assuming chemical equilibrium between mantle
    and core. P, T at CMB and xi = xiFe, xiH, xiS, xiSi, xiO in core.
    """
    xiFe, xiH, xiS, xiSi, xiO = xi
    xiSiO2 = xi_SiO2(P, T, xi)
    xiFeO = xi_FeO(P, T, xi)
    SiMg = xiSiO2 / (1 - xiFeO - xiSiO2)
    result = SiMg / (1 + SiMg)
    # print ('P, T, xi =', P, T, xi)
    # print ('xiSiO2 =', xiSiO2)
    # print ('xiFeO =', xiFeO)
    # print ('Si/Mg =', SiMg)
    # result = min(result, Si_number_max(1-xiFeO))
    # result = max(result, Si_number_min(1-xiFeO))

    return result


def xi_Fe_mantle_min(MgSi):
    """
    Compute the min. Fe content in the silicates for given Mg/Si ratio.

    Returns
    -------
    Float.

    """

    # Calculation performed for (Mg,Fe)2SiO4 because for pyroxene Mg/Si
    return 1.0 - MgSi / 2.0


def xi_mantle(T, P, ll, x_core=0.0, x_Fe_core=0.0, x_O_core=0.0):
    """
    Computes the composition of the silicates from the core composition
    according to the single stage core-segregation model.

    Parameters:
    T (float): Temperature in K.
    P (float): Pressure in Pa.
    x_core (float):
    x_Fe_core (float): Molar Fe concentration in metals.
    x_O_core (float): Molar O concentration in metals.

    Returns:
    tuple: Molar abundance of oxygen and iron in the silicates.
    """
    KD = partition_coefficient(T, P, ll)

    # # Oxygen
    # if ll == 1:
    #     return xi_core / KD

    # else:

    #     # Compute FeO content in mantle first
    #     KD_FeO = partition_coefficient(T, P, 1)
    #     xi_FeO_mantle = xi_Fe_core * xi_O_core / KD_FeO
    #     # print("log KD_FeO =", np.log10(KD_FeO))
    #     # print("xi_FeO_mantle =", xi_FeO_mantle)
    #     return xi_core / xi_Fe_core * xi_FeO_mantle / KD

    # Compute FeO content in mantle first
    KD_FeO = partition_coefficient(T, P, 1)
    xi_FeO_mantle = x_Fe_core * x_O_core / KD_FeO
    return x_core / KD, x_core / x_Fe_core * xi_FeO_mantle / KD