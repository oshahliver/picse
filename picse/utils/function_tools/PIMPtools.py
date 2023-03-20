# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 19:59:16 2018

@author: os18o068
"""

import matplotlib.ticker
from matplotlib import pyplot as plt
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
    mMg,
    mFe,
    mSi,
    mO,
    mS,
    q_list,
    G,
    m_earth,
    r_earth,
    n_list,
    aT_list,
    bT_list,
    cT_list,
    gamma0_list,
    thetaD0_list,
)
from picse.runparams import (
    dens_min,
    dens_max,
    eps_P,
    bisectionstep_limit,
    pres_round,
    temp_round,
    dens_round,
    mass_round,
    radius_round,
)

# import PIMPeos
import numpy as np
import time
from scipy.integrate import quad
import sys


def C2(x2, y2):
    return mFe + (1.0 - y2) * mMg + x2 * mSi + (2.0 * x2 + 1.0) * mO


def C3(x3, y3):
    return (
        (x3 * y3 + 2.0 * y3) * mFe
        + (2.0 - 2.0 * y3) * mMg
        + (2.0 - x3 - 4.0 * x3 * y3) * mSi
        + (6.0 - 2.0 * x3) * mO
    )


def MFe_mantle(M2, M3, x2, x3, y2, y3):
    return (M2 / C2(x2, y2) + (x3 * y3 + 2.0 * y3) * M3 / C3(x3, y3)) * mFe


def MFe_core(M1):
    return mFe * M1 / (0.87 * mFe + 0.13 * (mFe + mS))


def MSi_mantle(M2, M3, x2, x3, y2, y3):
    return (M2 / C2(x2, y2) + (2.0 - x3) * M3 / C3(x3, y3)) * mSi


def FeSi_tot(M1, M2, M3, C2, C3, x3, y3):
    return (
        (MFe_core(M1) + MFe_mantle(M2, M3, x2, x3, y2, y3))
        / MSi_mantle(M2, M3, x2, x3, y2, y3)
        * mSi
        / mFe
    )


def dydx(
    x,
    y,
    derivative_type,
    whichEOS=None,
    lay=None,
    materials=None,
    density_type=None,
    Ttype="isothermal ambient",
    matdens=None,
    laymatfrac_list=None,
):
    """the mean density is used to reconstruct the different density contributions
    in the given layer. These contributions are then in turn used to compute
    the weighted average density gradient which is returned by the function
    output formag: [dP/dr, dm/dr, dT/dr, drho/dr]"""
    temperature_type = Ttype
    if derivative_type == "structure":
        matdensRK = density(
            whichEOS,
            lay,
            y,
            matdens,
            dens_min,
            dens_max,
            eps_P,
            materials,
            laymatfrac_list,
            Ttype=Ttype,
            identifier="Runge Kutta",
        )[1]

        if density_type == "constant":
            dgrad = 0.0
        elif density_type == "integrate":
            dgrad = meandens_grad(
                x, y, matdensRK, whichEOS, lay, materials, laymatfrac_list
            )[0]

        if temperature_type == "isothermal":
            Tgrad = 0.0
        elif (
            temperature_type == "isothermal ambient"
            or temperature_type == "ambient isothermal"
        ):
            Tgrad = 0.0
            y[2] = 300.0

        elif temperature_type == "adiabatic":
            md = [rho0_list[m] for m in materials]
            if lay == 3:
                gamma_bulk = 1.36 * (rho0_list[9] / matdensRK[0]) ** q_list[9]
                Tgrad = (
                    -gamma_bulk
                    * y[2]
                    / (
                        y[3]
                        * mean_dEOSdrho(
                            y,
                            matdensRK,
                            whichEOS,
                            lay,
                            materials,
                            laymatfrac_list,
                            Ttype=temperature_type,
                        )
                    )
                    * G
                    * y[1]
                    * y[3]
                    / x ** 2
                )
            elif lay == 2:
                gamma_bulk = (
                    1.96
                    * (mean_dens(md, lay, materials, laymatfrac_list) / y[3]) ** 2.5
                )
                Tgrad = (
                    -gamma_bulk
                    * y[2]
                    / (
                        y[3]
                        * mean_dEOSdrho(
                            y,
                            matdensRK,
                            whichEOS,
                            lay,
                            materials,
                            laymatfrac_list,
                            Ttype=temperature_type,
                        )
                    )
                    * G
                    * y[1]
                    * y[3]
                    / x ** 2
                )
            elif lay == 1:
                gamma_bulk = (
                    1.26
                    * (mean_dens(md, lay, materials, laymatfrac_list) / y[3]) ** 2.9
                )
                Tgrad = (
                    -gamma_bulk
                    * y[2]
                    / (
                        y[3]
                        * mean_dEOSdrho(
                            y,
                            matdensRK,
                            whichEOS,
                            lay,
                            materials,
                            laymatfrac_list,
                            Ttype=temperature_type,
                        )
                    )
                    * G
                    * y[1]
                    * y[3]
                    / x ** 2
                )
            elif lay == 0:
                Tgrad = 0.0

            else:
                print("WARNING: invalid layer identification!\nlayer =", lay)
        try:
            return (
                [-G * y[1] * y[3] / x ** 2, 4.0 * np.pi * x ** 2 * y[3], Tgrad, dgrad],
                Tgrad,
                dgrad,
            )
        except UnboundLocalError:
            print("WARNING!")
            print("layer = ", lay)
            print("Tgrad = ", Tgrad)
            sys.exit()

    elif derivative_type == "Eth":
        return [x ** 3 / (np.exp(x) - 1.0)], None, None

    elif derivative_type == "something else":
        # pass a coupled derivative vector you want to be integrated numerically
        pass


def RungeKutta(
    x,
    y,
    h,
    derivType,
    whichEOS=None,
    lay=None,
    materials=None,
    order=4,
    density_type=None,
    Ttype="isothermal",
    materialdensities=None,
    laymatfrac_list=None,
):
    """This is the Runge Kutta solver for the structure equations. Possible is a 4th order and a 5th order scheme. The 6th order scheme is equiped with an intrinsic
    integration step controler. Note that the relative deviation between RK4 and RK5 on the final mass and radius of a planet can amount to as much as 5%! The disadvantage
    of the intrinsic step controler in the 5th order scheme is that it requires considerably more integration steps. If no step controler is used, RK4 and RK5 require a
    comparable amount of time roughly the same amount of integration steps."""
    temperature_type = Ttype
    lmfl = laymatfrac_list

    if order == 4:
        # Butcher-table for the 4th order solver
        alist = [
            [
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            [0.5, 0.0, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ]
        blist = [1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0]
        clist = [0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0]

    k_list = []
    for i in range(order):
        ki = [y[s] for s in range(len(y))]
        for j in range(i):
            for s in range(len(y)):
                ki[s] += h * k_list[j][s] * alist[i][j]

        newdydx, Tgrad, dgrad = dydx(
            x + h * clist[i],
            ki,
            derivType,
            whichEOS,
            lay,
            materials,
            density_type,
            Ttype=temperature_type,
            matdens=materialdensities,
            laymatfrac_list=lmfl,
        )
        k_list.append(newdydx)

    for i in range(len(y)):
        for j in range(len(blist)):
            y[i] += h * k_list[j][i] * blist[j]

    if derivType == "structure":
        # mass and pressure gradient are computed and outputed to allow for
        # new integration step estimation outside of RK but otherwise they
        # are not needed
        Pgrad = pres_grad(x, y)
        mgrad = mass_grad(x, y)
        nextmd = density(
            whichEOS,
            lay,
            y,
            materialdensities,
            dens_min,
            dens_max,
            eps_P,
            materials,
            laymatfrac_list,
            "nextmatdens",
            Ttype=temperature_type,
        )[1]
        # dgrad, matdgrad = meandens_grad(x + dx, y, materialdensities, whichEOS, lay, materials, Ttype = temperature_type)
        return y, nextmd, h, Pgrad, mgrad, Tgrad, dgrad

    elif derivType == "Eth":
        return y, h


def alphaT(T, aT, bT, cT):
    """Compute thermal expansion coefficient at zero pressure"""
    return aT + bT * T - cT / T ** 2


def rho0T(T, T0, d0, ll, aT, bT, cT):
    # This is for BME only
    try:
        result = d0 / np.exp(
            aT * T
            + 0.5 * bT * T ** 2
            + 2.0 * cT / T ** 3
            - aT * T0
            - 0.5 * bT * T0 ** 2
            - 2.0 * cT / T0 ** 3
        )

        return result

    except TypeError:
        return d0


def gamma(q, gamma0, rho0, d):
    try:
        return gamma0 * (rho0 / d) ** q

    except ZeroDivisionError:
        print("WARNING: zero density passed to PIMPtools.gamma()!")
        print("returning:", gamma0)
        return gamma0


def thetaD(thetaD0, rho0, gamma, d):
    try:
        return (
            thetaD0 * (d / rho0) ** gamma
        )  # np.exp((gamma0_list[ll] - gamma(d, ll))/q_list[ll])
    except ZeroDivisionError:
        print("WARNING: zero density passed to PIMPtools.thetaD()!")
        print("returning:", thetaD0)
        return thetaD0


def Debye_integral(z):
    def integrand(x):
        return x ** 3 / (np.exp(x) - 1.0)

    def expint(a, b):
        return quad(integrand, a, b)[0]

    return 3 / z ** 3 * expint(0.0, z)


def Eth(d, T, thetaD0, rho0, gamma0, q, molmass, nparticles, ll, **kwargs):
    """Computes thermal energy contribution to the MGD EOS at given temperature
    and density. The density is relevant to compute the debye temperature
    """
    # using scipy is like soooo much faster than the RK scheme I coded for
    # the structure integration and that was previously also employed for
    # computing Eth. But at least a comparison showed, that my method and the
    # scipy result agreed well within << 0.01 % relative error, so my stuff
    # was accurate and robust, but very slow
    def integrand(t):
        return t ** 3 / (np.exp(t) - 1.0)

    def expint(a, b):
        """integrate given integrand between a and b using the scipy library"""
        return quad(integrand, a, b)[0]

    gam = gamma(q, gamma0, rho0, d)
    Tdebye = thetaD(thetaD0, rho0, gam, d)
    xendEth = Tdebye / T
    yendEth = expint(1.0e-12, xendEth)

    try:
        return (T / Tdebye) ** 3 * 9.0 * n_list[ll] * Rgas * T * yendEth * d / molmass

    except OverflowError:
        return 0e0


def EOS(
    d,
    T,
    which,
    ll,
    temperature_type="isothermal ambient",
    MGD_acc=1.0e-4,
    MGD_eps_x=0.5,
):

    if (
        temperature_type == "isothermal ambient"
        or temperature_type == "ambient isothermal"
    ):
        T = T0_list[ll]

    K0prime = K0prime_list[ll]
    d0T = rho0T(T, ll)
    d0 = rho0_list[ll]
    K0 = K0_list[ll]
    aP = aP_list[ll]
    T0 = T0_list[ll]

    # Birch-Murnaghan (BM)
    if which == 0:
        try:
            K = K0 + aP * (T - T0)
        except TypeError:
            K = K0
        eta = d / d0T
        return (
            3.0
            / 2.0
            * K
            * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
            * (1.0 + 3.0 / 4.0 * (K0prime - 4.0) * (eta ** (2.0 / 3.0) - 1.0))
        )

    # Mie-Grüneisen-Debye (MGD)
    elif which == 1:
        eta = d / d0
        if (
            temperature_type == "isothermal ambient"
            or temperature_type == "ambient isothermal"
        ):
            deltaP = 0.0
        else:
            acc = MGD_acc
            eps_x = MGD_eps_x
            E1 = Eth(d, T, acc, eps_x, ll)
            E2 = Eth(d, 300.0, acc, eps_x, ll)
            deltaP = (E1 - E2) * gamma(d, ll)
        return (
            3.0
            / 2.0
            * K0
            * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
            * (1.0 + 3.0 / 4.0 * (K0prime - 4.0) * (eta ** (2.0 / 3.0) - 1.0))
            + deltaP
        )

    # Vinet (Vinet)
    elif which == 2:
        eta = d / d0
        return (
            3.0
            * K0
            * eta ** (2.0 / 3.0)
            * (1.0 - eta ** (-1.0 / 3.0))
            * np.exp(3.0 / 2.0 * (K0prime - 1.0) * (1.0 - eta ** (-1.0 / 3.0)))
        )

    # Belonoshko (Bel)
    elif which == 3:  # for pure iron
        eta = d / d0
        if (
            temperature_type == "isothermal ambient"
            or temperature_type == "ambient isothermal"
        ):
            deltaP = 0.0
        else:
            deltaP = (
                3.0
                * Rgas
                * gamma(d, ll)
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

    # ideal gas law (IGL)
    elif which == 4:  # for atmospheres
        return d * kB * T / mH


def dEOSdrho(d, T, which, ll, temperature_type="isothermal ambient"):

    if (
        temperature_type == "isothermal ambient"
        or temperature_type == "ambient isothermal"
    ):
        T = T0_list[ll]

    # Birch-Murnaghan (BM)
    K0prime = K0prime_list[ll]
    d0T = rho0T(T, ll)
    d0 = rho0_list[ll]
    K0 = K0_list[ll]
    aP = aP_list[ll]
    T0 = T0_list[ll]

    if which == 0:
        K = K0 + aP * (T - T0)
        eta = d / d0T
        return (
            3.0
            / 2.0
            * K
            / d0T
            * (
                (7.0 / 3.0 * eta ** (4.0 / 3.0) - 5.0 / 3.0 * eta ** (2.0 / 3.0))
                * (1.0 - 3.0 / 4.0 * (4.0 - K0prime) * (eta ** (2.0 / 3.0) - 1.0))
                + 1.0
                / 2.0
                * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0))
                * (4.0 - K0prime)
                * eta ** (-1.0 / 3.0)
            )
        )

    # Mie-Grüneisen-Debye (MGD)
    elif which == 1:
        eta = d / d0
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
        )

    # Vinet (Vinet)
    elif which == 2:
        eta = d / d0
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
    elif which == 3:  # for pure iron
        eta = d / d0
        if (
            temperature_type == "isothermal ambient"
            or temperature_type == "ambient isothermal"
        ):
            delta = 0.0
        else:
            delta = 3.0 * Rgas * gamma(d, ll) * (T - T0_list[ll]) / molar_mass_list[ll]
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

    # ideal gas law (IGL)
    elif which == 4:  # for atmospheres
        return kB * T / mH


def simplemass(r, r0, d):
    return 4.0 / 3.0 * np.pi * d * (r - r0) ** 2


def mean_dens(md, lay, ll, laymatfrac_list):
    return 1.0 / sum([laymatfrac_list[lay][i] / md[i] for i in range(len(ll))])


def dens_grad(r, y, md, EOS, ll, Ttype="isothermal ambient"):
    return -G * y[1] * y[3] / (r ** 2 * PIMPeos.dEOSdrho(md, y[2], EOS, ll, Ttype))


def meandens_grad(
    r, y, materialdensities, EOS, lay, ll, laymatfrac_list, Ttype="isothermal ambient"
):
    """This computes the mean density gradient as weighted average over all the
    density contributions of the different materials present in the layer at hand"""
    dgrad_list = []
    for l in range(len(ll)):
        try:
            dgrad_list.append(
                dens_grad(r, y, materialdensities[l], EOS[ll[l]], ll[l], Ttype)
            )
        except IndexError:
            print("WARNING: IndexError in function meandens_grad()")
            if len(dgrad_list) != len(laymatfrac_list[lay]):
                print(
                    "REASON: density gradient and layer fractions have not the same dimensions"
                )
            print(
                "parameters that lead to falure are:\nlayer =",
                lay,
                ", l =",
                l,
                ", ll =",
                ll,
            )
            print(
                "layer fractions =", laymatfrac_list[lay], ", dgrad_list =", dgrad_list
            )
            sys.exit()

    return (
        1.0 / sum([laymatfrac_list[lay][i] / dgrad_list[i] for i in range(len(ll))]),
        dgrad_list,
    )


def mean_dEOSdrho(
    y, matdens, whichEOS, lay, ll, laymatfrac_list, Ttype="isothermal ambient"
):
    return sum(
        laymatfrac_list[lay][i]
        * dEOSdrho(matdens[i], y[2], whichEOS[ll[i]], ll[i], temperature_type=Ttype)
        for i in range(len(ll))
    )


def pres_grad(r, y):
    return -G * y[1] * y[3] / r ** 2


def mass_grad(r, y):
    return 4.0 * np.pi * r ** 2 * y[3]


def temp_grad(r, y):
    return None


def bisec(
    which,
    y,
    a,
    b,
    eps,
    ll,
    identifier,
    Ttype="isothermal ambient",
    materialdensity=None,
    firstlimit=10,
    secondlimit=20,
    thirdlimit=45,
    forthlimit=200,
    firstepsmin=0.995,
    firstepsmax=1.005,
    secondepsmin=0.9,
    secondepsmax=1.1,
    thirdepsmin=0.5,
    thirdepsmax=1.5,
    densmin=5.0e2,
    densmax=2.0e4,
):
    """compute the density at a given pressure for a specific EOS by numerically inverting
    the P(rho) relation to a decired precission eps"""
    a = a * firstepsmin
    b = b * firstepsmax
    a0 = a
    b0 = b
    P = y[0]
    T = y[2]
    limit = bisectionstep_limit
    if not materialdensity == None:
        d = materialdensity
    else:
        d = y[3]
    counter = 0
    limit = firstlimit
    exceeded = False
    doubleexceeded = False
    tripleexceeded = False
    if P <= 0.0:
        print("WARNING: non-positive pressure given")
        return None
    elif P > 0.1:
        c0 = (a + b) / 2.0
        c = c0
        EOSc = EOS(c0, T, which, ll, Ttype)
        if EOSc == P:
            pass
        else:
            while abs(EOSc - P) / P > eps:
                counter += 1
                c = (a + b) / 2.0
                EOSc = EOS(c, T, which, ll, Ttype)
                EOSb = EOS(b, T, which, ll, Ttype)
                if not P == EOSc:
                    if (EOSc - P) / abs(EOSc - P) == (EOSb - P) / abs(EOSb - P):
                        b = c
                    else:
                        a = c
                else:
                    pass
                    # print ('NOTE: P = EOSc')

                if counter > limit and not exceeded:
                    exceeded = True
                    limit = secondlimit
                    counter = 0
                    a = a0 * secondepsmin
                    b = b0 * secondepsmax
                elif counter > limit and exceeded and not doubleexceeded:
                    doubleexceeded = True
                    limit = thirdlimit
                    counter = 0
                    a = a0 * thirdepsmin
                    b = b0 * thirdepsmax
                elif counter > limit and doubleexceeded and not tripleexceeded:
                    tripleexceeded = True
                    limit = forthlimit
                    counter = 0
                    a = densmin
                    b = densmax
                    print(
                        "NOTE: setting maximal density range for bisection to [",
                        a,
                        b,
                        "]",
                    )
                elif counter > limit and tripleexceeded:
                    print(
                        "WARNING: bisection step limit exceeded in <", identifier, "> !"
                    )
                    print("-> parameters that lead to failure are:")
                    print("counter =", counter, "\nlimit =", limit)
                    print(
                        "a =",
                        a,
                        "\nb =",
                        b,
                        "\nc =",
                        c,
                        "\na0 =",
                        a0,
                        "\nb0 =",
                        b0,
                        "\nc0 =",
                        c0,
                        "\ndens =",
                        d,
                        "\npres = ",
                        P * 1.0e-9,
                        "GPa",
                        "\ntemp =",
                        T,
                        "[K] \nTtype =",
                        Ttype,
                        "\nmatdens =",
                        materialdensity,
                        "\nM/Mearth =",
                        y[1] / m_earth,
                    )
                    print("EOS(c) =", EOSc, "\nEOS type =", which, "\nmaterial =", ll)
                    eps = 10.0
        return c


def density(
    whichEOS,
    lay,
    param_list,
    matdens_list,
    a,
    b,
    eps,
    ll,
    laymatfrac_list,
    identifier,
    Ttype="isothermal ambient",
    layvic=False,
):
    """Here the mass density is computed as a weighted average over all
    contributions of the different materials in each layer. Note that the
    materialdensity in this function is only used to estimate the bisection
    interval to find the density to the corresponding pressure. It does therefore
    not matter if the old or updated materialdensities are used in terms of the
    computed value for the new density and materialdensities"""
    matdens_dummy_list = []
    for l in range(len(ll)):
        try:
            matdens = matdens_list[l]
        except TypeError:
            matdens = None
        except IndexError:
            pass
        # compute mols per volume
        if layvic:
            a = dens_min
            b = dens_max
        else:
            try:
                a = matdens_list[l]
                b = matdens_list[l]
            except TypeError:
                a = dens_min
                b = dens_max
        matdens = bisec(
            whichEOS[ll[l]],
            param_list,
            a,
            b,
            eps,
            ll[l],
            identifier,
            Ttype=Ttype,
            materialdensity=matdens,
        )
        matdens_dummy_list.append(matdens)
    try:
        return (
            1.0
            / sum(
                [
                    laymatfrac_list[lay][i] / matdens_dummy_list[i]
                    for i in range(len(ll))
                ]
            ),
            matdens_dummy_list,
        )
    except TypeError:
        print(
            "laymatfrac =", laymatfrac_list[lay], "matdens dummy =", matdens_dummy_list
        )
        print(
            "M/Mearth =",
            param_list[1] / m_earth,
            "pres =",
            param_list[0],
            "temp =",
            param_list[2],
            "dens =",
            param_list[3],
        )
        sys.exit()


def layer_transition_notification(x, y, Mtot, newlayer):
    print("-----------------------------------------")
    print(
        "-> transition to layer",
        newlayer,
        "at: ",
        round(x * 1.0e-3, radius_round),
        "km",
        ",",
        round(y[0] * 1.0e-9, pres_round),
        "GPa",
        ", M/Mtot = ",
        round(y[1] / Mtot, mass_round),
        "T =",
        round(y[2], temp_round),
        "K",
    )
