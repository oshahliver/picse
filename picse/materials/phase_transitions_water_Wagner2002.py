# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 15:19:16 2018

@author: os18o068

This module computes the melting and sublimation curve for pure water following
the the IAPWS formulation (Wagner 2011)
"""

from matplotlib import pyplot as plt
import matplotlib.ticker
from matplotlib import rc
import matplotlib.transforms as transforms

import numpy as np
from scipy.interpolate import interp1d
from picse.utils.function_tools import functionTools as ftool
import iapws
from picse.runparams import eos_pres_trans, plot_params, color_list
import matplotlib
import seafreeze

# -----------------
# physical constants
Rgas = 8.3144598  # universal gas constant
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

# ------------------------
# EOS specific parameters

T_triple = 273.16  # triple point temperature
P_triple = 611.657  # triple point pressure in Pa
T_critical = 647.096
P_critical = 22.064e6
rho_critical = 322.0
P0 = 1.01325e5  # normal pressure
T0 = 273.152519  # normal melting point
T0_melt = 273.152519  # normal melting temperature
s0 = 189.13  # j kg-1 K-1


# complex constants [re, im]
t1 = complex(3.68017112855051e-02, 5.10878114959572e-02)
r1 = complex(44.7050716285388, 65.6876847463481)
t2 = complex(0.337315741065416, 0.335449415919309)
r2_list = [
    complex(-72.597457432922, -78.100842711287),
    complex(-5.57107698030123e-05, 4.64578634580806e-05),
    complex(2.34801409215913e-11, -2.85651142904972e-11),
]
t_list = [t1, t2]

# real constants
g0_list = [
    -632020.233449497,
    0.655022213658955,
    -1.89369929326131e-08,
    3.39746123271053e-15,
    -5.56464869058991e-22,
]
amelt_list = [0.119539337e7, 0.808183159e5, 0.333826860e4]
bmelt_list = [0.3000000e1, 0.257500e2, 0.103750e3]

# coefficients to compute sublimation curve
asubl_list = [-0.212144006e2, 0.273203819e2, -0.610598130e1]
bsubl_list = [0.333333333e-2, 0.120666667e1, 0.170333333e1]

# coefficients to compute vapor-liquid transition curve (Wagner 2002)
aevap_list = [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502]

# convention: Ih/l/v, Ih/III/l, III/V/l, V/VI/l, VI/VII/l, Ih/III/II, III/V/II,
# VI/V/II, VII/VI/VIII, Ih/II/IX
# in K
t_triple_list = [
    T_triple,
    251.165,
    256.164,
    273.31,
    355.0,
    238.46,
    248.86,
    218.16,
    268.16,
    73.4,
    130.0,
]
# in Pa
p_triple_list = [
    P_triple,
    208.566e6,
    350.1e6,
    632.4e6,
    2216.0e6,
    212.9e6,
    344.3e6,
    620.0e6,
    2100.0e6,
    89.6e6,
    828.5e6,
]

# high pressure transition Schwager 2008
# [P(GPa), T(K)]
data_schwager = [
    [1.59004e1, 7.50567e2],
    [1.79379e1, 9.25928e2],
    [1.91294e1, 1.01505e3],
    [2.11794e1, 1.13581e3],
    [2.37529e1, 1.23935e3],
    [2.68479e1, 1.33429e3],
    [2.97783e1, 1.39187e3],
    [3.30557e1, 1.44659e3],
    [3.65110e1, 1.48119e3],
    [3.92795e1, 1.48991e3],
    [4.11624e1, 1.58481e3],
    [4.38178e1, 1.70847e3],
    [4.73337e1, 1.85514e3],
    [5.05093e1, 1.97595e3],
    [5.47225e1, 2.10254e3],
    [5.96303e1, 2.22053e3],
    [6.32514e1, 2.28675e3],
    [6.79986e1, 2.35014e3],
    [7.25787e1, 2.38766e3],
    [7.76829e1, 2.40509e3],
    [8.28744e1, 2.41964e3],
    [8.74617e1, 2.42555e3],
    [8.96262e1, 2.42563e3],
]


data_dunaeva = [
    [889.529817905081, 2387.7099716545163],
    [768.0059999779412, 2369.9232723409104],
    [636.9572170690546, 2239.798750739887],
    [597.3683920279705, 2209.662832121941],
    [577.5574354506048, 2177.030599152209],
    [547.3413700684923, 2097.641185142702],
    [530.3605501450361, 2069.6706997400743],
    [531.2076058543901, 1968.961511170915],
    [511.3988551512677, 1938.671181355951],
]

data_French_2016 = [
    [15.272150521285202, 726.07339616628],
    [31.468989889748443, 1047.8645904475707],
    [58.80718699340868, 1420.422984109639],
    [107.4907126067832, 1875.290260893321],
    [185.46815943936951, 2089.730455766516],
    [304.5948978849551, 2063.732185949096],
    [440.2698867475116, 1683.4790329190005],
    [492.6901436682674, 1328.894049082548],
    [528.7639637094252, 961.7083276419426],
    [537.0886644278124, 683.2058461958741],
]

# Table 2
LtoX_coeffs_Dunaeva2010 = [0.2524, 0.0019, 0.2795, 0.5, 1.1675]

# From fit to extracted data from figure 2
LtoVII_coeffs_Dunaeva2010 = [-2.29294290e-01, 4.11453047e01, 5.46863218e02]

# Fit solidus T(P), T output in K P input in GPa
LtoVII_coeffs_Sun = [1.14666190e-04, -7.94260914e-02, 1.83673796e01, 5.52107003e02]

data_Jiming_Sun = [
    [51.214508580343235, 2522.1029641185646],
    [68.28315132605304, 2742.0218408736346],
    [94.23946957878319, 3005.7098283931355],
    [117.69656786271453, 3181.2605304212166],
    [142.43369734789394, 3343.525741029641],
    [167.8307332293292, 3461.647425897035],
    [199.59984399375975, 3566.2839313572536],
]

data_dunaeva_LtoVII = [
    [97.22, 2384.0],
    [87.2637087599545, 2373.174061433446],
    [75.14994311717861, 2359.3515358361765],
    [62.289874857792945, 2230.8532423208185],
    [58.265301478953354, 2201.4675767918084],
]

# 4th order polynomial fit for solidus French 2016 Figure 4
# Input is P(GPa) output is T(K)
coeffs_French2016 = [
    -2.26536372e-07,
    2.64778698e-04,
    -1.19895845e-01,
    2.37743503e01,
    3.96728823e02,
]

# -------------
# functions


def bisec(a, b, y, f, *args, eps=1.0e-8, plot=False, axis=None):
    limit = 100
    counter = 0
    c0 = (a + b) / 2.0
    c = c0
    fc = f(c0, *args)
    if fc == y:
        print("NOTE: fc = y in bisection")

    else:
        while abs(fc - y) / y > eps:
            if plot:
                axis.plot(
                    fc,
                    c,
                    linestyle="None",
                    marker="x",
                    color="k",
                    zorder=100,
                    markersize=6,
                )
            counter += 1
            c = (a + b) / 2.0
            fc = f(c, *args)
            fb = f(b, *args)
            if not y == fc:
                if (fc - y) / abs(fc - y) == (fb - y) / abs(fb - y):
                    b = c
                else:
                    a = c
            else:
                break

            if counter > limit:
                break
                print("WARNING: bisection step limit exceeded!")
        if plot:
            axis.scatter(
                fc, c, marker="o", color="r", s=100, zorder=1000, facecolor="none"
            )
    return c


def iceI(t):
    # pressure in MPa
    theta = t / 273.16
    return (
        1.0 + sum([amelt_list[i] * (1.0 - theta ** bmelt_list[i]) for i in range(3)])
    ) * 611.657e-6


def iceIII(t):
    # pressure in MPa
    theta = t / 251.165
    return (1.0 - 0.299948 * (1.0 - theta**60)) * 208.566


def iceV(t):
    # pressure in MPa
    theta = t / 256.164
    return (1.0 - 1.18721 * (1.0 - theta**8)) * 350.1


def iceVI(t):
    # pressure in MPa
    theta = t / 273.31
    return (1.0 - 1.07476 * (1.0 - theta**4.6)) * 632.4


def iceVII(t):
    # pressure in MPa
    return iapws._Melting_Pressure(t)

    """np.exp(1.73683 * (1. - 1./theta) - 0.0544606 * (1. - theta**5) \
                  + 0.806106e-7 * (1. - theta**22)) * 2216.0"""


def sublPres(t):
    """compute sublimation pressure in MPa at given temperature with
    analytical expressin from IAPWS
    """

    theta = t / 273.16

    return (
        np.exp(
            1.0
            / theta
            * sum([asubl_list[i] * theta ** bsubl_list[i] for i in range(3)])
        )
        * 611.657e-6
    )


def meltPres(T=None):
    pass


def meltTemp(p):
    """Compute metling temperature in Kelvin with analytical expression from
    IAPWS
    """
    tmin = 251.165
    tmax = 715.0

    if p < 208.566e6:
        fct = iceI

    elif p >= 208.566e6 and p < 350.1e6:
        fct = iceIII

    elif p >= 350.1e6 and p < 632.4e6:
        fct = iceV

    elif p >= 632.4e6 and p < 2216e6:
        fct = iceVI

    elif p >= 2216e6:
        fct = iceVII

    # fct gives pressure in MPa but input to meltTemp is Pa
    return bisec(tmin, tmax, p * 1.0e-6, fct)


def evapPres(t):
    """Copmute evaporation pressure in Pa at given temperature"""
    theta = 1.0 - t / T_critical
    a1, a2, a3, a4, a5, a6 = aevap_list
    res = P_critical * np.exp(
        T_critical
        / t
        * (
            a1 * theta
            + a2 * theta**1.5
            + a3 * theta**3
            + a4 * theta**3.5
            + a5 * theta**4
            + a6 * theta**7.5
        )
    )
    return res


def evapTemp(p):
    tmin = T_triple
    tmax = T_critical
    return bisec(tmin, tmax, p, evapPres)


def sublTemp(p):
    """Compute sublimation temperature in Kelvin"""
    tmin = 10.0
    tmax = T_critical
    # subl gives pressure in MPa but input to sublTemp is Pa, so factor 1.0e-6
    return bisec(tmin, tmax, p * 1.0e-6, sublPres)


def solidusTempHigh(P=None):

    if P < 1.0e9:
        # For T < 2380K fit to high T data from Figure 2
        a, b, c = LtoVII_coeffs_Dunaeva2010
        P = P * 1.0e-9

        return c + P * b + P**2 * a

    else:
        a, b, c, d, e = coeffs_French2016
        P = P * 1.0e-9

        return e + d * P + c * P**2 + b * P**3 + a * P**4
    """
    elif P < 7.0e10:
        #Pressure input to fit is in bar and outbput in K
        #Convert Pa to bar
        P *=1.0e-5
        
        a,b,c,d,e = LtoX_coeffs_Dunaeva2010
        
        #Return solodus temperature in K
        return a + b*P + c*np.log(P) + d/P + e*np.sqrt(P)
    """


def solidusTemp(P=None, order=1):
    # high pressure region
    if P >= 5.0e9:
        # fit = fitSchwager(order=2, whicharg='P')
        # temp = fit(p*1.0e-9)
        temp = solidusTempHigh(P=P)

    # interpolation region
    # elif p >= 1.0e10 and p < 1.59e10:
    #   inter = interpolateWagnerSchwager(order=order, whicharg='P')
    #  temp = inter(p*1.0e-9)

    # low pressure region
    else:
        if P < P_triple:
            temp = sublTemp(P)

        else:
            temp = meltTemp(P)

    return temp


def solidusPres(T=None, order=2, acc=1e-4):
    P = ftool.bisec(
        f=solidusTemp,
        whicharg="P",
        y=T,
        identity="solidusPres",
        a=1e3,
        b=1e11,
        noisy=False,
        order=2,
        acc=acc,
    )
    return P
    """
    if T >= 750.:
        fit = fitSchwager(order=2, whicharg='T')
        return fit(T)
    
    elif T >= 600. and T < 750.:
        inter = interpolateWagnerSchwager(order=order, whicharg='T')
        return inter(T)
    
    else:
        if T < T_triple:
            return sublPres(T)
        
        else:
            return meltPres(T)
    """


def fitSchwager(order=2, whicharg="P"):
    # rearange data
    # pressure in GPa
    x = [d[0] for d in data_schwager]
    # temperature in K
    y = [d[1] for d in data_schwager]

    if whicharg == "P":
        z = np.polyfit(x, y, order)

    elif whicharg == "T":
        z = np.polyfit(y, x, order)

    poly = np.poly1d(z)
    return poly


def interpolateWagnerSchwager(N_Schwager=3, N_Wagner=3, order=3, whicharg="T"):
    """Interpolates solidus line between Wagner and Schwager data around T_crit
    Returns polyfit instance with x=T (K) and y=P (GPa) (whicharg = 'T')
    or x=P (GPa) and y=T (K) (whicharg = 'P') of solidus line
    """
    xx = []
    yy = []

    if whicharg == "T":
        for i in range(N_Wagner):
            p = 8.0e9 + i * 0.5e9
            xx.append(solidusTemp(p))
            yy.append(p * 1.0e-9)

        for i in range(N_Schwager):
            xx.append(data_schwager[i][1])
            yy.append(data_schwager[i][0])

    elif whicharg == "P":
        for i in range(N_Wagner):
            p = 8.0e9 + i * 0.5e9
            xx.append(p * 1.0e-9)
            yy.append(solidusTemp(p))

        for i in range(N_Schwager):
            xx.append(data_schwager[i][0])
            yy.append(data_schwager[i][1])

    z = np.polyfit(xx, yy, order)
    poly = np.poly1d(z)
    return poly


def solidLine(pmin=1.0e1, pmax=1.0e10, log=False, N=200, **kwargs):
    """Returns array containing a set of N P-T pairs on the phase transition line
    between pmin and pmax on a logarithmic or linear scale for the transition
    of water ices to liquid or vapor.
    """
    if log:
        pres = np.logspace(np.log10(pmin), np.log10(pmax), N)

    else:
        pres = np.linspace(pmin, pmax, N)

    return np.array([solidusTemp(p) for p in pres]), pres


def vaporLine(pmin=P_triple, pmax=P_critical, log=False, N=10, **kwargs):
    """Returns array containing the temperatures corresponding to the v-l transition
    for the given pressure range
    """
    if log:
        pres = np.logspace(np.log10(pmin), np.log10(pmax), N)

    else:
        pres = np.linspace(pmin, pmax, N)

    return np.array([evapTemp(p) for p in pres]), pres


def temp_6_to_7(P):
    """ice VI to ice VII transision from Haldemann et al. 2020"""
    x1, x2, x3, x4 = -1.4699e5, 6.10791e-6, 8.1529e3, -8.8439e-1

    return x1 + x2 * P + x3 * np.log(P) + x4 * np.sqrt(P)


def phaseCheck(p, t, double_check=False, dens=None, noisy=False, **kwargs):
    # Compute solidus temperature at given pressure
    temp_solidus = solidusTemp(p)
    temp_67 = temp_6_to_7(p)

    # ice Ih to ice IX at approx 200 MPa
    if p > 2.0e8:
        if t > temp_solidus:
            ph = "Mazevet"

        else:
            if t < temp_67:
                ph = "French"

            else:
                ph = "Journaux"

                PT = np.empty((1,), np.object)
                PT[0] = (p * 1.0e-6, t)

                phh = seafreeze.whichphase(PT)
                ph = seafreeze.phasenum2phase[phh[0]]

    else:
        if t < temp_solidus:
            ph = "Ih"

        elif t >= temp_solidus:
            if p < P_triple:
                if t < T_critical:
                    ph = "vapor"

                else:
                    ph = "gas"

            elif p >= P_triple:
                if p < P_critical:
                    if t < T_critical:
                        pres_vapor = evapPres(t)
                        if p > pres_vapor:
                            ph = "liquid"

                        else:
                            ph = "vapor"
                    else:
                        ph = "gas"
                else:
                    if t < T_critical:
                        ph = "liquid"
                    else:
                        ph = "supercritical"
            else:
                print("weird 2")

        else:
            print("ERROR in water phase check: T/P/d =", t, p, dens)

    # actually compute density and double check if it roughly matches the
    # density regime in the predicted phase
    if double_check:
        if ph == "vapor" or ph == "gas":
            try:
                if dens > 200.0:
                    if p < P_triple:
                        if noisy:
                            print("correcting", ph, "to solid")
                        ph = "ice Ih"

                    elif p >= P_triple and p < P_critical:
                        if noisy:
                            print("correcting", ph, "to liquid")
                        ph = "liquid"

                    else:
                        if noisy:
                            print("correcting", ph, "to supercritical")
                        ph = "supercritical"

            except TypeError:
                print("WARNING: cannot determine phase!")
                pass

        elif ph == "ice Ih" or ph == "liquid":
            try:
                if dens < 200.0:
                    if p < P_critical:
                        if noisy:
                            print("correcting", ph, "to vapor")
                        ph = "vapor"

                    else:
                        if noisy:
                            print("correcting", ph, "to supercritical")
                        ph = "supercritical"

            except TypeError:
                print("WARNING: cannot determine phase!")
                pass

    return ph


def phase(P=1e5, T=300.0):
    # Compute solidus temperature at given pressure
    temp_solidus = solidusTemp(P)

    if T > temp_solidus:
        if T > T_critical:
            if P > P_critical:
                ph = "supercritical"

            else:
                ph = "gas"

        else:
            pres_vapor = evapPres(T)
            if P < pres_vapor:
                ph = "vapor"
            else:
                ph = "liquid"

    else:
        ph = "solid"

    return ph


def plot_phases(N=16):

    temps = np.linspace(200, 300, N)
    pres = np.logspace(8, 10, N)
    phases = np.empty([N, N]).astype(int)
    phase_strings = [
        "Ih",
        "Mazevet",
        "French",
        "V",
        "VI",
        "II",
        "III",
        "liquid",
        "vapor",
        "gas",
    ]
    fig, ax = plt.subplots()
    cols = ["r", "g", "b", "k", "grey", "orange", "yellow", "purple", "magenta"]
    for i in range(N):
        for j in range(N):
            ph = phaseCheck(pres[j], temps[i])
            try:
                phases[i][j] = phase_strings.index(ph)
            except ValueError:
                phases[i][j] = 0

    print("plotting ...")
    ax.imshow(phases.T, origin="lower", cmap="Paired")

    # ax.set_yscale('log')


def plot_the_shit(
    xlims=[100.0, 1000.0],
    ylims=[1e5, 1e12],
    minor=[True, True],
    axis_scales=[1.0, 1.0],
    showOnlyStates=True,
):
    matplotlib.rc("text", usetex=True)
    matplotlib.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath, amssymb}"]

    # reformat Schwager 2008 data for plotting
    pres_schwager = np.empty([len(data_schwager)])
    temp_schwager = np.empty([len(data_schwager)])

    for i in range(len(data_schwager)):
        pres_schwager[i] = data_schwager[i][0]
        temp_schwager[i] = data_schwager[i][1]

    fnts1 = 18
    fnts2 = 18

    lwdth1 = 1
    lwdth2 = 2

    p = 1.0e8
    t = 900
    fig, ax = plt.subplots(figsize=(8, 8))
    # ax.set_facecolor(plot_params['backgroundcol'])

    ax.tick_params(
        which="major",
        labelsize=fnts2,
        direction="out",
        size=8,
        top=True,
        right=True,
        pad=10,
    )
    ax.tick_params(
        which="minor", labelsize=fnts2, direction="out", size=4, top=True, right=True
    )

    plot_color = color_list[2]
    phase_color = color_list[2]

    font1 = {
        "family": "sans",
        "color": "k",
        "weight": "bold",
        "size": fnts1,
    }

    font2 = {
        "family": "sans",
        "color": "k",
        "weight": "bold",
        "size": fnts2,
    }

    font3 = {
        "family": "sans",
        "color": "k",
        "weight": "normal",
        "size": fnts1,
    }

    font4 = {
        "family": "sans",
        "color": "k",
        "weight": "normal",
        "size": fnts1,
    }

    font5 = {
        "family": "sans",
        "color": "w",
        "weight": "normal",
        "size": 84,
    }

    font6 = {
        "family": "sans",
        "color": color_list[2],
        "weight": "normal",
        "size": fnts1,
    }

    plot_params["backgroundcol"] = "white"
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transAxes)
    if not showOnlyStates:
        ax.text(0.05, 0.1, s="solid", fontdict=font1, transform=trans)
    ax.text(0.6, 0.9, s="solid", fontdict=font1, transform=trans)
    ax.text(0.3, 0.5, s="liquid", fontdict=font1, transform=trans)
    ax.text(0.6, 0.1, s="vapor", fontdict=font1, transform=trans)
    # ax.text(.7, .1, s = 'gas', fontdict = font1, transform = trans)
    ax.text(0.7, 0.5, s="supercritical", fontdict=font1, transform=trans)
    ax.text(T_critical * 1.01, P_critical * 1.1, r"Critical  point", fontdict=font1)
    """
    ax.text(300., 1.0e8, s = 'IAPWS', 
            fontdict = font2)
    ax.text(1100., 5.0e9, s = 'Mazevet 2018', fontdict = font2)
    ax.text(100., 3.0e8, s = 'Journaux 2020', fontdict = font2)
    ax.text(330., 1.0e11, s = 'French 2015', fontdict = font2)
    ax.text(100., 1.0e6, s = 'Feistel 2006', fontdict = font2, rotation =0)
    
    #ax.text(350., 2.0e9, s = 'French 2016', fontdict = font4, rotation =5,
     #       bbox=dict(facecolor=plot_params['backgroundcol'], edgecolor='none', 
                 #     alpha=.75))
    ax.text(190., 1.5e9, s = 'Haldemann 2020', fontdict = font4, rotation =2,
            bbox=dict(facecolor=plot_params['backgroundcol'], edgecolor='none', 
                   
                      alpha=.75))
    ax.text(270., 1.0e5, s = 'Wagner 2002', fontdict=font4, rotation =90,
            bbox=dict(facecolor=plot_params['backgroundcol'], edgecolor='none', 
                      alpha=.75))
    #ax.text(280., 7.0e8, s = 'Wagner 2002', fontdict=font4, rotation =15,
     #       bbox=dict(facecolor=plot_params['backgroundcol'], edgecolor='none', 
      #                alpha=.75))
    ax.text(330., 5.0e4, s = 'Wagner 2002', fontdict=font6, rotation = 70,
            bbox=dict(facecolor=plot_params['backgroundcol'], edgecolor='none', 
                      alpha=.75))
    
    ax.text(1500., 1.6e8, s = r'$0.2  \ \rm GPa$', fontdict=font4, rotation =0,
            bbox=dict(facecolor=plot_params['backgroundcol'], edgecolor='none', 
                      alpha=1.))
    """
    ax.text(
        0.1,
        0.8,
        s=r"$\rm H_2 O$",
        fontdict=font5,
        rotation=0,
        # bbox=dict(facecolor='white', edgecolor='none', alpha=1.),
        alpha=1.0,
        transform=trans,
        ha="left",
    )
    """
    ax.text(242, 1.8e8, 'III', fontdict=font1)
    ax.text(250., 3.8e8, 'V', fontdict=font1)
    ax.text(220., 1.0e7, 'Ih', fontdict=font1)
    ax.text(220., 2.5e8, 'II', fontdict=font1)
    ax.text(230., 8.0e8, 'VI', fontdict=font1)
    ax.text(300, 4.0e9, 'VII', fontdict=font1)
    """
    # ax.text(T_triple + 20, P_triple, s = 'triple point', fontdict = font3,
    #       color = phase_color)

    # ax.text(T_critical + 20, P_critical*2, s = 'critical point',
    #           fontdict = font3, color = color_list[0])
    """
    #Plot gas-vapor interface
    ax.plot([T_critical, T_critical], [1.0, P_critical],
            color = (.8, .8, .8), linestyle = '--', zorder = 1)
    
    #Plot scf-liq interface
    ax.plot([T_critical, T_critical], [P_critical, 1.0e10], 
            [T_critical, 2500.], [P_critical, P_critical],
            color = color_list[0], linestyle = '--', zorder = 1)
    """
    inter = interpolateWagnerSchwager(order=1, whicharg="T", N_Schwager=1, N_Wagner=1)

    # Compute solidus line in entire pressure range
    pres_list = np.logspace(1, np.log10(5.0e11), 100)
    temp_list = [solidusTemp(p) for p in pres_list]

    # plot solidus
    ax.plot(temp_list, pres_list, color="k", linewidth=lwdth2)

    # T_inter = np.linspace(600, 750, 10)
    # ax.plot(T_inter, inter(T_inter)*1.0e9, color = plot_color)

    # T_dunaeva_list = [solidusTempHigh(P=p) for p in pres_list]

    # ax.plot(T_dunaeva_list, pres_list, color='r')

    # ax.plot(temp_schwager, pres_schwager*1.0e9, color = 'g')
    # convention: Ih/l/v, Ih/III/l, III/V/l, V/VI/l, VI/VII/l, Ih/III/II, III/V/II,
    # VI/V/II, VII/VI/VIII, Ih/II/IX, II/VI/XV
    """
    ax.plot([t_triple_list[1], t_triple_list[5]], [p_triple_list[1], 
             p_triple_list[5]], 
            [t_triple_list[6], t_triple_list[7]], [p_triple_list[6], 
             p_triple_list[7]],
            [t_triple_list[2], t_triple_list[6]], [p_triple_list[2], 
             p_triple_list[6]],
            [t_triple_list[6], t_triple_list[5]], [p_triple_list[6], 
             p_triple_list[5]],
            [t_triple_list[7], t_triple_list[3]], [p_triple_list[7], 
             p_triple_list[3]],
            [t_triple_list[8], t_triple_list[4]], [p_triple_list[8], 
             p_triple_list[4]],
            [t_triple_list[5], t_triple_list[9]], [p_triple_list[5], 
             p_triple_list[9]],
            [t_triple_list[7], t_triple_list[10]], [p_triple_list[7], 
             p_triple_list[10]],
            color = plot_color, linestyle = ':', zorder = 1)
    
    
    for i in range(len(t_triple_list)):
        ax.scatter(t_triple_list[i], p_triple_list[i], color='blue')
    """
    # ax.plot([301, 2500], [1.0e9, 1.0e9], color='k', linestyle='--',
    #       linewidth=lwdth2)
    phase_list = solidLine(log=True, N=50)
    x, vapor = vaporLine(log=True, N=50)

    # Plot vapor line
    # ax.semilogy(phase_list[0], phase_list[1], zorder = 1, color = phase_color)
    ax.plot(x, vapor, zorder=2, color="k", linestyle="--")

    ax.scatter(T_critical, P_critical, color="k", zorder=2, s=20)
    ax.scatter(t_triple_list[0], p_triple_list[0], color=phase_color, zorder=2, s=20)

    if not showOnlyStates:
        # Compute ice VI > ice VII transition line from Haldemann et al. 2020
        # up to the triple point VI-VII-liq at 2.216 GPa
        pres_list = np.logspace(9, np.log10(2.216e9), 10)
        temp_list = temp_6_to_7(pres_list)
        ax.plot(temp_list, pres_list, color="k", linewidth=lwdth2, linestyle="-")

    temp_fill_list_scf = np.linspace(T_critical, xlims[1])
    pres_fill_list_scf = np.array([solidusPres(t) for t in temp_fill_list_scf])

    temp_fill_list_liq = np.linspace(T_triple, T_critical)
    pres_fill_list_liq = np.array([solidusPres(t) for t in temp_fill_list_liq])

    pres_fill_list_liq = np.logspace(1, np.log10(solidusPres(T_critical)))
    temp_fill_list_liq = [solidusTemp(p) for p in pres_fill_list_liq]

    # solid, liquid, supercritical, vapor, gas
    phase_colors = [
        (0.7, 0.7, 0.9),
        (0.25, 0.25, 0.75),
        (0.75, 0.25, 0.75),
        (0.75, 0.75, 0.75),
        (0.75, 0.75, 0.75),
    ]
    # Fill gas region
    ax.fill_between(
        [T_critical, xlims[1]],
        [P_critical, P_critical],  # , [1.0, P_critical],
        color=phase_colors[4],
        zorder=1,
    )
    # Fill scf region
    ax.fill_between(
        temp_fill_list_scf,
        pres_fill_list_scf,  # , [1.0, P_critical],
        P_critical,
        color=phase_colors[2],
        zorder=0,
    )

    # Fill liquid region
    ax.fill_between(
        temp_fill_list_liq,
        pres_fill_list_liq,  # , [1.0, P_critical],
        color=phase_colors[1],
        zorder=0,
    )
    # Fill vapor region
    ax.fill_between(x, vapor, color=phase_colors[3], zorder=0)  # , [1.0, P_critical],

    ax.set_facecolor(phase_colors[0])

    ax.set_axisbelow(True)
    # ax.grid(which = 'both', color = (.9, .9, .9), zorder = 0, linewidth = 2)
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\rm Temperature \ [K]$", fontsize=fnts2, labelpad=20.0)
    ax.set_ylabel(r"$\rm Pressure \ [bar]$", fontsize=fnts2, labelpad=20.0)
    xticks = matplotlib.ticker.FuncFormatter(
        lambda x, pos: r"${0:g}$".format(x / axis_scales[0])
    )
    ax.xaxis.set_major_formatter(xticks)
    yticks = matplotlib.ticker.FuncFormatter(
        lambda y, pos: r"$10^{{{0:d}}}$".format(int(np.log10(y / axis_scales[1])))
    )
    ax.yaxis.set_major_formatter(yticks)

    plt.savefig(
        "/home/os18o068/Documents/PHD/Abbildungen/Figures_paper_2/water_phase_diagram.pdf",
        format="pdf",
        bbox_inches="tight",
        dpi=320,
    )
    plt.close(fig)
    # return fig, ax
    # plt.close('all')
