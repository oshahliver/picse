# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 15:19:16 2018

@author: os18o068

This module incorporates the equation of state (EOS) for water ice Ih in 
conjunction with the IAPWS formulation (Feistel 2006)
Range over which the EOS is constrained by experimental data:
    teperature: 0-273.16 K
    pressure: 0-210 MPa
    
The results of this module have been verified by comparison with table 10 in 
Feistel 2006
"""

from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
import numpy as np
import time
from pics.utils import functionTools as ftool
from pics.materials import phase_transitions_water_Wagner2002 as phase

# ------------------------------------------------------------------------------
# physical constants

Rgas = 8.3144598  # universal gas constant
hbar = 1.0545718e-34
EHartree = 4.359745e-18  # Hartree energy in joule
kB = 1.38065e-23  # Joule K-1
NA = 6.0221409e23
mFe = 55.845e-3  # molar mass in kg
mO = 15.999e-3
mH = 1.00794e-3
mMg = 24.305e-3
mS = 32.065e-3
mSi = 28.0855e-3
mH2O = 2 * mH + mO

# ------------------------------------------------------------------------------
# EOS specific parameters

T_triple = 273.16  # triple point temperature
P_triple = 611.657  # triple point pressure
P0 = 1.01325e5  # normal pressure
T0 = 273.152519  # normal melting point
T0_melt = 273.152519  # normal melting temperature
s0 = -3327.337564  # 189.13 #Joule kg-1 K-1


# complex constants [re, im]
t1 = np.complex(3.68017112855051e-02, 5.10878114959572e-02)
r1 = np.complex(44.7050716285388, 65.6876847463481)
t2 = np.complex(0.337315741065416, 0.335449415919309)
r2_list = [
    np.complex(-72.597457432922, -78.100842711287),
    np.complex(-5.57107698030123e-05, 4.64578634580806e-05),
    np.complex(2.34801409215913e-11, -2.85651142904972e-11),
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


# set bisection pressure limits
P_min_bisec = 1.0e-4
P_max_bisec = 1.0e9

# set bisection temperature limits
T_min_bisec = 1.0e0
T_max_bisec = 1000.0

# set interpolation range for extremely low pressures
P_min_interpol = 1.0e3
P_max_interpol = 1.0e4
# ------------------------------------------------------------------------------
# functions


def g0(p):
    pi0 = P0 / P_triple
    pi = p / P_triple
    return sum([g0_list[k] * (pi - pi0) ** k for k in range(len(g0_list))])


def g0P(P):
    pi0 = P0 / P_triple
    pi = P / P_triple

    cc = 0.0
    for k in range(4):
        kk = k + 1
        cc += g0_list[kk] * kk / P_triple * (pi - pi0) ** (kk - 1)

    return cc


def g0PP(P):
    pi0 = P0 / P_triple
    pi = P / P_triple

    cc = 0.0
    for k in range(3):
        kk = k + 2
        cc += g0_list[kk] * kk * (kk - 1) / P_triple**2 * (pi - pi0) ** (kk - 2)

    return cc


def r2(P):
    pi0 = np.complex(P0 / P_triple, 0.0)
    pi = np.complex(P / P_triple, 0.0)
    return sum([r2_list[k] * (pi - pi0) ** k for k in range(len(r2_list))])


def r2P(P):
    pi0 = np.complex(P0 / P_triple, 0.0)
    pi = np.complex(P / P_triple, 0.0)

    cc = 0.0
    for k in range(2):
        kk = k + 1
        cc += r2_list[kk] * kk / P_triple * (pi - pi0) ** (kk - 1)

    return cc


def r2PP(P):
    return r2_list[2] * 2.0 / P_triple**2


def complex_number(p, t):
    tau = t / T_triple
    r_list = [r1, r2(p)]
    return sum(
        [
            r_list[k]
            * (
                (t_list[k] - tau) * np.log(t_list[k] - tau)
                + (t_list[k] + tau) * np.log(t_list[k] + tau)
                - 2.0 * t_list[k] * np.log(t_list[k])
                - tau**2 / t_list[k]
            )
            for k in range(len(t_list))
        ]
    )


def g(P, T):
    """This function computes the specific Gibbs free energy for pure water
    ice Ih accoring to Feistel 2009
    """
    tau = T / T_triple
    cc = complex_number(P, T)
    return g0(P) - s0 * T_triple * tau + T_triple * cc.real


def gT(P=None, T=None):
    tau = T / T_triple

    # compute complex number
    r_list = [r1, r2(P)]
    cc = np.complex(0.0, 0.0)
    for k in range(2):
        cc += r_list[k] * (
            -np.log(t_list[k] - tau) + np.log(t_list[k] + tau) - 2.0 * tau / t_list[k]
        )

    return -s0 + np.real(cc)


def gTT(T=None, P=None):
    tau = T / T_triple

    # compute complex number
    r_list = [r1, r2(P)]
    cc = np.complex(0.0, 0.0)
    for k in range(2):
        cc += r_list[k] * (
            1.0 / (t_list[k] - tau) + 1.0 / (t_list[k] + tau) - 2.0 / t_list[k]
        )

    return 1 / T_triple * np.real(cc)


def gP(T=None, P=None):
    tau = T / T_triple
    t2 = t_list[1]
    # compute complex number
    cc = r2P(P) * (
        (t2 - tau) * np.log(t2 - tau)
        + (t2 + tau) * np.log(t2 + tau)
        - 2 * t2 * np.log(t2)
        - tau**2 / t2
    )

    return g0P(P) + T_triple * np.real(cc)
    # return ftool.deriv(f=g, whicharg='p', x0=P, t=T, acc=1.0e-4)


def gPP(T=None, P=None):
    tau = T / T_triple
    t2 = t_list[1]
    # compute complex number
    cc = r2PP(P) * (
        (t2 - tau) * np.log(t2 - tau)
        + (t2 + tau) * np.log(t2 + tau)
        - 2 * t2 * np.log(t2)
        - tau**2 / t2
    )

    return g0PP(P) + T_triple * np.real(cc)


def gTP(T=None, P=None):
    tau = T / T_triple
    t2 = t_list[1]
    cc = r2P(P) * (-np.log(t2 - tau) + np.log(t2 + tau) - 2 * tau / t2)

    return np.real(cc)
    # return ftool.deriv(f=gP, whicharg='T', x0=T, P=P, acc=1.0e-4)


def beta(T=None, P=None, d=None):
    return -gTP(T=T, P=P) / gPP(T=T, P=P)


def u_spec(T=None, P=None, d=None):
    """Compute specific internal energy"""

    return g(P, T) - T * gT(T=T, P=P) - P * gP(T=T, P=P)


def density(P=None, T=None, **kwargs):
    return 1.0 / gP(T=T, P=P)


def pressure(d=None, T=None, **kwargs):
    temp = T
    return ftool.bisec(
        a=P_min_bisec,
        b=P_max_bisec,
        y=d,
        whicharg="P",
        f=density,
        T=temp,
        identity="Feistel pres",
        **kwargs
    )


def dPdrho_T(T=None, P=None, d=None):
    return -gP(T=T, P=P) ** 2 / gPP(T=T, P=P)


def temperature(d=None, p=None, **kwargs):
    pres = p
    return ftool.bisec(
        a=T_min_bisec,
        b=T_max_bisec,
        y=d,
        whicharg="T",
        f=density,
        P=pres,
        identity="Feistel temp",
        **kwargs
    )


def alpha_th_p(T=None, P=None, d=None):
    """Compute isobaric thermal expansion coefficientin 1/K"""
    return gTP(T=T, P=P) / gP(T=T, P=P)


def s_spec(T=None, P=None, d=None):
    """Compute specific entropy"""
    """
    if d==None and not T==None:
        s = -ftool.deriv(f=g, whicharg='t', x0=T, p=P,
                    identity='Feistel entropy', acc=1.0e-6)

    elif T==None and not d==None:
        T = temperature(d=d, p=P)
        s = -ftool.deriv(f=g, whicharg='t', x0=T, p=P,
                    identity='Feistel entropy')     
    """
    s = -gT(T=T, P=P)

    return s


def f(T=None, P=None, d=None):
    """Compute specific Helmholty free energy"""
    return g(T=T, P=P) - P * gP(T=T, P=P)


def cp_spec(T=None, P=None, d=None):
    """Compute specific isobaric heat capacity"""
    # dSdT_P = ftool.deriv(f=s_spec, whicharg='T', x0=T, P=P)

    return -T * gTT(T=T, P=P)


def cV_spec(T=None, P=None, d=None):
    """Compute specific isochoric heat capacity"""

    dSdT_V = ftool.deriv(f=s_spec, whicharg="T", x0=d, P=P)


def vS(T=None):
    """Compute speed of sound of P-waves according to Vogt 2008"""
    T0 = 273.15
    return 3837.9 - 2.812 * (T - T0)


def dTdP_S(T=None, P=None, d=None, type=1):
    """Compute adiabatic gradient. type==0 does't really work for some reason"""
    dTdP_S = None
    if d == None:
        d = density(P=P, T=T)

    if type == 0:
        # compute density derivative at given pressure
        dsdrho = ftool.deriv(f=s_spec, whicharg="d", x0=d, P=P)

        # compute dT/dP_S from Maxwell relation
        # dT/dP_S = dV/dS_P= = -1/rho*2 * (dS_spec/drho)**(-1)
        dTdP_S = -1.0 / (d**2 * dsdrho)

    elif type == 1:
        a = alpha_th_p(T=T, P=P, d=d)
        cp = cp_spec(T=T, P=P)
        dTdP_S = T * a / (d * cp)

    return dTdP_S


def together(T=None, P=None):

    P = min(5.0e8, P)

    K1 = gP(T=T, P=P)
    K2 = gT(T=T, P=P)
    K3 = gTT(T=T, P=P)
    K4 = gPP(T=T, P=P)
    K5 = gTP(T=T, P=P)
    K6 = g(T=T, P=P)

    d = 1.0 / K1
    u = K6 - T * K2 - P * K1
    s = -K2
    cp = -T * K3
    alpha = K5 / K1
    dPdrho = -(K1**2) / K4
    dTdP = T * alpha / (d * cp)

    return d, dTdP, dPdrho, alpha, cp, s, u


def test(T=None, P=None):
    t0 = time.time()
    a, b, c, d, e, f, g = together(P=P, T=T)
    t = time.time()
    print("elapsed time for simultaneous evaluation:", t - t0)

    t0 = time.time()
    dens = density(T=T, P=P)
    dPdrho = dPdrho_T(T=T, P=P, d=dens)
    dTdP = dTdP_S(T=T, P=P, d=dens)
    cp = cp_spec(T=T, P=P, d=dens)
    alpha = alpha_th_p(T=T, P=P, d=dens)
    u = u_spec(T=T, P=P, d=dens)
    s = s_spec(T=T, P=P, d=dens)
    t = time.time()
    print("elapsed time for seprate evaluation:", t - t0)

    print("d:", a, dens)
    print("dTdP_S:", b, dTdP)
    print("dPdhro_T:", c, dPdrho)
    print("alpha:", d, alpha)
    print("cp:", e, cp)
    print("s:", f, s)
    print("u:", g, u)


# this is the function that should be called from outside of the module
# it computes pressure, density or temperature as function of the other two
# parameters
def compute(what="dens", **kwargs):
    if what == "dens":
        return density(**kwargs)

    elif what == "pres":
        return pressure(**kwargs)

    elif what == "temp":
        return temperature(**kwargs)


def compare_s_spec(N=5):
    # compute profile as function of P and T
    pres = np.logspace(2, 8, N)
    temp = np.linspace(100.0, 150.0, 2)
    dens = []
    s_normal = []
    s_not_normal = []
    for t in range(len(temp)):
        print("\nT=", temp[t])
        s_normal.append([])
        s_not_normal.append([])
        dens.append([])
        for i in range(N):
            s_normal[t].append(s_spec(T=temp[t], P=pres[i]))

            d = density(t=temp[t], p=pres[i])
            dens[t].append(d)
            print("d=", d)
            try:
                s_not_normal[t].append(s_spec(P=pres[i], d=d))
            except RuntimeError:
                s_not_normal[t].append(None)

    fig, ax = plt.subplots(1, 2)
    for t in range(len(temp)):
        ax[0].semilogx(pres, s_normal[t], label=str(temp[t]))
        ax[0].scatter(pres, s_not_normal[t])
        ax[1].plot(dens[t], s_not_normal[t])
    ax[0].legend()


"""
#plot some stuff
fig, ax = plt.subplots(1, 2, sharey=True)
ax1, ax2 = ax
fig.subplots_adjust(wspace = .1)


a = compute('dens', P = 1.0e8, T = 273.)
print ('a =', a)
b = compute('pres', d = a, T = 273.)
print ('b =', b)

N=45

P_min=1.0e-3
P_max=1.0e9

temp_list = [50., 100., 200., 250, 260, 270, 273.16]
pres_list = np.logspace(np.log10(P_min), np.log10(P_max), N)

scatter_list = []
for temp in temp_list:
    dens_list = []
    for pres in pres_list:
        check = phase.phaseCheck(pres, temp)
        if check == 'solid':
            dens_list.append(compute('dens', p = pres, t = temp, noisy=False))
        else:
            dens_list.append(None)
    ax1.semilogy(dens_list, pres_list, label = 'T = '+ str(int(temp)) + ' K',
                 linewidth = 2, zorder = 2)


dens_list = [930, 931, 932, 933, 934, 935, 940, 960, 980]

temp_list = np.linspace(0., 350., N)

for dens in dens_list:
    pres_list = []
    for temp in temp_list:
        p = (compute('pres', d = dens, T = temp))
        check = phase.phaseCheck(p, temp)
        if check == 'solid':
            pres_list.append(p)
        else:
            pres_list.append(None)
    ax2.semilogy(temp_list, pres_list, label = r'$\rho = $'+str(round(dens,2)),
                 linewidth = 2, zorder = 2)


tosolid_list = phase.solidLine(log = True)
tovapor_list = phase.vaporLine(log = True)

ax2.semilogx(tosolid_list[0], tosolid_list[1], color = 'grey', linestyle = '--', 
             zorder = 1)

ax2.plot(tovapor_list[0], tovapor_list[1], color = 'grey', linestyle = '--', 
             zorder = 1)
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

ax1.set_title(r'$Water \ ice \ isotherms$')
ax1.set_xlabel(r'$\rm Density \ [kg \ m^{-3}]$')
ax1.set_ylabel(r'$\rm Pressure \ [Pa]$')
ax1.set_ylim(P_min, P_max)
ax1.set_xlim(900, 1000)

ax1.tick_params(which = 'both', axis = 'both', direction = 'in', top = 'on',
			left = 'on', right='on')


ax1.tick_params(which ='major', axis='both', length=10)
ax1.tick_params(which ='minor', axis='both', length=5)

ax2.set_title(r'$Water \ ice \ isochors$')
ax2.set_xlabel(r'$\rmTemperature \ [K]$')
ax2.yaxis.set_label_position('right')
ax2.yaxis.tick_right()
ax2.set_xlim(0, 1.0e3)

ax2.tick_params(which = 'both', axis = 'both', direction = 'in', top = 'on',
			left = 'on', right='on')

ax2.yaxis.set_tick_params(which='both', labelright=True)

ax2.tick_params(which = 'major', axis = 'y', length = 10)
ax2.tick_params(which = 'minor', axis = 'y', length = 5)
ax2.tick_params(which = 'major', axis = 'x', length = 10)
ax2.tick_params(which = 'minor', axis = 'x', length = 5)

ax2.set_axisbelow(True)
ax1.grid(which = 'both', color = (.8, .8, .8))
ax2.grid(which = 'both', color = (.8, .8, .8))


ax1.legend()
ax2.legend()
plt.show()
"""
