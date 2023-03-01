# this is a simple python program that uses a fortran subroutine
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from pics.materials import Mazevet2018EOS
from pics.materials import Wagner2002EOS
from pics.materials import French2015EOS
from pics.materials import Feistel2006EOS
from pics.materials import phase_transitions_water_Wagner2002 as phaseTrans
import iapws
import seafreeze
import time
from pics.utils.function_tools import functionTools as ftool
import sys
from pics.utils.function_tools import analytical_functions as anfct
import matplotlib.ticker
from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pics.physicalparams import (
    T_triple,
    T_critical,
    P_triple,
    P_critical,
    mO,
    mH,
    kB,
    NA,
    Rgas,
)

from pics.runparams import eos_pres_trans
from pics.materials.phase_transitions_water_Wagner2002 import data_schwager

T_triple = 273.16  # triple point temperature
P_triple = 611.657  # triple point pressure in Pa
T_critical = 647.096  # temperature at critical point
P_critical = 22.064e6  # pressure at critical point

# -----------------------------------------------------------------------------
# main
# -----------------------------------------------------------------------------
def IAPWS_density(P=None, T=None, **kwargs):

    # Ideal gas law as IAPWS gives spurious behaviour at low P
    if P < 1.0e3:
        d = P * (2 * mH + mO) / (kB * T * NA)

    else:
        d = iapws.IAPWS95(T=T, P=P * 1.0e-6).rho

    return d


def phase(P=None, T=None, d=None, **kwargs):
    return phaseTrans.phaseCheck(p=P, t=T, d=d, **kwargs)


def get_phase(P, T):

    phh = phaseTrans.phaseCheck(p=P, t=T)
    ph = phh
    phase_strings = ["Ih", "Mazevet", "French", "V", "VI", "II", "III", "Wagner"]
    if type(ph) is str:

        if ph == "Ih":
            ph = 0

        elif ph == "Mazevet":
            ph = 1

        elif ph == "French":
            ph = 2

        elif ph == "V":
            ph = 3

        elif ph == "VI":
            ph = 4

        elif ph == "II":
            ph = 5

        elif ph == "III":
            ph = 6

        else:
            ph = 7

    else:
        phh = phase_strings[ph]
    return ph, phh


def together(T=None, P=None, ph=None):
    """Computes all relevant parameters for eos table"""
    if ph == None:
        # phh = phaseTrans.phaseCheck(p=P, t=T)
        # ph = phh
        ph, phh = get_phase(P, T)

    # phase_strings = ["Ih", "Mazevet", "French", "V", "VI", "II", "III", "Wagner"]
    # if type(ph) is str:
    #     phh = ph
    #     if ph == "Ih":
    #         ph = 0

    #     elif ph == "Mazevet":
    #         ph = 1

    #     elif ph == "French":
    #         ph = 2

    #     elif ph == "V":
    #         ph = 3

    #     elif ph == "VI":
    #         ph = 4

    #     elif ph == "II":
    #         ph = 5

    #     elif ph == "III":
    #         ph = 6

    #     else:
    #         ph = 7

    # else:
    #     phh = phase_strings[ph]

    # Ice Ih Feistel 2006
    if ph == 0:
        a, b, c, d, e, f, g = Feistel2006EOS.together(T=T, P=P)

    # high pressure Mazevet 2018
    elif ph == 1:
        a, b, c, d, e, f, g = Mazevet2018EOS.together(T=T, P=P)

    # high pressure French 2015
    elif ph == 2:
        a, b, c, d, e, f, g = French2015EOS.together(T=T, P=P)

    # high pressure Journaux 2020
    elif ph == 3 or ph == 4 or ph == 5 or ph == 6:
        PT = np.empty((1,), np.object)
        PT[0] = (P * 1.0e-6, T)

        out = seafreeze.seafreeze(PT, phh)

        (a,) = out.rho
        (b,) = T * out.alpha / (out.Cp * out.rho)
        (c,) = out.Kt / out.rho * 1e6
        (d,) = out.alpha
        (e,) = out.Cp
        (f,) = out.S
        (g,) = out.U

    # IAPWS Wagner 2009 & 2016
    elif ph == 7:
        a, b, c, d, e, f, g = Wagner2002EOS.together(T=T, P=P)

    # rho, dTdP_S, dPdrho_T, alpha, cp_spec, s_spec, u_spec
    return [a, b, c, d, e, f, g, ph]


# def Plot(
#     P_min=10.0,
#     P_max=1.0e13,
#     T_min=100.0,
#     T_max=1.0e4,
#     log=True,
#     res=3,
#     data_log=False,
#     Nx_ticks=10,
#     Ny_ticks=10,
#     params=[0],
# ):
#     labels = [
#         r"$\rho \ \rm [kg \ m^{-3}]$",
#         r"$\log\left(dT/dP_S / \rm \left(K \ Pa^{-1}\right)\right)$",
#         r"$dP/d\rho_T$",
#         r"$\alpha_{\rm th} \ \rm [10^{-5} \ K^{-1}]$",
#         r"$c_P$",
#         r"$s$",
#         r"$u$",
#     ]
#     factors = ([1.0, 200.0, 500.0],)

#     scales = ["lin", "log", "log", "lin", "lin", "lin", "lin"]

#     N = 2**res

#     if log:
#         x = np.logspace(np.log10(T_min), np.log10(T_max), N)
#         y = np.logspace(np.log10(P_min), np.log10(P_max), N)
#         # x = np.logspace(np.log10(T_min), np.log10(T_max), N)

#         # Nx_ticks = int(np.log10(T_max)) - int(np.log10(T_min)) + 1
#         Ny_ticks = int(np.log10(P_max)) - int(np.log10(P_min)) + 1

#         x_ticks = np.linspace(0, N, Nx_ticks)
#         y_ticks = np.linspace(0, N, Ny_ticks)

#         # x_tick_labels = np.logspace(np.log10(T_min), np.log10(T_max), Nx_ticks)
#         dT = (T_max - T_min) / (Nx_ticks - 1)
#         x_tick_labels = [int(T_min + dT * i) for i in range(Nx_ticks)]
#         x_tick_labels = np.linspace(np.log10(T_min), np.log10(T_max), Nx_ticks)
#         x_tick_labels = [
#             r"$10^{}$".format("{" + str(int(xx)) + "}") for xx in x_tick_labels
#         ]
#         y_tick_labels = np.linspace(np.log10(P_min), np.log10(P_max), Ny_ticks)
#         y_tick_labels = [
#             r"$10^{}$".format("{" + str(int(yy)) + "}") for yy in y_tick_labels
#         ]

#         print(y_tick_labels)
#     else:
#         y = np.linspace(P_min, P_max, N)
#         x = np.linspace(T_min, T_max, N)

#         x_ticks = np.linspace(0, N, Nx_ticks)
#         y_ticks = np.linspace(0, N, Ny_ticks)

#         dT = (T_max - T_min) / (Nx_ticks - 1)
#         dP = (P_max - P_min) / (Ny_ticks - 1)

#         x_tick_labels = [int(T_min + dT * i) for i in range(Nx_ticks)]
#         y_tick_labels = [int(1.0e-9 * (P_min + dP * j)) for j in range(Ny_ticks)]

#     zz = np.zeros([N, N])
#     data = np.empty([N, N, len(params)])
#     for i in range(N):
#         for j in range(N):
#             dat = together(T=x[i], P=y[j])
#             data[i][j] = np.array([dat[p] for p in params])

#     for i in range(len(params)):
#         p = params[i]
#         fig, ax = plt.subplots()
#         if scales[p] == "lin":
#             plot_dat = data.T[i] * factors[p]
#         elif scales[p] == "log":
#             plot_dat = np.log10(data.T[i] * factors[p])
#         im = ax.imshow(plot_dat, origin="lower", cmap="terrain")
#         divider = make_axes_locatable(ax)
#         cax = divider.append_axes("right", size="5%", pad=0.5)
#         cbar = fig.colorbar(im, cax=cax)
#         cbar.set_label(labels[params[i]])

#         ax.set_xticks(x_ticks)
#         ax.set_yticks(y_ticks)

#         ax.set_xticklabels(x_tick_labels)
#         ax.set_yticklabels(y_tick_labels)

#         ax.set_xlabel(r"$\rm Temperature \ [K]$")
#         ax.set_ylabel(r"$\rm Pressure \ [Pa]$")
#     return fig, ax


# def analyze(
#     P_min=10.0,
#     P_max=1.0e12,
#     T_min=100.0,
#     T_max=1.0e5,
#     P_log=True,
#     T_log=True,
#     analytical=None,
#     order=2,
#     resolution=16,
#     param=2,
#     **kwargs
# ):

#     import eosTables

#     watab = eosTables.Table()
#     watab.load("eos_0.tab")

#     # compute number of grid points along each axis
#     N_T = resolution
#     N_P = resolution

#     # initiated pressure and temperature arrays in log or linear scales
#     if T_log:
#         nnT = int(np.log10(T_max)) - int(np.log10(T_min)) + 1
#         T_list = np.logspace(np.log10(T_min), np.log10(T_max), N_T)
#         major_ticks_x = np.linspace(0, 1, nnT)
#         major_tick_labels_x = np.linspace(
#             np.log10(T_min), np.log10(T_max), nnT, dtype=int
#         )

#     else:
#         T_list = np.linspace(T_min, T_max, N_T)
#         major_ticks_x = np.linspace(0, 1, 11)
#         major_tick_labels_x = np.linspace(int(T_min), int(T_max), 11)

#     if P_log:
#         nnP = int(np.log10(P_max)) - int(np.log10(P_min)) + 1
#         P_list = np.logspace(np.log10(P_min), np.log10(P_max), N_P)
#         major_ticks_y = np.linspace(0, 1, nnP)
#         major_tick_labels_y = np.linspace(
#             np.log10(P_min), np.log10(P_max), nnP, dtype=int
#         )

#     else:
#         P_list = np.linspace(P_min, P_max, N_P)
#         major_ticks_y = np.linspace(0, 1, 6)
#         major_tick_labels_y = np.linspace(int(P_min * 1.0e-9), int(P_max * 1.0e-9), 6)

#     try:
#         if analytical == None:
#             dat_analytical = np.zeros([len(T_list), len(P_list)])

#     except ValueError:
#         dat_analytical = analytical

#     dat_table = np.zeros([len(T_list), len(P_list)])
#     dat_residual = np.zeros([len(T_list), len(P_list)])

#     for i in range(len(T_list)):
#         for j in range(len(P_list)):
#             dat_table[i][j] = watab.interpolate(
#                 order=order, x=T_list[i], y=P_list[j], param=param
#             )

#             try:
#                 if analytical == None:
#                     if param == 2:
#                         dat_analytical[i][j] = density(T=T_list[i], P=P_list[j])

#                     elif param == 4:
#                         dat_analytical[i][j] = dPdrho_S(t=T_list[i], p=P_list[j])

#             except ValueError:
#                 pass

#             dat_residual[i][j] = (
#                 dat_analytical[i][j] - dat_table[i][j]
#             ) / dat_analytical[i][j]

#             """
#             #print (dat_residual[i][j])
#             if abs(dat_residual[i][j]) > 1.:
#                 print ('weird value encountered at T=', str(round(T_list[i],2))+' P='+
#                        str(round(P_list[j]/(10**int(np.log10(P_list[j]))),2))+'e'+
#                        str(int(np.log10(P_list[j]))),
#                        ': table=', round(dat_table[i][j],2),
#                        'analytical=', round(dat_analytical[i][j],2))
#             """

#     print("maximal deviation before scaling=", np.nanmax(abs(dat_residual)) * 100, "%")
#     # transform residuals for logarithmic plotting with negative and positive
#     # values
#     def_zero = -5

#     with np.nditer(dat_residual, op_flags=["readwrite"]) as it:
#         for x in it:
#             if np.isnan(x):
#                 pass

#             else:
#                 if abs(x) < 10**def_zero:
#                     x[...] = 10**def_zero * x / abs(x)

#                 elif abs(x) > 1.0:
#                     x[...] = x / abs(x)

#     dat_max = np.nanmax(dat_analytical)
#     dat_min = np.nanmin(dat_analytical)
#     cbar_max_log = int(np.log10(dat_max))
#     cbar_min_log = int(np.log10(dat_min))
#     cbar_max_lin = 10 ** (int(np.log10(dat_max)))
#     print("maximal deviation=", np.nanmax(abs(dat_residual)) * 100, "%")

#     dat_residual = (
#         (np.log10(abs(dat_residual)) - def_zero) * dat_residual / abs(dat_residual)
#     )

#     # plot data
#     fig, axes = plt.subplots(2, 2)

#     # set axis labels
#     axes[1][0].set_xlabel(r"$log_{10}(T) \ [K]$")
#     axes[0][0].set_ylabel(r"$log_{10}(P) \ [Pa]$")
#     axes[1][1].set_xlabel(r"$log_{10}(T) \ [K]$")
#     axes[1][0].set_ylabel(r"$log_{10}(P) \ [Pa]$")

#     im00 = axes[0][0].imshow(
#         dat_analytical.T,
#         cmap="rainbow",
#         extent=[0, 1, 0, 1],
#         norm=None,
#         aspect="equal",
#         origin="lower",
#         vmin=0.0,
#         vmax=cbar_max_lin,
#     )

#     im01 = axes[0][1].imshow(
#         np.log10(dat_analytical.T),
#         cmap="rainbow",
#         norm=None,
#         aspect="equal",
#         origin="lower",
#         extent=[0, 1, 0, 1],
#         vmin=cbar_min_log,
#         vmax=cbar_max_log,
#     )

#     im10 = axes[1][0].imshow(
#         dat_table.T,
#         cmap="rainbow",
#         extent=[0, 1, 0, 1],
#         norm=None,
#         aspect="equal",
#         origin="lower",
#         vmin=0.0,
#         vmax=cbar_max_lin,
#     )

#     im11 = axes[1][1].imshow(
#         np.log10(dat_table.T),
#         cmap="rainbow",
#         norm=None,
#         origin="lower",
#         vmin=cbar_min_log,
#         extent=[0, 1, 0, 1],
#         vmax=cbar_max_log,
#     )

#     minor_ticks_x = []
#     minor_ticks_y = []
#     T_range = int(np.log10(T_max)) - int(np.log10(T_min))

#     # compute locations of minor ticks on logscale
#     for decs in range(T_range):
#         for ticks in range(8):
#             minor_ticks_x.append((np.log10(ticks + 2) + decs) / (T_range))

#     P_range = int(np.log10(P_max)) - int(np.log10(P_min))
#     for decs in range(P_range):
#         for ticks in range(8):
#             minor_ticks_y.append((np.log10(ticks + 2) + decs) / (P_range))

#     im = [[im00, im01], [im10, im11]]
#     i = 0
#     for axx in axes:
#         j = 0
#         for ax in axx:
#             ax.set_xticks([])
#             ax.set_yticks([])
#             ax.set_xticks(minor_ticks_x, minor=True)
#             ax.set_yticks(minor_ticks_y, minor=True)

#             ax.set_xticks([0.0, 0.5, 1.0])
#             ax.set_yticks(major_ticks_y)
#             ax.set_xticks(major_ticks_x)

#             ax.set_xticklabels(major_tick_labels_x)
#             ax.set_yticklabels(major_tick_labels_y)

#             ax.tick_params(top=True, right=True, which="both")
#             ax.tick_params(size=10, which="major")
#             ax.tick_params(size=5, which="minor")

#             divider = make_axes_locatable(ax)
#             cax = divider.append_axes("right", size="5%", pad=0.5)
#             cbar = fig.colorbar(im[i][j], cax=cax)
#             if j == 0:
#                 cbar.set_label(r"$\rho \ [kg \ m^{-3}]$")
#             else:
#                 cbar.set_label(r"$log_{10}(\rho) \ [kg \ m^{-3}]$")
#             j += 1
#         i += 1

#     # plot residuals
#     fig, ax = plt.subplots()

#     # compute colorbar tick positions
#     cbar_ticks = np.linspace(def_zero, -def_zero, 2 * abs(def_zero) + 1)

#     # set colorscheme for residuals
#     cmap = mpl.colors.ListedColormap(
#         [
#             (1.0, 0.5, 0.5),
#             (1.0, 0.0, 0.0),
#             (0.6, 0.0, 0.0),
#             (0.2, 0.0, 0.0),
#             (0.0, 0.0, 0.0),
#             (0.0, 0.0, 0.0),
#             (0.0, 0.2, 0.0),
#             (0.0, 0.6, 0.0),
#             (0.0, 1.0, 0.0),
#             (0.5, 1.0, 0.5),
#         ]
#     )

#     bounds = cbar_ticks
#     norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

#     im = ax.imshow(
#         dat_residual.T,
#         cmap=cmap,
#         extent=(0, 1, 0, 1),
#         norm=norm,
#         aspect="equal",
#         origin="lower",
#         vmin=def_zero,
#         vmax=-def_zero,
#     )

#     ax.set_xticks([])
#     ax.set_yticks([])
#     ax.set_xticks(minor_ticks_x, minor=True)
#     ax.set_yticks(minor_ticks_y, minor=True)

#     ax.set_xticks(major_ticks_x)
#     ax.set_yticks(major_ticks_y)

#     ax.set_xticklabels(major_tick_labels_x)
#     ax.set_yticklabels(major_tick_labels_y)

#     ax.tick_params(top=True, right=True, which="both")
#     ax.tick_params(size=10, which="major")
#     ax.tick_params(size=5, which="minor")

#     cbar_ticks = np.linspace(def_zero, -def_zero, 2 * abs(def_zero) + 1)
#     cbar_ticklabels = []
#     for c in range(abs(def_zero)):
#         cbar_ticklabels.append(10 ** (-c) * 100)

#     cbar_ticklabels.append(0)
#     for c in range(abs(def_zero)):
#         cbar_ticklabels.append(10 ** (def_zero + c + 1) * 100)

#     ax.set_xlabel(r"$log_{10}(T) \ [K]$")
#     ax.set_ylabel(r"$log_{10}(P) \ [Pa]$")

#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="5%", pad=0.5)
#     cbar = fig.colorbar(im, cax=cax, ticks=cbar_ticks)
#     cbar.set_label(r"$residuals \ [\%]$")
#     cbar.ax.set_yticklabels(
#         cbar_ticklabels,
#     )
#     return dat_analytical, dat_table, dat_residual


# def replot(
#     dat_analytical,
#     dat_table,
#     dat_residual,
#     P_min=10.0,
#     P_max=1.0e12,
#     T_min=100.0,
#     T_max=1.0e5,
#     P_log=True,
#     T_log=True,
#     color_scheme="rainbow",
#     **kwargs
# ):

#     rho_max = np.nanmax(dat_analytical)
#     rho_min = np.nanmin(dat_analytical)
#     cbar_max_log = int(np.log10(rho_max))
#     cbar_min_log = int(np.log10(rho_min))
#     cbar_max_lin = 10 ** (int(np.log10(rho_max)))

#     if T_log:
#         major_ticks_x = matplotlib.ticker.LinearLocator(3)
#         major_tick_labels_x = np.linspace(
#             np.log10(T_min), np.log10(T_max), 3, dtype=int
#         )

#     else:
#         major_ticks_x = np.linspace(0, 1, 3)
#         major_tick_labels_x = np.linspace(int(T_min), int(T_max), 3)

#     if P_log:
#         major_ticks_y = np.linspace(0, 1, 12)
#         major_tick_labels_y = np.linspace(
#             np.log10(P_min), np.log10(P_max), 12, dtype=int
#         )

#     else:
#         major_ticks_y = np.linspace(0, 1, 12)
#         major_tick_labels_y = np.linspace(int(P_min * 1.0e-9), int(P_max * 1.0e-9), 12)

#     fig, axes = plt.subplots(2, 2)
#     im00 = axes[0][0].imshow(
#         dat_analytical.T,
#         cmap="rainbow",
#         extent=(0, 1, 0, 1),
#         norm=None,
#         aspect="equal",
#         origin="lower",
#         vmin=0.0,
#         vmax=cbar_max_lin,
#     )

#     im01 = axes[0][1].imshow(
#         np.log10(dat_analytical.T),
#         cmap="rainbow",
#         extent=(0, 1, 0, 1),
#         norm=None,
#         aspect="equal",
#         origin="lower",
#         vmin=cbar_min_log,
#         vmax=cbar_max_log,
#     )

#     im10 = axes[1][0].imshow(
#         dat_table.T,
#         cmap="rainbow",
#         extent=(0, 1, 0, 1),
#         norm=None,
#         aspect="equal",
#         origin="lower",
#         vmin=0.0,
#         vmax=cbar_max_lin,
#     )

#     im11 = axes[1][1].imshow(
#         np.log10(dat_table.T),
#         cmap="rainbow",
#         extent=(0, 1, 0, 1),
#         norm=None,
#         aspect="equal",
#         origin="lower",
#         vmin=cbar_min_log,
#         vmax=cbar_max_log,
#     )

#     im = [[im00, im01], [im10, im11]]
#     i = 0
#     for axx in axes:
#         j = 0
#         for ax in axx:
#             # ax.set_xticks(major_ticks_x, minor=True)
#             ax.xaxis.set_major_locator(major_ticks_x)
#             ax.set_yticks(major_ticks_y)

#             ax.set_xticklabels(major_tick_labels_x)
#             ax.set_yticklabels(major_tick_labels_y)

#             ax.set_xlabel(r"$log_{10}(T) \ [K]$")
#             ax.set_ylabel(r"$log_{10}(P) \ [Pa]$")

#             divider = make_axes_locatable(ax)
#             cax = divider.append_axes("right", size="5%", pad=0.5)
#             cbar = fig.colorbar(im[i][j], cax=cax)

#             if j == 0:
#                 cbar.set_label(r"$\rho \ [kg \ m^{-3}]$")
#             else:
#                 cbar.set_label(r"$log_{10}(\rho) \ [kg \ m^{-3}]$")
#             j += 1
#         i += 1

#     fig, ax = plt.subplots()
#     im = ax.imshow(
#         dat_residual.T,
#         cmap="gist_stern",
#         extent=(0, 1, 0, 1),
#         norm=None,
#         aspect="equal",
#         origin="lower",
#         vmin=-1,
#         vmax=1,
#     )

#     ax.xaxis.set_major_locator(major_ticks_x)
#     ax.set_yticks(major_ticks_y)

#     ax.set_xticklabels(major_tick_labels_x)
#     ax.set_yticklabels(major_tick_labels_y)

#     cbar_ticks = np.linspace(-5, 5, 11)
#     cbar_ticklabels = [
#         -1,
#         -0.1,
#         -0.01,
#         -0.001,
#         -0.0001,
#         -0.00001,
#         0.0,
#         0.00001,
#         0.0001,
#         0.001,
#         0.01,
#         1,
#     ]
#     print(cbar_ticks)
#     ax.set_xlabel(r"$log_{10}(T) \ [K]$")
#     ax.set_ylabel(r"$log_{10}(P) \ [Pa]$")
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="5%", pad=0.5)
#     cbar = fig.colorbar(im, cax=cax, ticks=cbar_ticks)
#     cbar.set_label(r"$residual$")
#     cbar.ax.set_yticklabels(
#         cbar_ticklabels,
#     )


# #plot comparison between iapws, ideal gas and van der waals
# temp=np.linspace(410.5, 1000, 50)
# pres_igl=1.0e3
# pres_vdw=1.0e4
# dens_iapws_high=np.array([IAPWS_density(P = pres_vdw, T = t) for t in temp])
# dens_iapws_low=np.array([IAPWS_density(P = pres_igl, T = t) for t in temp])
# dens_vdw=np.array([anfct.rho_VDW(T=t, P=pres_vdw, ll=5) for t in temp])
# dens_igl=np.array([anfct.rho_IGL(T=t, P=pres_igl, ll=5) for t in temp])

# fig, ax=plt.subplots(1,2)
# ax[0].plot(temp, dens_igl)
# ax[0].plot(temp, dens_vdw)
# ax[0].plot(temp, dens_iapws_high)
# ax[0].plot(temp, dens_iapws_low)
# ax[0].set_ylim(0., .5)

# res_vdw=(dens_iapws_high-dens_vdw)/dens_iapws_high
# res_igl=(dens_iapws_low-dens_igl)/dens_iapws_low

# ax[1].semilogy(temp, res_vdw)
# ax[1].semilogy(temp, res_igl)

# print ('maximal deviation vdw:', round(max(abs(res_vdw*100)),3),'%')
# print ('maximal deviation igl:', round(max(abs(res_igl*100)),3),'%')

# #plot eos over entire PT range specified
# grid_color =(.9,.9,.9)
# plot_color = (.4,.5,.7)
# phase_color = (.2,.3,.5)
# color_list = [plot_color, (.2, .8, .8), (.2, .8, .2), (.5, .5, 0.), \
#               (.8, .4, 0.), 'grey', 'k']

# font1 = {'family': 'sans',
#         'color':  phase_color,
#         'weight': 'bold',
#         'size': 10,
#         }

# font2= {'family': 'sans',
#         'color':  'k',
#         'weight': 'normal',
#         'size': 14,
#         }

# font3= {'family': 'sans',
#         'color':  phase_color,
#         'weight': 'bold',
#         'size': 6,
#         'backgroundcolor': 'None'
#         }


# N = 150
# fig, ax = plt.subplots(1, 2)
# ax1, ax2 = ax
# fig.subplots_adjust(wspace = .1)

# lwdth1 = 1.5
# lwdth = .9

# #plot density as function of pres for different isotherms on the lhs
# temp_list = [300., 500., 550., 600., 1000., 5000., 10000.]
# pres_list = np.logspace(np.log10(1.0e6), np.log10(1.0e12), N)
# ax1label_list = [[], []]
# ax1plot_list = [[], []]

# satcurve_list = []

# #high pressure water according to Mazevet 2018
# for temp in temp_list:
#     ii = temp_list.index(temp)
#     dens_list = []
#     for pres in pres_list:
#         dens = Mazevet2018EOS.density(P = pres, T = temp)
#         dens_list.append(dens)

#     plot, = ax1.loglog(dens_list, pres_list, linewidth = lwdth1, zorder = 2,
#                  color = color_list[ii], linestyle = '--')
#     #ax1plot_list[0].append(plot)
#     #ax1label_list[0].append('T = '+str(int(temp))+' K')

# ax1plot_list[1].append(plot)

# #------------------------------------------------------------------------------
# #water according to IAPWS 2016

# dsat_list = []
# psat_list = []
# for temp in temp_list:
#     ii = temp_list.index(temp)
#     dens_list = []
#     try:
#         psat = iapws.iapws97._PSat_T(temp)*1.0e6
#         dsat = iapws.IAPWS95(T = temp, P = psat*1.0e-6).rho
#         dsat_list.append(dsat)
#         psat_list.append(psat)

#         if temp > 400:
#             satcurve_list.append([temp, psat])

#     except NotImplementedError:
#         pass

#     for pres in pres_list:
#         try:
#             if pres < 1.0e10:
#                 dens = iapws.IAPWS95(T = temp, P = pres*1.0e-6).rho

#             else:
#                 dens = None
#             dens_list.append(dens)

#         except OverflowError:
#             pass

#         except TypeError:
#             pass
#     #np.sort(pres_list)
#     #np.sort(np.asarray(dens_list))
#     plot, = ax1.plot(dens_list, pres_list, linewidth = lwdth1, zorder = 2,
#                  color = color_list[ii], linestyle = '-', markersize = 10,
#                  marker = '', markevery=2)
#     ax1plot_list[0].append(plot)
#     ax1label_list[0].append('T = '+str(int(temp))+' K')

# ax1plot_list[1].append(plot)

# sc, = ax1.plot(dsat_list, psat_list, zorder = 3, color = 'k', marker = 'x', \
#                linestyle = 'none', markersize=8)

# ax1plot_list[1].append(sc)

# ax1label_list[1] = ['ab initio (Mazevet 2018)', 'IAPWS95 (Wagner 2009)', 'saturation points']


# #plot isochores on the rhs
# ax2label_list = []
# ax2plot_list = []

# #------------------------------------------------------------------------------
# #ice Ih according to Feistel 2006

# dens_list = [930, 933, 933.5, 934, 934.5, 935, 940, 945]
# posx_list = [50., 80., 150., 20., 20., 20., 20., 20., 20.]
# posy_list = [2.0e2, 2.0e2, 2.0e2, 2.0e6, 6.5e6, 1.4e7, 6.0e7, 1.0e8]

# temp_list = np.logspace(1, np.log10(300), N)
# pres_list = []
# for ii in range(len(dens_list)):
#     dens = dens_list[ii]
#     pres_list.append([])
#     for temp in temp_list:
#         pres = Feistel2006EOS.compute(what='pres', d = dens, T = temp, type = 'tangent', \
#                                acc = 1.0e-6)

#         state = phaseTrans.phaseCheck(pres, temp)

#         #print ('dens:', dens, 'pres:',  round(pres*1.0e-9,4), 'temp=', temp, state)
#         if state == 'ice Ih' or state == 'high pressure French':
#             pass

#         else:
#             pres = None

#         pres_list[ii].append(pres)

#     plot, = ax2.plot(temp_list, pres_list[ii], linestyle = ':',
#                      color = plot_color, linewidth = lwdth)

#     if ii == 0:
#         ax2plot_list.append(plot)
#         ax2label_list.append('Feistel 2006')

# for i in range(len(posy_list)):
#     dens = dens_list[i]
#     dens_str = str(dens)
#     ax2.text(posx_list[i], posy_list[i], s = dens_str, fontdict = font3)

# #------------------------------------------------------------------------------
# #liquid and vapor according to Wagner 2009

# dens_list = [.001, .01, .1, 1, 10, 100, 600, 800, 1000, 1200, 1500, 2000]
# posx_list = [2.0e4, 2.0e4, 2.0e4, 2.0e4, 2.0e4, 2.0e4, 9.0e2, 6.0e2, 3.6e2,
#              4.0e2, 1.0e3, 2.0e3]
# posy_list = [5.0e3, 5.0e4, 5.0e5, 5.0e6, 5.0e7, 5.0e8, 1.0e8, 2.0e8, 3.0e8,
#              8.0e8, 0.6e10, 2.0e10]

# temp_list = np.logspace(np.log10(273.16), 5, N)
# pres_list = []
# for dens in dens_list:
#     ii = dens_list.index(dens)
#     pres_list.append([])
#     for temp in temp_list:
#         pres = iapws.IAPWS95(T = temp, rho = dens).P*1.0e6
#         state = phaseTrans.phaseCheck(pres, temp)

#         if state == 'ice Ih' or state == 'high pressure French':
#             pres = None

#         pres_list[ii].append(pres)

#     plot, = ax2.plot(temp_list, pres_list[ii], linestyle = '-', label = str(dens),
#     color = plot_color, linewidth = lwdth)
#     if ii == 0:
#         ax2plot_list.append(plot)
#         ax2label_list.append('Wagner 2009')

# for i in range(len(posy_list)):
#     dens = dens_list[i]
#     dens_str = str(dens)
#     ax2.text(posx_list[i], posy_list[i], s = dens_str, fontdict = font3)

# #------------------------------------------------------------------------------
# #high pressure regime according to Mazevet 2018

# dens_list = [1, 10, 100, 600, 800, 1000, 1200, 1500, 2000, 3000, 4000, 6000, 8000,
#              10000, 20000]
# posx_list = [None, None, None, None, None, None, None, None, None, 3.0e3, 3.0e3,
#              3.0e3, 3.0e3, 3.0e3, 3.0e3]
# posy_list = [None, None, None, None, None, None, None, None, None, 1.0e11, 2.0e11,
#              7.1e11, 1.6e12, 3.0e12, 1.5e13]

# temp_list = np.logspace(np.log10(300), 5, N)
# pres_list = []
# for dens in dens_list:
#     ii = dens_list.index(dens)
#     pres_list.append([])
#     for temp in temp_list:
#         pres = Mazevet2018EOS.compute('pres', d = dens, T = temp)
#         state = phaseTrans.phaseCheck(pres, temp)

#         if not state == 'high pressure Mazevet':
#             pres = None

#         pres_list[ii].append(pres)

#     plot, = ax2.plot(temp_list, pres_list[ii], linestyle = '--',
#                      color = plot_color, linewidth = lwdth)
#     if ii == 0:
#         ax2plot_list.append(plot)
#         ax2label_list.append('Mazevet 2018')

# #add isochor annotations
# for i in range(len(posy_list)):
#     dens = dens_list[i]
#     dens_str = str(dens)
#     if not posx_list[i] == None:
#         #print (dens_str, 'at:', posx_list[i], posy_list[i])
#         ax2.text(posx_list[i], posy_list[i], s = dens_str, fontdict = font3)


# #------------------------------------------------------------------------------
# #French 2015

# dens_list = [1420, 1430, 1440, 1450, 1500, 1600, 1800, 2000, 2500, 3000, 4000]
# posx_list = [150., 100, 50, 30, 20, 20, 20, 20, 300, 300, 300]
# posy_list = [3.0e8, 3.2e8, 4.0e8, 5.0e8, 1.5e9, 3.0e9, 6.0e9, 1.0e10, 3.5e10,
#              7.0e10, 2.1e11]

# temp_list = np.logspace(1, np.log10(3000.), N)
# pres_list = []

# for dens in dens_list:
#     ii = dens_list.index(dens)
#     pres_list.append([])
#     for temp in temp_list:
#         pres = French2015EOS.Pressure(d = dens, T = temp)
#         state = phaseTrans.phaseCheck(pres, temp)

#         if not state == 'high pressure French':
#             pres = None

#         pres_list[ii].append(pres)

#     plot, = ax2.plot(temp_list, pres_list[ii], linestyle = '-.',
#                      color = plot_color, linewidth = lwdth)
#     if ii == 0:
#         ax2plot_list.append(plot)
#         ax2label_list.append('French 2015')

# #add isochor annotations
# for i in range(len(posy_list)):
#     dens = dens_list[i]
#     dens_str = str(dens)
#     if not posx_list[i] == None:
#         #print (dens_str, 'at:', posx_list[i], posy_list[i])
#         ax2.text(posx_list[i], posy_list[i], s = dens_str, fontdict = font3)

# #------------------------------------------------------------------------------
# #VDW eos for comparison

# dens_list = [0.001, .01, .1, 1, 10, 100]
# temp_list = np.logspace(2, 5, 20)
# pres_list = []
# for dens in dens_list:
#     ii = dens_list.index(dens)
#     pres_list.append([])
#     for temp in temp_list:
#         pres = anfct.P_VDW(d = dens, T = temp, ll=5)
#         state = phaseTrans.phaseCheck(pres, temp)
#         if state == 'ice Ih':
#             pres = None

#         pres_list[ii].append(pres)

#     plot, = ax2.plot(temp_list, pres_list[ii], linestyle = 'none',
#                      color = plot_color, linewidth = lwdth, marker='s',
#                      markersize=3, markevery=1, markerfacecolor='None')
#     if ii == 0:
#         ax2plot_list.append(plot)
#         ax2label_list.append('van der Waals')

# #------------------------------------------------------------------------------
# #ideal gas for comparison

# temp_list = np.logspace(2, 5, 21)
# pres_list = []
# for dens in dens_list:
#     ii = dens_list.index(dens)
#     pres_list.append([])
#     for temp in temp_list:
#         pres = anfct.P_IGL(d = dens, T = temp, ll=5)
#         state = phaseTrans.phaseCheck(pres, temp)
#         if state == 'ice Ih':
#             pres = None

#         pres_list[ii].append(pres)

#     plot, = ax2.plot(temp_list, pres_list[ii], linestyle = 'none',
#                      color = plot_color, linewidth = lwdth, marker='v',
#                      markersize=3, markevery=1)
#     if ii == 0:
#         ax2plot_list.append(plot)
#         ax2label_list.append('ideal gas law')

# #------------------------------------------------------------------------------
# #phase boundaries

# tosolid_list = phaseTrans.solidLine(log = True, N = 100)
# tovapor_list = phaseTrans.vaporLine(log = True, N = 100)

# #reformat Schwager 2008 data for plotting
# pres_schwager = np.empty([len(data_schwager)])
# temp_schwager = np.empty([len(data_schwager)])


# #Compute solidus line over entire pressure range
# pres_solidus = np.logspace(1, np.log10(5.0e11), 200)
# temp_solidus = [phaseTrans.solidusTemp(p) for p in pres_solidus]


# for i in range(len(data_schwager)):
#     pres_schwager[i]=data_schwager[i][0]
#     temp_schwager[i]=data_schwager[i][1]

# #inter = phaseTrans.interpolateWagnerSchwager()
# #T_inter = np.linspace(600, 750, 10)
# #ax2.plot(T_inter, inter(T_inter)*1.0e9, color = plot_color)
# #ax2.plot(temp_schwager, pres_schwager*1.0e9, color = plot_color)

# for s in range(len(satcurve_list)):
#     sc, = ax2.plot(satcurve_list[s][0], satcurve_list[s][1], color = 'k', \
#                 marker = 'x', markersize=8, zorder = 4, linestyle = 'None')
#     if s == 0:
#         ax2plot_list.append(sc)
#         ax2label_list.append('saturation points')

# ax2.text(50., 1.0e8, s = 'solid (ice Ih)', fontdict = font1)
# ax2.text(50., 1.0e10, s = 'solid (ices VII & X)', fontdict = font1)
# ax2.text(300., 1.0e8, s = 'liquid', fontdict = font1)
# ax2.text(5000., 1.0e8, s = 'supercritical', fontdict = font1)
# ax2.text(5000., 1.0e6, s = 'gas', fontdict = font1)
# ax2.text(300., 3.0e2, s = 'vapor', fontdict = font1)

# ax2.loglog(temp_solidus, pres_solidus, color = phase_color)


# ax2.plot([0., 251.], [2.08e8, 2.086e8], color=phase_color)

# ax2.plot([T_critical, T_critical], [1.0, P_critical],
#         [0., T_critical], [P_critical, P_critical],
#         color = (.8, .8, .8), linestyle = '--', zorder = 1)

# ax2.plot([T_critical, T_critical], [P_critical, 1.2e10],
#         [T_critical, 1.0e5], [P_critical, P_critical],
#         color = (1, .6, .6), linestyle = '--', zorder = 1)

# legend1a = ax1.legend(ax1plot_list[0], ax1label_list[0], loc = 4,
#                       frameon = False)

# ax1.add_artist(legend1a)

# legend1b = ax1.legend(ax1plot_list[1], ax1label_list[1], loc = 2, fontsize=10,\
#                      fancybox=True, shadow=True)

# ax1.add_artist(legend1b)

# legend2 = ax2.legend(ax2plot_list, ax2label_list, loc = 2, fontsize=10,\
#                      fancybox=True, shadow=True)

# for legobj in legend2.legendHandles:
#     legobj.set_linewidth(1.5)

# ax2.add_artist(legend2)

# #ax2.loglog(tosolid_list[0], tosolid_list[1], color = phase_color,
#  #          linestyle = '-',
#   #       zorder = 3)

# ax2.plot(tovapor_list[0], tovapor_list[1], color = phase_color,
#          linestyle = '-',
#          zorder = 3)

# ax2.scatter(phaseTrans.T_critical, phaseTrans.P_critical, color = phase_color,
#             zorder = 3)

# ax2.scatter(phaseTrans.t_triple_list, phaseTrans.p_triple_list,
#             color = phase_color)

# ax1.set_title(r'$\rm H_2 O\ isotherms$', fontdict = font2)
# ax1.set_xlabel(r'$\rm Density \ [kg \ m^{-3}]$', fontdict = font2)
# ax1.set_ylabel(r'$\rm Pressure \ [Pa]$', fontdict = font2)
# ax1.set_ylim(1.0e6, 2.0e12)
# ax1.set_xlim(10., 5000.)

# ax2.set_title(r'$\rm H_2 O \ isochors \ [kg \ m^{-3}]$', fontdict = font2)
# ax2.set_xlabel(r'$\rmTemperature \ [K]$', fontdict = font2)
# ax2.set_ylabel(r'$\rm Pressure \ [Pa]$', fontdict = font2)
# ax2.yaxis.set_label_position('right')
# ax2.yaxis.tick_right()
# ax2.set_xlim(10, 1.0e5)
# ax2.set_ylim(1.0e2, 1.0e14)

# ax1.tick_params(which = 'both', axis = 'both', direction = 'in', top = 'on',
# 			right = 'on')
# ax1.tick_params(which = 'major', axis = 'y', length = 10, labelsize = 14, pad=10)
# ax1.tick_params(which = 'minor', axis = 'y', length = 5, labelsize = 14)
# ax1.tick_params(which = 'major', axis = 'x', length = 5, labelsize = 14, pad=10)

# ax2.tick_params(which = 'both', axis = 'both', direction = 'in', top = 'on',
# 			left = 'on')
# ax2.tick_params(which = 'major', axis = 'y', length = 10, labelsize = 14, pad=10)
# ax2.tick_params(which = 'minor', axis = 'y', length = 5, labelsize = 14)
# ax2.tick_params(which = 'major', axis = 'x', length = 5, labelsize = 14, pad=10)

# ax1.set_axisbelow(True)
# ax2.set_axisbelow(True)
# ax1.grid(which = 'both', color = grid_color)
# ax2.grid(which = 'both', color = grid_color)
# ax1.grid(which = 'major', linewidth=1, color = (.85, .85, .85))
# ax2.grid(which = 'major', linewidth =1, color = (.85, .85, .85))

# plt.show()
# fig.savefig('/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/test.png')
