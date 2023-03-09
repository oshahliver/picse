#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 14:38:11 2020

@author: os18o068
"""

import numpy as np
from matplotlib import pyplot as plt
from pics.runparams import color_list, pure_curves_color_list, material_plot_list_fort
from pics.physicalparams import m_earth, r_earth, mH, mFe
import matplotlib
import astropy.table
from astropy.io import ascii
import pandas as pd

matplotlib.rc("text", usetex=True)
matplotlib.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath}"]


suff = "double_check"

file_dir = "/home/os18o068/Documents/PHD/Abbildungen/Figures_paper_1/"
filename = "test"
format = "pdf"

delta_R_labels = ["Type-1/2", "Type-2/3"]
delta_R_label_pos = [[2.0, 0.065], [2.0, 850.0]]
pures = [1, 2, 4]

pure_curves = []

pure_names = ["pure_curve_" + str(p) + ".atab" for p in pures]

pure_positions = [[0.2, 1.1], [1.0, 0.65], [0.4, 0.9]]
pure_rotations = [30, 30, 30]

# Load pre-generated curves for pure substances
for p in pure_names:
    pure_data = ascii.read(p)

    df = pure_data.to_pandas()
    pure_curves.append(df)

data_hyd = np.load("planet_output_hyd_" + suff + ".npy")
data_dry = np.load("planet_output_dry_" + suff + ".npy")
data_oce = np.load("planet_output_oce_" + suff + ".npy")

type_names = ["Type-1", "Type-2", "Type-3"]

acc_require = np.array([1.0e-3, 5.0e-3, 1.0e-2, 3.0e-2, 5.0e-2, 1.0e-2])

data = [data_dry, data_hyd, data_oce]

temp_label_color = "grey"

legend_alpha = 0.3

mark = "s"
fnts0 = 12
fnts1 = 9
fnts2 = 7
fnts3 = 5
lwdth = 0.75

yaxes_labels = [
    [r"$\rm tot. \ M_{\rm H_2O}/M$", r"$M_{\rm H_2O}/M$", r"$\delta R_1/R_1$"],
    [r"$\rm Ocean \ depth \ [km]$", r"$M_{\rm H_2O}/M$", r"$\delta R_2/R_2$"],
]
temp_labels = [
    [
        r"${\Delta T_{\rm TBL} = 200 \ \rm K}$",
        r"${\Delta T_{\rm TBL} = 700 \ \rm K}$",
        r"${\Delta T_{\rm TBL} = 1200 \ \rm K}$",
        r"${\Delta T_{\rm TBL} = 1700 \ \rm K}$",
    ],
    [
        r"${\Delta T_{\rm TBL} = 200 \ \rm K}$",
        r"${\Delta T_{\rm TBL} = 700 \ \rm K}$",
        r"${\Delta T_{\rm TBL} = 1200 \ \rm K}$",
        r"${\Delta T_{\rm TBL} = 1700 \ \rm K}$",
    ],
]

pres_labels = [r"$P_{\rm TBL} = $", r"$P_{S} = $"]

axes_lims = [
    [[0.0, 0.065], [0.0, 0.065], [0.0, 0.03]],
    [[0.0, 850.0], [0.0, 0.065], [-0.03, 0.06]],
]

temp_label_ypos = [[5.0, sum(axes_lims[0][2]) / 2.0], [5.0, sum(axes_lims[1][2]) / 2.0]]

N_temps = len(data[0][0][0])
N_Mg = len(data[0][0][0][0])
N_Si = len(data_hyd[0])

mean_T_surface = np.empty([3, N_Si, N_temps])
mean_Mg_number = np.empty([3, N_Si, N_temps, N_Mg])
mean_P_surface = np.empty([3, N_Si, N_temps])

# Check if planets have converged properly
for l in range(3):
    # loop over Si#
    for s in range(len(data_hyd[0])):
        # loop over temp
        for i in range(len(data[l][0][0])):
            all_T_surface = []
            all_P_surface = []
            all_Mg_number = []
            # loop over Mg#
            for j in range(len(data[l][0][0][0])):
                all_Mg_number.append([])

                # loop over planets
                for p in range(len(data[l][0][s][i][j])):
                    pl = [data[l][d][s][i][j][p] for d in range(len(data[l]))]
                    reldevs = np.zeros([6])

                    # total mass
                    reldevs[0] = abs(pl[0] - pl[17]) / pl[17]

                    # Mg#
                    reldevs[1] = abs(pl[2] - pl[16]) / pl[16]

                    # T_surface
                    reldevs[2] = abs(pl[5] - pl[18]) / pl[18]

                    # xi_H_core
                    if pl[5] < 1000.0:
                        reldevs[3] = abs(pl[14] - pl[15])

                    else:
                        reldevs[3] = abs(pl[14] - pl[15]) / pl[15]

                    # ocean fraction
                    reldevs[4] = abs((10 ** pl[20] - 10 ** pl[19])) / 10 ** pl[19]

                    match = True

                    if l == 1:
                        w = 5
                        if reldevs[w] > acc_require[w]:
                            print("not converged")
                            print(l, s, i, j, p)
                            print("reldev =", reldevs[w])
                            print(
                                "mass should =",
                                pl[17],
                                "radius =",
                                pl[1],
                                "Mg# should =",
                                pl[16],
                                "T_surface_is =",
                                pl[5],
                            )
                            print(
                                "T should =",
                                data[l][18][s][i][j][p],
                                "OF should =",
                                data[l][19][s][i][j][p],
                                "OF is =",
                                data[l][20][s][i][j][p],
                            )

                    for r in range(len(reldevs)):
                        if reldevs[r] > acc_require[r]:
                            match = False
                            for d in range(len(pl)):
                                data[l][d][s][i][j][p] = None

                    if match == False:
                        print("not converged")
                        print(l, s, i, j, p)
                        print("reldev =", reldevs[w])
                        print(
                            "mass should =",
                            pl[17],
                            "radius =",
                            pl[1],
                            "Mg# should =",
                            pl[16],
                            "T_surface_is =",
                            pl[5],
                        )
                        print(
                            "T should =",
                            data[l][18][s][i][j][p],
                            "OF should =",
                            data[l][19][s][i][j][p],
                            "OF is =",
                            data[l][20][s][i][j][p],
                        )

                    if match:
                        all_T_surface.append(data[l][5][s][i][j][p])
                        all_Mg_number[j].append(data[l][2][s][i][j][p])
                        all_P_surface.append(data[l][6][s][i][j][p])

                mean_Mg_number[l][s][i][j] = np.mean(np.asarray(all_Mg_number[j]))
            mean_T_surface[l][s][i] = np.mean(np.asarray(all_T_surface[i]))
            mean_P_surface[l][s][i] = np.mean(np.asarray(all_P_surface[i]))


# Compute mean values and standard deviations for BC
mean_T1 = np.mean((mean_T_surface[0][0][0] + mean_T_surface[1][0][0]) / 2.0)
mean_T2 = np.mean(mean_T_surface[1][0][0])


# Plot T_C, P_C, P_H2, xi_H and M-R
# loop over Si#
for s in range(len(data_hyd[0])):
    fig, ax = plt.subplots(max(len(data[1][0][0]), 2), 2, sharex=True)
    fig.subplots_adjust(hspace=0.2, wspace=0.2)

    figs = []
    axes = []

    ffigs = []
    aaxes = []

    # loop over different planet types
    for k in range(len(data)):
        f, a = plt.subplots(max(len(data[1][0][0]), 2), 3, sharex=True)
        f.subplots_adjust(hspace=0.2, wspace=0.3)

        figs.append(f)
        axes.append(a)

        a[-1][0].set_xlabel(r"$M/M_{\oplus}$", fontsize=fnts2)
        a[-1][1].set_xlabel(r"$M/M_{\oplus}$", fontsize=fnts2)
        a[-1][2].set_xlabel(r"$M/M_{\oplus}$", fontsize=fnts2)

        for kk in range(len(a)):
            a[kk][0].set_ylabel(r"$T_{\rm C} \ \rm [10^3 \ K]$", fontsize=fnts2)
            a[kk][1].set_ylabel(r"$P_{\rm C} \ \rm [GPa]$", fontsize=fnts2)
            a[kk][2].set_ylabel(r"$R/R_{\oplus}$", fontsize=fnts2)

        # for ocean planets plot P_Ocean and delta_T ocean
        if k == 2:
            ff, aa = plt.subplots(max(len(data[1][0][0]), 2), 2, sharex=True)
            ff.subplots_adjust(hspace=0.2, wspace=0.3)

            ffigs.append(ff)
            aaxes.append(aa)

            aa[-1][0].set_xlabel(r"$M/M_\oplus$", fontsize=fnts2)
            aa[-1][1].set_xlabel(r"$M/M_\oplus$", fontsize=fnts2)

            bbox_props = dict(edgecolor="k", facecolor="white", linewidth=0.5)

            aa[0][0].text(1.5, 100.0, "Ocean Properties", ha="center", bbox=bbox_props)

            for kk in range(len(aa)):
                aa[kk][0].set_ylabel(r"$P_{\rm Ocean} \ \rm [GPa]$", fontsize=fnts2)
                aa[kk][1].set_ylabel(
                    r"$\Delta T_{\rm Ocean} \ \rm [K]$", fontsize=fnts2
                )

                aa[kk][0].set_ylim(0.005, 100.0)
                aa[kk][1].set_ylim(0, 500)

                aa[kk][0].set_xlim(0.1, 3.0)
                aa[kk][1].set_xlim(0.1, 3.0)

                for kkk in range(len(aa[kk])):
                    aa[kk][kkk].tick_params(
                        direction="in",
                        top=True,
                        right=True,
                        which="both",
                        labelsize=fnts3,
                    )

    ax[-1][0].set_xlabel(r"$M/M_{\oplus}$", fontsize=fnts2)
    ax[-1][1].set_xlabel(r"$M/M_{\oplus}$", fontsize=fnts2)

    bbox_props = dict(edgecolor="k", facecolor="white", linewidth=0.5)

    ax[0][0].text(1.0, 1040, "Hydrogen in the Core", ha="center", bbox=bbox_props)

    # loop over temp
    for i in range(len(data[1][0][0])):
        plot_list = []
        legend_list = []

        for k in range(len(ax[i])):
            ax[i][k].tick_params(
                direction="in", top=True, right=True, which="both", labelsize=fnts3
            )

        for k in range(len(axes)):
            axes[k][i][-1].text(
                5.0,
                0.9,
                temp_labels[0][i],
                rotation=90,
                fontsize=fnts2,
                va="center",
                color=temp_label_color,
                weight="bold",
            )

            axes[k][i][-1].set_ylim(0.4, 1.4)

        ax[i][0].set_ylabel(r"$P_{\rm H_2, CMB} \ \rm [MPa]$", fontsize=fnts2)
        ax[i][1].set_ylabel(r"$X_{\rm H, Core} \ \rm [wt \%]$", fontsize=fnts2)
        ax[i][0].set_xlim(0.1, 3.0)
        ax[i][1].set_xlim(0.1, 3.0)

        ax[i][-1].text(
            5.0,
            0.6,
            temp_labels[0][i],
            rotation=90,
            fontsize=fnts2,
            color=temp_label_color,
            va="center",
            weight="bold",
        )

        # loop over Mg#
        for j in range(len(data[1][0][0][0])):
            cc = color_list[j]

            # plot P_H2
            (pl,) = ax[i][0].semilogx(
                data[1][0][s][i][j],
                data[1][21][s][i][j] * 1.0e-6,
                color=cc,
                linewidth=lwdth,
            )

            xi = data[1][14][s][i][j]
            X = mH * xi / (mFe * (1.0 - xi) + xi * mH)

            # plot X_H
            ax[i][1].semilogx(data[1][0][s][i][j], X * 100.0, color=cc, linewidth=lwdth)

            ax[i][0].set_ylim(0.0, 1200.0)
            ax[i][1].set_ylim(0.0, 1.5)

            # Loop over planet types
            for k in range(len(data)):

                # Plot T_C
                axes[k][i][0].semilogx(
                    data[k][0][s][i][j],
                    data[k][11][s][i][j] * 1.0e-3,
                    color=cc,
                    linewidth=lwdth,
                )

                # Plot P_C
                axes[k][i][1].loglog(
                    data[k][0][s][i][j],
                    data[k][12][s][i][j] * 1.0e-9,
                    color=cc,
                    linewidth=lwdth,
                )

                # Plot M-R
                axes[k][i][2].semilogx(
                    data[k][0][s][i][j], data[k][1][s][i][j], color=cc, linewidth=lwdth
                )

                for p in range(len(pures)):
                    axes[k][i][2].semilogx(
                        pure_curves[p]["col0"] / m_earth,
                        pure_curves[p]["col1"] / r_earth,
                        color=pure_curves_color_list[pures[p] - 1],
                        linewidth=1.5,
                    )

                    axes[k][i][2].text(
                        pure_positions[p][0],
                        pure_positions[p][1],
                        material_plot_list_fort[pures[p] - 1],
                        color=pure_curves_color_list[pures[p] - 1],
                        fontsize=fnts2,
                        rotation=pure_rotations[p],
                    )

                if k == 2:
                    # P_Ocean
                    aaxes[0][i][0].loglog(
                        data[k][0][s][i][j],
                        data[k][22][s][i][j] * 1.0e-9,
                        color=cc,
                        linewidth=lwdth,
                    )

                    # Delta_T_Ocean
                    aaxes[0][i][1].semilogx(
                        data[k][0][s][i][j],
                        data[k][23][s][i][j] - data[0][5][s][i][j],
                        color=cc,
                        linewidth=lwdth,
                    )

                    aaxes[0][i][1].text(
                        5.0,
                        250.0,
                        temp_labels[0][i],
                        rotation=90,
                        fontsize=fnts2,
                        va="center",
                        color=temp_label_color,
                        weight="bold",
                    )

                axes[k][i][0].set_ylim(0.0, 20.0)
                axes[k][i][1].set_ylim(10, 2000)
                axes[k][i][0].set_xlim(0.1, 3.0)

                for l in range(len(axes[k][i])):
                    axes[k][i][l].tick_params(
                        direction="in",
                        top=True,
                        right=True,
                        which="both",
                        labelsize=fnts3,
                    )

            plot_list.append(pl)
            legend_list.append(
                r"$\rm Mg\# = \ $" + str(round(mean_Mg_number[0][s][i][j], 1))
            )

            if i == 0:
                legend = ax[i][0].legend(
                    plot_list,
                    legend_list,
                    loc=2,
                    fontsize=fnts3,
                    framealpha=legend_alpha,
                )

                ax[i][0].add_artist(legend)

                for k in range(3):
                    legend1 = axes[k][0][0].legend(
                        plot_list,
                        legend_list,
                        loc=2,
                        fontsize=fnts3,
                        framealpha=legend_alpha,
                    )

                    axes[k][0][0].add_artist(legend1)

                    if k == 2:
                        legend2 = aaxes[0][0][0].legend(
                            plot_list,
                            legend_list,
                            loc=2,
                            fontsize=fnts3,
                            framealpha=legend_alpha,
                        )

                        aaxes[0][0][0].add_artist(legend2)

    fig.align_labels()
    fig.savefig(file_dir + "core_content" + "." + format, bbox_inches="tight")

    for k in range(len(data)):
        bbox_props = dict(edgecolor="k", facecolor="white", linewidth=0.5)

        axes[k][0][0].text(
            3.5,
            20.0,
            type_names[k],
            fontsize=fnts0,
            ha="right",
            va="bottom",
            bbox=bbox_props,
        )

        figs[k].align_labels()
        figs[k].savefig(
            file_dir + "internals_type_" + str(k + 1) + "." + format,
            bbox_inches="tight",
        )

        if k == 2:
            ffigs[0].align_labels()
            ffigs[0].savefig(file_dir + "ocean_props." + format, bbox_inches="tight")


# Plot core and mantle reservoirs and delta R
for l in range(2):
    # loop over Si#
    for s in range(len(data_hyd[0])):
        plot_list1 = []
        plot_list2 = []

        fig, ax = plt.subplots(max(len(data[l + 1][0][0]), 2), 3, sharex=True)
        fig.subplots_adjust(hspace=0.2, wspace=0.4)

        bbox_props = dict(edgecolor="k", facecolor="white", linewidth=0.5)

        ax[0][0].text(
            delta_R_label_pos[l][0],
            delta_R_label_pos[l][1],
            delta_R_labels[l],
            ha="center",
            bbox=bbox_props,
        )

        # loop over temp
        for i in range(len(data[l + 1][0][0])):
            legend_list1 = []
            legend_list2 = []

            # loop over Mg#
            for j in range(len(data[l + 1][0][0][0])):
                cc = color_list[j]

                # Compute mean temperature and standarddeviation
                T_mean = np.nanmean(data[1][5][s][i])
                Mg_mean = np.nanmean(data[1][2][s][i][j])
                Si_mean = np.nanmean(data[1][3][s][i][j])
                P_mean = np.nanmean(data[1][6][s][i][j])

                markers = []
                for p in range(len(data[l + 1][0][s][i][j])):
                    Mg = data[l + 1][2][s][i][j][p]
                    T = data[l + 1][5][s][i][j][p]
                    M = data[l + 1][0][s][i][j][p]

                if l == 0:
                    (pl,) = ax[i][0].semilogx(
                        data[l + 1][0][s][i][j],
                        data[l + 1][8][s][i][j],
                        color=cc,
                        linewidth=lwdth,
                    )

                    (pl1,) = ax[i][1].semilogx(
                        data[l + 1][0][s][i][j],
                        data[l + 1][10][s][i][j],
                        color=cc,
                        linestyle="--",
                        linewidth=lwdth,
                    )

                    (pl2,) = ax[i][1].semilogx(
                        data[l + 1][0][s][i][j],
                        data[l + 1][8][s][i][j] - data[l + 1][10][s][i][j],
                        color=cc,
                        linewidth=lwdth,
                    )

                elif l == 1:
                    (pl,) = ax[i][0].semilogx(
                        data[l + 1][0][s][i][j],
                        data[l + 1][13][s][i][j] * r_earth / 1000.0,
                        color=cc,
                        linewidth=lwdth,
                    )

                    (pl1,) = ax[i][1].semilogx(
                        data[l + 1][0][s][i][j],
                        data[l + 1][8][s][i][j],
                        color=cc,
                        linewidth=lwdth,
                    )

                if j == 0:
                    plot_list2 = [pl1, pl2]
                    legend_list2 = ["Core reservoir", "Mantle reservoir"]

                ax[i][2].semilogx(
                    data[l + 1][0][s][i][j],
                    (data[l + 1][1][s][i][j] - data[l][1][s][i][j])
                    / data[l][1][s][i][j],
                    color=cc,
                    linewidth=lwdth,
                )

                ax[i][2].plot(
                    [0.0, 10.0], [0.0, 0.0], color="k", linestyle="--", linewidth=0.5
                )

                legend_list1.append(r"$\rm Mg \# = \ $" + str(round(Mg_mean, 2)))
                plot_list1.append(pl)

            # loop over columns
            for k in range(3):
                ax[i][k].set_ylim(axes_lims[l][k])
                ax[i][k].set_ylabel(yaxes_labels[l][k], fontsize=fnts2)
                ax[i][k].tick_params(
                    direction="in", top=True, right=True, which="both", labelsize=fnts3
                )
                ax[i][k].set_xlim(0.1, 3.0)

                ax[-1][k].set_xlabel(r"$M/M_{\oplus}$", fontsize=fnts2)

            if i == 0:
                legend = ax[i][0].legend(
                    plot_list1,
                    legend_list1,
                    loc=2,
                    fontsize=fnts3,
                    framealpha=legend_alpha,
                )
                ax[i][0].add_artist(legend)

            if l == 0 and i == 0:
                legend = ax[i][1].legend(
                    plot_list2,
                    legend_list2,
                    loc=2,
                    fontsize=fnts3,
                    framealpha=legend_alpha,
                )

                ax[i][1].add_artist(legend)

            ax[i][-1].text(
                temp_label_ypos[l][0],
                temp_label_ypos[l][1],
                temp_labels[l][i],
                rotation=90,
                fontsize=fnts2,
                va="center",
                color=temp_label_color,
                weight="bold",
            )

        fig.align_labels()
        fig.savefig(
            file_dir + "delta_R" + str(l + 1) + "." + format, bbox_inches="tight"
        )

plt.close("all")
