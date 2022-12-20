#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 18:05:17 2019

@author: oshah
"""
import numpy as np
from matplotlib import pyplot as plt
import copy
from PIMPrunparams import color_list, grid_color
from matplotlib.ticker import AutoMinorLocator

# define tripple points of Mg(OH)2

T_triple1 = 794.5
T_triple2 = 1087.0

P_triple1 = 29.67e9
P_triple2 = 33.65e9


coeffs = [
    np.array(
        [-7.19993352e-02, 8.20956584e00, -3.46793056e02, 6.48875298e03, -4.50699821e04]
    ),
    np.array(
        [-1.01962812e01, 1.26707689e03, -5.90347936e04, 1.22208342e06, -9.48250541e06]
    ),
    np.array([-14.93208847, 1589.76095947]),
    np.array(
        [-8.71692725e-04, 1.19579527e-01, -5.57672252e00, 8.55574995e01, 1.08400218e03]
    ),
    np.array([-31.11306139, 252.46893265, 863.14809939]),
]

data = np.array(
    [
        [
            [19.38, 20.12, 22.27, 24.94, 27.1, 29.67, 31.01, 32.66, 33.65],
            [0.0, 214.8, 383.7, 535.9, 671.6, 794.5, 891.7, 1013.0, 1087.0],
            "Hermann 2016 (ab initio)",
        ],
        [
            [29.67, 31.28, 32.66, 33.45, 33.45, 30.41],
            [794.5, 618.8, 418.3, 196.7, 0.0, 725.2],
            "Hermann 2016",
        ],
        [[27.23, 33.65, 39.97], [1198.0, 1087.0, 1008.0], "Hermann 2016"],
        [
            [
                1.0,
                1.5,
                2.0,
                3.0,
                4.0,
                4.0,
                6.0,
                6.0,
                8.0,
                8.0,
                8.0,
                10.0,
                11.5,
                11.5,
                13.0,
                15.0,
                15.0,
                13.0,
                11.5,
                4.0,
            ],
            [
                1056,
                1204.0,
                1258.0,
                1349.0,
                1374.0,
                1393.0,
                1414.0,
                1441.0,
                1457.0,
                1464.0,
                1442.0,
                1459.0,
                1443.0,
                1498.0,
                1521.0,
                1550.0,
                1521.0,
                1497.0,
                1392.0,
                1293.0,
            ],
            "Johnson & Walker 1993 (measured)",
        ],
        [[1.67, 2.7, 3.29], [1198.0, 1318.0, 1357.0], "Irving et al. 1977 (measured)"],
    ]
)


def fit_phases(data=None, order=[4, 4, 1, 4, 2, 2]):
    coeffs = []
    for d in range(len(data)):

        dat = copy.deepcopy(data[d])
        o = order[d]

        # for the upper most curve use experimental data in the low pressure
        # regime and the theoretical studies by Hermann 2016 for intermediate
        # pressures

        if d == 3:
            dat = [[], []]
            for i in range(len(data[2][0])):
                dat[0].append(data[2][0][i])
                dat[1].append(data[2][1][i])

            for i in range(len(data[3][0])):
                dat[0].append(data[3][0][i])
                dat[1].append(data[3][1][i])

            for i in range(len(data[4][0])):
                dat[0].append(data[4][0][i])
                dat[1].append(data[4][1][i])

        # make sure that fits go through triple points
        if d == 3 or d == 0 or d == 2:
            for i in range(100):
                dat[0].append(P_triple2 * 1.0e-9)
                dat[1].append(T_triple2)

        if d == 1 or d == 0:
            for i in range(100):
                dat[0].append(P_triple1 * 1.0e-9)
                dat[1].append(T_triple1)

        z = np.polyfit(dat[0], dat[1], o)
        coeffs.append(z)

    return coeffs


def T_trans(P=None, i=0):
    T = [coeffs[i][c] * P ** (len(coeffs[i]) - c - 1) for c in range(len(coeffs[i]))]

    for c in range(len(coeffs[i])):
        print(P ** (len(coeffs[i]) - c - 1))

    # print ('T_Trans =', sum(T))
    return sum(T)


def phaseCheck(T=None, P=None):
    """Determines the phase region in the PT-plane of Brucite at given pressure
    in Pa and temperature in K. Returns integer value for the phase at hand.

    0: alpha (Br-P3)
    1: beta (Br-P41212)
    2: gamma1 & gamma2 (Per + H2O)
    """
    try:
        pres = P * 1.0e-9

    except TypeError:
        pass

    if P < P_triple2:

        if T_trans(P=pres, i=3) < T:
            return 2

        else:
            if T_trans(P=pres, i=0) < T:
                return 0

            else:
                if P < P_triple1:
                    return 1

                else:
                    if T_trans(P=pres, i=1) < T:
                        return 2

                    else:
                        return 1

    else:
        return 2


def plot(pmax=60.0, log=False):

    pres = np.array(
        [
            np.linspace(0.0, P_triple2 * 1.0e-9, 100),
            np.linspace(P_triple1 * 1.0e-9, pmax, 100),
            np.linspace(P_triple2 * 1.0e-9, pmax, 100),
            np.linspace(0.0, P_triple2 * 1.0e-9, 100),
            np.linspace(0.0, pmax, 100),
        ]
    )

    fig, ax = plt.subplots()

    ax.set_ylim(0, 2000)
    ax.set_xlim(0, pmax)

    ax.set_xlabel(r"$Pressure \ [GPa]$")
    ax.set_ylabel(r"$Temperature \ [K]$")

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="both", top="on", right="on", direction="in")
    ax.grid(color=grid_color, linewidth=2, which="major")
    ax.grid(color=grid_color, linewidth=1, which="minor")
    ax.text(10, 600, r"$Br-P \bar3$")
    ax.text(22, 200, r"$Br-P4_1 2_1 2$")
    ax.text(40, 600, r"$Per \ + \ H_2 O(s)$")
    ax.text(40, 1500, r"$Per \ + \ H_2 O(f)$")

    legend_list = []
    plot_list = []
    legend_list = []

    markers = ["o", "o", "o", "s", "v", "s"]

    for d in range(len(data)):

        dat = data[d]
        facecol = color_list[d]
        edgecol = "k"

        if d < 3:
            facecol = color_list[d]

        else:
            facecol = color_list[d]

        sc = ax.scatter(
            dat[0],
            dat[1],
            zorder=10,
            marker=markers[d],
            s=40,
            facecolor=facecol,
            edgecolor=edgecol,
        )

        if d < 4:
            ax.plot(
                pres[d],
                T_trans(P=pres[d], i=d),
                color="grey",
                linestyle="-",
                linewidth=2,
            )

        if d == 0 or d == 3 or d == 4:
            plot_list.append(sc)
            legend_list.append(data[d][2])

    legend = ax.legend(plot_list, legend_list, loc=2, frameon=True)

    ax.add_artist(legend)

    if log:
        ax.set_xscale("log")
