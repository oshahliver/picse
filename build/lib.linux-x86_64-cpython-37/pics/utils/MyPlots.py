#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:06:11 2019

@author: oshah
"""

from matplotlib import pyplot as plt
import numpy as np
from PIMPrunparams import color_list, grid_color, gridalpha
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes


class plot:
    def __init__(self, dim=[], data=None):
        if len(dim) == 0:
            self.a, self.b = None, None

        else:
            self.a, self.b = dim

        self.data = data

    def display(self):
        # data must be in the shape of ixjxp where i is number of plot rows,
        # j number of plot columns and p number of curves for each ij pair
        data = self.data
        if not self.a == None and not self.b == None:
            self.fig, self.ax = plt.subplots(self.a, self.b)

        else:
            self.fig, self.ax = plt.subplots()

        # iterate over all plot fields and gather corresponding data
        self.plots = []

        for i in range(self.b):
            self.plots.append([])

            if self.a > 1:
                print("weird")
                for j in range(self.a):
                    dat = data[i][j]
                    self.plots[i].append([])
                    for p in range(len(dat)):
                        (pl,) = self.ax[i][j].plot(dat[p][0], dat[p][1])
                        self.plots[i][j].append(pl)

            else:
                print("not weird")
                dat = data[i]

                for p in range(6):

                    (pl,) = self.ax[i].plot(
                        dat[0][p], dat[1][p], marker="s", linestyle="-"
                    )
                    self.plots[i].append(pl)

    def adjust(self, labels=None):
        fnts = 12
        ticks = 8

        # first adjust curve colors
        for i in range(len(self.plots)):
            for j in range(len(self.plots[i])):
                self.plots[i][j].set_color(color_list[j * 2])
                self.plots[i][j].set_marker("s")
                self.plots[i][j].set_markersize(6)

        for i in range(len(self.ax)):
            if self.a > 1:
                for j in range(len(self.ax[i])):
                    axx = self.ax[i][j]
                    axx.set_xlabel(r"$x$")
                    axx.set_ylabel(r"$y$")

                    if j == 0:
                        axx.yaxis.set_label_position("left")

                    elif j == self.a - 1:
                        axx.yaxis.set_label_position("right")

            else:
                axx = self.ax[i]
                axx.set_xlabel(r"$M/M_{\oplus}$", fontsize=fnts)
                # general adjustments
                axx.tick_params(
                    which="both",
                    top=True,
                    right=True,
                    left=True,
                    direction="in",
                    size=ticks,
                )

                axx.set_xscale("log")
                axx.grid(color=(0.95, 0.95, 0.95), which="both", linewidth=1)
                axx.grid(color=(0.95, 0.95, 0.95), which="major", linewidth=2)

                if i == 0:
                    axx.yaxis.set_label_position("left")

                elif i == self.b - 1:
                    axx.yaxis.set_label_position("right")
                    axx.tick_params(right="on", left="off")
                    axx.yaxis.tick_right()

            self.ax[0].set_ylim(0.25, 2.2)
            self.ax[0].set_xlim(0.1, 10)
            self.plots[0][-1].set_marker("")
            self.plots[0][-2].set_marker("")
            self.plots[0][-1].set_color(color_list[7])
            self.plots[0][-2].set_color(color_list[10])
            self.ax[0].set_ylabel(r"$R/R_{\oplus}$", fontsize=fnts)
            self.ax[1].set_ylabel(r"$M_{Core}/M$", fontsize=fnts)

            self.ax[0].text(
                0.15,
                2.05,
                r"$\rm Mg$" + "#" + r"$=0.76 \ ($" + "~" + r"$\rmsolar)$",
                fontsize=14,
            )

            self.ax[0].text(0.15, 1.92, r"$\rm T = 300 \ K$", fontsize=14)

            self.ax[0].text(
                0.6,
                0.4,
                r"$\rm (R_{hyd} - R_{dry})/R_{dry}$" + " ~" + r"$5-10$" + "%",
                color=(0.4, 0.4, 0.4),
                fontsize=14,
            )

            try:
                # add plot labels
                for p in range(len(self.plots[0])):
                    col = color_list[p * 2]
                    if p == 4:
                        col = color_list[7]
                    elif p == 5:
                        col = color_list[10]

                    self.ax[0].text(0.15, 1.76 - 0.1 * p, labels[p], color=col)

            except TypeError:
                pass

            errorsx = [0.369, 0.43]
            errorsy = [0.108, 0.093]
            names = ["HD219134 b", "HD 219134 c"]

            self.ax[0].errorbar(
                [4.736, 4.357],
                [1.602, 1.511],
                zorder=10,
                xerr=errorsx,
                yerr=errorsy,
                color="k",
                marker="o",
                linestyle="None",
            )

            self.ax[0].text(4.7, 1.76, names[0], fontsize=8)
            self.ax[0].text(4.3, 1.35, names[1], fontsize=8)
