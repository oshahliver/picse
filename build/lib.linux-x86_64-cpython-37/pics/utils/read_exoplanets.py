#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 19:12:30 2021

@author: os18o068
"""


import pandas as pd
import numpy as np
from astropy.io import ascii
from matplotlib import pyplot as plt
import matplotlib as mpl
from adjustText import adjust_text
from tabulate import tabulate
import matplotlib.transforms as transforms
from pics.utils import functionTools as ftool

from pics.materials.phase_transitions_water_Wagner2002 import T_critical, P_critical
from pics.physicalparams import M_solar, R_solar, T_solar, names_solar, r_earth, m_earth

mpl.rc("text", usetex=True)
mpl.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath, amssymb}"]


def _filter(M, T, **kwargs):
    a, b, c, d = 0.1468, 0.2694, 0.05546, -0.03538
    return 10 ** (a + b * np.log10(M) + c * T / 1000 + d * np.log10(M) * T / 1000)


class Archive:
    def __init__(self, M_max=10.0, T_max=T_critical, dM_max=0.15):
        self.data = ascii.read(
            "/home/os18o068/Documents/PHD/Projects/Planets/Code/PS_2021.12.02_02.23.23.csv"
        )
        self.pureDir = "/home/os18o068/Documents/PHD/Projects/Planets/Data/pure_curves/"
        self.df = self.data.to_pandas()

        # self.df.sort_values(by = 'pl_eqt', na_position='arlast')
        self.planets = []
        self.names = []
        self.clean_names = []
        self.hostnames = []
        self.filtered_data = []
        self.distances = []
        self.T_max = T_max
        self.M_max = M_max
        self.dM_max = dM_max

        # self.clean()
        self.add_max_water_radius()
        self.readPureCurves()

    def readPureCurves(self):
        self.pureCurveNames = ["pure_iron", "pure_iron_rock_50_50_300_K", "pure_rock"]
        self.pureCurvesPlanets = []
        for i in range(len(self.pureCurveNames)):
            self.pureCurvesPlanets.append([])
            planets = ftool.load_objects(
                self.pureDir + self.pureCurveNames[i] + "_cali.pkl"
            )
            # Only consider planets which converged
            for j in range(len(planets)):
                planets[j].check_convergence()
                if planets[j].status == "constructed":
                    if not planets[j].converged:
                        pass

                    else:
                        self.pureCurvesPlanets[i].append(planets[j])

    def clean(self, **kwargs):
        """Drops all data points for which one of the needed parameters is NaN"""
        try:
            params = kwargs["params"]

        except KeyError:
            params = [
                "pl_bmasse",
                "pl_bmasseerr1",
                "pl_bmasseerr2",
                "pl_rade",
                "pl_radeerr1",
                "pl_radeerr2",
                "pl_eqt",
                "pl_eqterr1",
                "pl_eqterr2",
                "pl_name",
                "hostname",
                "sy_dist",
            ]
        print("params =", params)
        for param in params:
            self.df = self.df.dropna(subset=[param])

        # Filter out planets for which upper mass limit is too large
        self.df = self.df.drop(
            self.df[self.df["pl_bmasse"] + self.df["pl_bmasseerr1"] > self.M_max].index
        )

        # Filter out planets for which lower mass error is too large
        self.df = self.df.drop(
            self.df[
                abs(self.df["pl_bmasseerr1"] / self.df["pl_bmasse"]) > self.dM_max
            ].index
        )

        # Filter out planets for which upper mass error is too large
        self.df = self.df.drop(
            self.df[
                abs(self.df["pl_bmasseerr2"] / self.df["pl_bmasse"]) > self.dM_max
            ].index
        )

        # Filter out planets for which equilibrium temperature is too high
        self.df = self.df.drop(self.df[self.df["pl_eqt"] > self.T_max].index)

        # Add series with clean planet names for file labeling
        if "pl_cleanname" in self.df.columns:
            self.df = self.df.drop(columns=["pl_cleanname"])

        clean_names = self.df["pl_name"].copy(deep=True)
        for ind in self.df.index:
            clean_names[ind] = clean_names[ind].replace(" ", "_")
            clean_names[ind] = clean_names[ind].replace("-", "_")

        self.df.insert(1, "pl_cleanname", clean_names, True)

    def reset(self):
        self.df = self.data.to_pandas()

    def unfilter(self):
        self.df = self.data.to_pandas()
        self.clean()
        self.add_max_water_radius()

    def add_max_water_radius(self, fct=_filter, **kwargs):
        """Computes the maximum radius of a pure water sphere within the
        uncertainties of the equilibrium temperature and total mass
        """
        # Upper limit for equilibrium temperature
        # Compute the max radius for the lower mass limit

        for i in range(2):
            for j in range(2):
                str1 = "pl_bmasseerr" + str(i + 1)
                str2 = "pl_eqterr" + str(j + 1)
                str3 = "pl_maxrade" + str(i + 1) + str(j + 1)
                rad = [
                    fct(
                        self.df["pl_bmasse"][i] + self.df[str1][i],
                        self.df["pl_eqt"][i] + self.df[str2][i],
                        **kwargs,
                    )
                    for i in self.df.index
                ]

                if str3 in self.df.columns:
                    self.df[str3] = rad

                else:
                    self.df.insert(1, str3, rad, True)

    def prt(self):
        dat = []
        for i in range(len(self.filtered_names)):
            dat.append([])
            dat[i].append(self.filtered_names[i])
            for j in range(len(self.filtered_data[i])):
                dat[i].append(self.filtered_data[i][j])

            dat[i].append(
                int(self.filtered_data[i][1] / self.filtered_data[i][0] * 100)
            )
            dat[i].append(
                int(self.filtered_data[i][2] / self.filtered_data[i][0] * 100)
            )

        tab = tabulate(
            dat,
            headers=[
                "Name",
                "Mass",
                "+err",
                "-err",
                "Radius",
                "+err",
                "-err",
                "T_eq",
                "+err",
                "-err",
                "+%",
                "-%",
            ],
        )

        print()
        print(f"{tab}")
        print()

    def filter(self, fct=filter, max_mass_err=0.1, **kwargs):
        """Apply filter function to filter out unwanted data points"""
        self.clean()
        self.add_max_water_radius()

        # Remove stupid guy that can only be explained by pure iron
        self.df = self.df.drop(self.df[self.df["pl_name"] == "HD 219134 f"].index)

        # Remove planets for which a pure water sphere is not possible
        self.df = self.df.drop(
            self.df[
                (self.df["pl_rade"] + self.df["pl_radeerr1"] > self.df["pl_maxrade11"])
                & (
                    self.df["pl_rade"] + self.df["pl_radeerr1"]
                    > self.df["pl_maxrade12"]
                )
                & (
                    self.df["pl_rade"] + self.df["pl_radeerr1"]
                    > self.df["pl_maxrade21"]
                )
                & (
                    self.df["pl_rade"] + self.df["pl_radeerr1"]
                    > self.df["pl_maxrade22"]
                )
                & (
                    self.df["pl_rade"] + self.df["pl_radeerr2"]
                    > self.df["pl_maxrade11"]
                )
                & (
                    self.df["pl_rade"] + self.df["pl_radeerr2"]
                    > self.df["pl_maxrade12"]
                )
                & (
                    self.df["pl_rade"] + self.df["pl_radeerr2"]
                    > self.df["pl_maxrade21"]
                )
                & (
                    self.df["pl_rade"] + self.df["pl_radeerr2"]
                    > self.df["pl_maxrade22"]
                )
            ].index
        )

    def plot_filter(self, ax=None, vmin=0, vmax=2500, cmap="jet"):
        x = np.logspace(-4, 1)
        t = np.linspace(vmin, vmax, 6)

        data = np.empty([len(t), len(x)])

        for i in range(len(t)):
            data[i] = _filter(x, t[i])

        if ax == None:
            fig, ax = plt.subplots()

        cm = mpl.cm.get_cmap(cmap)
        for i in range(len(t)):
            c = (t[i] - vmin) / (vmax - vmin)
            ax.plot(x, data[i], color=cm(c), zorder=0, linewidth=2)

    def plot(
        self,
        ax=None,
        vmin=0,
        vmax=1000,
        cmap="jet",
        n_ticks=5,
        fnts1=14,
        fnts2=12,
        cbar_dims=(0.05, 0.7, 0.03, 0.25),
        annotate=False,
        xy_slices=["pl_bmasse", "pl_rade"],
        z="pl_eqt",
        save=True,
        ylims=[0, 4],
        xlims=[0, 10],
        plotSolar=False,
        stage=0,
    ):
        # Check if ax is given
        try:
            ax.set_xlabel("test")
        except AttributeError:
            fig, ax = plt.subplots(figsize=(12, 8))

        if stage > 0:
            # Plot the MR data for the archive data
            x, y = xy_slices
            if x == "pl_bmasse" and y == "pl_rade":
                ax.errorbar(
                    x=self.df[x],
                    y=self.df[y],
                    xerr=[self.df["pl_bmasseerr1"], -self.df["pl_bmasseerr2"]],
                    yerr=[self.df["pl_radeerr1"], -self.df["pl_radeerr2"]],
                    color="k",
                    linestyle="",
                    marker="",
                    zorder=10,
                )

            sc = ax.scatter(
                self.df[x],
                self.df[y],
                c=self.df[z],
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                edgecolor="k",
                zorder=20,
                marker="s",
            )

            trafo = transforms.blended_transform_factory(ax.transAxes, ax.transAxes)

            if stage > 1:
                print("stage {}".format(stage))
                # Plot pure curves
                pureCols = [(0.0, 0.0, 0.0), (0.25, 0.25, 0.25), (0.5, 0.5, 0.5)]
                pureStyles = ["-", "--", ":"]
                labels = [
                    r"$\rm Iron$",
                    r"$\rm Rock \ + \ Iron \ (Fe/Mg = 1)$",
                    r"$\rm Rock \ (Si/Mg = 1)$",
                ]
                pos = {"x": [0.45, 0.6, 0.5], "y": [0.3, 0.35, 0.45], "rot": [5, 5, 5]}

                order = [0, 2, 1]

                for i in range(stage - 1):
                    if i < 2:

                        ii = order[i]
                        curve = self.pureCurvesPlanets[ii]
                        ax.plot(
                            [pl.M_surface_is / m_earth for pl in curve],
                            [pl.R_surface_is / r_earth for pl in curve],
                            color="k",  # pureCols[ii],
                            linestyle=pureStyles[ii],
                            zorder=0,
                            linewidth=2,
                            label=labels[ii],
                        )

                    elif i == 2:
                        # Plot the filter function for T between vmin and vmax
                        self.plot_filter(ax=ax, vmin=vmin, vmax=vmax, cmap=cmap)
                        ax.text(
                            0.2,
                            0.6,
                            r"$\rm Water$",
                            # rotation = 20,
                            transform=trafo,
                            fontsize=fnts1,
                        )

                    else:
                        print("i =", i)
                        ii = order[i - 1]
                        curve = self.pureCurvesPlanets[ii]
                        ax.plot(
                            [pl.M_surface_is / m_earth for pl in curve],
                            [pl.R_surface_is / r_earth for pl in curve],
                            color="k",  # pureCols[ii],
                            linestyle=pureStyles[ii],
                            zorder=0,
                            linewidth=2,
                            label=labels[ii],
                        )
            # Plot solar system planets
            if plotSolar:
                names = ["earth", "venus", "mars", "mercury"]
                annots = []
                for name in names:
                    ax.scatter(
                        M_solar[name.lower()],
                        R_solar[name.lower()],
                        c=T_solar[name.lower()],
                        cmap=cmap,
                        vmin=vmin,
                        vmax=vmax,
                    )

                    annots.append(
                        ax.text(
                            M_solar[name.lower()] + 0.05,
                            R_solar[name.lower()] - 0.1,
                            f"${name.capitalize()}$",
                            fontsize=fnts2,
                        )
                    )
                # adjust_text(annots)#, expand_text = (1.1, 1.5))

            if annotate:

                pos = {
                    "x": [1.66, 8.5, 2.5, 2.2, 6.6, 2.6, 6.4, 5.0],
                    "y": [0.9, 2.75, 1.385, 1.521, 1.635, 1.242, 2.355, 2.133],
                }
                print(self.df["pl_name"].values)
                annots = []
                for i, txt in enumerate(self.df["pl_name"]):
                    if self.df[x].values[i] < xlims[1]:
                        annots.append(
                            ax.text(pos["x"][i], pos["y"][i], f"{txt}", fontsize=fnts2)
                        )

                adjust_text(annots)  # , expand_text = (1.1, 1.5))

            inset_ax = ax.inset_axes(cbar_dims)
            cbar = plt.colorbar(
                sc, cax=inset_ax, ticks=np.linspace(vmin, vmax, n_ticks)
            )

            cbar.ax.set_ylabel(r"$\rm Temperature \ (K)$", fontsize=fnts2)
            cbar.ax.tick_params(labelsize=fnts2)

            ticklabels = np.linspace(vmin, vmax, n_ticks).astype(int)
            cbar.ax.set_yticklabels(ticklabels)

        ax.set_xlabel(r"${\rm Mass} \ (M_\oplus)$", fontsize=fnts1)
        ax.set_ylabel(r"${\rm Radius} \ (R_\oplus)$", fontsize=fnts1)
        ax.tick_params(labelsize=fnts1)
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        # ax.legend(loc = 4, fontsize = fnts2, frameon = False)

        # plt.xticks(rotation = 45)
        if save:
            plt.savefig(
                "/home/os18o068/Documents/PHD/Projects/Planets/Data/arch_stage_{}.pdf".format(
                    stage
                ),
                format="pdf",
                bbox_inches="tight",
            )
            plt.savefig(
                "/home/os18o068/Documents/PHD/Projects/Planets/Data/arch_stage_{}.png".format(
                    stage
                ),
                format="png",
                dpi=300,
                bbox_inches="tight",
                transparent=True,
            )
            plt.close()


def MR_Grasset_2009(M, X_w):
    """Empirical M-R relation for ocean planets up to 100 M_E from
    Grasset et al. 2009.

    M: Planet mass in M_earth
    X_w: Water content in wt%
    """

    all_xis = np.array(
        [
            [1.010, 2.859e-1, -5.518e-4, 2.096e-6],
            [5.714e-3, -2.116e-4, -4.221e-6, 3.085e-8],
            [-1.704e-5, -9.015e-7, 5.397e-8, -3.166e-10],
        ]
    )

    def xi(x, xis):
        dummy = [xis[i] * x ** i for i in range(3)]
        return sum(dummy)

    coefs = [xi(X_w, all_xis.T[i]) for i in range(4)]

    a, b, c, d = coefs

    return 10 ** (np.log10(a) + (b + c * M + d * M ** 2) * np.log10(M))


def plot_MR(N_mass=20, N_oc=3, oc=[0.0, 100.0]):
    masses = np.linspace(1.0, 10.0, N_mass)
    waters = np.linspace(oc[0], oc[1], N_oc)

    data = np.empty([len(waters), len(masses)])
    fig, ax = plt.subplots()
    ax.set_xlabel(r"Mass [$M_\oplus$]")
    ax.set_ylabel(r"Radius [$R_\oplus$]")
    ax.set_xlim(1.0, 10.0)
    ax.set_ylim(0.75, 2.5)

    for i in range(len(waters)):
        for j in range(len(masses)):
            data[i][j] = MR_Grasset_2009(masses[j], waters[i])

    ticklabels = np.linspace(0, 100, 6)
    cmap = plt.cm.get_cmap("viridis")
    sm = plt.cm.ScalarMappable(cmap=cmap)
    cbar = plt.colorbar(sm, ticks=np.linspace(0, 1, 6))
    cbar.set_label(r"Water mass fraction [%]")
    cbar.ax.set_yticklabels(ticklabels.astype(int))

    c1 = 0
    c2 = 100
    for i in range(len(waters)):
        c = (waters[i] - c1) / (c2 - c1)
        print("c =", c)
        ax.plot(masses, data[i], color=cmap(c))
