# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 14:05:05 2018

@author: os18o068
"""

import matplotlib.ticker
from matplotlib import pyplot as plt
from matplotlib.ticker import (
    MultipleLocator,
    FormatStrFormatter,
    AutoMinorLocator,
    LogLocator,
    FixedLocator,
)
from pics.physicalparams import (
    T0_list,
    K0_list,
    K0prime_list,
    rho0_list,
    aP_list,
    molar_mass_list,
    Rgas,
    mH,
    mO,
    kB,
    G,
    EOS_type,
    material_list,
    EOSstring_list,
    NA,
    P_earth,
    m_earth,
    r_earth,
    material_plot_list,
    solar_system,
    abbrevations_solar,
    P_zero,
    T_zero,
    mFe,
    mMg,
    mSi,
    mH2O,
    Mg_number_solar,
    temperature_jumps,
    mS,
    material_YMg,
    material_YSi,
    material_YO,
    material_YH,
    material_YS,
)
from pics.runparams import (
    eps_Psurf,
    eps_Mtot,
    eps_layer,
    param_colors,
    plot_units,
    suffix,
    eps_Mg_number,
    plot_params,
    eps_Tsurf,
    color_list,
)
import numpy as np
import time
from pics.utils.function_tools import functionTools as ftool
import sys
from tabulate import tabulate
from decimal import Decimal
from pics.atmos import Atmosphere
import astropy.table
from astropy.io import ascii
import warnings
from pics.materials import Material

# import eosfort
from pics.utils.plot_tools.plot_tools import plotTools
from pics.utils import fortfunctions
import matplotlib

matplotlib.rc("text", usetex=True)
matplotlib.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath, amssymb}"]


warnings.filterwarnings("ignore")

materials = [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 2]
# eosfort.initiate(materials=materials, ntables=len(materials))


def convergence_test(
    eps_list=np.array([1.0, 0.5, 0.1]), contents=[], T_center=None, P_center=None
):
    """This methode investigates the convergence of the integration scheme
    with increasing refinement level (i.e. decreasing eps_r) and compares
    integration using stepwise constant density and integrated density
    """
    yaxis_labels = [
        [r"$R_{tot} \ [R_{\oplus}]$", r"$P_{surf} \ [GPa]$"],
        [r"$T_{surf} \ [K]$", r"$\rho_{surf} \ [kg/m^3]$"],
    ]

    data = np.zeros([2, 2, 2, len(eps_list)])

    for d in range(len(data)):
        if d == 0:
            rhoType = "constant"
        elif d == 1:
            rhoType = "integrate"

        for e in range(len(eps_list)):
            pl = Planet(
                contents=contents,
                T_center=T_center,
                P_center=P_center,
                majorConstraint="M_tot",
                tempType="adiabatic",
                differentiated=True,
                layermasses=[0.3, 1.0],
                eps_r=eps_list[e],
                rhoType=rhoType,
            )

            pl.construct()

            data[d][0][0][e] = pl.R_tot_is / r_earth
            data[d][0][1][e] = pl.P_surface_is * 1.0e-9
            data[d][1][0][e] = pl.T_surface_is
            data[d][1][1][e] = pl.layers[-1].dens

    fig, ax = plt.subplots(2, 2)

    for dat in data:
        for i in range(len(dat)):
            for j in range(len(dat[i])):
                ax[i][j].semilogx(
                    eps_list, dat[i][j], marker="s", markerfacecolor="white"
                )

                ax[i][j].set_ylabel(yaxis_labels[i][j])

                if i == 1:
                    ax[i][j].set_xlabel(r"$ \epsilon_r $")

                if j == 1:
                    ax[i][j].tick_params(
                        labelright="on", labelleft="off", right="on", top="on"
                    )
                    ax[i][j].yaxis.set_label_position("right")

                if j == 0:
                    ax[i][j].tick_params(
                        labelright="off", labelleft="on", right="on", top="on"
                    )


def ComputeCoreMass(
    contents=None,
    Mg_number=None,
    M_surface=1.0,
    Mg_number_mantle=None,
    SiMg=None,
    M_ocean=0.0,
    xi_H_core=0.0,
):
    """Computes the core mass of a planet at given total mass, composition and
    value for Mg#
    """

    Mg_number_mantle = min(Mg_number_mantle, 0.9999999)
    Si_number_mantle = SiMg / (SiMg + 1.0 / Mg_number_mantle - 1)
    O_number_mantle = (1.0 + 2.0 * SiMg * Mg_number_mantle) / (
        2.0 + (2.0 * SiMg - 1.0) * Mg_number_mantle
    )

    fractions = fortfunctions.functions.compute_abundance_vector(
        simg=SiMg,
        femg=1.0 / Mg_number_mantle - 1.0,
        n_mats=len(contents[2]),
        ymgi=[material_YMg[i - 1] for i in contents[2]],
        ysii=[material_YSi[i - 1] for i in contents[2]],
        xih2oi=[0.0 for i in contents[2]],
        xifei=[1.0 - Mg_number_mantle for i in contents[2]],
        xialsii=[0.0 for i in contents[2]],
        xialmgi=[0.0 for i in contents[2]],
        contents=contents[2],
        additional=[],
    )

    Q1 = sum(
        [
            fractions[i] * Mg_number_mantle * material_YMg[contents[2][i] - 1]
            for i in range(len(contents[2]))
        ]
    )
    Q2 = (
        sum(
            [
                fractions[i] * Mg_number_mantle * material_YMg[contents[2][i] - 1]
                for i in range(len(contents[2]))
            ]
        )
        * mMg
        + sum(
            [
                fractions[i]
                * (1.0 - Mg_number_mantle)
                * material_YMg[contents[2][i] - 1]
                for i in range(len(contents[2]))
            ]
        )
        * mFe
        + sum(
            [
                fractions[i] * material_YSi[contents[2][i] - 1]
                for i in range(len(contents[2]))
            ]
        )
        * mSi
        + sum(
            [
                fractions[i] * material_YO[contents[2][i] - 1]
                for i in range(len(contents[2]))
            ]
        )
        * mO
    )
    Q3 = sum(
        [
            fractions[i] * (1.0 - Mg_number_mantle) * material_YMg[contents[2][i] - 1]
            for i in range(len(contents[2]))
        ]
    )
    Q4 = 1.0 + xi_H_core / 2.0
    Q5 = xi_H_core / 2.0 * (mFe + mO) + mFe + xi_H_core * mH

    core_frac = (
        (1.0 - M_ocean / M_surface)
        * (Q1 / Q2 - Mg_number * (Q1 / Q2 + Q3 / Q2))
        / (Mg_number * (Q4 / Q5 - Q1 / Q2 - Q3 / Q2) + Q1 / Q2)
    )

    Q = (
        Mg_number_mantle / (1.0 - Mg_number_mantle) * mMg
        + Si_number_mantle / (1.0 - Si_number_mantle) * mSi
        + O_number_mantle / (1.0 - O_number_mantle) * mO
        + mFe
    )

    Q1 = (M_surface - M_ocean) * (Mg_number_mantle / Mg_number - 1.0)
    Q2 = (1.0 - Mg_number_mantle) / mFe * Q + Mg_number_mantle / Mg_number - 1
    return Q1 / Q2


def PlotCoreMasses(
    M_ocean=np.linspace(0, 0.1, 3),
    SiMg=np.linspace(0, 1.0, 32),
    Mg_numbers=np.linspace(0.5, 1.0, 32),
    dim=1,
    figsize=(10, 5),
):

    if dim == 1:
        linestyles = ["-", "--", ":"]

        Mg_numbers_mantle = np.linspace(0.5, 1.0)
        Mg_numbers = np.linspace(0.2, 0.7, 6)

        core_masses = np.zeros([len(M_ocean), len(Mg_numbers), len(Mg_numbers_mantle)])
        mantle_masses = np.zeros(
            [len(M_ocean), len(Mg_numbers), len(Mg_numbers_mantle)]
        )

        for o in range(len(M_ocean)):
            for i in range(len(Mg_numbers)):
                for j in range(len(Mg_numbers_mantle)):
                    core_masses[o][i][j] = ComputeCoreMass(
                        Mg_number=Mg_numbers[i],
                        Mg_number_mantle=Mg_numbers_mantle[j],
                        M_ocean=M_ocean[o],
                        M_surface=1.0,
                        SiMg=0.1,
                        contents=[[2], [2], [6, 7], [4, 5], [1]],
                    )

                    mantle_masses[o][i][j] = 1.0 - core_masses[o][i][j] - M_ocean[o]

        plots = plotTools.Plot(
            col=2,
            row=1,
            axis_limits=[[[[0.5, 1.0], [0.0, 1.0]], [[0.5, 1.0], [0.0, 1.0]]]],
            majorlocatorx=[[0.1, 0.1]],
            majorlocatory=[[0.1, 0.1]],
            minorlocatorx=[[0.05, 0.05]],
            minorlocatory=[[0.05, 0.05]],
            axis_labels=[
                [
                    [r"$1-\xi_{\rm Fe}$", r"$M_{\rm Core}/M$"],
                    [r"$1-\xi_{\rm Fe}$", r"$M_{\rm Mantle}/M$"],
                ]
            ],
            sharex=False,
            sharey=False,
            hspace=0.1,
            wspace=0.1,
            majorlabelsize=12,
            axislabelsize=14,
            titlefontsize=12,
            logx=[[False, False]],
            logy=[[False, False]],
            title_pad=50,
            figsize=figsize,
            grid=False,
            right=False,
            top=False,
            labelright=False,
            labeltop=False,
        )

        legend_list1 = []
        legend_list2 = []
        plot_list1 = []
        plot_list2 = []

        ax = plots.ax[0]
        ax[1].tick_params(labelright=True)

        for o in range(len(M_ocean)):

            for i in range(len(Mg_numbers)):

                (pl,) = ax[0].plot(
                    Mg_numbers_mantle,
                    core_masses[o][i],
                    linestyle=linestyles[o],
                    color=color_list[i],
                    zorder=10,
                )

                ax[1].plot(
                    Mg_numbers_mantle,
                    mantle_masses[o][i],
                    linestyle=linestyles[o],
                    color=color_list[i],
                    zorder=10,
                )

                if o == 0:
                    legend_list1.append(
                        r"$\rm Mg \# = \ $" + str(round(Mg_numbers[i], 3))
                    )
                    plot_list1.append(pl)

                if i == 0:
                    legend_list2.append(
                        r"$M_{\rm Ocean}/M = \ $" + str(round(M_ocean[o], 3))
                    )
                    plot_list2.append(pl)

        legend1 = ax[0].legend(
            plot_list1,
            legend_list1,
            loc=0,
            framealpha=0.75,
            labelspacing=0,
            fancybox=True,
        )

        legend3 = ax[0].legend(
            plot_list2,
            legend_list2,
            loc=1,
            framealpha=0.75,
            labelspacing=0,
            fancybox=True,
        )
        # legend2 = ax[1].legend(plot_list1, legend_list1)

        legend1.set_zorder(100)
        legend3.set_zorder(100)

        ax[0].add_artist(legend1)
        ax[0].add_artist(legend3)
        ax[0].tick_params(labelright=False)
        ax[1].tick_params(labelleft=False)

        for axx in ax:
            axx.tick_params(direction="in", right=True, top=True)
            # axx.grid()

        plots.fig.savefig(
            "/home/os18o068/Documents/PHD/Abbildungen/Figures_paper_1/core_masses.pdf",
            format="pdf",
            bbox_inches="tight",
        )

        plt.close("all")

    elif dim == 2:

        xx, yy = np.meshgrid(Mg_numbers, SiMg)

        data = np.zeros([len(Mg_numbers), len(SiMg)])

        for i in range(len(Mg_numbers)):
            for j in range(len(SiMg)):
                data[i][j] = ComputeCoreMass(
                    Mg_number_mantle=Mg_numbers[i],
                    SiMg=SiMg[j],
                    M_ocean=0.1,
                    M_surface=1.0,
                    Mg_number=0.5,
                )

        fig, ax = plt.subplots()

        ax.imshow(data.T)

        fig.savefig(
            "/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/core_masses.pdf",
            format="pdf",
        )


class Shell:
    def __init__(
        self,
        radius=0.1,
        T=None,
        m=None,
        P=None,
        d=None,
        layer=None,
        tempType=1,
        status="bare",
        contents=[],
        FeMg=None,
        SiMg=None,
        fractions=[],
        Fe_number=[],
        X_H2O=1.0,
        gamma=1.36,
        eps_r=0.5,
        rhoType="integrate",
        d0=None,
        adiabatType=1,
        q=1.0,
        eps_T_zero=0.0,
        saturation=False,
        xiBr=1.0,
        **kwargs,
    ):
        """A shell is defined as a the volume enclosed from r to r+dr. It
        consists of a 'Mixture' instance containing all the relevant material
        properties. The unit cell is defined as the volume normalized shell
        size. In this sense, each shell, regardless of the actual size,
        consists of a unit cell. Note that the enclosed mass and the radius of
        the shell are shell properties and not material properties. A mixture
        instance does not have these properties.
        """

        # after integration the all parameters correspond to the top of the layer
        # at r+dr, at initiation radius is the inner radius of the shell
        self.radius = radius
        self.mass = m  # mass contained inside of the shell at given radius
        self.pres = P
        self.temp = T
        self.layer = layer  # this cannot change in the life time of the shell
        self.SiMg = SiMg
        self.FeMg = FeMg
        self.N_Fe = 0.0
        self.N_Si = 0.0
        self.N_Mg = 0.0
        self.N_H2O = 0.0
        self.N_tot = 0.0
        self.adiabatType = adiabatType
        self.eps_T_zero = eps_T_zero
        self.saturation = saturation

        # compute gravitanional acceleration
        try:
            self.gravity = G * self.mass / self.radius**2

        except TypeError:
            self.gravity = None

        # compute escape velocity at upper edge of the shell
        try:
            self.v_esc = np.sqrt(2.0 * G * self.mass / self.radius)

        except TypeError:
            self.v_esc = None

        # compute thermal velocity of water vapor at current conditions
        try:
            self.v_th_H2O = np.sqrt(6 * kB * self.temp * NA / (mO + 2 * mH))

        except TypeError:
            self.v_th_H2O = None

        self.indigenous_mass = 0.0
        self.contents = contents
        self.fractions = fractions
        self.Fe_number = Fe_number

        # We assume that in a given shell all material components have the
        # same hydration level (betwen 0=dry and 1=sat)
        self.X_H2O = X_H2O
        self.dPdr = None
        self.t_fortran = 0.0
        self.t_python = 0.0
        self.t_update = 0.0
        self.t_mixcompute = 0.0
        self.xiH2O_Ws = 0.0
        self.xiH2O_Ol = 0.0
        self.xiH2O_Perov = 0.0
        self.xiEn = None
        self.xiOl = None
        self.xiWs = None
        self.xiPerov = None
        self.xiBr = xiBr
        self.xiFe_mantle = FeMg / (1.0 + FeMg)

        if not sum(self.fractions) == 1.0:
            pass
            """
            print ('WARNING: invalid material fractions given')
            print ('got:', self.fractions)
            """
        # Initiate shell material as mixture using the correct abundances
        self.mix = Material.Mixture(
            T=self.temp,
            P=self.pres,
            contents=self.contents,
            fractions=self.fractions,
            Fe_number=self.Fe_number,
            X_H2O=self.X_H2O,
            saturation=self.saturation,
            **kwargs,
        )

        # Compute materials and thermal properties of the mixture
        self.mix.Compute(**kwargs)

        # Compute material fractions for current layer and water content
        # The water content is selfconsistently obtained from the phases
        # in the layer at the current temperature and pressure
        self.getAbundances()

        # Update the mean parameters of the mixture using the updated
        # material fractions
        self.mix.UpdateFractions(newfractions=self.fractions)

        # if initial density value is given set shell density accordingly
        # if not, use the temperature and pressure to compute it via the EOS
        # for the given materials and the specified mixing law. Note, that during
        # layer construction, the PT profile and the corresponding densities
        # obtained during integration will differ due to numerical errors
        if not d == None:
            self.dens = d

        else:
            self.dens = self.mix.dens

        self.dPdrho = self.mix.dPdrho
        self.dPdr = None
        self.dTdP = None
        self.drhodr = None
        self.dTdr = None
        self.dmdr = None
        self.d0 = d0
        self.gamma = gamma  # default is value for core from Sotin 2007
        self.status = status
        self.force_bisection = False
        self.eps_r = eps_r
        self.q = q

        self.tempType = tempType
        self.rhoType = rhoType

        params = [self.pres, self.mass, self.temp, self.dens]
        self.gradients = eosfort.gradients(
            r=self.radius,
            y=params,
            fractions=self.fractions,
            nmat=len(self.fractions),
            ll=self.contents,
            gammag0=self.gamma,
            temptype=self.tempType,
            q=self.q,
            d0=self.d0,
            adiabattype=self.adiabatType,
            eps_t_zero=self.eps_T_zero,
        )

        # store all initial inputs for the shell in a dictionary which can later
        # be used torestore the shells initial state at any time if needed
        self.initials = {
            "radius": self.radius,
            "T": self.temp,
            "m": self.mass,
            "P": self.pres,
            "d": self.dens,
            "layer": self.layer,
            "tempType": self.tempType,
            "dPdrho": self.dPdrho,
            "mix": self.mix,
            "gradients": self.gradients,
        }

    def prt(self, digits=3, **kwargs):
        formatstr = "{:." + str(digits) + "E}"

        print("=======================================")
        print(
            "Shell properties:\nradius:",
            self.radius * 1.0e-3,
            "km \nenclosing mass:",
            self.mass,
            "kg",
        )
        print("\nstatus:", self.status)
        print("\ntempType:", self.tempType)
        print("\nThe gradients are:")
        print(
            "\ndP/dr [Pa/m]:",
            formatstr.format(Decimal(str(self.gradients[0]))),
            "\ndm/dr [kg/m]:",
            formatstr.format(Decimal(str(self.gradients[1]))),
            "\ndT/dr [K/m]:",
            formatstr.format(Decimal(str(self.gradients[2]))),
            "\ndrho/dr [kg/m4]:",
            round(self.gradients[3], digits),
        )

        self.mix.prt(digits=digits, **kwargs)

    def getAbundances(self, **kwargs):
        """Computes the material fractions for the phases in the given layer
        selfconsistently as a function of iron content and water content. The
        water content in the saturated case is itself dependant on the location
        in the PT-plane as the saturation depends on P and T for Olivine and
        Periclase. In the latter case it is either hydrous or anhydrous. This is
        important to be done for each shell individually as it is initiated
        because in order to maintain fixed Fe/Mg and Si/Mg the material
        fractions need to be adjusted for varying water content.
        """
        # Compute material fractions for hydration model
        if self.layer < 2:
            pass

        # in layer 1 the water content in Brucite is constant and there
        # are no other hydrated phases
        elif self.layer == 2:
            # Compute water content in Periclase as weight fraction
            # phPer = eosfort.interpolate(x=self.temp, y=self.pres,
            #                           params=[10], ll=11,
            #                          nparams=1)

            phPer = self.mix.mix[0].phase
            phPv = self.mix.mix[2].phase

            if phPer == 2:
                self.xiH2O_Per = 0.0

            else:

                # Compute water content in Periclase as mole fraction
                self.xiH2O_Per = 0.5 * self.xiBr

            if phPv == 1:
                self.XH2O_Perov = 0.03 * self.X_H2O[2]

            else:
                self.XH2O_Perov = 0.0 * self.X_H2O[2]

            self.XH2O_Per = Material.eta_H2O_Ws(
                xi_H2O_Per=self.xiH2O_Per, FeMg=self.FeMg
            )  # self.mix.mix[0].X_H2O*self.X_H2O[0]

            # Compute water content in Perovskite as mole fraction
            self.xiH2O_Perov = Material.xi_H2O_Perov(
                X_H2O_Perov=self.XH2O_Perov, FeMg=self.FeMg, lay=self.layer
            )

            # Compute water content in Wuestite as mole fraction
            self.xiH2O_Ws = Material.xi_H2O_Ws(
                X_H2O_Per=self.XH2O_Per, FeMg=self.FeMg, lay=self.layer
            )

            # Compute mole fraction of Wuestite
            self.xiWs = Material.xi_Ws(
                SiMg=self.SiMg, FeMg=self.FeMg, xi_H2O_Ws=self.xiH2O_Ws, lay=self.layer
            )

            # Compute partitioning between Fe and Mg in Wuestite
            xiWs_Fe = self.xiFe_mantle * self.xiWs
            xiWs_Mg = self.xiWs - xiWs_Fe

            # Compute mole fraction of Perovskite
            self.xiPerov = 1.0 - self.xiWs

            # Compute partitioning between Fe and Mg Perovskite
            xiPerov_Fe = self.xiFe_mantle * self.xiPerov
            xiPerov_Mg = self.xiPerov - xiPerov_Fe

            self.fractions = [
                xiWs_Mg * self.xiBr,
                xiWs_Mg * (1.0 - self.xiBr),
                xiWs_Fe,
                self.xiPerov,
            ]

        # in layer 2 Olivine can take up a variable amount of water depending
        # on the pressure and temperature. In order to keep Fe/Mg and Si/Mg
        # at fixed values in the shell the molar abundances have to be
        # adjusted to the present water content
        elif self.layer == 3:
            # Get water content in Olivine as weight fraction
            # Note that for the current compositions xi_Fe = Fe#
            self.XH2O_Ol = self.mix.mix[0].X_H2O * self.X_H2O[0]

            # Compute water content in Olivine as mole fraction
            self.xiH2O_Ol = Material.xi_H2O_Ol(
                X_H2O_Ol=self.XH2O_Ol, FeMg=self.FeMg, lay=self.layer
            )

            # Compute mole fraction of Olivine
            self.xiOl = Material.xi_Ol(
                SiMg=self.SiMg, FeMg=self.FeMg, xi_H2O_Ol=self.xiH2O_Ol, lay=self.layer
            )

            # Compute mole fraction of Enstatite
            self.xiEn = 1.0 - self.xiOl

            # Compute partitioning between Fe and Mg in Enstatite
            xiEn_Fe = self.xiFe_mantle * self.xiEn
            xiEn_Mg = self.xiEn - xiEn_Fe

            # Update material fractions in the shell
            self.fractions = [self.xiOl, xiEn_Mg, xiEn_Fe]

    def getContents(self, **kwargs):
        """Computes the total amount of particles and the molar abundances of
        Mg, Fe, Si and H2O at given Fe/Mg, Si/Mg and water content
        """
        # composition is pure Fe (inner core) or Fe + FeS mixture (outer core)
        # but the inner core is not always taken into account and both inner
        # and outer core are Fe + FeS
        if self.layer < 2:

            if self.contents == [9] or self.contents == [9, 9]:
                self.N_tot = self.indigenous_mass / (mFe)

            elif self.contents == [9, 10]:
                self.N_tot = self.indigenous_mass / (mFe + 0.13 * mS)

            self.N_Fe = self.N_tot

        # composition is (Mg,Fe)O + (Mg,Fe)SiO3
        elif self.layer == 2:
            # Compute particle mass at given iron content
            m_wuestite = Material.m_Ws(xi_Fe=self.xiFe_mantle)
            m_perovskite = Material.m_Perov(xi_Fe=self.xiFe_mantle)

            # Compute total amount of moles in the shell
            self.N_tot = self.indigenous_mass / (
                self.xiPerov * (1.0 - self.xiH2O_Perov) * m_perovskite
                + self.xiWs * (1.0 - self.xiH2O_Ws) * m_wuestite
                + (self.xiWs * self.xiH2O_Ws + self.xiPerov * self.xiH2O_Perov) * mH2O
            )

            # Comptue molar abundances of constituents
            self.N_Mg = self.N_tot * (1.0 - self.xiFe_mantle)
            self.N_Fe = self.N_tot * self.xiFe_mantle
            self.N_Si = self.N_tot * self.xiPerov
            self.N_H2O = self.N_tot * (
                self.xiWs * self.xiH2O_Ws + self.xiPerov * self.xiH2O_Perov
            )

        # composition is (Mg,Fe)2SiO4 + (Mg,Fe)2Si2O6
        elif self.layer == 3:
            # Compute particle mass at given iron content
            m_enstatite = Material.m_En(xi_Fe=self.xiFe_mantle)
            m_olivine = Material.m_Ol(xi_Fe=self.xiFe_mantle)

            # Compute total amount of moles in the shell
            self.N_tot = self.indigenous_mass / (
                self.xiEn * m_enstatite
                + self.xiOl * (1 - self.xiH2O_Ol) * m_olivine
                + self.xiOl * self.xiH2O_Ol * mH2O
            )

            # Compute molar abundances of constituents
            self.N_Mg = (
                self.N_tot
                * 2
                * (1 - self.xiFe_mantle)
                * (1 - self.xiOl * self.xiH2O_Ol)
            )
            self.N_Fe = (
                self.N_tot * 2 * self.xiFe_mantle * (1 - self.xiOl * self.xiH2O_Ol)
            )
            self.N_Si = self.N_tot * (2 - self.xiOl * (1 + self.xiH2O_Ol))
            self.N_H2O = self.N_tot * (self.xiOl * self.xiH2O_Ol)

        """   
        print ('tot: ',self.N_tot)
        print ('Mg:', self.N_Mg)
        print ('Fe:', self.N_Fe)
        print ('Si:', self.N_Si)
        print ('H2O:', self.N_H2O)
        """

    def Update(self, update_grads=True, T=None, P=None, **kwargs):
        """This function is meant to update all relevant parameters
        after the Construction method has been employed
        """

        # t0=time.time()
        # first update material parameters

        if not P == None:
            self.pres = P

        if not T == None:
            self.temp = T

        # t0mix=time.time()
        self.mix.Update(P=self.pres, T=self.temp, **kwargs)
        # tmix=time.time()
        self.dens = self.mix.dens

        # Update gradients
        if update_grads:

            params = [self.pres, self.mass, self.temp, self.dens]

            self.gradients = eosfort.gradients(
                r=self.radius,
                y=params,
                fractions=self.fractions,
                nmat=len(self.fractions),
                ll=self.contents,
                gammag0=self.gamma,
                temptype=self.tempType,
                q=self.q,
                d0=self.d0,
                adiabattype=self.adiabatType,
                eps_t_zero=self.eps_T_zero,
            )

        # if a unit detects a problem, the force_bisection call is here
        # passed to the next level (from mixture to shell)
        if self.mix.force_bisection:
            self.force_bisection = True

        # update gravity and escape velocity
        self.gravity = G * self.mass / self.radius**2
        self.v_esc = np.sqrt(2.0 * G * self.mass / self.radius)

        # update thermal velocity of water vapor at these conditions assuming
        # 3+3 degrees of freedom (vibrational modes neglected)
        self.v_th_H2O = np.sqrt(6 * kB * self.temp * NA / (mO + 2 * mH))

        # Update abundances
        self.getContents()

        # t=time.time()
        # self.t_update += t-t0
        # print ('\nshell.update =', t-t0)
        # print ('mix.update=', tmix-t0mix)

    def Reset(self, **kwargs):
        """If the shell has already been constructed, reset all parameters to
        the initial values and marks the shell as bare
        """

        (
            self.radius,
            self.pres,
            self.mass,
            self.temp,
            self.dens,
            self.dPdrho,
            self.mix,
            self.gradients,
        ) = (
            self.initials["radius"],
            self.initials["P"],
            self.initials["m"],
            self.initials["T"],
            self.initials["d"],
            self.initials["dPdrho"],
            self.initials["mix"],
            self.initials["gradients"],
        )

        self.indigenous_mass = 0.0
        # self.mix = Material.Mixture(T=self.temp, P=self.pres, contents=self.contents,\
        #                  fractions = self.fractions,**kwargs)
        # t0=time.time()
        # self.mix.Compute(**kwargs)
        self.mix.Update(T=self.temp, P=self.pres)
        # t=time.time()
        # self.t_mixcompute += t-t0
        # Shell.Update(self, **kwargs)

        # if a unit detects a problem, the force_bisection call is here
        # passed to the next level (from mixture to shell)
        if self.mix.force_bisection:
            self.force_bisection = True

        self.status = "bare"

    def Construct(self, dr=1.0, overconstruct=False, fortran=False, **kwargs):
        """Here the shells are integrated from r to r+dr where dr is determined
        by the typical length scales of the shell paremeters. After integration
        the shell type of the currently processed cell is updated to 'constructed'
        and a new, blank shell is initiated. The next time the 'Construct'
        function is called, this new shell will be integrated and so forth.
        """
        if self.status == "constructed" and not overconstruct:
            print("NOTE: this shell has already been constructed!")
            print(
                """pass overconstruct=True to the Construct method to ignore
            this message"""
            )

        else:
            # generate y-vector for integration methode
            params = [self.pres, self.mass, self.temp, self.dens]
            old_mass = self.mass

            # call the integrate methode from ftool to use 4th order RK scheme
            """
            a, b = ftool.integrate(dydx=Material.gradients, \
                        start = self.radius, whicharg = 'r', y=params, h=dr, \
                        dPdrho = self.dPdrho, tempType = self.tempType, \
                        identity = 'construct shell', gamma = self.gamma, \
                        materials = self.mix.contents, N=1,\
                        fractions = self.mix.fractions, eps=self.eps_r,\
                        oldgrad=self.gradients, rhoType=self.rhoType, **kwargs)
            """
            """
            print('input to integrate:', params, self.radius, dr, len(self.contents),
                   self.fractions, self.tempType, self.contents, self.gamma, self.q,
                   self.d0, self.adiabatType, self.eps_T_zero)
            """
            # call the integrate methode from eosfort to use 4th order RK scheme
            self.radius, newy, newgrads = eosfort.integrate(
                r_start=self.radius,
                h=dr,
                y_in=params,
                gammag0=self.gamma,
                fractions=self.mix.fractions,
                temptype=self.tempType,
                nmat=len(self.contents),
                ll=self.contents,
                q=self.q,
                d0=self.d0,
                adiabattype=self.adiabatType,
                eps_t_zero=self.eps_T_zero,
            )

            self.pres, self.mass, self.temp, self.dens = newy
            self.gradients = newgrads

            # compute indigenous mass of the isolated shell
            self.indigenous_mass = self.mass - old_mass

            # print ('T/P/d after1=', self.temp, self.pres, self.dens)
            # use new pressure, temperature, radius and enclosed mass to update
            # all other parameters
            self.Update(update_grads=False)

            # print ('->T/P/d after2=', self.temp, self.pres, self.dens)

            # indicate that the integration has been performed for this shell
            self.status = "constructed"

            # if a unit detects a problem, the force_bisection call is here
            # passed to the next level (from mixture to shell)
            if self.mix.force_bisection:
                self.force_bisection = True


class Layer:
    def __init__(
        self,
        r_in=None,
        r_out=None,
        T=None,
        P=None,
        contents=[],
        fractions=[],
        m=None,
        eps_r=0.5,
        tempType=1,
        Fe_number=[],
        X_H2O=1.0,
        rhoType="integrate",
        mass_should=None,
        FeMg=None,
        SiMg=None,
        layer=None,
        gamma=1.36,
        d0=None,
        adiabatType=1,
        q=1.0,
        eps_T_zero=0.0,
        saturation=False,
        xiBr=1.0,
        **kwargs,
    ):
        """Initiates a new Layer. Temperature, pressure and enclosed mass that
        are passed __init__ refer to the layer bottom properties and serve
        as initial conditions for the construction of the layer as a multitude
        of individual shells.
        """
        self.layer = layer
        self.tempType = tempType
        self.rhoType = rhoType
        self.r_in = r_in
        self.r_out = r_out
        self.eps_r = eps_r
        self.radius = self.r_in
        self.FeMg = FeMg
        self.SiMg = SiMg
        self.Fe_number = Fe_number
        self.X_H2O = X_H2O
        self.N_tot = 0.0
        self.N_Mg = 0.0
        self.N_Si = 0.0
        self.N_H2O = 0.0
        self.xiBr = xiBr
        self.d0 = d0
        self.q = q
        self.eps_T_zero = (
            eps_T_zero  # parameter that defines how low temperatures are handeled
        )

        # set initial conditions
        self.pres, self.mass, self.temp = P, m, T
        self.mass_should = mass_should
        self.indigenous_mass = 0.0  # mass contained within the layer
        self.contents = contents
        self.fractions = fractions
        self.saturation = saturation
        self.materials = [material_list[ll] for ll in self.contents]

        # integration step which will constantly be updated as the structure
        # is computed
        self.dr = None
        self.radius_list = []
        self.plot_list = None
        self.bisec = False
        self.shellIteration = True
        self.shellIterationCount = 0
        self.changebisec = True
        # this list will be filled with all the shells in the layer
        self.shells = []
        self.force_bisection = False
        self.t_construct = 0.0
        self.t_others = 0.0
        self.t_python = 0.0
        self.t_fortran = 0.0
        self.t_update = 0.0
        self.t_mixcompute = 0.0
        self.gamma = gamma
        self.adiabatType = adiabatType

        # if no fractions have been specified, distribute the components evenly
        # over the mixture
        if len(fractions) == 0:
            print("setting up uniform fractions")
            frac = 1.0 / len(self.contents)
            self.fractions.append(frac)

        if len(Fe_number) == 0:
            print("setting Fe# = 0.0")
            self.Fe_number = [0.0 for i in range(len(self.contents))]

        # generate boundary shell
        # this shell will NOT be constructed to preserve the initial state
        # at the layer boundary

        try:
            shell_boundary = Shell(
                T=self.temp,
                P=self.pres,
                contents=self.contents,
                fractions=self.fractions,
                radius=self.radius,
                m=self.mass,
                tempType=self.tempType,
                layer=self.layer,
                Fe_number=self.Fe_number,
                X_H2O=self.X_H2O,
                rhoType=self.rhoType,
                FeMg=self.FeMg,
                SiMg=self.SiMg,
                gamma=self.gamma,
                d0=self.d0,
                adiabatType=self.adiabatType,
                q=self.q,
                eps_T_zero=self.eps_T_zero,
                saturation=self.saturation,
                xiBr=self.xiBr,
                **kwargs,
            )

        except TypeError:
            print("Type Error in initiating layer!")
            print(self.temp, "K", self.pres * 1.0e-9, "GPa")

        # generate initial shell for construction
        shell = Shell(
            T=self.temp,
            P=self.pres,
            contents=self.contents,
            fractions=self.fractions,
            radius=self.radius,
            m=self.mass,
            tempType=self.tempType,
            Fe_number=self.Fe_number,
            X_H2O=self.X_H2O,
            rhoType=self.rhoType,
            FeMg=self.FeMg,
            SiMg=self.SiMg,
            layer=self.layer,
            gamma=self.gamma,
            d0=self.d0,
            adiabatType=self.adiabatType,
            q=self.q,
            eps_T_zero=self.eps_T_zero,
            saturation=self.saturation,
            xiBr=self.xiBr,
            **kwargs,
        )

        self.fractions = shell.mix.fractions

        # the consturction methode targets the last shell in self.shells and
        # thus the boundary shell will remain bare
        self.shells.append(shell_boundary)
        self.shells.append(shell)
        self.dens = shell.mix.dens
        self.params = [self.pres, self.mass, self.temp, self.dens]
        self.gradients = shell.gradients

        # compute typical length scales
        lengthscales_list = []
        for i in range(len(self.gradients)):
            try:
                lengthscale = abs(self.params[i] / self.gradients[i])

            # if a gradient is zero (e.g. for isothermal) the typical length scale
            # is technically infinity -> set it to a very large value
            except ZeroDivisionError:
                lengthscale = 1.0e10

            lengthscales_list.append(lengthscale)

        # determine step size for next integration
        self.dr = min(lengthscales_list) * self.eps_r

    def Update(self, **kwargs):
        """Here all the gradients are updated after the current shell of the
        layer has been constructed and the next integration step size is
        computed as the minimum of the typical length scales of all the
        relevant structure parameters multiplied by a specified scaling factor.
        """
        shell = self.shells[-1]
        self.radius, self.pres, self.mass, self.temp, self.dens = (
            shell.radius,
            shell.pres,
            shell.mass,
            shell.temp,
            shell.dens,
        )

        # self.shells[-1].Update(**kwargs)
        # if a unit detects a problem, the force_bisection call is here
        # passed to the next level (from shell to layer)
        if self.shells[-1].force_bisection:
            self.force_bisection = True

        self.params = [self.pres, self.mass, self.temp, self.dens]
        self.gradients = shell.gradients
        self.indigenous_mass = sum([s.indigenous_mass for s in self.shells])

        if not self.bisec:
            # compute typical length scales as Q/(dQ/dr) where Q denotes a parameter
            lengthscales_list = []
            for i in range(len(self.gradients)):
                try:
                    lengthscale = abs(self.params[i] / self.gradients[i])

                # allow for constant parameters (e.g. if isothermal)
                # and set length scale to large number so that it does not effect
                # the integratino step size
                except ZeroDivisionError:
                    lengthscale = 1.0e10

                lengthscales_list.append(lengthscale)

            # determine step size for next integration as scaled min. of typ. len.
            self.dr = min(lengthscales_list) * self.eps_r

    def major_bisection(self, param_is, param_should, eps, direction):
        """Performs major bisection step.
        Param_should is exceeded or negative -> perform bisection step
        (all parameters valid as majorConstraints must physically always be
        positive, e.g. total mass, surface pressure, surface temperature).
        The direction argument defines the deviation sign for exceeded and
        not reached respectively. If direction=[-1, -1] (e.g. P_surf) that means
        that param_is decreases with increasing iteration variable (in this case
        the iteration variable is the integration variable r). If direction
        = [1, -1] (e.g. M_tot) that means that param_is increases with
        increasing iteration variable (again, here with r)
        """

        if (
            direction[0] * (param_should - param_is) / param_is < direction[1] * eps
            or param_is < 0.0
        ):
            self.shells.pop(-1)
            self.shells[-1].reset()
            self.Update()
            self.dr = self.dr * 0.5
            self.bisec = True
            self.changebisec = False

        elif abs(param_should - param_is) / param_is < eps:
            self.majorComplete = True
            # the shell.Construct() methode automatically adds a new
            # shell to the layer after construction. If the surface
            # of the planet is already reached this shell will not be
            # needed for further integration -> remove it
            self.shells.pop(-1)
            self.Update()
            self.shellIteration = False

        else:
            if self.changebisec:
                self.bisec = False

    def Construct(self, **kwargs):
        """The Construct method integrates the current shell from r to r+dr
        and subsequently updates all the corresponding parameters. To construct
        the entire layer, the construction method needs to be called several
        times until the desired layer conditions are met.
        """
        # t0=time.time()
        # Perform integration step
        # print ('bisec =', self.bisec)
        # print ('dr =', self.dr)
        # print (self.shells[-1].radius, self.shells[-1].mass, self.shells[-1].pres)
        # print (self.shells[-1].gradients)
        self.shells[-1].Construct(dr=self.dr, **kwargs)

        self.shellIterationCount += 1
        # t=time.time()
        # self.t_construct += t-t0
        # self.t_fortran += self.shells[-1].t_fortran
        # self.t_python += self.shells[-1].t_python

        # if a unit detects a problem, the force_bisection call is here
        # passed to the next level (from shell to layer)
        if self.shells[-1].force_bisection:
            self.force_bisection = True

        self.indigenous_mass = sum([s.indigenous_mass for s in self.shells])

        # update all parameters using the freshly comptued shell properties
        self.Update()

        self.r_out = self.radius

        # if bisection is enforced, the old shell will be re-constructed
        # using a smaller integration step and no new shell needs to be
        # initiated at this stage

        if not self.force_bisection:
            # initiate new shell
            shell = Shell(
                T=self.temp,
                P=self.pres,
                contents=self.contents,
                fractions=self.fractions,
                radius=self.radius,
                d=self.dens,
                m=self.mass,
                tempType=self.tempType,
                eps_r=self.eps_r,
                Fe_number=self.Fe_number,
                X_H2O=self.X_H2O,
                rhoType=self.rhoType,
                FeMg=self.FeMg,
                SiMg=self.SiMg,
                layer=self.layer,
                gamma=self.gamma,
                d0=self.d0,
                adiabatType=self.adiabatType,
                q=self.q,
                eps_T_zero=self.eps_T_zero,
                saturation=self.saturation,
                xiBr=self.xiBr,
                **kwargs,
            )

            # shell.Update()
            # if the density is held constant in each integration step, the
            # gradients are computed with the old density during integration. If
            # the

            # in the end a new shell is generated and added to the layer, ready
            # for the next integration step
            self.shells.append(shell)

        # print ('grads =', self.gradients)
        # print ('r=', self.radius, 'm=', self.mass,  'rho=', self.dens, 'temp=', self.temp, 'pres =', self.pres)

    def construct_all_shells(self):
        """Equivalent to Planet.construct method but only the layer at hand
        is integrated. Integration is performed from the initial layer
        properties up to a fixed mass. The normal bisection algorithm is
        employed to match the given total layer mass at the top by adjusting
        dr iteratively. This method was included to allow for more flexibility
        in terms of the strategy to create a planet. For instance, if some
        layer properties depend on the properties of the layer bellow it, a
        planet can first be created and constructed. Then the properties can
        be extracted in a second step and in a third step this method can be
        used to add an additional layer on top of the constructed planet. This
        is useful e.g. for the Mg(OH)2 dissociation model. But if no such
        specific constraints are required, the Planet.consturct() method is
        in most cases more straight forward to use to construct the entire
        planet in one go.
        """

        self.shellIteration = True
        self.shellIterationCount = 0
        while self.shellIteration:

            self.construct()
            self.shellIteration += 1

            # if outer pressure or temperature drop bellow
            # threshold value force cancellation of iteration
            if self.pres < P_zero:
                print(
                    "Zero pressure reached at",
                    self.pres,
                    "Pa after",
                    self.shellIterationCount,
                    "iterations",
                )

                # remove last integrated shell and the already prepared
                # bare shell for the next integration step
                self.shells.pop(-1)
                self.shells.pop(-1)
                self.Update()
                self.shellIteration = False

            if self.temp < T_zero:
                print(
                    "Zero pressure reached at",
                    self.temp,
                    "K after",
                    self.shellIterationCount,
                    "iterations",
                )

                # remove last integrated shell and the already prepared
                # bare shell for the next integration step
                self.shells.pop(-1)
                self.shells.pop(-1)
                self.Update()
                self.shellIteration = False

            # Check for water vapor in pure H2O layer
            if self.contents == [0]:
                if self.shells[-1].mix.mix[0].phase == "vapor":
                    print("phase=", self.shells[-1].mix.mix[0].phase)
                    print(
                        "Water vapor phase reached at",
                        round(self.temp, 3),
                        "K, ",
                        round(self.pres * 1.0e-5, 3),
                        " bar",
                        "after",
                        self.shellIterationCount,
                        "iterations",
                    )

                    # remove last integrated shell and the already prepared
                    # bare shell for the next integration step
                    self.shells.pop(-1)
                    self.shells.pop(-1)
                    self.Update()
                    self.shellIteration = False

            self.major_bisection(
                param_is=self.mass / m_earth,
                param_should=self.mass_should,
                eps=eps_Mtot,
                direction=[1, -1],
            )

    def Plot(
        self,
        pres_unit="GPa",
        temp_unit="K",
        dens_unit="kg/m3",
        mass_unit="m_earth",
        log=True,
        **kwargs,
    ):

        title_list = [
            "Pressure [" + pres_unit + "]",
            "Mass [" + mass_unit + "]",
            "Density [" + dens_unit + "]",
            "Temperature [" + temp_unit + "]",
        ]

        radius_list = np.array([shell.radius for shell in self.shells])
        plot_list = np.array(
            [
                [shell.pres for shell in self.shells],
                [shell.mass / m_earth for shell in self.shells],
                [shell.dens for shell in self.shells],
                [shell.temp for shell in self.shells],
            ]
        )

        fig, ax = plt.subplots(2, 2, sharex=True)
        plt.subplots_adjust(hspace=0.5, wspace=0.5)
        axes = [ax[0][0], ax[0][1], ax[1][0], ax[1][1]]
        axes[2].set_xlabel("Radius [km]")
        axes[3].set_xlabel("Radius [km]")

        for i in range(len(axes)):
            if log:
                axes[i].semilogy(radius_list / 1000, plot_list[i])
            else:
                axes[i].plot(radius_list / 1000, plot_list[i])

            axes[i].set_ylabel(title_list[i])


class Planet:
    def __init__(
        self,
        P_center=1.0e12,
        T_center=3000,
        R_seed=0.1,
        contents=[[9]],
        tempType=1,
        layerType="homogeneous",
        fractions=[[1]],
        P_surface_should=P_earth,
        majorConstraint="P_surface",
        layerradii=[],
        layerConstraint=["none", "none", "none", "none"],
        add_atmosphere=False,
        eps_r=0.5,
        layermasses=[],
        layerpres=[],
        M_surface_should=1.0,
        T_surface_should=300,
        R_surface_should=1.0,
        P_dissoc=5.0e11,
        f_dissoc=1.0,
        Fe_number=[],
        X_H2O=[],
        rhoType="integrated",
        echo=False,
        brucite_dissoc=False,
        vapor_stop=False,
        Mg_number_should=Mg_number_solar,
        Si_number_should=0.0,
        Mg_number_layer=False,
        subphase_res=32,
        modelType=0,
        SiMg=0.0,
        FeMg_mantle=0.0,
        temp_jumps=[0.0, 0.0, 0.0, 0.0, 0.0],
        gammas_layer=[2.43, 2.43, 1.96, 1.26, 1.0],
        layertemps=[0, 0, 0, 0, 0],
        T_zero=[10, 10, 10.0, 10.0, 10.0],
        inner_core_frac=None,
        M_ocean_should=0.0,
        M_core_should=None,
        ocean_frac_should=0.0,
        adiabatType=1,
        q=[0.489, 0.489, 2.5, 2.9, 1.0],
        silence=False,
        measureTime=False,
        eps_T_zero=0.0,
        xiBr=1.0,
        **kwargs,
    ):

        self.measureTime = measureTime

        # print ('initializing planet ...')
        self.status = "shadow of the future"
        self.M_H2 = 0.0

        # pressure level at which molecular dissociation of water occurs
        self.P_dissoc = P_dissoc

        # dissociation fraction
        self.f_dissoc = f_dissoc

        # dissociation mechnism for Mg(OH)2 -> MgO + H2O in phase 2
        self.brucite_dissoc = brucite_dissoc

        # if set to True construction is stopped if vapor phase is reached
        self.vapor_stop = vapor_stop
        self.M_tot_should = None
        self.R_tot_should = None
        self.M_tot_is = None
        self.R_tot_is = None
        self.R_seed = R_seed
        self.rho_center = None
        self.P_center = P_center
        self.T_center = T_center
        self.P_surface_should = P_surface_should
        self.M_surface_should = M_surface_should * m_earth
        self.P_space_should = P_surface_should
        self.P_space_is = None
        self.T_surface_should = T_surface_should
        self.R_surface_should = R_surface_should * r_earth
        self.contents = contents
        self.temp_jumps = temp_jumps
        self.materials = [[material_list[i] for i in con] for con in contents]
        self.fractions = fractions
        self.Fe_number = Fe_number
        self.X_H2O = X_H2O
        self.v_esc_surface = 0.0
        self.gravity = 0.0
        self.Mg_number_is = 0.0
        self.Mg_number_should = Mg_number_should
        self.Si_number_is = 0.0
        self.Si_number_should = Material.Si_number(SiMg, Mg_number_should)
        self.SiMg = SiMg
        self.FeMg_mantle = FeMg_mantle
        self.echo = echo
        self.M_H2O_should = 0.0
        self.M_H2O_is = 0.0
        self.M_H2O_hidden_is = 0.0
        self.delta_M_H2O = 0.0
        self.delta_M_H2O_hidden = 0.0
        self.H2O_count = 0.0
        self.Fe_count = 0.0
        self.density_inversion = False
        self.match = True
        self.vapor_reached = False
        self.N_shells = 0
        self.test_variable = 0.0
        self.Mg_count = 0.0
        self.delta_Mg_count = 0.0
        self.Mg_number_layer = Mg_number_layer
        self.subphase_refined = False
        self.subphase_res = subphase_res
        self.subphase_refine_inhibit = False
        self.gammas_layer = gammas_layer
        self.adiabatType = adiabatType
        self.xiBr = xiBr
        self.eps_T_zero = eps_T_zero

        self.M_outer_mantle = 0.0
        self.M_ocean_is = 0.0
        self.M_DWR_is = 0.0
        self.M_ocean_should = M_ocean_should
        self.M_core_should = M_core_should
        self.M_core_is = 0.0
        self.ocean_frac_should = ocean_frac_should
        self.ocean_frac_is = -10
        self.T_zero = T_zero

        # defines the physical model that is invoked for the planet
        # e.g Brucite dissociation model without silicates: modelType = 0
        # " " " with hydrated OlivFalseine: modelType = 1
        self.modelType = modelType

        self.layermasses = layermasses
        self.layerradii = layerradii
        self.layertemps = layertemps
        self.layerpres = layerpres
        self.d0 = []
        self.q = q

        # It is important that the inner core fraction is accounted for properly
        # If a BC are to met via toolkit.iterate(), the inner core fraction
        # must be updated each time the core mass is updated otherwise the
        # iteration will in many cases crash and no solution can be found
        if not inner_core_frac == None:
            self.inner_core_frac = inner_core_frac
            self.layermasses[0] = (
                self.layermasses[1] * inner_core_frac / (1.0 - inner_core_frac)
            )
            # self.layermasses[1] = self.M_core_should*(1.-self.inner_core_frac)

        else:
            try:
                M_core = self.layermasses[0] + self.layermasses[1]
                self.inner_core_frac = self.layermasses[0] / M_core

            except IndexError:
                self.inner_core_frac = inner_core_frac

        self.eps_r = eps_r

        if len(self.contents) > 0:
            self.differentiated = True

        else:
            self.differentiated = False

        self.shellIterationCount = 0
        self.layerIterationCount = 0
        self.changebisec = True
        self.shellIteration = True
        self.layerIteration = True

        self.tempType = tempType
        self.layerType = layerType
        self.rhoType = rhoType
        self.majorConstraint = majorConstraint
        self.layerConstraint = layerConstraint
        self.layers = []
        self.lay = 0
        self.test = False
        self.atmosphere = None
        self.add_atmosphere = add_atmosphere
        self.have_atmosphere = False
        self.majorComplete = False
        self.layerComplete = False
        self.layer_properties = [
            {
                "P_outer": None,
                "T_outer": None,
                "rho_outer": None,
                "indigenous_mass": 0.0,
            }
            for i in range(len(self.contents))
        ]

        self.silence = silence

        # if no fractions have been specified, distribute the components evenly
        # over the mixture
        if len(fractions) == 0:
            if not self.silence:
                print("setting up uniform fractions")
            for c in range(len(self.contents)):
                frac = 1.0 / len(self.contents[c])
                self.fractions.append([frac for i in self.contents[c]])

        # if modelType = 1 compute material fractions in the mantle from
        # the given values for Mg# and Si#
        if self.modelType == 1:
            xi_per = Material.xi_per(self.SiMg, lay=1)
            xi_ol = Material.xi_ol(self.SiMg, lay=1)
            self.fractions[1] = [xi_per, xi_ol, xi_ol]

            xi_per = Material.xi_per(self.SiMg, lay=2)
            xi_ol = Material.xi_ol(self.SiMg, lay=2)
            self.fractions[2] = [xi_per, xi_ol]

            self.contents[1] = [11, 5, 6]
            self.contents[2] = [11, 1]
            # print ('setting up mantle fractions for modelType 1:')
            # print (self.fractions)

        if not len(Fe_number) == len(contents):
            if not self.silence:
                print("setting Fe# = 0.0")

            self.Fe_number = [
                [0.0 for j in range(len(self.contents[i]))]
                for i in range(len(self.contents))
            ]

        if not len(X_H2O) == len(contents):
            if not self.silence:
                print("setting X_H2O = 0.0")
            self.X_H2O = [
                [0.0 for j in range(len(self.contents[i]))]
                for i in range(len(self.contents))
            ]

        if len(self.layermasses) != len(self.contents):
            if not self.silence:
                print("\nWARNING: number of layers and layermasses does not match")
                print("given layermasses:", len(self.layermasses))
                print("given layers: ", len(self.contents))

        self.saturation = []

        for l in range(len(self.contents)):
            # Check if material is water saturated
            self.saturation.append([])
            for ll in range(len(self.contents[l])):
                if self.X_H2O[l][ll] == 1.0:

                    self.saturation[l].append(True)

                else:
                    self.saturation[l].append(False)

            mix = Material.Mixture(
                contents=self.contents[l],
                Fe_number=[FeMg_mantle for ll in range(len(self.contents[l]))],
                saturation=self.saturation[l],
            )

            # Compute properties at ambient conditions
            mix.Compute(T=300, P=1.0e4)

            self.d0.append(mix.dens)
            # print (self.d0)
        # print ('-> completed')
        # initalize core seed for integration
        # print ('generating seed...')
        self.seed_material = Material.Mixture(
            contents=self.contents[0],
            fractions=self.fractions[0],
            P=self.P_center,
            T=self.T_center,
            saturation=self.saturation[0],
        )
        self.seed_material.Compute()
        self.rho_center = self.seed_material.dens
        self.M_seed = 4.0 / 3.0 * np.pi * self.rho_center * R_seed**3
        self.seed = Layer(
            m=self.M_seed,
            T=self.T_center,
            P=self.P_center,
            tempType=self.tempType,
            r_in=self.R_seed,
            contents=self.contents[0],
            fractions=self.fractions[0],
            eps_r=self.eps_r,
            Fe_number=self.Fe_number[0],
            X_H2O=self.X_H2O[0],
            rhoType=self.rhoType,
            layer=self.lay,
            FeMg=self.FeMg_mantle,
            SiMg=self.SiMg,
            gamma=self.gammas_layer[0],
            d0=self.d0[0],
            q=self.q[0],
            adiabatType=self.adiabatType,
            eps_T_zero=self.eps_T_zero,
            saturation=self.saturation[0],
            xiBr=self.xiBr,
            **kwargs,
        )

        self.layers.append(self.seed)
        self.status = "seed"

        self.M_surface_is = self.M_seed
        self.T_surface_is = T_center
        self.P_surface_is = P_center
        self.R_surface_is = R_seed
        self.exit_code = 0
        self.direction = None
        self.constraintValue_is = None
        self.constraintValue_should = None

        # gather all initial properties in a dictionary
        self.initials = {
            "P_center": self.P_center,
            "T_center": self.T_center,
            "R_seed": self.R_seed,
            "contents": self.contents,
            "tempType": self.tempType,
            "layerType": self.layerType,
            "fractions": self.fractions,
            "P_surface_should": self.P_surface_should,
            "majorConstraint": self.majorConstraint,
            "layerConstraint": self.layerConstraint,
            "add_atmosphere": self.add_atmosphere,
            "eps_r": self.eps_r,
            "layermasses": self.layermasses,
            "layerradii": self.layerradii,
            "differentiated": self.differentiated,
            "M_surface_should": self.M_surface_should / m_earth,
            "T_surface_should": self.T_surface_should,
            "R_surface_should": self.R_surface_should / r_earth,
            "Fe_number": self.Fe_number,
            "X_H2O": self.X_H2O,
            "rhoType": self.rhoType,
            "echo": self.echo,
            "brucite_dissoc": self.brucite_dissoc,
            "Mg_number_should": self.Mg_number_should,
            "Si_number_should": self.Si_number_should,
            "SiMg": self.SiMg,
            "FeMg_mantle": self.FeMg_mantle,
            "vapor_stop": self.vapor_stop,
            "Mg_number_layer": self.Mg_number_layer,
            "modelType": self.modelType,
            "temp_jumps": self.temp_jumps,
            "gammas_layer": self.gammas_layer,
            "inner_core_frac": self.inner_core_frac,
            "M_ocean_should": self.M_ocean_should,
            "ocean_frac_should": self.ocean_frac_should,
            "layertemps": self.layertemps,
            "T_zero": self.T_zero,
            "M_core_should": self.M_core_should,
            "subphase_res": self.subphase_res,
            "layerpres": self.layerpres,
            "xiBr": self.xiBr,
        }

        self.finals = {
            "P_surface_is": self.P_surface_is,
            "M_surface_is": self.M_surface_is,
            "T_surface_is": self.T_surface_is,
            "R_surface_is": self.R_surface_is,
            "Mg_number_is": self.Mg_number_is,
            "Si_number_is": self.Si_number_is,
            "M_H2O_is": self.M_H2O_is,
            "M_H2O_should": self.M_H2O_should,
            "M_H2O_hidden_is": self.M_H2O_hidden_is,
            "density_inversion": self.density_inversion,
            "match": self.match,
            "vapor_reached": self.vapor_reached,
            "layer_properties": self.layer_properties,
            "N_shells": self.N_shells,
            "M_ocean_is": self.M_ocean_is,
            "ocean_frac_is": self.ocean_frac_is,
            "M_core_is": self.M_core_is,
            "M_DWR_is": self.M_DWR_is,
        }

    def prt(self, digits=4, **kwargs):
        print("=======================================")
        print("Planet properties:")
        print("\n--------------------------------")
        print("Layer details:")
        print("--------------------------------")
        for i in range(len(self.layers)):
            print("\nlayer: ", i)
            first_shell = self.layers[i].shells[0]
            for j in range(len(self.contents[i])):
                print(
                    round(first_shell.fractions[j] * 100, digits),
                    "%",
                    material_list[self.contents[i][j]],
                    " ph =",
                    self.layers[i].shells[-1].mix.mix[j].phase,
                )

            print(
                "layer mass [M_earth]: ",
                round(self.layers[i].indigenous_mass / m_earth, digits),
            )
            print(
                "outer radius [R_earth]:",
                round(self.layers[i].radius / r_earth, digits),
            )
            print("outer P [GPa]: ", round(self.layers[i].pres * 1.0e-9, digits))
            print("outer T [K]: ", round(self.layers[i].temp, digits))

        try:
            print("\n--------------------------------")
            print("Major parameters:")
            print("--------------------------------")
            print(
                "R_surface_is [R_earth]:",
                round(self.R_surface_is / r_earth, digits),
                "\nM_surface_is [M_earth]:",
                round(self.M_surface_is / m_earth, digits),
                "\nT_surface_is [K]:",
                round(self.T_surface_is, digits),
                "\nT_center [K]:",
                round(self.T_center, digits),
                "\nP_surface_is [bar]:",
                round(self.P_surface_is * 1.0e-5, digits),
                "\nP_center [GPa]:",
                round(self.P_center * 1.0e-9, digits),
                "\nMg_number_is:",
                round(self.Mg_number_is, digits),
                "\n\nM_surface_should [M_earth]:",
                round(self.M_surface_should / m_earth, digits),
                "\nM_H2O [wt%]:",
                round(self.H2O_count * mH2O / self.M_surface_is * 100, digits),
                "\nSi_number_should:",
                round(self.Si_number_should, digits),
                "\nSi/Mg:",
                self.SiMg,
            )

        except TypeError:
            print("WARNING: Type Error in Planet.prt()")

        print("\n--------------------------------")
        print("Layer overview:")
        print("--------------------------------")

        dat = []
        for i in range(len(self.layers)):
            lay = self.layers[i]
            first_shell = self.layers[i].shells[0]
            material_str = ""
            for j in range(len(self.materials[i])):
                material_str += str(round(first_shell.fractions[j] * 100, 1)) + "% "
                material_str += self.materials[i][j]

                if j < len(self.materials[i]) - 1 and len(self.materials[i]) > 1:
                    material_str += ", "

            dat.append(
                [
                    i,
                    material_str,
                    ftool.scinot(lay.radius / r_earth, digits=digits),
                    ftool.scinot(lay.indigenous_mass / m_earth, digits=digits),
                    ftool.scinot(lay.pres * 1.0e-9, digits=digits),
                    ftool.scinot(lay.temp, digits=digits),
                    ftool.scinot(lay.dens, digits=digits),
                ]
            )

        tabl = tabulate(
            dat,
            headers=[
                "layer",
                "contents",
                "R [R_e]",
                "m [M_e]",
                "P [GPa]",
                "T [K]",
                "rho [kg m-3]",
            ],
        )

        print()
        print(f"{tabl}")
        print()

    def Update_finals(self, **kwargs):
        # gather final properties in a dictionary
        self.finals = {
            "P_surface_is": self.P_surface_is,
            "M_surface_is": self.M_surface_is / m_earth,
            "T_surface_is": self.T_surface_is,
            "R_surface_is": self.R_surface_is / r_earth,
            "Mg_number_is": self.Mg_number_is,
            "Si_number_is": self.Si_number_is,
            "M_H2O_is": self.M_H2O_is,
            "M_H2O_should": self.M_H2O_should,
            "M_H2O_hidden_is": self.M_H2O_hidden_is,
            "density_inversion": self.density_inversion,
            "match": self.match,
            "vapor_reached": self.vapor_reached,
            "layer_properties": self.layer_properties,
            "N_shells": self.N_shells,
            "M_ocean_is": self.M_ocean_is,
            "ocean_frac_is": self.ocean_frac_is,
            "M_core_is": self.M_core_is,
            "M_DWR_is": self.M_DWR_is,
        }

    def Update_initials(self, **kwargs):
        """ """
        if self.M_surface_should == float("inf"):
            print("HERE infinity in Planet.update_initials()")

        # gather all initial properties in a dictionary
        self.initials = {
            "P_center": self.P_center,
            "T_center": self.T_center,
            "R_seed": self.R_seed,
            "contents": self.contents,
            "tempType": self.tempType,
            "layerType": self.layerType,
            "fractions": self.fractions,
            "P_surface_should": self.P_surface_should,
            "majorConstraint": self.majorConstraint,
            "layerConstraint": self.layerConstraint,
            "add_atmosphere": self.add_atmosphere,
            "eps_r": self.eps_r,
            "layermasses": self.layermasses,
            "layerradii": self.layerradii,
            "differentiated": self.differentiated,
            "M_surface_should": self.M_surface_should / m_earth,
            "T_surface_should": self.T_surface_should,
            "R_surface_should": self.R_surface_should / r_earth,
            "Fe_number": self.Fe_number,
            "X_H2O": self.X_H2O,
            "rhoType": self.rhoType,
            "echo": self.echo,
            "brucite_dissoc": self.brucite_dissoc,
            "Mg_number_should": self.Mg_number_should,
            "Si_number_should": self.Si_number_should,
            "SiMg": self.SiMg,
            "FeMg_mantle": self.FeMg_mantle,
            "vapor_stop": self.vapor_stop,
            "Mg_number_layer": self.Mg_number_layer,
            "modelType": self.modelType,
            "temp_jumps": self.temp_jumps,
            "gammas_layer": self.gammas_layer,
            "inner_core_frac": self.inner_core_frac,
            "M_ocean_should": self.M_ocean_should,
            "ocean_frac_should": self.ocean_frac_should,
            "layertemps": self.layertemps,
            "T_zero": self.T_zero,
            "M_core_should": self.M_core_should,
            "subphase_res": self.subphase_res,
            "layerpres": self.layerpres,
            "xiBr": self.xiBr,
        }

    def Reset(self, **kwargs):
        """Resets all planetary properties to the initial state. This allows
        to use the exact same specifications for re-constructing the planet.
        """

        # delete old properties to release memory
        for lay in self.layers:
            del lay.shells
            del lay

        self.__init__(**self.initials)

        if self.initials["M_surface_should"] > 1.0e30:
            print("HERE large after Planet.reset() in Planet")
            print("large value=", self.initials["M_surface_should"])

        if self.M_surface_should == float("inf"):
            print("HERE infinity in Planet.reset() after")

    def AddStuff(self):
        """If a integration is performed, the newly added molar abundances
        need to be added to the counts
        """
        self.Mg_count += self.layers[self.lay].shells[-2].N_Mg
        self.Fe_count += self.layers[self.lay].shells[-2].N_Fe
        self.H2O_count += self.layers[self.lay].shells[-2].N_H2O

        try:
            self.Mg_number_is = self.Mg_count / (self.Mg_count + self.Fe_count)

        except ZeroDivisionError:
            self.Mg_number_is = 0.0

        self.M_H2O_is = self.H2O_count * mH2O / m_earth
        self.M_DWR_is = self.H2O_count * mH2O / m_earth

    def RemoveStuff(self):
        """If an integration is undone, the previously added molar abundances
        need to be subtracted from the counts again
        """
        self.Mg_count -= self.layers[self.lay].shells[-2].N_Mg
        self.Fe_count -= self.layers[self.lay].shells[-2].N_Fe
        self.H2O_count -= self.layers[self.lay].shells[-2].N_H2O

        try:
            self.Mg_number_is = self.Mg_count / (self.Mg_count + self.Fe_count)

        except ZeroDivisionError:
            self.Mg_number_is = 0.0

        self.M_H2O_is = self.H2O_count * mH2O / m_earth
        self.M_DWR_is = self.H2O_count * mH2O / m_earth

    def GetLayerConstraints(self, lay=None, skip=0, **kwarg):
        """Extract layer constraint parameter and compute the relative deviation
        with respect to the desired value.
        """
        if self.layerConstraint[lay + skip] == "mass":
            self.constraintValue_is = self.layers[lay].indigenous_mass

            # if the indigenous mass is probed in the next layer, it is currently
            # zero if it is checked from the bisection of the previous layer.
            # this has to be set manually as the next layer does not exist yet.
            if skip > 0:
                self.constraintValue_is = 0.0

            self.constraintValue_should = self.layermasses[lay + skip] * m_earth

            self.direction = 1

        elif self.layerConstraint[lay + skip] == "enclosed mass":
            self.constraintValue_is = self.layers[lay].mass

            self.constraintValue_should = (
                self.layermasses[0] + self.layermasses[1] + self.layermasses[3]
            ) * m_earth

            self.direction = 1

        # extract layer radius in m
        elif self.layerConstraint[lay + skip] == "radius":
            self.constraintValue_is = self.layers[lay].radius

            self.constraintValue_should = self.layerradii[lay + skip] * r_earth

            self.direction = 1

        elif self.layerConstraint[lay + skip] == "pres":
            self.constraintValue_is = self.layers[lay].pres

            self.constraintValue_should = self.layerpres[lay + skip]

            self.direction = -1

        elif self.layerConstraint[lay + skip] == "temp":
            self.constraintValue_is = self.layers[lay].temp

            self.constraintValue_should = self.layertemps[lay + skip]

            self.direction = -1

    def minor_bisection(self, **kwargs):

        if not self.layerConstraint[self.lay] == "none":
            # Compute constraint values
            self.GetLayerConstraints(lay=self.lay)

            # Compute relative deviation
            reldev = (
                self.constraintValue_should - self.constraintValue_is
            ) / self.constraintValue_should

            # Constraint exceeded, reset and re-integrate old shell
            # note that layermasses are in units of earth masses
            if reldev * self.direction < -eps_layer:

                # Undo the Mg, Fe and H2O counting
                self.RemoveStuff()
                # remove bare shell that was prepared for next
                # integration step
                self.layers[-1].shells.pop(-1)
                # reset the just integrated shell to re-integrate
                # it using a smaller dr
                self.layers[-1].shells[-1].Reset()
                self.layers[-1].Update()
                self.layers[-1].dr = self.layers[-1].dr * 0.5
                self.layers[-1].bisec = True
                self.subphase_refine_inhibit = True
                self.changebisec = False

                if self.layers[-1].dr < 1.0e-10:
                    # sys.exit('WARNING: dr reached zero')
                    self.layerIteration = False
                    self.ShellIteration = False
                    self.exit_code = 3

            # Check if current layer constraint has approached
            # the desired value to specified accuracy
            elif abs(reldev) < eps_layer:
                self.layerComplete = True
                self.shellIteration = False
                self.changebisec = True
                skip = 0

                # if this was the last layer, then
                # terminate internal structure integration
                if self.lay == len(self.contents) - 1:
                    if self.echo:
                        print("terminating layer iteration")

                    self.layerIteration = False

                # initiate next layer
                else:

                    # Check for the next layer if the condition is already met and
                    # in this case skip the next layer. This can occur for instance
                    # if the current layer was constrained by mass or radius and
                    # the next one via the pressure but the pressure for the next
                    # layer is already exceeded at this point. In this case the next
                    # layer does not exist.

                    if not self.lay == len(self.contents) - 1:
                        self.GetLayerConstraints(lay=self.lay, skip=1)

                        # Compute relative deviation
                        reldev = (
                            self.constraintValue_should - self.constraintValue_is
                        ) / self.constraintValue_should

                        # Check if next layer constraint has approached
                        # the desired value to specified accuracy or is exceeded already

                        if (
                            abs(reldev) < eps_layer
                            or reldev * self.direction < -eps_layer
                        ):
                            skip = 1
                            # if this was the last layer, then
                            # terminate internal structure integration
                            if self.lay + 1 == len(self.contents) - 1:
                                if self.echo:
                                    print("terminating layer iteration")

                                self.layerIteration = False

                    for i in range(1 + skip):
                        radius = self.layers[-1].radius
                        mass = self.layers[-1].mass
                        pres = self.layers[-1].pres
                        temp = self.layers[-1].temp
                        # print ('lay=', self.lay)
                        deltaT = self.temp_jumps[self.lay]
                        new_T = max(temp - deltaT, 150)
                        nextlayer = Layer(
                            m=mass,
                            T=new_T,
                            P=pres,
                            tempType=self.tempType,
                            r_in=radius,
                            contents=self.contents[self.lay + 1],
                            fractions=self.fractions[self.lay + 1],
                            eps_r=self.eps_r,
                            X_H2O=self.X_H2O[self.lay + 1],
                            Fe_number=self.Fe_number[self.lay + 1],
                            rhoType=self.rhoType,
                            layer=self.lay + 1,
                            FeMg=self.FeMg_mantle,
                            SiMg=self.SiMg,
                            gamma=self.gammas_layer[self.lay + 1],
                            d0=self.d0[self.lay + 1],
                            adiabatType=self.adiabatType,
                            q=self.q[self.lay + 1],
                            eps_T_zero=self.eps_T_zero,
                            saturation=self.saturation[self.lay + 1],
                            xiBr=self.xiBr,
                            **kwargs,
                        )

                        self.lay += 1
                        self.layers.append(nextlayer)

            # if the condition is neither met nor exceeded
            # (re-)deactivate bisection proceedure and perform
            # normal shell integration in the next iteration
            else:
                # if surface pressure bisection is ongoing, don't allow
                # to change bisection to False
                if self.changebisec:

                    self.layers[-1].bisec = False
                else:
                    pass

    def major_bisection(self, **kwargs):
        """Performs major bisection step.
        Param_should is exceeded or negative -> perform bisection step
        (all parameters valid as majorConstraints must physically always be
        positive, e.g. total mass, surface pressure, surface temperature).
        The direction argument defines the deviation sign for exceeded and
        not reached respectively. If direction=[-1, -1] (e.g. P_surf) that means
        that param_is decreases with increasing iteration variable (in this case
        the iteration variable is the integration variable r). If direction
        = [1, -1] (e.g. M_tot) that means that param_is increases with
        increasing iteration variable (again, here with r)
        """

        if self.majorConstraint == "P_surface":
            param_is = self.layers[-1].pres
            param_should = self.P_surface_should
            eps = eps_Psurf
            direction = [-1, -1]

        elif self.majorConstraint == "T_surface":
            param_is = self.layers[-1].temp
            param_should = self.T_surface_should
            eps = eps_Tsurf
            direction = [-1, -1]

        elif self.majorConstraint == "M_surface":
            param_is = self.layers[-1].mass
            param_should = self.M_surface_should
            eps = eps_Mtot
            direction = [1, -1]

        reldev = (param_should - param_is) / param_should

        # Overshoot
        if direction[0] * reldev < direction[1] * eps or param_is < 0.0:
            # Undo Mg, Fe and H2O count
            self.RemoveStuff()

            self.layers[-1].shells.pop(-1)
            self.layers[-1].shells[-1].Reset()
            self.layers[-1].Update()
            self.layers[-1].dr = self.layers[-1].dr * 0.5
            self.layers[-1].bisec = True
            self.changebisec = False

        # Surface reached
        elif abs(param_should - param_is) / param_should < eps:
            self.majorComplete = True
            # the shell.Construct() methode automatically adds a new
            # shell to the layer after construction. If the surface
            # of the planet is already reached this shell will not be
            # needed for further integration -> remove it
            self.layers[-1].shells.pop(-1)
            self.layers[-1].Update()
            self.shellIteration = False
            self.layerIteration = False
            self.layers[-1].temp -= self.temp_jumps[self.lay]

        else:
            if self.changebisec:
                self.layers[-1].bisec = False

    def Construct(self, print_time=False, fortran=True, echo=False, **kwargs):
        """N sets the number of shells in the case that no surface constrain
        type has been specified. Otherwise N will be ignored.
        """
        self.t0 = time.time()
        if self.status == "constructed":
            print(
                "This planet has already been constructed. Use the",
                "reconstruct() methode to reconstruct the planet from seed",
                "with different specifications if you wish.",
            )

        else:
            self.status = "constructing"
            # iterate over all layers
            self.lay = 0
            self.layerIterationCount = 0
            self.shellIterationCount = 0
            self.layerIteration = True

            while self.layerIteration:
                self.layerIterationCount += 1
                self.shellIteration = True
                self.changebisec = True
                self.layerComplete = False
                self.subphase_refine_inhibit = False

                if self.echo:
                    # print (f'\n> processing layer {lay}')
                    print("\nlayer properties:")

                    for m in range(len(self.layers[-1].materials)):
                        print(
                            self.layers[-1].fractions[m] * 100,
                            "%",
                            self.layers[-1].materials[m],
                        )

                    tabl = tabulate(
                        [],
                        headers=[
                            "R [R_e]",
                            "P [MPa]",
                            "m [M_e]",
                            "T [K]",
                            "rho [kg m-3]",
                            "dr [km]",
                            "layer",
                            "contents",
                        ],
                    )
                    print()
                    print(f"{tabl}")
                    print()
                    print("check")

                # if an empty layer is given, skip the integration
                if self.layermasses[self.lay] == 0.0:
                    self.shellIteration = False
                    self.layerIteration = False

                # iterate over shells that build up the current layer
                while self.shellIteration:

                    if self.initials["M_surface_should"] > 1.0e30:
                        print("HERE large in Planet.construct()")
                        print("large value=", self.initials["M_surface_should"])
                        sys.exit()

                    self.shellIterationCount += 1

                    sys.stdout.flush()
                    olddens = self.layers[-1].dens
                    # perform next integration step
                    try:
                        # print ()
                        # print ('Constructing shell', len(self.layers[-1].shells)-1)
                        self.layers[-1].Construct()
                        # print ('mass =', self.layers[-1].mass)
                        # print ('layermass =', self.layers[-1].indigenous_mass)
                        # print ('shell mass =', self.layers[-1].shells[-2].indigenous_mass)

                        # Update Mg, Fe and H2O contents
                        self.AddStuff()

                    except ValueError:
                        print(self.initials)
                        self.prt()
                        print(self.layers[-1].pres, self.layers[-1].temp)
                        sys.exit()

                    # check if brucite dissociation is reached and if so
                    # re-integrate with smaller integration step to increase
                    # dissociation resolution
                    last_shell = self.layers[-1].shells[-3]
                    shell = self.layers[-1].shells[-2]

                    # it is important, that Mg(OH)2 is the first
                    # component of self.contents because the index
                    # zero is targeted here, assuming that this corresponds
                    # to the desired Mg(OH)2 component of the mixture
                    # Now also (Mg,Fe)2SiO4 subphase refinement is performed
                    ph0 = shell.mix.mix[0].phase
                    ph1 = last_shell.mix.mix[0].phase

                    # if first component does not have phase transition and the
                    # current layer is the lower mantle, check for Pv->PPv
                    if ph0 == ph1 and self.lay == 2:
                        ph0 = shell.mix.mix[2].phase
                        ph1 = last_shell.mix.mix[2].phase

                    if ph0 != ph1 and not self.subphase_refine_inhibit:

                        # At the beginning, no additional refinement has been
                        # performed and it is False. If the Brucite dissociation
                        # is reached, it will then set it to True and perform
                        # the refinement. Because it is now set to True it
                        # cannot be done twice by mistake if another phase
                        # transition were to occur (which can't happen for
                        # Brucite, but that's a safety switch that might come in
                        # handy at some point in the future). Also this inhibits
                        # triggering additional refinement if phases change in
                        # a actual layer transition where the materials change
                        # Undo the Mg, Fe and H2O counting
                        self.RemoveStuff()
                        self.layers[-1].shells.pop(-1)
                        self.layers[-1].shells[-1].Reset()
                        self.layers[-1].Update()

                        # divide current shell into several smaller shells
                        # to better resolve the subphase transition
                        # the number of refined shells is given by subphase_res
                        new_dr = self.layers[-1].dr / self.subphase_res
                        self.layers[-1].dr = new_dr

                        # refine the shell structure n_refine times
                        for i in range(self.subphase_res):
                            self.layers[-1].Construct()
                            self.layers[-1].dr = new_dr

                            # Update Mg, Fe and H2O contents
                            self.AddStuff()

                            # Check if major bisectio needs to kick in to avoid
                            # that the subphase refine misses the major constraint
                            # and then the iteration just contunious infinitely
                            self.major_bisection()

                            # furthermore in each refinement step it needs
                            # to be checked if the Olivine dissociation
                            # pressure and therefore the lower to upper
                            # mantle transition is reached. If so, the previous
                            # refinement step needs to be re-done with a
                            # smaller integration step
                            layer_before = self.lay

                            # Probe layer transition
                            self.minor_bisection()

                            # If a layer transition has occured,
                            # the refinement has to be interupted
                            if layer_before != self.lay:
                                break

                    # check for spurious numerical behaviour by making sure
                    # that the density decreases towards the surface
                    newdens = self.layers[-1].dens
                    if (newdens - olddens) / olddens > 1.0e-3:
                        # Undo the Mg, Fe and H2O counting
                        self.RemoveStuff()

                        # use constant density for shell integration for this
                        # shell if density inversion is encountered
                        self.layers[-1].shells.pop(-1)
                        self.layers[-1].shells[-1].Reset()
                        self.layers[-1].Update()

                        # perform next integration step
                        try:
                            self.layers[-1].shells[-1].rhoType = "constant"
                            self.layers[-1].Construct()

                            # Update Mg, Fe and H2O contents
                            self.AddStuff()

                        except TypeError:
                            print(
                                """WARNING: an error occured while constructing
                                   the layer after density inversion"""
                            )
                            print("layer =", self.lay)
                            print("T=", self.layers[-1].temp, "K")
                            print("P=", self.layers[-1].pres * 1.0e-5, "bar")
                            print("rho=", self.layers[-1].dens)

                            self.layerIteration = False
                            self.shellIteration = False

                        newdens = self.layers[-1].dens

                    # if outer pressure or temperature drop bellow
                    # threshold value force cancellation of iteration
                    if self.layers[-1].pres < P_zero:
                        print(
                            "Zero pressure reached at",
                            self.layers[-1].pres,
                            "Pa after",
                            self.shellIterationCount,
                            "iterations",
                        )

                        # remove last integrated shell and the already prepared
                        # bare shell for the next integration step
                        self.layers[-1].shells.pop(-1)
                        self.layers[-1].shells.pop(-1)
                        self.layers[-1].Update()
                        self.layerIteration = False
                        self.shellIteration = False
                        self.exit_code = 1

                    if self.layers[-1].temp < self.T_zero[self.lay]:
                        print(
                            "Zero temperature reached at",
                            round(self.layers[-1].temp, 2),
                            "K after",
                            self.shellIterationCount,
                            "iterations",
                        )

                        print("layer =", self.lay)
                        print("T_zero =", self.T_zero[self.lay])

                        print("P =", self.layers[-1].pres * 1.0e-9, "GPa")

                        # remove last integrated shell and the already prepared
                        # bare shell for the next integration step
                        self.layers[-1].shells.pop(-1)
                        # self.layers[-1].shells.pop(-1)
                        self.layers[-1].Update()
                        self.layerIteration = False
                        self.shellIteration = False
                        self.exit_code = 1

                    # Check for water vapor in pure H2O layer
                    if self.contents[-1] == [0] and self.vapor_stop:
                        if (
                            self.layers[-1].shells[-1].mix.mix[0].phase == "vapor"
                            or self.layers[-1].dens < 10
                        ):
                            print(
                                "Water vapor phase reached at",
                                round(self.layers[-1].temp, 3),
                                "K, ",
                                round(self.layers[-1].pres * 1.0e-5, 3),
                                " bar",
                                "after",
                                self.shellIterationCount,
                                "iterations",
                            )
                            print("density is:", self.layers[-1].dens)
                            self.vapor_reached = True

                            # remove last integrated shell and the already prepared
                            # bare shell for the next integration step
                            self.layers[-1].shells.pop(-1)
                            self.layers[-1].shells.pop(-1)

                            self.layers[-1].Update()
                            self.shellIteration = False
                            self.layerIteration = False
                            self.exit_code = 2

                    # Check for major constraint
                    self.major_bisection()

                    if self.layers[-1].dr < 1.0e-10:
                        self.layerIteration = False
                        self.shellIteration = False
                        self.exit_code = 3

                    # second: check layer constraint condition to probe transition
                    # to the next layer

                    if (
                        self.differentiated
                        and not self.vapor_reached
                        and not self.majorComplete
                    ):

                        # Probe layer transition
                        self.minor_bisection()

                    radius = self.layers[-1].radius
                    mass = self.layers[-1].mass
                    pres = self.layers[-1].pres
                    temp = self.layers[-1].temp
                    dens = self.layers[-1].dens

                    if self.echo:
                        digits = 4

                        try:
                            dat = [
                                [
                                    ftool.scinot(radius / r_earth, digits=digits),
                                    ftool.scinot(pres * 1.0e-9, digits=digits),
                                    ftool.scinot(mass / m_earth, digits=digits),
                                    ftool.scinot(temp, digits=digits),
                                    ftool.scinot(dens, digits=digits),
                                    ftool.scinot(
                                        self.layers[-1].dr / 1000.0, digits=digits
                                    ),
                                    self.lay,
                                    self.layers[-1].contents,
                                ]
                            ]

                        except UnboundLocalError:
                            pass

                        tabl = tabulate(
                            dat,
                            headers=[
                                "R [R_e]",
                                "P [GPa]",
                                "m [M_e]",
                                "T [K]",
                                "rho [kg m-3]",
                                "dr [km]",
                                "layer",
                                "contents",
                            ],
                        )

                        sys.stdout.write("\033[3A")

                        print(f"{tabl}")
                        print()

                        if self.layerComplete:
                            print(
                                "\nDesired pressicion for layer",
                                self.lay - 1,
                                "reached after",
                                self.shellIterationCount,
                                "iterations",
                            )

                if self.majorComplete and self.echo:
                    print(
                        "\nDesired pressicion for",
                        self.majorConstraint,
                        "reached after",
                        self.shellIterationCount,
                        "iterations",
                    )

                if self.M_surface_should == float("inf"):
                    print("HERE infinity in Planet.construct()")
                    sys.exit()

            # print ('\n-> completed')
            # print ('\n==============================================='+
            #          '=========')
            self.R_tot_is = self.layers[-1].radius
            self.M_tot_is = self.layers[-1].mass
            self.R_surface_is = self.R_tot_is
            self.M_surface_is = self.M_tot_is
            self.T_surface_is = self.layers[-1].temp
            self.P_surface_is = self.layers[-1].pres
            self.P_space_is = self.P_surface_is
            self.N_shells = sum([len(layer.shells) for layer in self.layers])

            try:
                self.M_core_is = self.layers[1].mass / m_earth

            except IndexError:
                self.M_core_is = self.layers[0].mass / m_earth

            # If no ocean exists the planet will have only 4 layers and
            # the ocean layer index cannot be targeted
            try:
                self.M_ocean_is = self.layers[4].indigenous_mass / m_earth

                try:
                    self.ocean_frac_is = np.log10(
                        self.layers[4].indigenous_mass / self.layers[4].mass
                    )

                except TypeError:
                    self.ocean_frac_is = -10

            except IndexError:
                self.M_ocean_is = 0.0
                self.ocean_frac_is = -10

            # construct atmosphere on top of naked planet
            if self.add_atmosphere:
                self.atmosphere = Atmosphere.Atmosphere(
                    r_in=self.R_tot,
                    T=self.T_surface_is,
                    P=self.P_surface_is,
                    contents=[0],
                    m=self.M_surface_is,
                    tempType="isothermal",
                )

                self.atmosphere.Construct()
                self.P_space_is = self.atmosphere.P_space_is
                self.layers.append(self.atmosphere.layers[0])
                self.have_atmosphere = True

            self.v_esc_surface = np.sqrt(2 * G * self.M_tot_is / self.R_tot_is)

            self.status = "constructed"

        for layer in self.layers:
            if layer.contents == [0]:
                self.M_H2O_is = layer.indigenous_mass / m_earth

        if print_time:
            self.t = time.time()
            print("total elapsed time=", round(self.t - self.t0, 4), "s")

        if self.echo:
            self.prt()

        self.layer_properties = [
            {
                "P_out": None,
                "T_out": None,
                "rho_out": None,
                "R_out": None,
                "indigenous_mass": 0.0,
            }
            for i in range(len(self.layers))
        ]

        for i in range(len(self.layers)):
            self.layer_properties[i]["P_out"] = self.layers[i].pres
            self.layer_properties[i]["T_out"] = self.layers[i].temp
            self.layer_properties[i]["R_out"] = self.layers[i].radius
            self.layer_properties[i]["rho_out"] = self.layers[i].dens
            self.layer_properties[i]["indigenous_mass"] = self.layers[i].indigenous_mass

        self.Update_finals()
        self.Update_initials()

    def construct_atmosphere(self, contents=[0], tempType="isothermal", **kwargs):
        """With this method an atmospheric layer is added on top of the naked
        planet if <add_atmosphere=False>. If <add_atmosphere=True>, the
        atmosphere is generated automatically and calling this method will
        do nothing expect telling you that the atmosphere already exists.
        """
        if self.have_atmosphere:
            print(
                """NOTE: atmosphere has already been constructed or will
                  \n be constructed upon calling the <Construct> method."""
            )

        else:
            # construct atmosphere on top of naked planet
            self.atmosphere = Atmosphere.Atmosphere(
                r_in=self.R_surface_is / r_earth,
                T=self.T_surface_is,
                P=self.P_surface_is,
                contents=contents,
                m=self.M_surface_is / m_earth,
                tempType=tempType,
                majorConstraint="M_atmos",
                **kwargs,
            )

            self.atmosphere.construct()
            self.layers.append(self.atmosphere.layers[0])
            self.have_atmosphere = True

            self.R_tot_is = self.layers[-1].radius
            self.M_tot_is = self.layers[-1].mass
            self.P_space_is = self.atmosphere.P_space_is

    def Strip_atmosphere(self):
        if self.have_atmosphere:
            self.layers.pop(-1)
            self.have_atmosphere = False
            self.P_space = self.P_surface_is
            self.R_tot_is = self.layers[-1].radius
            self.M_tot_is = self.layers[-1].mass

        else:
            print("Planet has no atmosphere to strip away")

    def Update(self, force=False):
        """Updates global parameters. This is useful if e.g. a layer is
        added after constructing a planet to update automatically the total
        mass, surface pressure, temperature, layer compositions etc. Note that
        this method should only be used if the planet is already constructed and
        subsequently modified.
        """

        if not self.status == "constructed":
            print(
                """WARNING: You are about to update a non-constructed planet.
                   This is likely to cause trouble and is by default ommited.
                   Use the keyword argument force=True to overwrite this
                   behaviour."""
            )

        if self.status == "constructed" or force == True:
            print("updating")
            self.contents = [lay.contents for lay in self.layers]
            self.fractions = [lay.fractions for lay in self.layers]
            self.materials = [lay.materials for lay in self.layers]
            self.layermasses = [lay.indigenous_mass / m_earth for lay in self.layers]
            self.X_H2O = [lay.X_H2O for lay in self.layers]
            self.Fe_number = [lay.Fe_number for lay in self.layers]

            self.Mg_number_is = Material.Mg_number(
                contents=self.contents, layers=self.layers, X_H2O=self.X_H2O
            )

            self.R_tot_is = self.layers[-1].radius
            self.M_tot_is = self.layers[-1].mass
            self.R_surface_is = self.R_tot_is
            self.M_surface_is = self.M_tot_is
            self.T_surface_is = self.layers[-1].temp
            self.P_surface_is = self.layers[-1].pres
            self.P_space_is = self.P_surface_is

        else:
            print("wtf")

    def computeEquilibriumTemperature(
        self, N=1, T_star=5800.0, a_mean=1.0, albedo=0.3, eps_emm=0.78, eta_cloud=0.0
    ):
        """Computes the surface equilibrium temperature of the planet based
        on a variety of free parameters:

            T_star (K): effective temperature of central star
            a_mean (au): mean orbital distance
            albedo: planetary albedo
            eps_emm: emissivity of atmospherif gases
            eta_cloud: cloud coverage
            N: refinement level of atmospheric strucutre

        By default, the parameters are set to represent the idealized, cloud-
        less greenhouse model for the earth (e.g. Wikipedia).
        """
        pass

    def Plot(
        self,
        scatter=False,
        layerColoring=False,
        axis=[],
        x_axis="radius",
        save=False,
        filename="planet_structure",
        suffix="png",
        path="./",
        **kwargs,
    ):

        if axis == []:
            fig, axis = plt.subplots(2, 3, sharex=True)
            fig.subplots_adjust(hspace=0, wspace=0.3)

        title_list = [
            r"$\rmPressure \ [GPa]$",
            r"$\rmMass \ [M_\oplus]$",
            r"$\rmDensity \ [10^3 \ kg/m^3]$",
            r"$\rmTemperature  \ [K]$",
            r"$\rmGravity \ [m/s^2]$",
            r"$\rmv_{esc} \ [km/s]$",
        ]

        radius_list_dummy = np.array(
            [
                [shell.radius / r_earth for shell in layer.shells]
                for layer in self.layers
            ]
        )

        # collect data for all layers
        plot_list_dummy = [
            [[shell.pres * 1.0e-9 for shell in layer.shells] for layer in self.layers],
            [[shell.mass / m_earth for shell in layer.shells] for layer in self.layers],
            [[shell.dens * 1.0e-3 for shell in layer.shells] for layer in self.layers],
            [[shell.temp for shell in layer.shells] for layer in self.layers],
            [[shell.gravity for shell in layer.shells] for layer in self.layers],
            [[shell.v_esc * 1.0e-3 for shell in layer.shells] for layer in self.layers],
        ]

        # collect thermal velocity data for water vapor to overplot it on
        # the graph for the escape velocity for comparison

        v_th_H2O_dummy = [
            [shell.v_th_H2O * 1.0e-3 for shell in layer.shells] for layer in self.layers
        ]

        # extract data of individual layers to one array for each parameter
        plot_list = []
        radius_list = []
        v_th_H2O_list = []

        for p in range(len(plot_list_dummy)):
            plot_list.append([])

        for p in range(len(plot_list_dummy)):
            for l in range(len(self.layers)):
                for i in range(len(plot_list_dummy[p][l])):
                    plot_list[p].append(plot_list_dummy[p][l][i])

        for l in range(len(self.layers)):
            for i in range(len(radius_list_dummy[l])):
                radius_list.append(radius_list_dummy[l][i])
                v_th_H2O_list.append(v_th_H2O_dummy[l][i])

        ax = [axis[0][0], axis[0][1], axis[1][0], axis[1][1], axis[0][2], axis[1][2]]

        ax[2].set_xlabel(r"$\rmRadius  \ [R_\oplus]$")
        ax[3].set_xlabel(r"$\rmRadius \ [R_\oplus]$")
        ax[5].set_xlabel(r"$\rmRadius \ [R_\oplus]$")
        # for axx in ax:
        #   axx.grid(which='both', axis='both')

        for i in range(len(ax)):
            if x_axis == "radius":
                ax[i].plot(radius_list, plot_list[i], color=param_colors[i], zorder=3)

            elif x_axis == "pres":
                ax[i].semilogx(
                    plot_list[0], plot_list[i], color=param_colors[i], zorder=3
                )

            ax[i].tick_params(right=True, top=True, direction="in", which="both")
            ax[i].tick_params(which="major", axis="both", length=8)
            ax[i].grid(
                True,
                which="both",
                zorder=0,
                color=plot_params["gridcolor"],
                alpha=plot_params["gridalpha"],
            )

            ax[i].set_facecolor(plot_params["backgroundcol"])
            ax[i].grid(which="major", axis="both", linewidth=2, zorder=0)
            ax[i].xaxis.set_minor_locator(AutoMinorLocator())
            ax[i].yaxis.set_minor_locator(AutoMinorLocator())
            """
            for j in range(len(radius_list)):
                if scatter:
                    ax[i].scatter(radius_list[j], plot_list[i][j])
                else:
                    ax[i].plot(radius_list[j], plot_list[i][j])
             """
            ax[i].set_ylabel(title_list[i])

        ax[-1].plot(radius_list, v_th_H2O_list, color=param_colors[-1], linestyle="--")

        # save figure as image file if turned on
        if save:
            fig.savefig(path + filename + "." + suffix)

        fig2, ax2 = plt.subplots()
        ax2.loglog(np.asarray(plot_list[3]), np.asarray(plot_list[0]) * 1.0e9)

    def write(self, out="planet", loc="./"):
        """Creates data table of the planetary structure and writes it to ascii
        file. By default, the data is organized as follows:

        R           M           T       P       rho         v_esc       g
        (r_earth)   (m_earth)   (K)     (GPa)   (kg/m3)     (km/s)      (m/s2)
        -----------------------------------------------------------------------
        R_center    M_center    ...
        ...         ...         ...
        R_tot       M_tot       ...

        The meta data contains the properties that are stored in Planet.initials
        and Planet.finas by default
        """

        rad, mass, temp, pres, dens, vesc, grav = [], [], [], [], [], [], []

        # gather structure data
        for lay in self.layers:
            for sh in lay.shells:
                rad.append(sh.radius / r_earth)
                mass.append(sh.mass / m_earth)
                temp.append(sh.temp)
                pres.append(sh.pres * 1.0e-9)
                dens.append(sh.dens)
                vesc.append(sh.v_esc / 1000.0)
                grav.append(sh.gravity)

        data = [rad, mass, temp, pres, dens, vesc, grav]

        # gather all planetary parameters as meta data
        meta = {"initials": self.initials, "finals": self.finals}

        data_table = astropy.table.Table(data, meta=meta)

        names = (
            "R (r_earth)",
            "M (m_earth)",
            "T (K)",
            "P (GPa)",
            "rho (kg/m3)",
            "v_esc (km/s)",
            "g (m/s2)",
        )

        ascii.write(
            data_table,
            loc + out + suffix["planet"],
            overwrite=True,
            format="ecsv",
            names=names,
        )

    def load(self, loc="./", file_name="planet"):
        """Loads a Planet.Planet() object into a session by reading in all
        the relevant data from a planet.out file that has been generated
        by the Planet.write() methode. Note, at this point only the structure
        data and the Planet.finals and Planet.initials dictonaries are loaded.
        The Planet.Shell and Planet.Layer objects are not recreated, so they
        will not be accessable for a Planet.Planet object after this methode
        has been invoked. This might change in the future as not beeing able
        to use the full Planet.Planet instance from file is rather restrictive.
        However, just reading in the important parameters is more efficient
        and will therefore still be available in the future as there are many
        cases in which the Layer and Shell instances are not needed anyways
        (e.g. for simple replotting or adapting some initial conditions and
        re-construct the Planet in which case the Layer and Shell instances
        have to be re-initialized anyways by the Planet.py module).
        """
        # read ascii file for eos table
        self.data_table = ascii.read(loc + file_name + suffix["planet"], format="ecsv")

        self.initials = self.data_table.meta["initials"]
        self.finals = self.data_table.meta["finals"]


class Population:
    def __init__(
        self,
        contents=[],
        P_center=1.0e12,
        itType="cont",
        fractions=[],
        layermasses=[],
        T_center=300.01,
        differentiated=False,
        tempType="isothermal",
        P_dissoc=5.0e11,
        f_dissoc=1.0,
        **kwargs,
    ):

        self.contents = contents
        self.fractions = fractions
        self.layermasses = layermasses
        self.differentiated = differentiated
        self.P_center = P_center
        self.T_center = T_center
        self.planets = []
        self.itType = itType  # iteration type
        self.tempType = tempType
        self.P_dissoc = P_dissoc
        self.f_dissoc = f_dissoc

        if type(P_center) is float:
            print("setting all central pressures to common value")
            self.P_center = [P_center for con in self.contents]

        # if no fractions have been specified, distribute the components evenly
        # over the mixture
        if len(self.fractions) == 0:
            # iterate over all planets in the population
            for c in range(len(self.contents)):
                self.fractions.append([])
                # iterate over all layers of a planet
                for cc in range(len(self.contents[c])):
                    frac = 1.0 / len(self.contents[c][cc])
                    self.fractions[c].append([frac for i in self.contents[c][cc]])
        else:
            self.fractions = fractions

        print("fractions in population.__init__() =", self.fractions)

        if not len(self.P_center) == len(self.contents):
            print("WARNING: central pressures and contents do not match!")

    def generate(self, **kwargs):
        # initiate empty planet list
        self.planets = []

        # iterate over different contents
        if self.itType == "cont":

            for p in range(len(self.contents)):
                print("\n####################################################")
                print("> processing planet", p)

                try:
                    lm = self.layermasses[p]

                except IndexError:
                    lm = []
                    print("empty layer masses given")

                pla = Planet(
                    P_center=self.P_center[p],
                    T_center=self.T_center,
                    tempType=self.tempType,
                    R_seed=2.0,
                    contents=self.contents[p],
                    layerType="homogeneous",
                    fractions=self.fractions[p],
                    layermasses=lm,
                    differentiated=self.differentiated,
                    majorConstraint="P_surf",
                    P_dissoc=self.P_dissoc,
                    f_dissoc=self.f_dissoc,
                    **kwargs,
                )

                pla.construct()
                self.planets.append(pla)

        # iterate over central pressures
        elif self.itType == "pres":
            for p in range(len(self.P_center)):
                print("\n####################################################")

                try:
                    lm = self.layermasses[p]

                except IndexError:
                    lm = []

                pla = Planet(
                    P_center=self.P_center[p],
                    T_center=self.T_center,
                    tempType=self.tempType,
                    R_seed=2.0,
                    contents=self.contents[p],
                    layerType="homogeneous",
                    fractions=self.fractions[p],
                    layermasses=lm,
                    differentiated=self.differentiated,
                    **kwargs,
                )

                pla.construct(N=400)
                self.planets.append(pla)

    def plot(self, **kwargs):
        fig, ax = plt.subplots(2, 2)
        for planet in self.planets:
            planet.Plot(axis=ax)


class MR_relation:
    def __init__(
        self,
        P_min=1.0e10,
        P_max=5.0e12,
        T=300,
        contents=[],
        fractions=[],
        layermasses=[],
        tempType="isothermal",
        T_center=300,
        differentiated=False,
        f_dissoc=[],
        P_dissoc=[],
        **kwargs,
    ):

        self.P_min = P_min
        self.P_max = P_max
        self.N = 5
        self.temp = T
        self.T_center = T_center
        self.tempType = tempType
        self.contents = contents
        self.layermasses = layermasses
        self.differentiated = differentiated
        self.material_list = []
        self.fractions = fractions
        self.M_tot = []
        self.R_tot = []
        self.populations = []
        self.label_list = None
        self.legend_list = None
        self.f_dissoc = f_dissoc
        self.P_dissoc = P_dissoc

        # if no fractions have been specified, distribute the components evenly
        # over the mixture
        if len(self.fractions) == 0:
            # iterate over all planets in the population
            for c in range(len(self.contents)):
                print("c:", c)
                self.fractions.append([])
                # iterate over all layers of a planet
                for cc in range(len(self.contents[c])):
                    frac = 1.0 / len(self.contents[c][cc])
                    self.fractions[c].append([frac for i in self.contents[c][cc]])

        else:
            self.fractions = fractions

        print("fractions=", self.fractions)
        print("contents=", self.contents)

    def generate(self, N=5, scaling="log", **kwargs):
        self.N = N

        if scaling == "lin":
            self.P_center = np.linspace(self.P_min, self.P_max, self.N, endpoint=True)

        elif scaling == "log":
            self.P_center = np.logspace(
                np.log10(self.P_min), np.log10(self.P_max), self.N, endpoint=True
            )

        for i in range(len(self.contents)):
            print("> processing population", i)
            con = self.contents[i]

            if self.differentiated:
                lm = self.layermasses[i]

            else:
                lm = self.layermasses

            contents = [con for j in range(len(self.P_center))]
            layermasses = [lm for j in range(len(self.P_center))]
            frac = self.fractions[i]
            fractions = [frac for j in range(len(self.P_center))]
            pop = Population(
                contents=contents,
                P_center=self.P_center,
                T_center=self.T_center,
                fractions=fractions,
                layermasses=layermasses,
                tempType=self.tempType,
                differentiated=self.differentiated,
                **kwargs,
            )

            pop.generate()
            self.populations.append(pop)
            mat = [material_plot_list[m] for m in contents[0][0]]
            self.material_list.append(mat)

    def plot(
        self,
        axis=None,
        replot=False,
        solar=False,
        log=True,
        colors=[],
        linestyles=[],
        legends=[],
        tick_params={},
        **kwargs,
    ):
        # replot means that the data has already been plotted and the
        # corresponding object attributes have previously been generated
        # this is useful for replotting with adjusted plot parameters
        # if replot is set to True, it is assumed, that all the plot lists
        # already exist or are manually added to the MR_relation instance
        if not replot:
            self.M_tot = []
            self.R_tot = []
            self.legend_list = []
            self.plot_list = []

            # gather planetary masses and radii for all planets in each population
            # iterate over populations
            for pop in self.populations:
                self.M_tot.append([])
                self.R_tot.append([])

                label_list = []
                # itearte over planets
                for pla in pop.planets:
                    self.M_tot[-1].append(pla.M_tot_is)
                    self.R_tot[-1].append(pla.R_tot_is)

                if not len(legends) == len(self.populations):

                    if len(legends) > 0:
                        print("WARNING: invalid input for <legends> given")

                    # generate content labels
                    # iterate over layers of the last planet in the population
                    for i in range(len(pla.contents)):
                        # iterate over layer compositions
                        for j in range(len(pla.contents[i])):
                            lab = " ".join(
                                [
                                    str(100 * pla.fractions[i][j]) + r"$\% \ $",
                                    material_plot_list[pla.contents[i][j]],
                                ]
                            )

                            label_list.append(lab)

                    label = ", ".join(label_list)
                    self.legend_list.append(label)

                else:
                    self.legend_list = legends

        else:
            pass

        # if an existing axis instance is passed, plot on this axis else
        # create a new one
        if axis == None:
            self.fig, self.ax = plt.subplots()

        else:
            print("adopting axis instance:", axis)
            self.ax = axis

        for m in range(len(self.M_tot)):
            self.ax.set_xlabel(r"$M/M_\oplus$")
            self.ax.set_ylabel(r"$R/R_\oplus$")

            if len(colors) == len(self.populations):
                col = colors[m]

            else:
                col = None

            if len(linestyles) == len(self.populations):
                lnst = linestyles[m]

            else:
                lnst = None

            if log:
                (plot,) = self.ax.loglog(
                    np.asarray(self.M_tot[m]) / m_earth,
                    np.asarray(self.R_tot[m]) / r_earth,
                    color=col,
                    linestyle=lnst,
                    marker="o",
                )

            else:
                (plot,) = self.ax.plot(
                    np.asarray(self.M_tot[m]) / m_earth,
                    np.asarray(self.R_tot[m]) / r_earth,
                    color=col,
                    linestyle=lnst,
                    marker="o",
                )

            self.plot_list.append(plot)

        # add selection of solar system objects for comparison
        if solar:
            self.ax.scatter(solar_system[0], solar_system[1], color="k", s=5, zorder=10)

            for i, txt in enumerate(abbrevations_solar):
                self.ax.annotate(txt, (solar_system[0][i], solar_system[1][i]))

        legend = self.ax.legend(self.plot_list, self.legend_list)
        self.ax.add_artist(legend)
        self.ax.tick_params(direction="in", top=True, right=True, which="both")
        self.ax.tick_params(which="major", length=8)
        self.ax.tick_params(which="minor", length=4)
        plt.show()
