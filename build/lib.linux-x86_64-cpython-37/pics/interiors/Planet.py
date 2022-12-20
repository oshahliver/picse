# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 14:05:05 2018

@author: os18o068
"""

import matplotlib.ticker
from matplotlib import pyplot as plt
from PIMPphysicalparams import (
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
    Mg_number_solar,
)

from PIMPrunparams import (
    eps_Psurf,
    eps_Mtot,
    eps_layer,
    param_colors,
    plot_units,
    suffix,
    eps_Mg_number,
)
import numpy as np
import time

# import eos
import pics.utils.functionTools as ftool
import sys
from tabulate import tabulate
from decimal import Decimal
import pics.interiors.Atmosphere as Atmosphere
import astropy.table
from astropy.io import ascii
import warnings
import pics.materials.Material as Material

# import eosfort

warnings.filterwarnings("ignore")

materials = [0, 1, 2, 5, 6, 9, 11]
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


class Shell:
    def __init__(
        self,
        radius=0.1,
        T=None,
        m=None,
        P=None,
        d=None,
        layer=0,
        tempType=0,
        status="bare",
        contents=[],
        fractions=[],
        Fe_number=[],
        X_H2O=[],
        gamma=1.36,
        eps_r=0.5,
        rhoType="integrate",
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

        # compute gravitanional acceleration
        try:
            self.gravity = G * self.mass / self.radius ** 2

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
        self.X_H2O = X_H2O
        self.dPdr = None
        self.t_fortran = 0.0
        self.t_python = 0.0
        self.t_update = 0.0
        self.t_mixcompute = 0.0

        if not sum(self.fractions) == 1.0:
            print("WARNING: invalid material fractions given")
            print("got:", self.fractions)

        self.mix = Material.Mixture(
            T=self.temp,
            P=self.pres,
            contents=self.contents,
            fractions=self.fractions,
            Fe_number=self.Fe_number,
            X_H2O=self.X_H2O,
            **kwargs,
        )
        # t0=time.time()
        self.mix.Compute(**kwargs)
        # t=time.time()
        # self.t_mixcompute += t-t0
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
        self.gamma = gamma  # default is value for core from Sotin 2007
        self.status = status
        self.layer = layer  # this cannot change in the life time of the shell
        self.force_bisection = False
        self.eps_r = eps_r

        self.tempType = tempType
        self.rhoType = rhoType

        self.dPdr = -G * self.mass * self.dens / self.radius ** 2

        if self.rhoType == "constant":
            self.drhodr = 0.0

        else:
            self.drhodr = self.dPdr / self.dPdrho

        try:
            self.dTdP = self.gamma * self.temp / (self.dens * self.dPdrho)

        except TypeError:
            print("WARNING: cannot compute dTdP!")
            print("dens=", self.dens, "\ndPdrho=", self.dPdrho)

        self.dTdr = self.dPdr * self.dTdP
        self.dmdr = 4.0 * np.pi * self.radius ** 2 * self.dens

        if self.tempType == 0:
            self.dTdr = 0.0
            self.dTdP = 0.0

        self.gradients = [self.dPdr, self.dmdr, self.dTdr, self.drhodr]
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
            formatstr.format(Decimal(str(self.dPdr))),
            "\ndm/dr [kg/m]:",
            formatstr.format(Decimal(str(self.dmdr))),
            "\ndT/dr [K/m]:",
            formatstr.format(Decimal(str(self.dTdr))),
            "\ndrho/dr [kg/m4]:",
            round(self.drhodr, digits),
            "\ndT/dP [K/Pa]:",
            formatstr.format(Decimal(str(self.dTdP))),
            "\ndP/drho [Pa m3/kg]:",
            formatstr.format(Decimal(str(self.dPdrho))),
        )

        self.mix.prt(digits=digits, **kwargs)

    def Update(self, **kwargs):
        """This function is meant to update all relevant parameters
        after the Construction method has been employed
        """

        # t0=time.time()
        # first update material parameters

        # t0mix=time.time()
        self.mix.Update(P=self.pres, T=self.temp, **kwargs)
        # tmix=time.time()
        self.dens = self.mix.dens

        # then update gradients
        if not self.mix.dPdrho == 0.0:
            self.dPdrho = self.mix.dPdrho

        self.dPdr = -G * self.mass * self.dens / self.radius ** 2

        try:
            if self.rhoType == "constant":
                self.drhodr = 0.0

            elif self.rhoType == "update":
                self.drhodr = self.dPdr / self.dPdrho

        # exit if the density derivative of the pressure is zero which should
        # not happen because the adiabatic gradient can then not be defined
        except ZeroDivisionError:
            print("ZeroDivisionError in Shell.Update(): dPdrho = 0")
            print("dPdr:", self.dPdr)
            print("dPdrho:", self.dPdrho)
            print("mass:", self.mass)
            print("dens:", self.dens)
            print("radius:", self.radius)
            print("mix:")
            self.mix.prt()
            sys.exit()

        self.dTdP = self.gamma * self.temp / (self.dens * self.dPdrho)
        self.dTdr = self.dPdr * self.dTdP
        self.dmdr = 4.0 * np.pi * self.radius ** 2 * self.dens

        # if isothermal set all temperature derivatives to zero
        if self.tempType == 0:
            self.dTdr = 0.0
            self.dTdP = 0.0

        self.gradients = [self.dPdr, self.dmdr, self.dTdr, self.drhodr]

        # if a unit detects a problem, the force_bisection call is here
        # passed to the next level (from mixture to shell)
        if self.mix.force_bisection:
            self.force_bisection = True

        # update gravity and escape velocity
        self.gravity = G * self.mass / self.radius ** 2
        self.v_esc = np.sqrt(2.0 * G * self.mass / self.radius)

        # update thermal velocity of water vapor at these conditions assuming
        # 3+3 degrees of freedom (vibrational modes neglected)
        self.v_th_H2O = np.sqrt(6 * kB * self.temp * NA / (mO + 2 * mH))

        # t=time.time()
        # self.t_update += t-t0
        # print ('\nshell.update =', t-t0)
        # print ('mix.update=', tmix-t0mix)

    def reset(self, **kwargs):
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

    def construct(self, dr=1.0, overconstruct=False, fortran=False, **kwargs):
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

            # call the integrate methode from eosfort to use 4th order RK scheme
            self.radius, newy = eosfort.integrate(
                r_start=self.radius,
                h=dr,
                y_in=params,
                gammag=self.gamma,
                fractions=self.mix.fractions,
                temptype=self.tempType,
                nmat=len(self.contents),
                ll=self.contents[0],
            )

            self.pres, self.mass, self.temp, self.dens = newy
            # print ('T/P/d after1=', self.temp, self.pres, self.dens)
            # use new pressure, temperature, radius and enclosed mass to update
            # all other parameters
            self.Update()
            # print ('->T/P/d after2=', self.temp, self.pres, self.dens)
            # compute indigenous mass of the isolated shell
            self.indigenous_mass = self.mass - old_mass

            # indicate that the integration has been performed for this shell
            self.status = "constructed"
            # self.Update(**kwargs)
            # if a unit detects a problem, the force_bisection call is here
            # passed to the next level (from mixture to shell)
            if self.mix.force_bisection:
                self.force_bisection = True

        # update gravity and escape velocity
        self.gravity = G * self.mass / self.radius ** 2
        self.v_esc = np.sqrt(2.0 * G * self.mass / self.radius)

        # update thermal velocity of water vapor at these conditions assuming
        # 3+3 degrees of freedom (vibrational modes neglected)
        self.v_th_H2O = np.sqrt(6 * kB * self.temp * NA / (mO + 2 * mH))


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
        tempType=0,
        Fe_number=[],
        X_H2O=[],
        rhoType="constant",
        mass_should=None,
        **kwargs,
    ):
        """Initiates a new Layer. Temperature, pressure and enclosed mass that
        are passed __init__ refer to the layer bottom properties and serve
        as initial conditions for the construction of the layer as a multitude
        of individual shells.
        """

        self.tempType = tempType
        self.rhoType = rhoType
        self.r_in = r_in
        self.r_out = r_out
        self.eps_r = eps_r
        self.radius = self.r_in

        self.Fe_number = Fe_number
        self.X_H2O = X_H2O

        # set initial conditions
        self.pres, self.mass, self.temp = P, m, T
        self.mass_should = mass_should
        self.indigenous_mass = 0.0  # mass contained within the layer
        self.contents = contents
        self.fractions = fractions
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

        # if no fractions have been specified, distribute the components evenly
        # over the mixture
        if len(fractions) == 0:
            print("setting up uniform fractions")
            frac = 1.0 / len(self.contents)
            self.fractions.append(frac)

        if len(Fe_number) == 0:
            print("setting Fe# = 0.0")
            self.Fe_number = [0.0 for i in range(len(self.contents))]

        if len(X_H2O) == 0:
            print("setting X_H2O = 0.0")
            self.X_H2O = [0.0 for i in range(len(self.contents))]

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
                Fe_number=self.Fe_number,
                X_H2O=self.X_H2O,
                rhoType=self.rhoType,
                **kwargs,
            )

        except TypeError:
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

    def construct(self, **kwargs):
        """The Construct method integrates the current shell from r to r+dr
        and subsequently updates all the corresponding parameters. To construct
        the entire layer, the construction method needs to be called several
        times until the desired layer conditions are met.
        """
        # t0=time.time()
        self.shells[-1].construct(dr=self.dr, **kwargs)
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
                **kwargs,
            )

            # shell.Update()
            # if the density is held constant in each integration step, the
            # gradients are computed with the old density during integration. If
            # the

            # in the end a new shell is generated and added to the layer, ready
            # for the next integration step
            self.shells.append(shell)

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

    def plot(
        self,
        pres_unit="GPa",
        temp_unit="K",
        dens_unit="kg/m3",
        mass_unit="m_earth",
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
            axes[i].plot(radius_list / 1000, plot_list[i])
            axes[i].set_ylabel(title_list[i])


class Planet:
    def __init__(
        self,
        P_center=1.0e12,
        T_center=300,
        R_seed=0.1,
        contents=[],
        tempType=0,
        layerType="homogeneous",
        fractions=[],
        P_surface_should=P_earth,
        majorConstraint="P_surface",
        layerConstraint="mass",
        add_atmosphere=False,
        eps_r=0.5,
        layermasses=[],
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
        n_refine=32,
        modelType=0,
        SiMg=0.0,
        **kwargs,
    ):

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
        self.echo = echo
        self.M_H2O_should = 0.0
        self.M_H2O_is = 0.0
        self.M_H2O_hidden_is = 0.0
        self.density_inversion = False
        self.match = True
        self.vapor_reached = False
        self.N_shells = 0
        self.test_variable = 0.0
        self.Mg_count = 0.0
        self.Mg_number_layer = Mg_number_layer
        self.additionally_refined = False
        self.n_refine = n_refine

        # defines the physical model that is invoked for the planet
        # e.g Brucite dissociation model without silicates: modelType = 0
        # " " " with hydrated Olivine: modelType = 1
        self.modelType = modelType

        self.layermasses = [laym for laym in layermasses]
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

        # if no fractions have been specified, distribute the components evenly
        # over the mixture
        if len(fractions) == 0:
            print("setting up uniform fractions")
            for c in range(len(self.contents)):
                frac = 1.0 / len(self.contents[c])
                self.fractions.append([frac for i in self.contents[c]])

        # if modelType = 1 compute material fractions in the mantle from
        # the given values for Mg# and Si#
        if self.modelType == 1:
            self.fractions[1] = [Material.xi_per(self.SiMg), Material.xi_ol(self.SiMg)]

            print("setting up mantle fractions for modelType 1")

        if len(Fe_number) == 0:
            print("setting Fe# = 0.0")
            for c in self.contents:
                self.Fe_number.append([0.0 for i in range(len(c))])

        if len(X_H2O) == 0:
            print("setting X_H2O = 0.0")
            for c in self.contents:
                self.X_H2O.append([0.0 for i in range(len(c))])

        if len(self.layermasses) != len(self.contents):
            print("\nWARNING: number of layers and layermasses does not match")
            print("given layermasses:", len(self.layermasses))
            print("given layers: ", len(self.contents))

        # print ('-> completed')

        # initalize core seed for integration
        # print ('generating seed...')
        self.seed_material = Material.Mixture(
            contents=self.contents[0], P=self.P_center, T=self.T_center
        )
        self.seed_material.Compute()
        self.rho_center = self.seed_material.dens
        self.M_seed = 4.0 / 3.0 * np.pi * self.rho_center * R_seed ** 3
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
            **kwargs,
        )

        self.layers.append(self.seed)
        self.status = "seed"

        self.M_surface_is = self.M_seed
        self.T_surface_is = T_center
        self.P_surface_is = P_center
        self.R_surface_is = R_seed

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
            "vapor_stop": self.vapor_stop,
            "Mg_number_layer": self.Mg_number_layer,
        }

        self.finals = {
            "P_surface_is": self.P_surface_is,
            "M_surface_is": self.M_surface_is,
            "T_surface_is": self.T_surface_is,
            "R_surface_is": self.R_surface_is,
            "Mg_number_is": self.Mg_number_is,
            "M_H2O_is": self.M_H2O_is,
            "M_H2O_should": self.M_H2O_should,
            "M_H2O_hidden_is": self.M_H2O_hidden_is,
            "density_inversion": self.density_inversion,
            "match": self.match,
            "vapor_reached": self.vapor_reached,
            "layer_properties": self.layer_properties,
            "N_shells": self.N_shells,
        }

    def prt(self, digits=5, **kwargs):
        print("=======================================")
        print("Planet properties:")
        print("\n--------------------------------")
        print("Composition:")
        print("--------------------------------")
        for i in range(len(self.layers)):
            print("\nlayer: ", i)
            for j in range(len(self.contents[i])):
                print(
                    round(self.fractions[i][j] * 100, 3),
                    "%",
                    material_list[self.contents[i][j]],
                )

            print(
                "mass [M_earth]: ",
                round(self.layers[i].indigenous_mass / m_earth, digits),
            )
            print("upper P [GPa]: ", round(self.layers[i].pres * 1.0e-9, digits))

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
                "\nM_surface_should [M_earth]:",
                round(self.M_surface_should / m_earth, digits),
                "\nM_H2O_should [M_earth]:",
                round(self.M_H2O_should, digits),
            )
            print("--------------------------------")

        except TypeError:
            print("WARNING: Type Error in Planet.prt()")

    def Update_finals(self, **kwargs):
        # gather final properties in a dictionary
        self.finals = {
            "P_surface_is": self.P_surface_is,
            "M_surface_is": self.M_surface_is / m_earth,
            "T_surface_is": self.T_surface_is,
            "R_surface_is": self.R_surface_is / r_earth,
            "Mg_number_is": self.Mg_number_is,
            "M_H2O_is": self.M_H2O_is,
            "M_H2O_should": self.M_H2O_should,
            "M_H2O_hidden_is": self.M_H2O_hidden_is,
            "density_inversion": self.density_inversion,
            "match": self.match,
            "vapor_reached": self.vapor_reached,
            "layer_properties": self.layer_properties,
            "N_shells": self.N_shells,
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
            "vapor_stop": self.vapor_stop,
            "Mg_number_layer": self.Mg_number_layer,
        }

    def Reset(self, **kwargs):
        """Resets all planetary properties to the initial state. This allows
        to use the exact same specifications for re-constructing the planet.
        """

        # delete old properties to release memory
        for lay in self.layers:
            del lay.shells
            del lay

        if self.M_surface_should == float("inf"):
            print("HERE infinity in Planet.reset() before")

        if self.initials["M_surface_should"] > 1.0e30:
            print("HERE large before Planet.reset() in Planet")
            print("large value=", self.initials["M_surface_should"])

        self.__init__(**self.initials)

        if self.initials["M_surface_should"] > 1.0e30:
            print("HERE large after Planet.reset() in Planet")
            print("large value=", self.initials["M_surface_should"])

        if self.M_surface_should == float("inf"):
            print("HERE infinity in Planet.reset() after")

    def CheckLayerTransition(self, transitionType="mass", **kwargs):
        """This methode is used to check in each integration step if a layer
        transition is reached by probing the according parameters. Layers can
        be defined by the <transitionType>. Default is 'mass' which probes
        the total layermasses which must be set in when initializing the planet
        """

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
            self.layers[-1].shells.pop(-1)
            self.layers[-1].shells[-1].reset()
            self.layers[-1].Update()
            self.layers[-1].dr = self.layers[-1].dr * 0.5
            self.layers[-1].bisec = True
            self.changebisec = False

        elif abs(param_should - param_is) / param_is < eps:
            self.majorComplete = True
            # the shell.Construct() methode automatically adds a new
            # shell to the layer after construction. If the surface
            # of the planet is already reached this shell will not be
            # needed for further integration -> remove it
            self.layers[-1].shells.pop(-1)
            self.layers[-1].Update()
            self.shellIteration = False
            self.layerIteration = False

        else:
            if self.changebisec:
                self.layers[-1].bisec = False

    def Construct(self, print_time=False, fortran=True, echo=False, **kwargs):
        """N sets the number of shells in the case that no surface constrain
        type has been specified. Otherwise N will be ignored.
        """
        self.echo = echo
        t0 = time.time()
        if self.status == "constructed":
            print(
                "This planet has already been constructed. Use the",
                "reconstruct() methode to reconstruct the planet from seed",
                "with different specifications if you wish.",
            )

        else:
            self.status = "constructing"
            # iterate over all layers
            lay = 0
            self.layerIterationCount = 0
            self.shellIterationCount = 0
            self.layerIteration = True

            while self.layerIteration:
                self.layerIterationCount += 1
                self.shellIteration = True
                self.changebisec = True

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

                self.layerComplete = False

                # if an empty layer is given, skip the integration
                if self.layermasses[lay] == 0.0:
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
                        # t0i = time.time()
                        self.layers[-1].construct()

                        # ti = time.time()
                        # tconstruct += (ti-t0i)

                        # check if brucite dissociation is reached and if so
                        # re-integrate with smaller integration step to increase
                        # dissociation resolution
                        last_shell = self.layers[-1].shells[-3]
                        shell = self.layers[-1].shells[-2]

                        if (
                            self.layers[-1].contents == [11]
                            and shell.mix.mix[0].phase != last_shell.mix.mix[0].phase
                        ):
                            if self.additionally_refined:
                                pass

                            else:
                                self.additionally_refined = True
                                self.layers[-1].shells.pop(-1)
                                self.layers[-1].shells[-1].reset()
                                self.layers[-1].Update()
                                new_dr = self.layers[-1].dr / self.n_refine

                                for i in range(self.n_refine - 1):

                                    self.layers[-1].dr = new_dr
                                    self.layers[-1].construct()

                                    # compute water contents as deltas so that it can be
                                    # removed again if a shell needs to be re-integrated
                                    delta_M_H2O_should = 0.0
                                    delta_M_H2O_hidden_is = 0.0

                                    # Compute the amount of water that is released by the
                                    # current shell if contains dissociated Mg(OH)2

                                    if (
                                        self.layers[lay].shells[-2].mix.mix[0].phase
                                        == 2
                                    ):

                                        delta_M_H2O_should = (
                                            self.layers[lay].shells[-2].indigenous_mass
                                            / (mMg + mO)
                                            * (mO + mH * 2.0)
                                            / m_earth
                                        )

                                    else:
                                        delta_M_H2O_hidden_is = (
                                            self.layers[lay].shells[-2].indigenous_mass
                                            / (mMg + 2.0 * mO + 2.0 * mH)
                                            * (mO + mH * 2.0)
                                            / m_earth
                                        )

                                    self.M_H2O_should += delta_M_H2O_should
                                    self.M_H2O_hidden_is += delta_M_H2O_hidden_is

                                    # compute amount of magnesium in current shell and
                                    # update the total Mg# accordingly
                                    delta_Mg = 0.0
                                    sh = self.layers[lay].shells[-2]
                                    if sh.mix.mix[0].phase == 2:
                                        delta_Mg = sh.indigenous_mass / (mMg + mO) * mMg

                                    else:
                                        delta_Mg = (
                                            sh.indigenous_mass
                                            / (mMg + 2 * mO + 2 * mH)
                                            * mMg
                                        )

                                    self.Mg_count += delta_Mg
                                    self.Mg_number_is = (
                                        self.Mg_count
                                        / mMg
                                        / (
                                            self.layermasses[0] * m_earth / mFe
                                            + self.Mg_count / mMg
                                        )
                                    )

                                    # if Mg# is exceeded, the mantle is too big
                                    # and the last shell is re-integrated with
                                    # smaller integration step. The water content
                                    # and Mg content must therefore be reset to
                                    # the previous integration step
                                    if (
                                        self.Mg_number_is - self.Mg_number_should
                                    ) / self.Mg_number_should > eps_Mg_number:
                                        # undo the Mg and H2O counting
                                        self.Mg_count -= delta_Mg
                                        self.M_H2O_should -= delta_M_H2O_should
                                        self.M_H2O_hidden_is -= delta_M_H2O_hidden_is
                                        self.layers[-1].shells.pop(-1)

                                        # reset the just integrated shell to re-integrate
                                        # it using a smaller dr
                                        self.layers[-1].shells[-1].reset()
                                        self.layers[-1].Update()
                                        break

                                # the last construction call shall be regualar
                                # in order to not compute the water mass twice
                                self.layers[-1].dr = new_dr
                                self.layers[-1].construct()

                    except ValueError:
                        print(self.initials)
                        self.prt()
                        print(self.layers[-1].pres, self.layers[-1].temp)
                        sys.exit()

                    except TypeError:
                        self.shellIteration = False
                        self.layerIteration = False

                    """
                    except TypeError:
                        print ('WARNING: an error occured while constructing'+\
                               'the layer')
                        print ('T=', self.layers[-1].temp, 'K')
                        print ('P=', self.layers[-1].pres*1.0e-5, 'bar')
                        print ('rho=', self.layers[-1].dens)
                        
                        self.layerIteration = False
                        self.shellIteration = False
                    """
                    newdens = self.layers[-1].dens

                    # check for spurious numerical behaviour by making sure
                    # that the density decreases towards the surface
                    if (newdens - olddens) / olddens > 1.0e-3:
                        #                          print ('WARNING: Density inversion encountered in'+\
                        #                                 ' layer', lay,', at:',
                        #                                 round(olddens, 4),'and', round(newdens, 4))
                        #                          print ('T/P =', self.layers[-1].temp,
                        #                                 self.layers[-1].pres)
                        # =============================================================================

                        # use constant density for shell integration for this
                        # shell if density inversion is encountered
                        self.layers[-1].shells.pop(-1)
                        self.layers[-1].shells[-1].reset()
                        self.layers[-1].Update()

                        # perform next integration step
                        try:
                            self.layers[-1].shells[-1].rhoType = "constant"
                            self.layers[-1].construct()

                        except TypeError:
                            print(
                                """WARNING: an error occured while constructing
                                   the layer"""
                            )
                            print("T=", self.layers[-1].temp, "K")
                            print("P=", self.layers[-1].pres * 1.0e-5, "bar")
                            print("rho=", self.layers[-1].dens)

                            self.layerIteration = False
                            self.shellIteration = False

                        newdens = self.layers[-1].dens

                        # sys.exit()
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

                    if self.layers[-1].temp < T_zero:
                        print(
                            "Zero temperature reached at",
                            self.layers[-1].temp,
                            "K after",
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

                    # Check for water vapor in pure H2O layer
                    if self.contents[lay] == [0] and self.vapor_stop:
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

                    # use surface pressure to define surface of planet
                    if self.majorConstraint == "P_surface":
                        # perform next integration step
                        # first check major constraint condition
                        # if pressure falls below specified surface pressure,
                        # reset and re-integrate old shell
                        pres = self.layers[-1].pres
                        self.major_bisection(
                            pres,
                            self.P_surface_should,
                            eps=eps_Psurf,
                            direction=[-1, -1],
                        )

                        # compute water contents as deltas so that it can be
                        # removed again if a shell needs to be re-integrated
                        delta_M_H2O_should = 0.0
                        delta_M_H2O_hidden_is = 0.0

                        # Compute the amount of water that is released by the
                        # current shell if contains dissociated Mg(OH)2
                        if self.brucite_dissoc:
                            if self.layers[lay].contents == [11]:

                                if self.layers[lay].shells[-2].mix.mix[0].phase == 2:

                                    delta_M_H2O_should = (
                                        self.layers[lay].shells[-2].indigenous_mass
                                        / (mMg + mO)
                                        * (mO + mH * 2.0)
                                        / m_earth
                                    )

                                else:
                                    delta_M_H2O_hidden_is = (
                                        self.layers[lay].shells[-2].indigenous_mass
                                        / (mMg + 2.0 * mO + 2.0 * mH)
                                        * (mO + mH * 2.0)
                                        / m_earth
                                    )

                                self.M_H2O_should += delta_M_H2O_should
                                self.M_H2O_hidden_is += delta_M_H2O_hidden_is

                    # use total mass to define surface of planet
                    elif self.majorConstraint == "M_surface":
                        mass = self.layers[-1].mass

                        self.major_bisection(
                            mass, self.M_surface_should, eps=eps_Mtot, direction=[1, -1]
                        )

                    # second: check layer constraint condition to probe transition
                    # to the next layer
                    if (
                        self.differentiated
                        and not self.vapor_reached
                        and not self.majorComplete
                    ):

                        # compute amount of magnesium in current shell and
                        # update the total Mg# accordingly
                        delta_Mg = 0.0

                        if (
                            self.layerConstraint == "mass"
                            and self.Mg_number_layer
                            and lay == 1
                        ):
                            if lay == 1:
                                sh = self.layers[lay].shells[-2]
                                if sh.mix.mix[0].phase == 2:
                                    delta_Mg = sh.indigenous_mass / (mMg + mO) * mMg

                                else:
                                    delta_Mg = (
                                        sh.indigenous_mass
                                        / (mMg + 2 * mO + 2 * mH)
                                        * mMg
                                    )

                                self.Mg_count += delta_Mg
                                self.Mg_number_is = (
                                    self.Mg_count
                                    / mMg
                                    / (
                                        self.layermasses[0] * m_earth / mFe
                                        + self.Mg_count / mMg
                                    )
                                )

                                # too much Mg
                                # note: if accuracy for Mg_number is > 1.0e-4
                                # i have found, that this error is enough for
                                # the corresponding fluctuation in the mantle
                                # mass to cause the total planetary mass to be
                                # too small and therefore it adds a small ocean
                                # on top even if it shouldn't. For acc = 1.0e-4
                                # this effect is still non-zero but negligible
                                if (
                                    self.Mg_number_is - self.Mg_number_should
                                ) / self.Mg_number_should > eps_Mg_number:
                                    # undo the Mg and H2O counting
                                    self.Mg_count -= delta_Mg
                                    self.M_H2O_should -= delta_M_H2O_should
                                    self.M_H2O_hidden_is -= delta_M_H2O_hidden_is

                                    self.layers[-1].shells.pop(-1)
                                    # reset the just integrated shell to re-integrate
                                    # it using a smaller dr
                                    self.layers[-1].shells[-1].reset()
                                    self.layers[-1].Update()
                                    self.layers[-1].dr = self.layers[-1].dr * 0.5
                                    self.layers[-1].bisec = True
                                    self.changebisec = False

                                # check if Mg number has approached
                                # the desired value
                                elif (
                                    abs(self.Mg_number_is - self.Mg_number_should)
                                    / self.Mg_number_should
                                    < eps_Mg_number
                                ):
                                    # update layermass for ocean
                                    self.layerComplete = True
                                    self.shellIteration = False
                                    self.changebisec = True

                                    # if this was the last layer, then
                                    # terminate internal structure integration
                                    if lay == len(self.contents) - 1:
                                        if self.echo:
                                            print("terminating layer iteration")

                                        self.layerIteration = False

                                    # initiate next layer
                                    else:
                                        radius = self.layers[-1].radius
                                        mass = self.layers[-1].mass
                                        pres = self.layers[-1].pres
                                        temp = self.layers[-1].temp

                                        nextlayer = Layer(
                                            m=mass,
                                            T=temp,
                                            P=pres,
                                            tempType=self.tempType,
                                            r_in=radius,
                                            contents=self.contents[lay + 1],
                                            fractions=self.fractions[lay + 1],
                                            eps_r=self.eps_r,
                                            X_H2O=self.X_H2O[lay + 1],
                                            Fe_number=self.Fe_number[lay + 1],
                                            rhoType=self.rhoType,
                                        )

                                        lay += 1
                                        self.layers.append(nextlayer)

                        # note to future Oliver: i am very sorry for this ugly
                        # crap, but i'm kinda in a hurry (stupid epsc!)
                        # best, present Oliver (3.9.2019)
                        checkLayer = False
                        if not self.Mg_number_layer:
                            checkLayer = True

                        elif self.Mg_number_layer and not lay == 1:
                            checkLayer = True

                        if (
                            self.layerConstraint == "mass"
                            and checkLayer
                            and not self.layerComplete
                        ):
                            laymass = self.layers[lay].indigenous_mass
                            # sometimes, for whatever reason, a shell is not
                            # constructed in the layer.construct method. In this
                            # case the indigenous mass will be zero and laymass
                            # can be zero. To avoid spurious behaviour manually
                            # construct the shell in question here
                            if laymass == 0.0:
                                # print ('WARNING: zero layermass in layer', lay)
                                # print (self.initials)
                                # for sh in self.layers[lay].shells:
                                #   sh.prt()

                                self.layers[lay].shells[-2].construct()
                                self.layers[lay].Update()
                                laymass = self.layers[lay].indigenous_mass

                            # mass exceeded, reset and re-integrate old shell
                            # note that layermasses are in units of earth masses
                            if (
                                self.layermasses[lay] * m_earth - laymass
                            ) / laymass < -eps_layer:
                                """
                                #Compute amount of dissociated or incorporated
                                #water in current shell. This then has to be
                                #substracted from the total amount of water
                                #in the mantle or the ocean because the shell
                                #is beeing reset and re-integrated
                                if self.brucite_dissoc:
                                    if self.layers[lay].contents == [11]:
                                        
                                        if self.layers[lay].shells[-2].mix.\
                                        mix[0].phase \
                                        == 2:
                                            self.M_H2O_should -= \
                                            self.layers[lay].shells[-2].\
                                            indigenous_mass/\
                                            (mMg+mO)*(mO+mH*2.)/m_earth
                                            
                                        else:
                                            self.M_H2O_hidden_is -= \
                                            self.layers[lay].shells[-2].\
                                            indigenous_mass/\
                                            (mMg+2.*mO+2.*mH)*(mO+mH*2.)/m_earth
                                """
                                # undo the Mg and H2O counting
                                self.Mg_count -= delta_Mg
                                self.M_H2O_should -= delta_M_H2O_should
                                self.M_H2O_hidden_is -= delta_M_H2O_hidden_is

                                # remove bare shell that was prepared for next
                                # integration step
                                self.layers[-1].shells.pop(-1)
                                # reset the just integrated shell to re-integrate
                                # it using a smaller dr
                                self.layers[-1].shells[-1].reset()
                                self.layers[-1].Update()
                                self.layers[-1].dr = self.layers[-1].dr * 0.5
                                self.layers[-1].bisec = True
                                self.changebisec = False

                            # check if current layer mass has approached
                            # the desired value for the total mass to
                            # specified accuracy
                            elif (
                                abs(self.layermasses[lay] * m_earth - laymass) / laymass
                                < eps_layer
                            ):

                                self.layerComplete = True
                                self.shellIteration = False
                                self.changebisec = True

                                if (
                                    len(self.layers) == len(self.contents)
                                    and self.layerComplete
                                ):
                                    sys.stdout.write("\033[1A")

                                # if this was the last layer, then
                                # terminate internal structure integration
                                if lay == len(self.contents) - 1:
                                    if self.echo:
                                        print("terminating layer iteration")

                                    self.layerIteration = False

                                # initiate next layer
                                else:

                                    radius = self.layers[-1].radius
                                    mass = self.layers[-1].mass
                                    pres = self.layers[-1].pres
                                    temp = self.layers[-1].temp

                                    next_eps_r = self.eps_r / 2.0

                                    nextlayer = Layer(
                                        m=mass,
                                        T=temp,
                                        P=pres,
                                        tempType=self.tempType,
                                        r_in=radius,
                                        contents=self.contents[lay + 1],
                                        fractions=self.fractions[lay + 1],
                                        eps_r=next_eps_r,
                                        X_H2O=self.X_H2O[lay + 1],
                                        Fe_number=self.Fe_number[lay + 1],
                                        rhoType=self.rhoType,
                                    )

                                    lay += 1
                                    self.layers.append(nextlayer)

                            # if the surface condition is neither met nor exceeded
                            # (re-)deactivate bisection proceedure and perform
                            # normal shell integration in the next iteration
                            else:
                                # if surface pressure bisection is ongoing, don't allow
                                # to change bisection to False
                                if self.changebisec:

                                    self.layers[-1].bisec = False
                                else:
                                    pass

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
                                    lay,
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
                                lay - 1,
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

        t = time.time()

        for layer in self.layers:
            if layer.contents == [0]:
                self.M_H2O_is = layer.indigenous_mass / m_earth

        if not self.Mg_number_layer:
            self.Mg_number_is = Material.Mg_number(
                contents=self.contents, layers=self.layers, X_H2O=self.X_H2O
            )

        if print_time:
            print("total elapsed time=", round(t - t0, 4), "s")

        if self.echo:
            self.prt()

        self.layer_properties = [
            {"P_out": None, "T_out": None, "rho_out": None, "indigenous_mass": 0.0}
            for i in range(len(self.layers))
        ]

        for i in range(len(self.layers)):
            self.layer_properties[i]["P_out"] = self.layers[i].pres
            self.layer_properties[i]["T_out"] = self.layers[i].temp
            self.layer_properties[i]["rho_out"] = self.layers[i].dens
            self.layer_properties[i]["indigenous_mass"] = self.layers[i].indigenous_mass

        self.Update_finals()

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
        save=False,
        filename="planet_structure",
        suffix="png",
        path="./",
        **kwargs,
    ):

        if axis == []:
            fig, axis = plt.subplots(2, 3, sharex=True)

        title_list = [
            r"$Pressure \ [GPa]$",
            r"$Mass \ [M_\oplus]$",
            r"$Density \ [kg/m^3]$",
            r"$Temperature  \ [K]$",
            r"$Gravity \ [m/s^2]$",
            r"$v_{esc} \ [km/s]$",
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
            [[shell.dens for shell in layer.shells] for layer in self.layers],
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

        plt.subplots_adjust(hspace=0.5, wspace=0.5)
        ax = [axis[0][0], axis[0][1], axis[1][0], axis[1][1], axis[0][2], axis[1][2]]

        ax[2].set_xlabel(r"$Radius  \ [R_\oplus]$")
        ax[3].set_xlabel(r"$Radius \ [R_\oplus]$")
        ax[5].set_xlabel(r"$Radius \ [R_\oplus]$")
        # for axx in ax:
        #   axx.grid(which='both', axis='both')

        for i in range(len(ax)):
            ax[i].plot(radius_list, plot_list[i], color=param_colors[i])
            ax[i].tick_params(right=True, top=True, direction="in")
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
