#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 16:59:25 2019

@author: oshah
"""

import matplotlib.ticker
from matplotlib import pyplot as plt
from pics.physicalparams import (
    T0_list,
    K0_list,
    K0prime_list,
    rho0_list,
    aP_list,
    molar_mass_list,
    Rgas,
    mH,
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
    gas_molar_mass_list,
    gas_EOS_type,
    T_zero,
    gas_material_list,
    gas_gamma_list,
    P_zero,
    rho_zero,
)

from pics.runparams import (
    eps_Psurf,
    eps_Pspace,
    eps_Mtot,
    param_colors,
    temp_round,
    plot_units,
)
import numpy as np

# import eos
from pics.utils import functionTools as ftool
import sys
from tabulate import tabulate
from decimal import Decimal


def P_barom(mu=None, r=None, g=None, R_planet=None, T=None, P0=None):
    pres = P0 * np.exp(-mu * g / (Rgas * T) * (r - R_planet))
    return pres


def rho_barom(mu=None, r=None, g=None, R_planet=None, T=None, rho0=None):
    return rho0 * np.exp(-mu * g / (Rgas * T) * (r - R_planet))


def m_barom(mu=None, r=None, g=None, R_planet=None, T=None, rho0=None, M0=None):
    H = mu * g / (Rgas * T)
    return (
        -4
        * np.pi
        * rho0
        / H**3
        * (
            (H**2 * r**2 + 2 * H * r + 2) * np.exp(-H * (r - R_planet))
            - H**2 * R_planet**2
            - 2 * H * R_planet
            - 2
        )
        + M0
    )


class Unit:
    def __init__(self, ll=None, P=1.0e5, T=300, d=0.0, **kwargs):
        self.which = ll
        self.pres = P
        self.temp = T
        self.dens = d
        self.dPdrho = 0.0
        self.fraction = 1
        self.material = gas_material_list[ll]
        self.whichEOS = gas_EOS_type[ll]
        self.phase = "gas"
        # specify here additional gas properties if you can't resist

    def prt(self, digits=20, level=0):
        """prints out selection of basic properties of the material object for
        more convenient handling
        """
        which = self.whichEOS
        print(
            "material:",
            self.material,
            "\nfraction:",
            self.fraction,
            "\nEOS:",
            EOSstring_list[which],
            "\nphase:",
            self.phase,
            "\n-----",
            "\ndens  [kg/m³]:",
            round(self.dens, digits),
            "\npres  [GPa]  :",
            round(self.pres * 1.0e-9, digits),
            "\ntemp  [K]    :",
            round(self.temp, digits),
            "\ndP/drho [GPa m3/kg]:",
            round(self.dPdrho * 1.0e-9, digits),
        )

    def Compute(self, what=None, d=None, P=None, T=None, **kwargs):
        """Here the material paremeters are computed. Normally values for P
        and T are passed and the density and pressure derivative are computed.
        Note: no argument for <what> should be passed if possible as only the
        passed parameter will be updated. This option is mainly for testing
        and debugging purposes"""

        # extract current gas paramterers
        if d == None:
            dens = self.dens
        else:
            dens = d
        if P == None:
            pres = self.pres
        else:
            pres = P
        if T == None:
            temp = self.temp
        else:
            temp = T

        dPdrho = self.dPdrho

        if what == "dens":
            dens = eos.rho_IGL(P=self.pres, T=self.temp, ll=self.which)

        elif what == "pres":
            pres = eos.P_IGL(d=self.dens, T=self.temp, ll=self.which)

        elif what == "dPdrho":
            dPdrho = eos.dPdrho_IGL(T=self.temp, ll=self.which)

        elif what == None:
            dens = eos.rho_IGL(P=self.pres, T=self.temp, ll=self.which)
            dPdrho = eos.dPdrho_IGL(T=self.temp, ll=self.which)

        # update all gas parameters
        self.pres = pres
        self.dens = dens
        self.temp = temp
        self.dPdrho = dPdrho


class Mixture:
    """This is a subclass of the 'Material' routine which accounts for
    multiple component material cells. Each component in the material is
    represented by a 'Unit' instance and associated to a certain fraction
    given of a component that is present in the unit cell as input to the class.
    If no fractions are specified, the unit cell is composed of the same amount
    of each material.

    Possible inputs:

        contents (kwarg, list, type = int, default = []):
            a python list containing the material id's (integers) that point to
            the associated EOS that has to be used to treat this material.

        fractions (kwarg, list, type = list, default = []):
            a python list containing the fractions for the materials specified
            in the contents list. Note, that there has to be one fraction for
            each material passed to the contents list and they have to add up
            to unity. If one of those conditions is not met, an error message
            will be promted.

        T (kwarg, float, default = 300 K):
            Temperature of the mixture in the unit cell in K

        P (kwarg, float, default = 1.0e5 Pa):
            Pressure of the mixture in the unit cell in Pa
    """

    def __init__(self, contents=[], fractions=[], T=300, P=1.0e5, d=None, **kwargs):
        self.contents = contents
        self.fractions = fractions
        self.materials = [gas_material_list[ll] for ll in self.contents]
        self.mix = []
        self.pres = P
        self.temp = T
        self.dens = d
        self.densities = None
        self.dPdrho = None
        # if no fractions have been specified, distribute the components evenly
        # over the mixture
        if len(fractions) == 0 and not len(contents) == 0:
            # print ('setting up uniform fractions')
            frac = 1.0 / len(self.contents)
            self.fractions = [frac for i in self.contents]

    def prt(self, individual=False, digits=20, **kwargs):
        print("-----------------------------------")
        print("Material properties of the mixture:")
        print("\ncontents:", self.materials, "\nfractions:", self.fractions)
        try:
            print(
                "\ndens  [kg/m³]:",
                round(self.dens, digits),
                "\npres  [GPa]  :",
                round(self.pres * 1.0e-9, digits),
                "\ntemp  [K]    :",
                round(self.temp, digits),
                "\ndP/drho [GPa m3/kg]:",
                round(self.dPdrho * 1.0e-9, digits),
            )
        except TypeError:
            print(
                "\ndens  [kg/m³]:",
                round(self.dens, digits),
                "\npres  [GPa]  :",
                round(self.pres * 1.0e-9, digits),
                "\ntemp  [K]    :",
                round(self.temp, digits),
                "\ndP/drho [GPa m3/kg]:",
                self.dPdrho,
            )
        print("-----------------------------------")

        if individual:
            print("\nThe Mixture contains the following components:")
            for m in range(len(self.mix)):
                print("\n#", m + 1)
                self.mix[m].prt(digits=digits, **kwargs)

    def Compute(self, **kwargs):
        """This method initializes a new material instance for each of the
        compontents in the mixture specified by the parameter 'contents' at
        given temperature and pressure and computes the density for each
        component individually and the resulting overall density of the mixture
        according to a specified mixing rule.
        """
        # in case a pressure and/or temperature value has been passed to the
        # construct method, update these properties for the unit cell
        try:
            self.pres = kwargs["P"]
        except KeyError:
            pass

        try:
            self.temp = kwargs["T"]
        except KeyError:
            pass

        # construct each component as a Unit object
        self.mix = [Unit(i, P=self.pres, T=self.temp, **kwargs) for i in self.contents]

        # update material fractions
        for m in range(len(self.mix)):
            self.mix[m].fraction = self.fractions[m]

        # compute density for different components individually
        for mm in self.mix:
            mm.Compute(**kwargs)

        # generate density list
        self.densities = [mm.dens for mm in self.mix]

        # compute mean density of the mixture using an adequate mixing rule
        self.dens = eos.mix1(self.densities, self.fractions)
        self.dPdrho = 0.0

        for m in range(len(self.materials)):
            self.dPdrho += self.mix[m].dPdrho * self.fractions[m]

    def Update(self, **kwargs):
        """This method updates the material instance of each component of the
        mixture individually and then computes the new mean density in the cell
        without re-initializing any component objects. If no new pressure or
        temperature is passed, nothing will be done.
        """
        # in case a pressure and/or temperature value has been passed to the
        # construct method, update these properties for the unit cell
        try:
            self.pres = kwargs["P"]
            found_pres = True

        except KeyError:
            found_pres = False

        try:
            self.temp = kwargs["T"]
            found_temp = True

        except KeyError:
            found_temp = False

        # update only if new pressure and|or temperature has been specified
        if found_pres == False and found_temp == False:
            pass

        else:
            # update properties of individual components
            for m in range(len(self.mix)):
                self.mix[m].pres = self.pres
                self.mix[m].temp = self.temp

            for mm in self.mix:
                mm.Compute(**kwargs)

            # update density list
            self.densities = [mm.dens for mm in self.mix]

            # update mean density of the mixture using an adequate mixing rule
            self.dens = eos.mix1(self.densities, self.fractions)

            # update mean pressure derivative using an adequate mixing rule
            self.dPdrho = 0.0
            for m in range(len(self.materials)):
                self.dPdrho += self.mix[m].dPdrho * self.fractions[m]


class Shell:
    def __init__(
        self,
        radius=0.1,
        T=None,
        m=None,
        P=None,
        d=None,
        tempType="adiabatic",
        status="bare",
        contents=[],
        fractions=[],
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
        self.contents = contents
        self.indigenous_mass = 0.0
        self.fractions = fractions
        self.v_esc = np.sqrt(2 * G * self.mass / (self.radius))
        self.v0 = []

        self.dPdr = None
        self.mix = Mixture(
            T=self.temp,
            P=self.pres,
            contents=self.contents,
            fractions=self.fractions,
            **kwargs,
        )
        self.mix.Compute(**kwargs)

        # compute mean particle velocity for each species
        for s in self.contents:
            self.v0.append(np.sqrt(2 * kB * self.temp / gas_molar_mass_list[s] * NA))

        self.dens = self.mix.dens
        self.dPdrho = self.mix.dPdrho
        self.dPdr = None
        self.dTdP = None
        self.drhodr = None
        self.dTdr = None
        self.dmdr = None
        self.gamma = sum([gas_gamma_list[l] for l in self.contents]) / len(
            [gas_gamma_list[l] for l in self.contents]
        )
        self.status = status
        self.layer = -1

        self.tempType = tempType

        self.dPdr = -G * self.mass * self.dens / self.radius**2
        self.drhodr = self.dPdr / self.dPdrho
        self.dTdP = self.gamma * self.temp / (self.dens * self.dPdrho)
        self.dTdr = self.dPdr * self.dTdP
        self.dmdr = 4.0 * np.pi * self.radius**2 * self.dens

        if self.tempType == "isothermal":
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
        }

    def prt(self, **kwargs):
        print(
            "Shell at r =", self.radius * 1.0e-3, "km, enclosing m =", self.mass, "kg"
        )
        print("\nstatus:", self.status)
        print("\ntempType:", self.tempType)
        print("\nThe gradients are:")
        print(
            "\ndP/dr:",
            self.dPdr,
            "\ndm/dr :",
            self.dmdr,
            "\ndT/dr :",
            self.dTdr,
            "\ndrho/dr:",
            self.drhodr,
            "\ndT/dP:",
            self.dTdP,
            "\ndP/drho :",
            self.dPdrho,
        )

        self.mix.prt(**kwargs)

    def Update(self, **kwargs):
        """This function is meant to update all relevant parameters
        after the Construction method has been employed
        """
        # first update gradients
        if not self.mix.dPdrho == 0.0:
            self.dPdrho = self.mix.dPdrho
        self.dPdr = -G * self.mass * self.dens / self.radius**2

        try:
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
        self.dmdr = 4.0 * np.pi * self.radius**2 * self.dens

        # if isothermal set all temperature derivatives to zero
        if self.tempType == "isothermal":
            self.dTdr = 0.0
            self.dTdP = 0.0

        self.gradients = [self.dPdr, self.dmdr, self.dTdr, self.drhodr]

        # update material parameters
        self.mix.Update(P=self.pres, T=self.temp, **kwargs)

        # compute escape velocity
        self.v_esc = np.sqrt(2 * G * self.mass / (self.radius))

        # compute mean particle velocity for each species
        self.v0 = []
        for s in self.contents:
            self.v0.append(np.sqrt(2 * kB * self.temp / gas_molar_mass_list[s] * NA))

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
        ) = (
            self.initials["radius"],
            self.initials["P"],
            self.initials["m"],
            self.initials["T"],
            self.initials["d"],
            self.initials["dPdrho"],
            self.initials["mix"],
        )

        self.indigenous_mass = 0.0
        self.mix = Mixture(
            T=self.temp,
            P=self.pres,
            contents=self.contents,
            fractions=self.fractions,
            **kwargs,
        )
        self.mix.Compute(**kwargs)
        Shell.Update(self, **kwargs)
        self.status = "bare"

    def Construct(self, dr=1.0, overconstruct=False, **kwargs):
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
            # call the integrate methode from ftool to use 4th order RK
            self.radius, newy = ftool.integrate(
                dydx=gradients,
                start=self.radius,
                whicharg="r",
                y=params,
                h=dr,
                dPdrho=self.dPdrho,
                tempType=self.tempType,
                identity="construct shell",
                gamma=self.gamma,
                materials=self.mix.contents,
                fractions=self.mix.fractions,
                **kwargs,
            )

            self.pres, self.mass, self.temp, self.dens = newy

            # compute indigenous mass of the isolated shell
            self.indigenous_mass = self.mass - old_mass

            # indicate that the integration has been performed for this shell
            self.status = "constructed"
            Shell.Update(self, **kwargs)


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
        tempType="adiabatic",
        **kwargs,
    ):
        self.tempType = tempType
        self.indigenous_mass = 0.0
        self.r_in = r_in
        self.r_out = r_out
        self.eps_r = eps_r  # integration step refinement parameter
        self.radius = self.r_in
        # set initial conditions
        self.pres, self.mass, self.temp = P, m, T
        self.contents = contents
        self.fractions = fractions
        self.materials = [gas_material_list[ll] for ll in self.contents]
        # integration step which will constantly be updated as the structure
        # is computed
        self.dr = None
        self.radius_list = None
        self.plot_list = None
        self.bisec = False
        # this list will be filled with all the shells in the layer
        self.shells = []
        # generate initial shell
        shell = Shell(
            T=self.temp,
            P=self.pres,
            contents=self.contents,
            fractions=self.fractions,
            radius=self.radius,
            m=self.mass,
            tempType=self.tempType,
            **kwargs,
        )

        self.fractions = shell.mix.fractions
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
        self.pres, self.mass, self.temp, self.dens = (
            shell.pres,
            shell.mass,
            shell.temp,
            shell.dens,
        )

        self.shells[-1].Update(**kwargs)

        self.params = [self.pres, self.mass, self.temp, self.dens]
        self.radius = self.radius + self.dr
        self.gradients = shell.gradients

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

    def Construct(self, **kwargs):
        """The Construct method integrates the current shell from r to r+dr
        and subsequently updates all the corresponding parameters
        """
        self.shells[-1].Construct(dr=self.dr, **kwargs)
        self.indigenous_mass = sum([s.indigenous_mass for s in self.shells])

        # update all the gradients before the next integration
        Layer.Update(self)

        # update outer radius of the atmoshpere
        self.r_out = self.radius
        """
        #if surface pressure is exceed, reset and re-integrate old shell
        if abs(self.pres - P_earth)/self.pres > .001 and self.pres < P_earth:
            self.dr = .5*self.dr
            self.shells[-1].Reset()
            self.bisec=True
        
        #create new shell at the top of the old shell
        else:            
            shell = Shell(T=self.temp, P=self.pres, contents=self.contents,\
                          fractions =self.fractions, radius=self.radius,\
                          m=self.mass, tempType = self.tempType, **kwargs)
            
            self.shells.append(shell)"""

        shell = Shell(
            T=self.temp,
            P=self.pres,
            contents=self.contents,
            fractions=self.fractions,
            radius=self.radius,
            m=self.mass,
            tempType=self.tempType,
            **kwargs,
        )

        self.shells.append(shell)

    def plot(
        self,
        pres_unit="bar",
        temp_unit="K",
        dens_unit="kg/m3",
        mass_unit="m_earth",
        rad_unit="r_earth",
        **kwargs,
    ):

        title_list = [
            "Pressure [" + pres_unit + "]",
            "Mass [" + mass_unit + "]",
            "Density [" + dens_unit + "]",
            "Temperature [" + temp_unit + "]",
        ]

        radius_list = np.array(
            [shell.radius / plot_units[rad_unit] for shell in self.shells]
        )
        plot_list = np.array(
            [
                [shell.pres / plot_units[pres_unit] for shell in self.shells],
                [shell.mass / plot_units[mass_unit] for shell in self.shells],
                [shell.dens / plot_units[dens_unit] for shell in self.shells],
                [shell.temp / plot_units[temp_unit] for shell in self.shells],
            ]
        )

        fig, ax = plt.subplots(2, 2, sharex=True)
        plt.subplots_adjust(hspace=0.5, wspace=0.5)

        self.axes = [ax[0][0], ax[0][1], ax[1][0], ax[1][1]]
        self.axes[2].set_xlabel(r"$R/R_\oplus$")
        self.axes[3].set_xlabel(r"$R/R_\oplus$")

        for i in range(len(self.axes)):
            self.axes[i].plot(radius_list, plot_list[i], color=param_colors[i])
            self.axes[i].set_ylabel(title_list[i])

            self.axes[i].tick_params(direction="in", right=True, top=True)


class Atmosphere:
    def __init__(
        self,
        r_in=1,
        r_out=None,
        T=300,
        P=P_earth,
        contents=[],
        fractions=[],
        m=1,
        eps_r=0.5,
        tempType="isothermal",
        P_space=P_zero,
        M_atmos=0.000001,
        majorConstraint="P_space",
        T_space=T_zero,
        rho_space=rho_zero,
        **kwargs,
    ):

        self.R_surface = r_in * r_earth
        self.P_surface = P
        self.T_surface = T
        self.height = 0.0
        self.temp = T
        self.pres = P
        self.P_space_is = P
        self.M_surface = m * m_earth
        self.M_atmos = M_atmos * m_earth
        self.contents = contents
        self.indigenous_mass = 0.0
        self.materials = [gas_material_list[l] for l in contents]

        # if no fractions have been specified, distribute the components evenly
        # over the mixture
        if len(fractions) == 0 and not len(contents) == 0:
            # print ('setting up uniform fractions')
            frac = 1.0 / len(self.contents)
            self.fractions = [frac for i in self.contents]

        else:
            self.fractions = fractions

        self.tempType = tempType
        self.eps_r = eps_r
        self.P_space_should = P_space
        self.T_space_should = None
        self.T_space_is = T
        self.majorConstraint = majorConstraint

        atm = Layer(
            r_in=self.R_surface,
            m=self.M_surface,
            P=self.P_surface,
            T=self.T_surface,
            contents=self.contents,
            fractions=self.fractions,
            tempType=self.tempType,
            eps_r=self.eps_r,
            **kwargs,
        )

        self.layers = [atm]

        # store surface gas properties
        self.dens_surface = self.layers[0].dens
        self.temp_surface = self.layers[0].temp
        self.pres_surface = self.layers[0].pres

    def construct(self, N=1, **kwargs):
        """N sets the number of shells in the case that no surface constrain
        type has been specified. Otherwise N will be ignored.
        """
        # iterate over all layers
        ll = 0
        shellIteration = True
        shellIterationCount = 0
        while shellIteration:
            shellIterationCount += 1
            sys.stdout.flush()
            print("\n===================================================", end="\n")
            print(f"> processing layer {ll+1}", end="\n")
            print("layer properties:")

            for m in range(len(self.layers[-1].materials)):
                print(
                    self.layers[-1].fractions[m] * 100,
                    "%",
                    self.layers[-1].materials[m],
                )

            # perform next integration step
            if not self.M_atmos == 0.0:
                self.layers[-1].Construct()

            # if outer pressure or temperature drop bellow
            # threshold value force cancellation of iteration
            if self.layers[-1].pres < P_zero:
                print(
                    "Zero pressure reached at",
                    self.layers[-1].pres,
                    "Pa after",
                    shellIterationCount,
                    "iterations",
                )
                shellIteration = False

            if self.layers[-1].temp < T_zero:
                print(
                    "Zero temperature reached at",
                    self.layers[-1].temp,
                    "K after",
                    shellIterationCount,
                    "iterations",
                )
                shellIteration = False

            if self.layers[-1].dens < rho_zero:
                print(
                    "Zero density reached at",
                    self.layers[-1].dens,
                    "kg m-3 after",
                    shellIterationCount,
                    "iterations",
                )
                shellIteration = False

            if self.majorConstraint == "P_space":

                # if surface pressure is exceeded, reset and re-integrate old shell
                pres = self.layers[-1].pres
                print(self.layers[-1].pres)
                if (pres - self.P_space_should) / pres < -eps_Pspace:
                    self.layers[-1].shells.pop(-1)
                    self.layers[-1].shells[-1].Reset()
                    self.layers[-1].Update()
                    self.layers[-1].dr = self.layers[-1].dr * 0.5
                    self.layers[-1].bisec = True

                elif abs(pres - self.P_space_should) / pres < eps_Pspace:
                    print(
                        "Desired precission for surface pressure reached",
                        "after",
                        shellIterationCount,
                        "iterations",
                    )
                    shellIteration = False

                else:
                    self.layers[-1].bisec = False

            elif self.majorConstraint == "M_atmos":

                # if atmosphere mass is exceeded, reset and re-integrate old shell
                mass = self.layers[-1].mass - self.M_surface
                # print ('mass:', mass/m_earth,'M_earth')
                # print ('M_atmos:', self.M_atmos/m_earth,'M_earth')
                # print ('M_surf:',self.M_surface/m_earth,'M_earth')
                # print ('layer mass:',self.layers[-1].mass/m_earth, 'M_earth')
                # print ('outer pressure:', self.layers[-1].pres)
                # print ('diff:', mass-self.M_atmos)

                if mass == 0.0 or self.M_atmos == 0.0:
                    dev = 0.0

                else:
                    dev = (mass - self.M_atmos) / self.M_atmos

                if dev > eps_Mtot:
                    print("enabeling bisection")
                    self.layers[-1].shells.pop(-1)
                    self.layers[-1].shells[-1].Reset()
                    self.layers[-1].Update()
                    self.layers[-1].dr = self.layers[-1].dr * 0.5
                    self.layers[-1].bisec = True

                elif abs(dev) < eps_Mtot:
                    print(
                        "Desired precission for atmosphere mass reached",
                        "after",
                        shellIterationCount,
                        "iterations",
                    )
                    shellIteration = False

                else:
                    self.layers[-1].bisec = False

            else:
                for j in range(N):
                    print(f"\r> constructing shell {j+1}", end="")
                    self.layers[-1].Construct()

            radius = self.layers[-1].radius
            dr = self.layers[-1].dr / 1000
            mass = self.layers[-1].mass
            pres = self.layers[-1].pres
            temp = self.layers[-1].temp
            dens = self.layers[-1].shells[-1].dens
            self.height = radius - self.R_surface
            self.indigenous_mass = self.layers[-1].indigenous_mass

            # generate iteration output
            digits = "{:.5E}"
            if pres > 1.0e9:
                pres_unit = "GPa"
                pres_scale = 1.0e-9
            else:
                pres_unit = "bar"
                pres_scale = 1.0e-5

            dat = [
                [
                    digits.format(Decimal(str(radius / r_earth))),
                    digits.format(Decimal(str(round(dr, 6)))),
                    digits.format(Decimal(str(pres * pres_scale))),
                    digits.format(Decimal(str(mass / m_earth))),
                    digits.format(Decimal(str(round(temp, temp_round)))),
                    digits.format(Decimal(str(dens))),
                ]
            ]

            tabl = tabulate(
                dat,
                headers=[
                    "R [R_e]",
                    "dr [km]",
                    "P [" + pres_unit + "]",
                    "m [M_e]",
                    "T [K]",
                    "rho [kg m-3]",
                ],
            )

            # print ('\nlayertop parameters:\n')
            print(f"\n\r{tabl}", end="")
            print("\n-> completed")

            # check if particle velocity exceeds escape velocity and stop
            # integration if so as gases are expected to be lost very quickly
            # to space beyond this point

            v_esc = self.layers[0].shells[-1].v_esc
            v0 = self.layers[0].shells[-1].v0[0]
            if abs((v_esc - v0)) / v_esc < 0.8:
                print(
                    "Escape velocity reached after", shellIterationCount, "iterations"
                )
                print(v_esc - v0)
                print(v_esc, v0)
                shellIteration = False

        self.R_tot = self.layers[-1].radius
        self.M_tot = self.layers[-1].mass
        self.P_space_is = self.layers[-1].pres
        self.T_space_is = self.layers[-1].temp

    def plot(
        self,
        check=False,
        tempType="isothermal",
        log=False,
        pres_unit="GPa",
        temp_unit="K",
        dens_unit="kg/m3",
        mass_unit="m_earth",
        rad_unit="r_earth",
        **kwargs,
    ):

        self.layers[-1].plot(
            pres_unit=pres_unit,
            temp_unit=temp_unit,
            dens_unit=dens_unit,
            mass_unit=mass_unit,
            rad_unit=rad_unit,
        )

        from pics.runparams import param_colors

        # plot barometric formula for linear or constant temperature profile
        # as reference
        if check:

            # compute mean molecular weight
            mu = sum(
                [
                    gas_molar_mass_list[self.contents[i]] * self.fractions[i]
                    for i in range(len(self.contents))
                ]
            )

            M_surf, R_surf, T_surf, P_surf = (
                self.layers[-1].shells[0].mass,
                self.layers[-1].shells[0].radius,
                self.layers[-1].shells[0].temp,
                self.layers[-1].shells[0].pres,
            )
            # compute surface gravity
            gravity = G * M_surf / R_surf**2

            # compute atmospheric density at surface level
            rho0 = self.layers[-1].shells[0].dens

            x = np.array([sh.radius for sh in self.layers[-1].shells])
            x = np.linspace(R_surf, self.R_tot, 20)

            temp = np.array([T_surf for xx in x])

            p = np.array(
                [
                    P_barom(
                        mu=mu, r=xx, P0=P_surf, T=T_surf, g=gravity, R_planet=R_surf
                    )
                    for xx in x
                ]
            )

            rho = np.array(
                [
                    rho_barom(
                        mu=mu, r=xx, rho0=rho0, T=T_surf, g=gravity, R_planet=R_surf
                    )
                    for xx in x
                ]
            )

            mass = np.array(
                [
                    m_barom(
                        mu=mu,
                        r=xx,
                        rho0=rho0,
                        T=T_surf,
                        g=gravity,
                        R_planet=R_surf,
                        M0=M_surf,
                    )
                    for xx in x
                ]
            )

            if len(x) <= 20:
                markevery = 1
            if len(x) > 20:
                markevery = 2
            if len(x) > 100:
                markevery = 3

            x = x / r_earth
            self.layers[0].axes[0].plot(
                x,
                p / plot_units[pres_unit],
                markersize=3,
                marker="o",
                markevery=markevery,
                linestyle="None",
                color=param_colors[0],
                label="analytical model",
            )
            self.layers[0].axes[2].plot(
                x,
                rho / plot_units[dens_unit],
                markersize=3,
                marker="o",
                markevery=markevery,
                linestyle="None",
                color=param_colors[2],
                label="analytical model",
            )
            self.layers[0].axes[1].plot(
                x,
                mass / plot_units[mass_unit],
                markersize=3,
                marker="o",
                markevery=markevery,
                linestyle="None",
                color=param_colors[1],
                label="analytical model",
            )
            self.layers[0].axes[3].plot(
                x,
                temp / plot_units[temp_unit],
                markersize=3,
                marker="o",
                markevery=markevery,
                linestyle="None",
                color=param_colors[3],
                label="analytical model",
            )

            self.layers[0].axes[0].scatter(
                np.amax(x) * 0.9725,
                np.amax(p) * 0.75 / plot_units[pres_unit],
                color=param_colors[0],
            )

            self.layers[0].axes[0].plot(
                [np.amax(x) * 0.97, np.amax(x) * 0.975],
                [
                    np.amax(p) * 0.9 / plot_units[pres_unit],
                    np.amax(p) * 0.9 / plot_units[pres_unit],
                ],
                color=param_colors[0],
                linestyle="-",
            )

            self.layers[0].axes[0].text(
                np.amax(x) * 0.9775,
                np.amax(p) * 0.874 / plot_units[pres_unit],
                "computed",
            )

            self.layers[0].axes[0].text(
                np.amax(x) * 0.9775,
                np.amax(p) * 0.72 / plot_units[pres_unit],
                "analytical",
            )
            """
            self.layers[0].axes[0].legend()
            self.layers[0].axes[2].legend()
            self.layers[0].axes[1].legend()
            self.layers[0].axes[3].legend()  
            """
