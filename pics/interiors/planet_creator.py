# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import itertools
from pics.utils.print_tools import print_planet
from pics.utils.plot_tools import plot_structure
from tabulate import tabulate

# from pics.utils import functionTools as ftool
import numpy as np
import copy
from pics.runparams import (
    fortplanet_input_keys,
    initial_predictor_keys,
    fortplanet_keys_translator,
    fortplanet_output_keys,
)

from pics.utils.initial_conditions import predict_initials
from pics.utils import fortplanet, fortfunctions
from pics.materials import Material
from pics.physicalparams import (
    material_list_fort,
    material_YMg,
    material_YSi,
    material_YO,
    material_YH,
    material_YS,
    material_YFe,
    mSi,
    mS,
    mO,
    mH,
    mFe,
    mMg,
    r_earth,
    m_earth,
    sigmaSB,
    G,
    mH2O,
    material_list,
)

# Load eos tables
fortplanet.wrapper.load_eos_tables(table_dir="{}/data/EoS_tables/".format("."))


class Parameters:
    def __init__(self, default_values):
        self.default_values = default_values
        self.allowed_keys = self.default_values.keys()

    def set_values(self, **kwargs):
        """Set planetary parameters to custom values."""
        for key, value in kwargs.items():
            if key in self.allowed_keys:
                setattr(self, key, value)

            else:
                raise KeyError("Invalid keyargument <{}> passed".format(key))

    def set_default_values(self):
        """Set planetary parameters to default values."""
        for key, value in self.default_values.items():
            setattr(self, key, value)


class PlanetaryOutputParams(Parameters):
    def __init__(self, **kwargs):
        self.default_values = {
            "S_count": 0,
            "Si_count": 0,
            "Fe_count": 0,
            "Mg_count": 0,
            "O_count": 0,
            "H_count": 0,
            "H2O_count": 0,
            "Mg_number_is": None,
            "Si_number_is": None,
            "M_surface_is": None,
            "P_surface_is": None,
            "T_surface_is": None,
            "moment_of_inertia_is": None,
            "ocean_fractions_is": None,
            "ocean_depth": None,
            "mean_density": None,
            "oxygen_fugacity_mantle": None,
            "luminosity_int_is": None,
            "luminosity_eff_is": None,
            "atomic_fractions_core": [],
            "inner_core_fraction": None,
            "layer_properties": [],
        }


class RunOutputParams(Parameters):
    def __init__(self, **kwargs):
        self.default_values = {
            "M_surface_is": None,
            "R_surface_is": None,
            "T_surface_is": None,
            "Mg_number_is": None,
            "Si_number_is": None,
            "Fe_count": None,
            "Si_count": None,
            "Mg_count": None,
            "O_count": None,
            "H2O_count": None,
            "H_count": None,
            "S_count": None,
            "ocean_fraction_is": None,
            "moment_of_inertia_is": None,
            "layer_properties": [],
            "profiles": [],
            "shell_count": 0,
            "layer_count": 0,
            "fractions_out": [],
            "x_Fe_mantle": None,
            "Si_number_mantle": None,
            "mantle_exists": False,
            "inner_core_exists": False,
            "outer_core_exists": False,
            "gravitational_energy": None,
            "internal_energy": None,
        }


class PlanetaryInputParams(Parameters):
    def __init__(self, type="telluric", **kwargs):

        if type == "telluric":
            self.default_values = {
                "Si_number_mantle": 0.4,
                "Fe_number_mantle": 0.0,
                "M_surface_should": 1.0,
                "Mg_number_should": 0.5,
                "T_surface_should": 300.0,
                "P_surface_should": 1e5,
                "R_surface_should": 1.0,
                "ocean_fraction_should": -10,
                "contents": [[2], [2, 9, 9, 9, 9], [4, 5], [6, 7]],
                "fractions": [[1.0], [1.0, 0.0, 0.0, 0.0, 0.0], [0.5, 0.5], [0.5, 0.5]],
                "layer_masses": [0.25, 0.25, 0.0, 100.0],
                "temperature_jumps": [0, 0, 0, 0],
                "grueneisen_gammas_layers": [1.36, 1.36, 1.96, 1.26, 1.0],
                "debye_exponents_layers": [0.91, 0.91, 2.5, 2.9, 1.0],
                "total_energy_should": 0,
                "luminosity_int_should": 0,
                "luminosity_eff_should": 0,
                "layer_pressures": [0, 0, 25.0e9, 0],
                "layer_radii": [0, 0, 0, 0],
                "layer_temperatures": [0, 0, 0, 0],
                "rotation_velocity": 0.0,
                "core_segregation_pressure": 0.0,
                "spin": 0.0,
                "x_all_core": [],
                "eta_all_core": [],
            }

        elif type == "aqua":
            self.default_values = {
                "Si_number_mantle": 0.4,
                "Fe_number_mantle": 0.0,
                "M_surface_should": 1.0,
                "Mg_number_should": 0.5,
                "T_surface_should": 300.0,
                "P_surface_should": 1e5,
                "ocean_fraction_should": 0.1,
                "contents": [[2], [2, 9, 9, 9, 9], [4, 5], [6, 7], [1]],
                "fractions": [
                    [1.0],
                    [1.0, 0.0, 0.0, 0.0, 0.0],
                    [0.5, 0.5],
                    [0.5, 0.5],
                    [1.0],
                ],
                "layermasses": [0.5, 0.5, 0.0, 0.0, 100.0],
            }

        elif type == "inferno":
            self.default_values = {
                "Si_number_mantle": 0.4,
                "Fe_number_mantle": 0.0,
                "M_surface_should": 1.0,
                "Mg_number_should": 0.5,
                "T_surface_should": 1000.0,
                "P_surface_should": 1e5,
                "ocean_fraction_should": 0.0,
                "contents": [[2], [2, 9, 9, 9, 9], [4, 5], [6, 7]],
                "fractions": [[1.0], [1.0, 0.0, 0.0, 0.0, 0.0], [0.5, 0.5], [0.5, 0.5]],
                "layermasses": [0.5, 0.5, 0.0, 100.0],
            }

        elif type == "ice":
            self.default_values = {
                "Si_number_mantle": 0.4,
                "Fe_number_mantle": 0.0,
                "M_surface_should": 1.0,
                "Mg_number_should": 0.5,
                "T_surface_should": 200.0,
                "P_surface_should": 1e5,
                "ocean_fraction_should": 0.1,
                "contents": [[2], [2, 9, 9, 9, 9], [4, 5], [6, 7]],
                "fractions": [[1.0], [1.0, 0.0, 0.0, 0.0, 0.0], [0.5, 0.5], [0.5, 0.5]],
                "layermasses": [0.5, 0.5, 0.0, 100.0],
            }

        else:

            if not "default_values" in kwargs:
                raise KeyError(
                    "If you don't specify the planet type you must provide all parameters!"
                )

            else:
                self.default_values = {
                    "Si_number_mantle": 0.4,
                    "Fe_number_mantle": 0.0,
                    "M_surface_should": 1.0,
                    "Mg_number_should": 0.5,
                    "T_surface_should": 300.0,
                    "P_surface_should": 1e5,
                    "ocean_fraction_should": 0.0,
                    "contents": [[2], [2, 9, 9, 9, 9], [4, 5], [6, 7]],
                    "fractions": [
                        [1.0],
                        [1.0, 0.0, 0.0, 0.0, 0.0],
                        [0.5, 0.5],
                        [0.5, 0.5],
                    ],
                    "layermasses": [0.5, 0.5, 0.0, 100.0],
                    "temperature_jumps": [0, 0, 0, 0],
                    "grueneisen_gammas_layers": [1.36, 1.36, 1.96, 1.26, 1.0],
                    "debye_exponents": [0.91, 0.91, 2.5, 2.9, 1.0],
                }

                self.default_values.update(kwargs["default_values"])

        Parameters.__init__(self, self.default_values)


class RunInputParams(Parameters):
    def __init__(self, **kwargs):
        self.default_values = {
            "adiabat_type": 0,
            "layer_constraints": [1, 1, 3, 1],
            "temperature_type": 1,
            "layer_type": 1,
            "major_constraint": "P_surface",
            "eps_r": 0.5,
            "density_type": 1,
            "subphase_resolution": 16,
            "core_segregation_type": 0,
            "inner_core_segregation_type": 0,
            "seed_radius": 0.1,
            "initial_predictor": 0,
        }

        Parameters.__init__(self, self.default_values)


class AllInputParams(RunInputParams, PlanetaryInputParams):
    def __init__(self):
        pass


class Planet:
    def __init__(self, label="A random planet", **kwargs):
        # self.default_values = {}
        # self.default_values.update(planetary_params.default_values)
        # self.default_values.update(run_params.default_values)
        self.default_values = {}
        self.initials = {}
        self.label = label
        self.vapor_reached = False
        self.status = "shadow of the future"

    def set_values(self, default=False, **kwargs):
        omit_keys = ["default_values", "allowed_keys", "label"]
        out_params = PlanetaryOutputParams()

        # initialize output parameters
        for key, value in out_params.default_values.items():
            setattr(self, key, value)

        # initialize input parameters
        if "planetary_params" in kwargs:
            # set attributes
            # print ("kwargs =", kwargs["planetary_params"].default_values.items())
            for key, value in kwargs["planetary_params"].default_values.items():
                if not key in omit_keys:
                    # print (key, value)
                    setattr(self, key, value)

        # print ("contents in set values =", self.contents)
        # print ("temp in set values =", self.T_surface_should)
        # initialize run parameters
        if "run_params" in kwargs:
            for key, value in kwargs["run_params"].__dict__.items():
                if not key in omit_keys:
                    setattr(self, key, value)

        # Set the given values as the default for this planetary object
        if default:
            self.default_values = copy.deepcopy(self.__dict__)

        # set all secondary parameters from the inputs
        self.layer_dims = [len(c) for c in self.contents]
        self.layer_properties = [
            {
                "P_outer": None,
                "T_outer": None,
                "rho_outer": None,
                "rho_inner": None,
                "indigenous_mass": 0.0,
                "R_outer": None,
                "n_shells": 0,
            }
            for i in range(len(self.contents))
        ]
        # predict central temperature, central pressure, and core mass
        # prediction for core mass is not relevant for basic models as it
        # is uniquely defined from the bulk composition in this case
        kwargs_pred = dict([key, self.__dict__[key]] for key in initial_predictor_keys)
        tc, pc, mc = predict_initials(**kwargs_pred)
        pc *= 1e9  # convert to Pa
        self.T_center = tc
        self.P_center = pc

        self.Si_number_layers = [
            0.0,
            0.0,
            self.Si_number_mantle,
            self.Si_number_mantle,
            0.0,
        ]
        self.Fe_number_layers = [
            1.0,
            1.0,
            self.Fe_number_mantle,
            self.Fe_number_mantle,
            0.0,
        ]
        self.M_ocean_should = self.M_surface_should * 10**self.ocean_fraction_should

        print("predicted central values are:", tc, pc, mc)

        try:
            self.x_all_core = Material.mat2at_core(xi=self.fractions[1], xiH=0.0)

            self.eta_all_core = Material.at2wt_core(self.x_all_core)

        except IndexError:
            self.x_all_core = None
            self.eta_all_core = None

        if self.initial_predictor == 0:
            # compute inner core mass from bulk composition
            # Note: inner core mass will be updated during the structure
            # integration from the solidus of iron-alloys

            # compute total core mass from bulk composition
            self.layer_masses[0] = self.compute_core_mass(M_IC=1.0)
            self.layer_masses[1] = self.compute_core_mass(M_IC=0.0)

        # Use more sophisticated multilinear regression models to predict the
        # core mass, central pressure, and central temperature.
        elif self.initial_predictor == 1:
            M_inner_core = self.M_core * self.inner_core_frac
            M_outer_core = self.M_core - self.M_inner_core
            self.layer_masses[0] = M_outer_core + M_inner_core
            self.layer_masses[1] = M_outer_core + M_inner_core

        self.M_core_should = self.layer_masses[1]
        print("the initial values are:", tc, pc, self.layer_masses)

    def update_initials(self):
        """Updates the initial planetary parameters from the current values"""
        self.intials = None

    def update_values(self):
        self.default_values = copy.deepcopy(self.__dict__)

    def compute_core_mass(self, n=3, M_IC=0.0, xi_all_core=[]):
        """Computes the core mass of a planet at given total mass, composition and
        value for Mg#
        """

        # print (contents, Mg_number, M_surface, Mg_number_mantle, SiMg, M_ocean, xi_H_core)
        Mg_number_mantle = min(1.0 - self.Fe_number_mantle, 0.9999999999)
        FeMg_mantle = (1.0 - Mg_number_mantle) / Mg_number_mantle
        FeMg = (1.0 - self.Mg_number_should) / self.Mg_number_should
        SiMg = self.Si_number_mantle / (1.0 - self.Si_number_mantle)

        # Compute the fractions in the mantle
        fractions = fortfunctions.functionspy.compute_abundance_vector(
            simg=SiMg,
            femg=FeMg_mantle,
            n_mats=len(self.contents[n]),
            ymgi=[material_YMg[i - 1] for i in self.contents[n]],
            ysii=[material_YSi[i - 1] for i in self.contents[n]],
            xih2oi=[0.0 for i in self.contents[n]],
            xifei=[1.0 - Mg_number_mantle for i in self.contents[n]],
            xialsii=[0.0 for i in self.contents[n]],
            xialmgi=[0.0 for i in self.contents[n]],
            contents=self.contents[n],
            additional=[],
        )
        # print ("fractions in compute core mass =", fractions)

        # Count
        # Note that the indices are shifted by one because the original definition
        # of the arrays comes from the fortran code.
        Q1 = sum(
            [
                fractions[i] * Mg_number_mantle * material_YMg[self.contents[n][i] - 1]
                for i in range(len(self.contents[n]))
            ]
        )

        # Compute total normalized mass in the mantle
        Q2 = (
            sum(
                [
                    fractions[i]
                    * Mg_number_mantle
                    * material_YMg[self.contents[n][i] - 1]
                    for i in range(len(self.contents[n]))
                ]
            )
            * mMg
            + sum(
                [
                    fractions[i]
                    * (1.0 - Mg_number_mantle)
                    * material_YMg[self.contents[n][i] - 1]
                    for i in range(len(self.contents[n]))
                ]
            )
            * mFe
            + sum(
                [
                    fractions[i] * material_YSi[self.contents[n][i] - 1]
                    for i in range(len(self.contents[n]))
                ]
            )
            * mSi
            + sum(
                [
                    fractions[i] * material_YO[self.contents[n][i] - 1]
                    for i in range(len(self.contents[n]))
                ]
            )
            * mO
            + sum(
                [
                    fractions[i] * material_YH[self.contents[n][i] - 1]
                    for i in range(len(self.contents[n]))
                ]
            )
            * mH
        )

        Q3 = sum(
            [
                fractions[i]
                * (1.0 - Mg_number_mantle)
                * material_YMg[self.contents[n][i] - 1]
                for i in range(len(self.contents[n]))
            ]
        )

        # Q4 = 1. #- xi_H_core
        # Q5 = (mFe + xi_S_core*mS)#*(1.-xi_H_core) + xi_H_core*mH
        # core_frac = (1.-M_ocean/self.M_surface_should)*(Q1/Q2-self.Mg_number_should*(Q1/Q2+Q3/Q2))/\
        #             (self.Mg_number_should*(Q4/Q5-Q1/Q2-Q3/Q2)+Q1/Q2)

        Q4 = self.x_all_core[0]
        m_core = [mFe, mH, mS, mSi, mO]
        Q5 = sum([m_core[i] * self.x_all_core[i] for i in range(len(self.x_all_core))])

        # print ('Q =', Q1, Q2, Q3, Q4, Q5)
        # core_frac = 1.0 - 10**self.ocean_fraction_should
        # core_frac *= Q1 / Q2 - self.Mg_number_should * (Q1 / Q2 + Q3 / Q2)
        # core_frac /= self.Mg_number_should * (Q4 / Q5 - Q1 / Q2 - Q3 / Q2) + Q1 / Q2
        # print ('core mass old =', core_frac * self.M_surface_should)

        core_frac = 1.0 - 10**self.ocean_fraction_should
        core_frac *= Q3 / Q2 - Q1 / Q2 * FeMg
        core_frac += M_IC / self.M_surface_should * (1.0 / mFe - Q4 / Q5)
        core_frac /= Q3 / Q2 - Q4 / Q5 - FeMg * Q1 / Q2
        # print ('core mass new =', core_frac * self.M_surface_should)
        # print ('')
        return core_frac * self.M_surface_should

    def update(self, default=False):
        """Computes all dependant planetary parameters"""
        self.update_composition()
        self.update_bulk()
        self.update_core()
        self.update_mantle()
        self.update_hydrosphere()

        # Set the given values as the default for this planetary object
        if default:
            self.default_values = copy.deepcopy(self.__dict__)

    def reset(self):
        """Resets all planetary parameters to the initial values"""
        for key in fortplanet_input_keys:
            setattr(self, key, self.initials[key])

    def print(self, style=0, digits=3):
        """Prints out a simple overview of all relevant planetary parameters."""
        print_planet(self, style=style, digits=digits)

    def update_layers(self, props):
        for i in range(len(props)):

            self.layer_properties[i]["P_outer"] = props[i][0]
            self.layer_properties[i]["T_outer"] = props[i][1]
            self.layer_properties[i]["rho_outer"] = props[i][2]
            self.layer_properties[i]["rho_inner"] = props[i][3]
            self.layer_properties[i]["indigenous_mass"] = props[i][4]
            self.layer_properties[i]["mass_fraction"] = (
                self.layer_properties_dummy[i][4] / self.M_surface_is
            )
            self.layer_properties[i]["R_outer"] = props[i][5]
            self.layer_properties[i]["n_shells"] = self.shell_count_layers[i]

            if i > 0:
                self.layer_properties[i]["P_inner"] = self.layer_properties[i - 1][
                    "P_outer"
                ]
                self.layer_properties[i][
                    "T_inner"
                ] = None  # self.layer_properties[i - 1]['T_outer']
                self.layer_properties[i]["R_inner"] = self.layer_properties[i - 1][
                    "R_outer"
                ]

            else:
                self.layer_properties[i]["P_inner"] = self.P_center
                self.layer_properties[i]["T_inner"] = self.T_center
                self.layer_properties[i]["R_inner"] = self.seed_radius

        self.M_core_is = (
            self.layer_properties[0]["indigenous_mass"]
            + self.layer_properties[1]["indigenous_mass"]
        ) / m_earth
        self.M_mantle_is = (
            self.layer_properties[2]["indigenous_mass"]
            + self.layer_properties[3]["indigenous_mass"]
        ) / m_earth

        try:
            self.M_ocean_is = self.layer_properties[4]["indigenous_mass"] / m_earth

        except IndexError:
            self.M_ocean_is = 0.0

    def check_convergence(self):
        self.converged = False
        pass

    def construct(self, echo=False):
        # Gather layer dims for fortplanet routine
        layer_dims = [len(item) for item in self.contents]
        # print ("contents =", self.contents)
        conts = list(itertools.chain.from_iterable(self.contents))
        fracs = list(itertools.chain.from_iterable(self.fractions))

        # generate intput parameters for the planet constructor
        kwargs = dict(
            (name, getattr(self, key))
            for name, key in zip(fortplanet_keys_translator, fortplanet_input_keys)
        )

        # convert contents and fractions nd lists to 1d lists
        kwargs["fractions"] = fracs
        kwargs["contents"] = conts

        # store initial values for subsequent use by the planet_iterator
        self.initials = dict([key, self.__dict__[key]] for key in fortplanet_input_keys)

        # fortran wrapper is called here
        # output = test_interface.interface.do_some_science_stuff(**kwargs)
        output = fortplanet.wrapper.create_planet(**kwargs)

        # update planetary output parameters
        for key, value in zip(fortplanet_output_keys, output):
            setattr(self, key, value)

        self.mean_density = self.M_surface_is / (
            4.0 / 3.0 * np.pi * self.R_surface_is**3
        )

        self.update_layers(self.layer_properties_dummy)
        self.trim_profiles()

        self.status = "very much alive"

    def trim_profiles(self):
        """The profiles output from the Fortran routine contains a sequence
        of zeros at the end. These elements are being removed here. The
        convention of the profiles is:

            Radius (m), Temperature (K), Pressure (Pa), Density (kg/m3),
            Mass (kg), Gravity (m/s2), Gravitational energy (J)
        """
        off = len(self.contents) - 1
        prof = np.empty([len(self.profiles), self.shell_count + off])

        i = 0
        while True:
            for p in range(len(self.profiles)):
                prof[p][i] = self.profiles[p][i]

            i += 1

            if i == self.shell_count + off:
                break

        self.profiles = np.array(prof)

        # normalize moi
        self.profiles[6] *= 1.0 / (self.M_surface_is * self.R_surface_is**2)

    def plot(
        self,
        **kwargs,
    ):

        plot_structure(self.profiles, **kwargs)

    def set_iterator_specs(self):
        pass


class TelluricPlanet(Planet):
    def __init__(self, planetary_params={}, run_params={}):

        Planet.__init__(self, label="telluric")

        pp = PlanetaryInputParams(type=self.label)
        rp = RunInputParams(type=self.label)

        pp.set_default_values()
        rp.set_default_values()

        # update planetary parameters if passed by user
        for key, val in planetary_params.items():
            pp.default_values.update({key: val})

        # update run parameters if passed by user
        for key, val in run_params.items():
            rp.default_values.update({key: val})

        Planet.set_values(self, planetary_params=pp, run_params=rp, default=True)


class AquaPlanet(Planet):
    def __init__(self):
        raise NotImplementedError("Aqua planets will be available soon!")
        Planet.__init__(self, label="aqua")
        pp = PlanetaryInputParams(type=self.label)
        rp = RunInputParams(type=self.label)
        pp.set_default_values()
        rp.set_default_values()
        Planet.set_values(self, planetary_params=pp, run_params=rp, default=True)


class InfernoPlanet(Planet):
    def __init__(self):
        raise NotImplementedError("Inferno planets will be available soon!")
        Planet.__init__(self, label="inferno")
        pp = PlanetaryInputParams(type=self.label)
        rp = RunInputParams(type=self.label)
        pp.set_default_values()
        rp.set_default_values()
        Planet.set_values(self, planetary_params=pp, run_params=rp, default=True)


class IcePlanet(Planet):
    def __init__(self):
        raise NotImplementedError("Ice planets will be available soon!")
        Planet.__init__(self, label="ice")
        pp = PlanetaryInputParams(type=self.label)
        rp = RunInputParams(type=self.label)
        pp.set_default_values()
        rp.set_default_values()
        Planet.set_values(self, planetary_params=pp, run_params=rp, default=True)


class CustomPlanet(Planet):
    def __init__(self, planetary_parameters={}, run_parameters={}):
        raise NotImplementedError("Custom type planets will be available soon!")
        default_values1 = {
            "ocean_fraction_should": 0.25,
            "Mg_number_should": 0.6,
            "T_surface_should": 320,
            "Fe_number_mantle": 0.15,
        }

        default_values2 = {"adiabat_type": 0}

        default_values1.update(planetary_parameters)
        Planet.__init__(self, label="custom")
        pp = PlanetaryInputParams(type=self.label, default_values=default_values1)
        rp = RunInputParams(type=self.label, default_values=default_values2)
        pp.set_default_values()
        rp.set_default_values()
        Planet.set_values(self, planetary_params=pp, run_params=rp, default=True)


# pp = PlanetaryParams()
# rp = RunParams()

# pp.set_values(**pp.default_values)
# rp.set_values(**rp.default_values)

# pl = Planet()
# pl.set_values(planetary_params = pp, run_params = rp, default = True)


# pl.Si_number_mantle = 6

# #print (pl.default_values)
# pl.reset()

# pl.update(default = True)
# pl.show()

# pl1 = AquaPlanet()
# pl2 = TelluricPlanet()
# pl3 = InfernoPlanet()
# pl4 = IcePlanet()
# pl5 = CustomPlanet(planetary_parameters={"Fe_number_mantle":0.2})

# pl1.show()
# pl2.show()
# pl3.show()
# pl4.show()
# pl5.show()
