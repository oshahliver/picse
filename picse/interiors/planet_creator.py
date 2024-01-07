# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import itertools
from picse.utils.print_tools.print_tools import print_planet
from picse.utils.plot_tools.plot_tools import plot_structure
from tabulate import tabulate
from picse.utils.file_tools.internal_data import get_eos_dir
from picse.interiors import core_creator
import random

# from picse.utils import functionTools as ftool
import numpy as np
import copy
from picse.runparams import (
    fortplanet_input_keys,
    initial_predictor_keys,
    fortplanet_keys_translator,
    fortplanet_output_keys,
    supported_base_types,
)

from picse.utils.calibration_tools.initial_conditions import predict_initials
from picse.utils.file_tools.internal_data import get_predictor_model
from picse.utils import fortplanet
from picse.materials import material
from picse.physicalparams import (
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


def load_eos_tables(**kwargs):
    """Loads the equation of state tables into memory for subsequent use during
    the structure integration.

    TODO. Handle internally!
    """
    fortplanet.wrapper.load_eos_tables(table_dir=get_eos_dir())


def load_predictors(**kwargs):
    """Loads the pre-calibrated predictor models for the different base types for
    predicting the initial conditions for the initial structure integration during
    the iteration process.

    TODO. Handle internally!
    """

    predictors = {}
    for bt in supported_base_types:
        mod = get_predictor_model("predictor_{}.pkl".format(bt))
        predictors.update({bt: mod})

    return predictors


predictors = load_predictors()


class Parameters:
    def __init__(self, default_values):
        # TODO. the default values make no sense. change the way the parameters are handled
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
            "ocean_fraction_is": None,
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
            "temperature_jumps": [0.0, 0.0, 0.0, 0.0],
            "grueneisen_gammas_layers": [1.36, 1.36, 1.96, 1.26],
            "debye_exponents_layers": [0.91, 0.91, 2.5, 2.9],
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

        if type == "telluric":
            self.default_values.update(
                dict(
                    Fe_number_layers=[
                        0.0,
                        0.0,
                        self.default_values["Fe_number_mantle"],
                        self.default_values["Fe_number_mantle"],
                    ],
                    Si_number_layers=[
                        0.0,
                        0.0,
                        self.default_values["Si_number_mantle"],
                        self.default_values["Si_number_mantle"],
                    ],
                )
            )

        elif type == "aqua":
            # Note. Add additional layer for the hydrosphere
            self.default_values.update(
                dict(
                    Fe_number_layers=[
                        0.0,
                        0.0,
                        self.default_values["Fe_number_mantle"],
                        self.default_values["Fe_number_mantle"],
                        0.0,
                    ],
                    Si_number_layers=[
                        0.0,
                        0.0,
                        self.default_values["Si_number_mantle"],
                        self.default_values["Si_number_mantle"],
                        0.0,
                    ],
                )
            )
            new_specs = dict(
                contents=[[2], [2, 9, 9, 9, 9], [4, 5], [6, 7], [1]],
                fractions=[
                    [1.0],
                    [1.0, 0.0, 0.0, 0.0, 0.0],
                    [0.75, 0.25],
                    [0.5, 0.5],
                    [1.0],
                ],
                ocean_fraction_should=np.log10(0.1),
                layer_pressures=[0.0, 0.0, 25.0e9, 0.0, 0.0],
                debye_exponents_layers=[0.91, 0.91, 2.5, 2.9, 1.0],
                grueneisen_gammas_layers=[1.36, 1.36, 1.96, 1.26, 1.0],
                temperature_jumps=[0.0, 0.0, 0.0, 0.0, 0.0],
                layer_radii=[0.0, 0.0, 0.0, 0.0, 0.0],
            )
            mantle_mass = self.default_values["M_surface_should"] * (
                1.0 - 10 ** new_specs["ocean_fraction_should"]
            )
            new_specs.update(
                {"layer_masses": [0.25, 0.25, mantle_mass, mantle_mass, 100.0]}
            )

            self.default_values.update(new_specs)

        else:

            if not "default_values" in kwargs:
                raise KeyError(
                    "If you don't specify the planet type you must provide all parameters!"
                )

            else:
                self.default_values.update(kwargs["default_values"])

        Parameters.__init__(self, self.default_values)


class RunInputParams(Parameters):
    def __init__(self, type="telluric", **kwargs):
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

        if type == "telluric":
            pass

        elif type == "aqua":
            # Add additional layer for the hydrosphere
            new_specs = dict(layer_constraints=[1, 1, 3, 1, 1])
            self.default_values.update(new_specs)

        Parameters.__init__(self, self.default_values)


class AllInputParams(RunInputParams, PlanetaryInputParams):
    def __init__(self):
        pass


class Planet:
    def __init__(
        self, predictor, predictor_type="reg", label="A random planet", **kwargs
    ):
        self.default_values = {}
        self.initials = {}
        self.label = label

        self.vapor_reached = False
        self.status = "shadow of the future"
        self.predictor_type = (
            predictor_type  # strategy for predicting initial conditions
        )
        self.iterator_specs = {}
        self.predictor = "predictor_{}.pkl".format(
            predictor
        )  # points to the predictor model for the iterator

    def set_values(
        self,
        default=False,
        pres_eps=[0.6, 0.1],
        planetary_params: PlanetaryInputParams = None,
        run_params: RunInputParams = None,
    ):
        """
        Sets up the planet and

        Parameters:
        default (bool, optional): Use default values for all properties. Defaults to False.
        pres_eps (list, optional): Relative pressure offsets for sampling. Defaults to [0.6, 01.]
        planetary_params (PlanetaryInputParams, optional): The physical parameters describing the planet. Defaults to {}.
        run_params (RunInputParams, optional): The run parameters for the simulation. Defaults to {}.
        """
        omit_keys = ["default_values", "allowed_keys", "label"]
        out_params = PlanetaryOutputParams()

        # initialize output parameters
        for key, value in out_params.default_values.items():
            setattr(self, key, value)

        # initialize input parameters
        if planetary_params is not None:
            # set attributes
            for key, value in planetary_params.default_values.items():
                if not key in omit_keys:
                    setattr(self, key, value)

        # initialize run parameters
        if run_params is not None:
            for key, value in run_params.__dict__.items():
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

        # Check if initial conditions are set manually or should be predicted
        if self.predictor_type == "reg":
            # predict central temperature, central pressure, and core mass
            # prediction for core mass is not relevant for basic models as it
            # is uniquely defined from the bulk composition in this case
            kwargs_pred = dict(
                [key, self.__dict__[key]] for key in initial_predictor_keys
            )
            # Values passed for pres_center and temp_center in the planetary_params
            # are overwritten here
            self.T_center, self.P_center, mc = predict_initials(
                predictors[self.label], **kwargs_pred
            )
            self.P_center *= 1e9  # convert to Pa

        elif self.predictor_type == "man":
            if planetary_params is not None:

                # Check if pressure and temperature values have been passed
                if (
                    "P_center" in planetary_params.default_values.keys()
                    and "T_center" in planetary_params.default_values.keys()
                ):
                    # values have already been set with the rest of the planetary params
                    # for manual predictor the core-to-mantle transition and the mantle-to-ocean
                    # transition are defined via a pressure value selected from within the range
                    # between the central pressure and the surface pressure

                    # terrestrial planets
                    eps1, eps2 = pres_eps

                    if self.label == "telluric":
                        self.layer_constraints[0] = 3
                        self.layer_constraints[1] = 3
                        self.layer_pressures[0] = 10 ** random.uniform(
                            np.log10(self.P_center * eps1),
                            np.log10(self.P_center * eps2),
                        )
                        self.layer_pressures[1] = self.layer_pressures[0]
                        print("updated layer pressures =", self.layer_pressures)

                    # aqua planets
                    elif self.label == "aqua":
                        self.layer_constraints[0] = 3
                        self.layer_constraints[1] = 3
                        self.layer_constraints[3] = 3

                        self.layer_pressures[0] = 10 ** random.uniform(
                            np.log10(self.P_center * eps1),
                            np.log10(self.P_center * eps2),
                        )
                        self.layer_pressures[1] = self.layer_pressures[0]
                        self.layer_pressures[3] = 10 ** random.uniform(
                            np.log10(self.layer_pressures[1] * eps1),
                            np.log10(self.P_surface_should * 1.1),
                        )

                    # self.P_center = kwargs["planetary_params"]["pres_center"]
                    # self.T_center = kwargs["planetary_params"]["temp_center"]

                else:
                    raise KeyError(
                        "If the initial conditions are set manually you need to pass values for temp_center and pres_center"
                    )

            else:
                raise KeyError(
                    "If predictor_type=man, values for the central pressure and temperature need to be passed to the planetary_params argument."
                )

        self.Si_number_layers = [
            0.0,
            0.0,
            self.Si_number_mantle,
            self.Si_number_mantle,
        ]
        self.Fe_number_layers = [
            1.0,
            1.0,
            self.Fe_number_mantle,
            self.Fe_number_mantle,
        ]

        self.M_ocean_should = self.M_surface_should * 10 ** self.ocean_fraction_should

        if self.label == "aqua":
            self.Si_number_layers.append(0.0)
            self.Fe_number_layers.append(0.0)

            # The iterator will automatically match the ocean fraction if the mantle
            # mass and core mass are reduced accordingly at a given total mass.
            self.layer_masses[2] = self.M_surface_should - self.M_ocean_should
            self.layer_masses[3] = self.M_surface_should - self.M_ocean_should
            # Outermost layer is defined via surface conditions and not layer mass
            self.layer_masses[4] = 100.0

        # try:
        self.x_all_core = material.mat2at_core(xi=self.fractions[1], xiH=0.0)
        self.eta_all_core = material.at2wt_core(self.x_all_core)

        # except IndexError:
        #     self.x_all_core = None
        #     self.eta_all_core = None

        if self.initial_predictor == 0:
            # compute inner core mass from bulk composition
            # Note: inner core mass will be updated during the structure
            # integration from the solidus of iron-alloys

            # Compute total core mass from bulk composition
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
        print(
            "the initial values are:", self.T_center, self.P_center, self.layer_masses
        )

    def update_initials(self):
        """Updates the initial planetary parameters from the current values"""
        self.intials = None

    def update_values(self):
        self.default_values = copy.deepcopy(self.__dict__)

    def get_hydrosphere_props(self, mode="basic"):
        """
        Extracts the structure and basic properties of the hydrosphere. These
        include: the total thickness, the layering into different phase regions,
        the liquid water mass fraction.

        Parameters:
        mode (str, optional): Specifies the amount of information to be extracted
        for the hydrosphere. Defaults to "basic".
        """
        raise NotImplementedError("Patience, young one!")

    def compute_ocean_depth(self):
        """
        Computes the depth thickness of the hydrosphere.
        """
        self.ocean_depth = self.R_surface_is - self.layer_properties[3]["R_outer"]

    def compute_oxide_fractions(self):
        """
        Computes the mole fractions of different oxides in the mantle from the atomic mantle composition. The planet needs to be constracted first.
        """
        # Compute mole fractions of oxides in the mantle
        if self.Fe_number_mantle > 0.0:
            self.xi_FeO_mantle = 1.0 - self.Si_number_mantle
            self.xi_FeO_mantle /= (
                (1.0 - self.Fe_number_mantle) / self.Fe_number_mantle
                + 1.0
                - self.Si_number_mantle
            )

        else:
            self.xi_FeO_mantle = 0.0

        self.xi_SiO2_mantle = (1.0 - self.xi_FeO_mantle) * self.Si_number_mantle
        self.xi_MgO_mantle = 1.0 - self.xi_FeO_mantle - self.xi_SiO2_mantle

    def compute_core_mass(self, n=3, M_IC=0.0):
        """Computes the core mass of a planet at given total mass, composition and
        value for Mg#
        """

        params = dict(
            M_surface_should=self.M_surface_should,
            Mg_number_should=self.Mg_number_should,
            contents=self.contents,
            Fe_number_mantle=self.Fe_number_mantle,
            Si_number_mantle=self.Si_number_mantle,
            ocean_fraction_should=self.ocean_fraction_should,
            x_all_core=self.x_all_core,
        )

        return core_creator.compute_core_mass(params, n=n, M_IC=M_IC)

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
        self.core_mass_fraction_is = self.M_core_is * m_earth / self.M_surface_is
        self.core_radius = self.layer_properties[1]["R_outer"]

        try:
            self.M_ocean_is = self.layer_properties[4]["indigenous_mass"] / m_earth

        except IndexError:
            self.M_ocean_is = 0.0

        self.mean_density = self.M_surface_is / (
            4.0 / 3.0 * np.pi * self.R_surface_is ** 3
        )
        self.compute_ocean_depth()

    def check_convergence(self, ocean=0.01):
        accs = self.iterator_specs["acc"]
        whats = self.iterator_specs["what"]

        checks = []

        # Check iterative parameters
        for acc, what in zip(accs, whats):
            if what == "M_surface":
                val_is = self.M_surface_is / m_earth
                val_should = self.M_surface_should

            elif what == "T_surface":
                val_is = self.T_surface_is
                val_should = self.T_surface_should

            reldev = abs(val_should - val_is) / val_should

            if reldev > acc:
                checks.append(1)
            else:
                checks.append(0)

        # Check passive parameters
        if self.label == "aqua":
            reldev = abs(self.M_ocean_should - self.M_ocean_is) / self.M_ocean_should

            if reldev > ocean:
                checks.append(1)
            else:
                checks.append(0)

        if sum(checks) > 0:
            self.converged = False
        else:
            self.converged = True

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
        self.profiles[6] *= 1.0 / (self.M_surface_is * self.R_surface_is ** 2)

    def plot(
        self,
        **kwargs,
    ):

        plot_structure(self.profiles, **kwargs)


class CustomPlanet(Planet):
    # TODO. Create interface to customize CustomPlanet types
    def __init__(self, label, planetary_params={}, run_params={}):
        predictor = "{}.pkl".format(label)

        Planet.__init__(self, predictor, label=label)

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


class BaseTypePlanet(Planet):
    def __init__(self, label, planetary_params={}, run_params={}, predictor_type="reg"):
        predictor = "predictor_{}.pkl".format(label)
        Planet.__init__(self, predictor, label=label, predictor_type=predictor_type)

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


class TelluricPlanet(BaseTypePlanet):
    def __init__(self, planetary_params={}, run_params={}, **kwargs):
        BaseTypePlanet.__init__(
            self,
            label="telluric",
            planetary_params=planetary_params,
            run_params=run_params,
            **kwargs,
        )


class AquaPlanet(BaseTypePlanet):
    def __init__(self, planetary_params={}, run_params={}, **kwargs):
        BaseTypePlanet.__init__(
            self,
            label="aqua",
            planetary_params=planetary_params,
            run_params=run_params,
            **kwargs,
        )


class PureSphere(BaseTypePlanet):
    def __init__(self, planetary_params={}, run_params={}):
        BaseTypePlanet.__init__(
            self, label="pure", planetary_params=planetary_params, run_params=run_params
        )


class YourPlanet(BaseTypePlanet):
    def __init__(self, planetary_params={}, run_params={}):
        BaseTypePlanet.__init__(
            self,
            label="custom",
            planetary_params=planetary_params,
            run_params=run_params,
        )

        # Create your own features and specifications here if you want ...
