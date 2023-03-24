""" TODO. Migrate all routines from the old PlanetFactory and perhaps PlanetInvestigator to this file.
This file should take care of higher level control of the modeeling workflow s.a. setting up
samples of planetary models, creating MR-relations and perform some post-processing tasks
"""

from picse.interiors import planet_creator, planet_iterator
from picse.materials import material
from picse.utils.file_tools import file_manager
from picse.utils.plot_tools.plot_tools import plot_mr
from picse.utils.file_tools import internal_data
from picse.physicalparams import r_earth, m_earth
from picse.materials.material import silicon_number_max, silicon_number_min
import sys
import numpy as np
import pandas as pd
import os
import copy
import random
from picse.physicalparams import xi_Fe_mantle_max, pres_core_seg_exponent

from alive_progress import alive_bar

# planet_creator.load_eos_tables()


class Toolkit:
    def __init__(self):
        self.iterator = planet_iterator.Toolkit()

    def sample_inputs(self, specs={}, n_planets = 10, iteration_limit = 10000):
        """Construct all input parameters for the most recent planetary model
        according to some pre-defined conditions. The input parameters are chosen
        from random sampling. Each call of this function generates a random planet
        within the allowed ranges.

        Convention is:

        mg#, cs pressure, FeS, FeSi, FeO, T TBL, T0 delta T core, mass
        """
        ranges = {
            "temp_surface":[200, 1500],
            "pres_surface":[1e-4, 1e8],
            "mass": [1e-2, 6.0],
            "Mg_number": [0.1, 0.9],
            "Fe_number_mantle":[0, 0.1],
            "Si_number_mantle":[.4, .6],
            "pres_core_seg": [3e9, 7e9],
            "x_FeS": [1e-6, 0.5],
            "x_FeSi": [1e-6, 0.3],
            "x_FeO": [1e-6, 0.3],
            "temp_TBL": [1400, 2000],
            "delta_temp_core": [1200, 1600],
            "oxygen_fug":[-10, 0],
        }
        sampling = {key: "lin" for key in ranges.keys()}
        
        # core composition sampled logarithmically by default
        # sampling["x_FeS"] = "log"
        sampling["x_FeSi"] = "log"
        sampling["x_FeO"] = "log"
        sampling["pres_surface"] = "log"

        try:
            specs.update({"core_segregation":specs["core_segregation"]})

        except KeyError:
            specs["core_segregation"] = "off"

        try:
            specs.update({"core_seg_type":specs["inverse"]})

        except KeyError:
            specs["core_seg_type"] = "inverse"

        if "ranges" in specs.keys():
            for key, val in specs["ranges"].items():
                if key in ranges.keys():
                    ranges.update({key: val})

        if "sampling" in specs.keys():
            for key, val in specs["sampling"].items():
                if key in sampling.keys():
                    sampling.update({key: val})
        
        all_inputs = {key:[] for key in ranges.keys()}
        all_inputs.update({"temp_core_seg":[]})

        # loop over number of planets in the sample
        with alive_bar(
            n_planets,
            title=f"Creating inputs",
            bar="bubbles",
            spinner="pulse",
            dual_line=True,
        ) as bar:
            for i in range(n_planets):
                bar.text = "Sampling parameters..."
                # First, create random values within the specified ranges
                inputs = {}
                for key, val in ranges.items():
                    p = random.random()

                    if sampling[key] == "lin":
                        value = val[0] + p * (val[1] - val[0])
                    elif sampling[key] == "log":
                        value = 10 ** (np.log10(val[0]) + p * np.log10(val[1] / val[0]))

                    inputs.update({key: value})

                # Impose additional constraints on the parameters and update
                # the sampled values if required.

                # The core segregation pressure scales with the mass
                inputs["pres_core_seg"] *= inputs["mass"] ** pres_core_seg_exponent

                # compute core segregation temperature
                temp_core_seg = material.temp_liquidus_pyrolite(inputs["pres_core_seg"])

                if specs["core_segregation"] == "on":
                    # Adjust the mantle composition until consistency is reached
                    if specs["core_seg_type"] == "forward":
                        raise NotImplementedError ("The 'forward' option for the core segregation model is not available yet.")
                        while True:
                            # if iteration exceeds a predefined limit resample mass and pressure as well
                            if iteration > iteration_limit:
                                for key in ["pres_core_seg", "mass"]:
                                    p = random.random()

                                    if sampling[key] == "lin":
                                        value = val[0] + p * (val[1] - val[0])
                                    elif sampling[key] == "log":
                                        value = 10 ** (np.log10(val[0]) + p * np.log10(val[1] / val[0]))

                                    inputs.update({key: value})

                                # The core segregation pressure scales with the mass
                                inputs["pres_core_seg"] *= inputs["mass"] ** pres_core_seg_exponent

                                # compute core segregation temperature
                                temp_core_seg = material.temp_liquidus_pyrolite(inputs["pres_core_seg"])
                                iteration = 0
                                bar.text = "Creating core composition (abort and resample!)"

                            # Si and Fe content in the mantle must be adjusted until a consistent
                            # set of values is obtained
                            adjust = ["Fe_number_mantle", "Si_number_mantle"]
                            
                            # Get max and min silicon content for the iron content
                            simmin = material.silicon_number_min(inputs["Fe_number_mantle"])
                            simmax = material.silicon_number_max(inputs["Fe_number_mantle"])
                            
                            # Probe silicon number in allowed range
                            p = random.random()
                            if sampling["Si_number_mantle"] == "lin":
                                    value = simmin + p * (simmax - simmin)
                            elif sampling["Si_number_mantle"] == "log":
                                value = 10 ** (np.log10(simmin) + p * np.log10(simmax / simmin))

                            inputs.update({"Si_number_mantle": value})                          

                            # Compute iron content in core from mantle composition and
                            # bulk magnesium number
                            x_Fe_core = random.random()

                            # TODO. Figure out smart way to invert partitioning model
                            # probably have to do it using brute force which sucks

                    # Adjust core composition until consistency is reached
                    elif specs["core_seg_type"] == "inverse":
                        # Adjust core composition until consistency is reached
                        # TODO. come up with a smart way to accelerate this process
                        iteration = 0
                        while True:
                            # if iteration exceeds a predefined limit resample mass and pressure as well
                            if iteration > iteration_limit:
                                
                                for key in ["pres_core_seg", "mass"]:
                                    p = random.random()
                                    val = ranges[key]
                                    if sampling[key] == "lin":
                                        value = val[0] + p * (val[1] - val[0])
                                    elif sampling[key] == "log":
                                        value = 10 ** (np.log10(val[0]) + p * np.log10(val[1] / val[0]))

                                    inputs.update({key: value})

                                # The core segregation pressure scales with the mass
                                inputs["pres_core_seg"] *= inputs["mass"] ** pres_core_seg_exponent
                                # compute core segregation temperature
                                temp_core_seg = material.temp_liquidus_pyrolite(inputs["pres_core_seg"])
                                iteration = 0
                                bar.text = f"Probing core composition (abort and resample!)"

                            # check if Si# and Fe# are in allowed range and resample if not
                            ocmf = [0.0, inputs["x_FeS"], inputs["x_FeSi"], inputs["x_FeO"]]

                            # convert to atomic abundances
                            xi = material.mat2at_core(ocmf, xiH=0.0)

                            # predict Fe and Si content in silicates according to core segregation
                            fem = material.Fe_number_mantle(inputs["pres_core_seg"], temp_core_seg, xi=xi)
                            sim = material.Si_number_mantle(inputs["pres_core_seg"], temp_core_seg, xi=xi)

                            # Get max and min silicon content for the iron content
                            simmin = material.silicon_number_min(fem)
                            simmax = material.silicon_number_max(fem)

                            margin = 1e-4
                            femmin = margin
                            femmax = xi_Fe_mantle_max * (1.0 + margin)

                            # check for inconsistencies
                            check = True
                            if sum(ocmf) > 1e0:
                                check = False

                            if fem < femmin or fem > femmax:
                                check = False

                            if sim < simmin * (1.0 + margin) or sim > simmax * (
                                1.0 - margin
                            ):
                                check = False

                            # No issues detected --> stop
                            if check:
                                break

                            # Inconsistency detected --> resample core composition
                            else:
                                # Si and O content in the core must be adjusted until a consistent
                                # set of values is obtained
                                for key in ["x_FeSi", "x_FeO"]:
                                    val = ranges[key]
                                    p = random.random()
                                    if sampling[key] == "lin":
                                        value = val[0] + p * (val[1] - val[0])
                                    elif sampling[key] == "log":
                                        value = 10 ** (np.log10(val[0]) + p * np.log10(val[1] / val[0]))

                                    inputs.update({key: value})

                            iteration += 1
                            bar.text = f"Probing core composition (Iteration {iteration})"
            
                    # get the oxide fractions in the silicates
                    oxides = material.compute_oxide_fractions(sim, fem)

                    # compute oxygen fugacity of the silicates
                    oxfug = material.logfO2(inputs["pres_core_seg"], xi[0], oxides[2])
                    inputs.update({"oxygen_fug":oxfug})
                    inputs.update({"Si_number_mantle":sim})
                    inputs.update({"Fe_number_mantle":fem})

                inputs.update({"temp_core_seg":temp_core_seg})
                
                
                for key in all_inputs.keys():
                    all_inputs[key].append(inputs[key])
                bar()

        return all_inputs

    def create_population(self):
        pass

    def create_sample(self):
        pass

    def create_mass_radius_relation(self):
        pass


class DataSet:
    def __init__(self, tag):
        self.tag = tag
        self.planets = []
        self.data = None
        self.meta = {
            "base_type": "None",
            "planetary_params": {},
            "planetary_params_ranges": {},
            "sampling_scales": {},
            "run_params": {},
        }

    def export_file(self, file_location, specs={}):

        if self.data == None:
            self.data = internal_data.create_data(self.planets)
        file_manager.write_to_csv(self.data, file_location, meta=self.meta, **specs)

    def import_file(self, file_location, overwrite=False, specs={}):
        self.data, self.meta = file_manager.read_csv(file_location, **specs)


class BlindSet(DataSet):
    """Creates a set of planets within specified ranges for the initial conditions
    without involing an iterator. This class can be used to create large data sets
    over a wide parameter space quickly as no iterations are required. However, it
    is not possible to specify boundary conditions and a uniform sampling of the
    parameter space is therefore not possible.
    """

    def __init__(self, tag="blind-1"):
        DataSet.__init__(self, tag)
        self.tag = tag
        self.ready = False

    def set_up(
        self, n, meta={}, sampling="uni",
    ):

        self.n_planets = n
        self.meta = meta

        # Screen for conflicts between planetary_params and planetary_params_ranges
        # A parameter must be either passed as a fixed parameter or a range, not both
        for key in self.meta["planetary_params"].keys():
            if key in self.meta["planetary_params_ranges"]:
                raise ValueError(
                    f"""Key '{key}' passed to planetary_params and planetary_params_ranges
                    simultaneously. Please resolve the conflict before proceeding."""
                )

        # adopt same sampling strategy for all parameters
        if len(list(self.meta["sampling_scales"].items())) == 0:
            # perform uniform sampling within ranges for all specified parameters
            if sampling == "uni":
                self.planetary_params_all = {
                    key: np.random.default_rng().uniform(val[0], val[1], self.n_planets)
                    for key, val in self.meta["planetary_params_ranges"].items()
                }

            # chose values with uniform spacings within given ranges
            elif sampling == "lin":
                self.planetary_params_all = {
                    key: np.linspace(val[0], val[1], self.n_planets)
                    for key, val in self.meta["planetary_params_ranges"].items()
                }

        # use parameter specific sampling strategies
        else:
            self.planetary_params_all = {}
            for key, val in self.meta["planetary_params_ranges"].items():
                if self.meta["sampling_scales"][key] == "lin":
                    if sampling == "uni":
                        self.planetary_params_all.update(
                            {
                                key: np.random.default_rng().uniform(
                                    val[0], val[1], self.n_planets
                                )
                            }
                        )
                    elif sampling == "lin":
                        self.planetary_params_all.update(
                            {key: np.linspace(val[0], val[1], self.n_planets)}
                        )

                elif self.meta["sampling_scales"][key] == "log":
                    if sampling == "uni":
                        self.planetary_params_all.update(
                            {
                                key: 10
                                ** np.random.default_rng().uniform(
                                    np.log10(val[0]), np.log10(val[1]), self.n_planets
                                )
                            }
                        )
                    elif sampling == "lin":
                        self.planetary_params_all.update(
                            {
                                key: np.logspace(
                                    np.log10(val[0]), np.log10(val[1]), self.n_planets
                                )
                            }
                        )

        # TODO. find a more elegant way to handle the Si# in the mantle ...
        # Sample the Si# in the mantle separately because the allowed ranges depend
        # on the Fe content of the mantle which might have been sampled as well
        Fe_number_mantle_all = self.planetary_params_all["Fe_number_mantle"]
        # Just use linear sampling by default
        # TODO. add option to chose different sampling strategies here
        Si_number_mantle_all = [
            np.random.uniform(silicon_number_min(fe), silicon_number_max(fe))
            for fe in Fe_number_mantle_all
        ]
        self.planetary_params_all.update({"Si_number_mantle": Si_number_mantle_all})
        self.ready = True

    def create(self, new=True, convergence_check=True, write_frequency=2):
        """
        Creates the planets and constructs them according to the specified properties.
        """

        write_frequency = min(write_frequency, self.n_planets)
        write_batches = int(self.n_planets / write_frequency)

        if not self.ready:
            raise AttributeError(
                "The population you are trying to create is not set up"
            )

        else:
            if new:
                self.planets = []

            planetary_params = copy.deepcopy(self.meta["planetary_params"])

            # Set up base type
            if self.meta["base_type"] == "telluric":
                planet_class = planet_creator.TelluricPlanet
            elif self.meta["base_type"] == "aqua":
                planet_class = planet_creator.AquaPlanet

            with alive_bar(
                self.n_planets,
                title=f"Creating population <{self.tag}>",
                bar="bubbles",
                spinner="pulse",
            ) as bar:

                for i in range(self.n_planets):
                    # export data batch wise according to specified frequency
                    if float(int(i / write_batches)) == i / write_batches:
                        pass

                    planetary_params.update(
                        {key: val[i] for key, val in self.planetary_params_all.items()}
                    )

                    # temporarely supress all prints
                    sys.stdout = open(os.devnull, "w")

                    pl = planet_class(
                        planetary_params=planetary_params,
                        run_params=self.meta["run_params"],
                        predictor_type="man",
                    )

                    pl.construct()
                    # pl.print()
                    # TODO. Check models and reject spurious ones
                    self.planets.append(pl)

                    sys.stdout = sys.__stdout__
                    bar()

            del planetary_params


class Population:
    """Simple class for handling populations of planets of specific base type"""

    def __init__(self, tag="pop-1", base_type="telluric"):
        self.tag = tag
        self.planets = []
        self.ready = False
        self.base_type = base_type

        if base_type == "telluric":
            self.planet_class = planet_creator.TelluricPlanet
        elif base_type == "aqua":
            self.planet_class = planet_creator.AquaPlanet
        else:
            raise ValueError("Passed unknown base tpye <{base_type}> to Population")

    def set_up(
        self,
        n,
        planetary_params_ranges={},
        run_params={},
        iterator_specs={},
        planetary_params={},
        sampling="uni",
    ):

        self.n_planets = n

        # perform uniform sampling within ranges for all specified parameters
        if sampling == "uni":
            self.planetary_params_all = {
                key: np.random.default_rng().uniform(val[0], val[1], self.n_planets)
                for key, val in planetary_params_ranges.items()
            }

        # chose values with uniform spacings within given ranges
        elif sampling == "lin":
            self.planetary_params_all = {
                key: np.linspace(val[0], val[1], self.n_planets)
                for key, val in planetary_params_ranges.items()
            }

        # chose values with logarithmic spacings within given ranges
        elif sampling == "log":
            self.planetary_params_all = {
                key: np.logspace(np.log10(val[0]), np.log10(val[1]), self.n_planets)
                for key, val in planetary_params_ranges.items()
            }

        self.run_params = run_params
        self.iterator_specs = iterator_specs
        self.planetary_params = planetary_params
        self.ready = True

    def create(self, iterator, new=True, convergence_check=True):
        if not self.ready:
            raise AttributeError(
                "The population you are trying to create is not set up"
            )

        else:
            if new:
                self.planets = []

            planetary_params = copy.deepcopy(self.planetary_params)

            with alive_bar(
                self.n_planets,
                title=f"Creating population <{self.tag}>",
                bar="bubbles",
                spinner="pulse",
            ) as bar:

                for i in range(self.n_planets):

                    planetary_params.update(
                        {key: val[i] for key, val in self.planetary_params_all.items()}
                    )

                    # temporarely supress all prints
                    sys.stdout = open(os.devnull, "w")

                    pl = self.planet_class(
                        planetary_params=planetary_params, run_params=self.run_params,
                    )

                    pl.construct()
                    iterator.iterate(planet=pl, iterator_specs=self.iterator_specs)

                    # Check convergence and reject spurious models
                    if convergence_check:
                        pl.check_convergence()
                        if pl.converged:
                            self.planets.append(pl)

                    else:
                        self.planets.append(pl)

                    sys.stdout = sys.__stdout__
                    bar()

            del planetary_params


class Sample(Population):
    """Creates populations of planets with additional constraints that cannot be
    specified by the initial or boundary conditions for the iterator. Samples can
    be used to create possible structure models for real objects for which e.g. ranges
    for the size or the MoI are known.
    """

    def __init__(self, tag="samp-1"):
        pass


class MassRadius:
    """Creates simple mass radius diagrams from populations of objects for given
    bulk compositions.
    """

    def __init__(self, tag="mrd-1", planetary_params={}):
        self.tag = tag
        self.ready = False

    def set_up(
        self,
        n,
        run_params=[],
        tags=[],
        base_types=[],
        iterator_specs=[],
        planetary_params=[],
        mass_range=[1.0, 2.0],
        sampling="lin",
    ):
        self.populations = []

        if len(tags) == 0:
            self.tags = [
                r"curve-{}".format(i + 1) for i in range(len(planetary_params))
            ]
        else:
            self.tags = tags

        # Set up the individual curves as populations
        for pps, base_type, tag in zip(planetary_params, base_types, self.tags):
            ppr = {"M_surface_should": mass_range}
            pop = Population(tag=self.tag, base_type=base_type)
            pop.set_up(
                n, planetary_params_ranges=ppr, planetary_params=pps, sampling=sampling
            )
            self.populations.append(pop)

        self.ready = True

    def add_population(self, **kwargs):
        if "tag" in kwargs:
            tag = kwargs[tag]
        else:
            tag = r"curve-{}".format(len(self.populations) + 1)

    def remove_population(self, tag):
        pass

    def create(self, iterator):
        if not self.ready:
            raise AttributeError(
                "The population you are trying to create is not set up"
            )

        else:
            for pop in self.populations:
                pop.create(iterator)

    def extract_data(self):
        self.data = {}
        for pop, tag in zip(self.populations, self.tags):
            dat = np.array([[pl.M_surface_is, pl.R_surface_is] for pl in pop.planets]).T
            self.data.update({tag: dat})

    def plot(self):
        plot_mr(self.data)


class SpecificObject(Sample):
    def __init__(self):
        pass
