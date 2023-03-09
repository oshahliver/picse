""" TODO. Migrate all routines from the old PlanetFactory and perhaps PlanetInvestigator to this file.
This file should take care of higher level control of the modeeling workflow s.a. setting up
samples of planetary models, creating MR-relations and perform some post-processing tasks
"""

from picse.interiors import planet_creator, planet_iterator
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

from alive_progress import alive_bar

planet_creator.load_eos_tables()


class Toolkit:
    def __init__(self):
        self.iterator = planet_iterator.Toolkit()

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
        self,
        n,
        meta={},
        sampling="uni",
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
                        planetary_params=planetary_params,
                        run_params=self.run_params,
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
