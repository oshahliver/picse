""" TODO. Migrate all routines from the old PlanetFactory and perhaps PlanetInvestigator to this file.
This file should take care of higher level control of the modeeling workflow s.a. setting up
samples of planetary models, creating MR-relations and perform some post-processing tasks
"""

from pics.interiors import planet_creator, planet_iterator
import sys
import numpy as np
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

    def create(self, iterator, new=True):
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
        

    def set_up(
        self,
        n,
        run_params=[],
        tags=[],
        base_types=[],
        iterator_specs=[],
        planetary_params=[],
        mass_range=[1.0, 2.0],
    ):

        self.populations = []
        if len(tags) == 0:
            tags = [r"curve-{}".format(i+1) for i in range(len(planetary_params))]

        # Set up the individual curves as populations
        for pps, base_type, tag in zip(planetary_params, base_types, tags):
            print ("check")
            ppr = {"M_surface_should":mass_range}
            pop = Population(tag=tag, base_type=base_type)
            pop.set_up(
                n, planetary_params_ranges=ppr, planetary_params=pps, sampling="lin"
            )
            self.populations.append(pop)

    def create(self, iterator):
        for pop in self.populations:
            pop.create(iterator)

class SpecificObject(Sample):
    def __init__(self):
        pass
