""" TODO. Migrate all routines from the old PlanetFactory and perhaps PlanetInvestigator to this file.
This file should take care of higher level control of the modeeling workflow s.a. setting up
samples of planetary models, creating MR-relations and perform some post-processing tasks
"""

from pics.interiors import planet_creator, planet_iterator
import sys
import numpy as np
import os
import copy
from progress.bar import Bar
from progressbar import progressbar
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
    def __init__(self, tag="", base_type="telluric"):
        self.tag = tag
        self.planets = []
        self.ready = False
        self.base_type = base_type

    def set_up(
        self,
        n,
        planetary_params_ranges={},
        run_params={},
        iterator_specs={},
        planetary_params={},
    ):
        self.ready = True
        self.n_planets = n
        self.planetary_params_all = {
            key: np.random.default_rng().uniform(val[0], val[1], self.n_planets)
            for key, val in planetary_params_ranges.items()
        }

        self.run_params = run_params
        self.iterator_specs = iterator_specs
        self.planetary_params = planetary_params

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
                title=f"Creating population {self.tag}",
                bar="bubbles",
                spinner="pulse",
            ) as bar:

                for i in range(self.n_planets):

                    planetary_params.update(
                        {key: val[i] for key, val in self.planetary_params_all.items()}
                    )

                    # temporarely supress all prints
                    sys.stdout = open(os.devnull, "w")

                    # TODO. bad --> let the workbench figure out the right planet class from the label
                    if self.base_type == "telluric":
                        pl = planet_creator.TelluricPlanet(
                            planetary_params=planetary_params,
                            run_params=self.run_params,
                        )

                    elif self.base_type == "aqua":
                        pl = planet_creator.AquaPlanet(
                            planetary_params=planetary_params,
                            run_params=self.run_params,
                        )

                    pl.construct()
                    iterator.iterate(planet=pl, iterator_specs=self.iterator_specs)
                    self.planets.append(pl)

                    sys.stdout = sys.__stdout__
                    bar()

            del planetary_params


class Sample:
    def __init__(self):
        pass


class MassRadius(Sample):
    def __init__(self):
        pass
