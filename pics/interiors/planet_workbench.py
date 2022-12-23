""" TODO. Migrate all routines from the old PlanetFactory and perhaps PlanetInvestigator to this file.
This file should take care of higher level control of the modeeling workflow s.a. setting up
samples of planetary models, creating MR-relations and perform some post-processing tasks
"""

from pics.interiors import planet_creator, planet_iterator
import sys
import os
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
    def __init__(self, label=""):
        self.label = label
        self.planets = []

    def create(self, n, iterator, new = True):
        if new:
            self.planets = []

        with alive_bar(
            n, title=f"Creating population {self.label}", bar="bubbles", spinner="pulse"
        ) as bar:
            for i in range(n):
                # temporarely supress all prints
                sys.stdout = open(os.devnull, "w")
                pl = planet_creator.TelluricPlanet()
                pl.construct()
                iterator.iterate(planet=pl)
                sys.stdout = sys.__stdout__
                bar()

                self.planets.append(pl)


class Sample:
    def __init__(self):
        pass


class MassRadius(Sample):
    def __init__(self):
        pass
