"""This is a test for the basic usage of the PlanetFactory Toolkits.

In this test we will create a planet using the new input parameter handler.
"""

from pics.interiors import planet_iterator
from pics.interiors import planet_creator

#planet_creator.load_eos_tables()

iterator = planet_iterator.Toolkit()

planetary_params = {
    "M_surface_should": 1.0,
    "T_surface_should": 400.0,
}

iterator_specs = dict(acc=[0.005, 0.005])


pl = planet_creator.TelluricPlanet(planetary_params=planetary_params)

pl.construct()

iterator.iterate(planet=pl, iterator_specs = iterator_specs)
pl.print()