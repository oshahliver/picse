"""This is a test for the basic usage of the PlanetFactory Toolkits.

In this test we will create a planet using the new input parameter handler.
"""

import numpy as np
from pics.physicalparams import m_earth
from pics.interiors import planet_iterator
from pics.interiors import creator

iterator = planet_iterator.Toolkit()

pl = creator.TelluricPlanet()

pl.construct()

print("Starting iteration to match boundary conditions.")

iterator_specs = {
    "whats": ["M_surface", "T_surface"],
    "hows": ["P_center", "T_center"],
    "vals_should": [1.0, 300.0],
    "predictors": ["linear", "linear"],  # --> no effect at this point
    "should_weights": ["log", "log"],
    "how_weights": ["exp", "exp"],
    "accs": [1e-3, 1e-2],
    "iterationLimit": 20,
    "deltaType": 0,
    "unpredictable": False,
}

iterator.iterate(planet=pl, **iterator_specs)

pl.show()
print("Exit Code 0")
