"""This is a test for the basic usage of the PlanetFactory Toolkits.

In this test we will create a planet using the new input parameter handler.
"""

from pics.physicalparams import m_earth
from pics.interiors import planet_iterator
from pics.interiors import planet_creator
from pics.physicalparams import m_earth

iterator = planet_iterator.Toolkit()

pl = planet_creator.TelluricPlanet()

pl.construct()

print("Starting iteration to match boundary conditions.")

iterator_specs = {
    "what": ["M_surface", "T_surface"],
    "how": ["P_center", "T_center"],
    "val_should": [m_earth, 300.0],
    "predictor": ["linear", "linear"],  # --> no effect at this point
    "all_val_should_weights": ["log", "log"],
    "all_howval_weights": ["exp", "exp"],
    "acc": [1e-3, 1e-2],
    "iterationLimit": 20,
    "deltaType": 0,
    "unpredictable": False,
}

iterator.iterate(planet=pl, **iterator_specs)

pl.show()
print("Exit Code 0")
