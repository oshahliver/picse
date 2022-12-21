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

whats = ["M_surface", "T_surface"]
hows = ["P_center", "T_center"]
vals_should = [m_earth, 300.0]
predictors = [
    "linear",
    "linear",

]  # Only linear predictors possible anyways, this parameter has no effect at this point
should_weights = ["log", "lin"]
how_weights = ["exp", "exp"]
accs = [1e-3, 1e-2]

iterator.iterate(
    planet=pl,
    what=whats,
    how=hows,
    val_should=vals_should,
    acc=accs,
    all_val_should_weights=should_weights,
    all_howval_weights=how_weights,
    unpredictable=True,
    deltaType=0,
    iterationLimit=20,
    update_val_should=False,
    write_learning_set=False,
)

pl.show()
print("Exit Code 0")
