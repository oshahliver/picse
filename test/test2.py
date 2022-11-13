"""This is a test for the basic usage of the PlanetFactory Toolkits.

In this test we will create a planet using the new input parameter handler.
"""

import fortplanet
import numpy as np
from PIMPphysicalparams import m_earth
#import planet_iterator
#import PlanetFort
import kit_planet

pl = kit_planet.TelluricPlanet()
pl.show()


print ("Exit Code 0")