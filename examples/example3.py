"""This is an example for the basic usage of the planet_workbench tools. A population of planets
is generated from uniform sampling within specified ranges for some of the planetary parameters
"""

from pics.interiors import planet_workbench
from pics.physicalparams import m_earth

# Create an instance of the workbench toolkit
# Note. The EoS tables are loaded internally by the workbench
workbench = planet_workbench.Toolkit()

#######################################################################
# Model specifications
#######################################################################

# define custom ranges for the planetary parameters
ppr = {
    "M_surface_should": [1.0, 2.0],
}

# define custom planetary parameters
pps = {"T_surface_should": 356.0}

#######################################################################
# Model creation and execution
#######################################################################

pop = planet_workbench.Population(tag = "example", base_type = "telluric")
pop.set_up(20, planetary_params_ranges=ppr, planetary_params=pps)
pop.create(workbench.iterator)

#######################################################################
# Model inspection
#######################################################################

for pl in pop.planets:
    print(
        "total mass / surface temp: {} / {}".format(
            round(pl.M_surface_is / m_earth, 3), round(pl.T_surface_is, 1)
        )
    )
