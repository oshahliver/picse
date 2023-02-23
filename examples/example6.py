"""This is an example for the basic usage of the planet_workbench tools. A blind set of planets
is generated from uniform sampling within specified ranges for the central pressure and temperature,
surface pressure, and bulk composition. This data set is then used to train a simple neural
network on to predict the initial conditions for given boundary conditions.
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
    "M_surface_should": [0.1, 5.],
    "Mg_number_should": [0., 1.],
    "T_center": [500, 1000],
    "P_center": [1e11, 1e12],
}

# define custom planetary parameters
pps = {}

#######################################################################
# Model creation and execution
#######################################################################

blind = planet_workbench.BlindSet(tag="example", base_type="telluric")
blind.set_up(100, planetary_params_ranges=ppr, planetary_params=pps, sampling="uni")
blind.create()

#######################################################################
# Model inspection
#######################################################################

for pl in blind.planets:
    print(
        "mass / temp surf / Mg#: {} /  {} / {}".format(
            round(pl.M_surface_is / m_earth, 3),
            round(pl.T_surface_is, 1),
            round(pl.Mg_number_is, 6),
        )
    )
