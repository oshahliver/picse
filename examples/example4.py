"""This is an example for the basic usage of the planet_workbench tools. A mass-radius diagram
is generated between 1 and 2 earth masses for different types of planets.
"""

from picse.interiors import planet_workbench
from picse.physicalparams import m_earth
import numpy as np

# Create an instance of the workbench toolkit
# Note. The EoS tables are loaded internally by the workbench
workbench = planet_workbench.Toolkit()

#######################################################################
# Model specifications
#######################################################################

# define the mass range in earth masses
mass_range = [0.5, 5.0]

# define custom planetary parameters
types = ["aqua", "aqua", "aqua", "aqua"]
pps = [
    {"T_surface_should": 300, "ocean_fraction_should": np.log10(0.2)},
    {"T_surface_should": 300, "ocean_fraction_should": np.log10(0.3)},
    {"T_surface_should": 300, "ocean_fraction_should": np.log10(0.4)},
    {"T_surface_should": 300, "ocean_fraction_should": np.log10(0.5)},
    # {"T_surface_should": 300, "Mg_number_should": 0.5},
    # {"T_surface_should": 300, "Mg_number_should": 0.25},
]


#######################################################################
# Model creation and execution
#######################################################################

mrd = planet_workbench.MassRadius(tag="example")
mrd.set_up(
    20, planetary_params=pps, base_types=types, mass_range=mass_range, sampling="log"
)
mrd.create(workbench.iterator)
mrd.extract_data()
mrd.plot()
#######################################################################
# Model inspection
#######################################################################
