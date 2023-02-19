"""This is an example for the basic usage of the planet_workbench tools. A mass-radius diagram
is generated between 1 and 2 earth masses for different types of planets.
"""

from pics.interiors import planet_workbench
from pics.physicalparams import m_earth

# Create an instance of the workbench toolkit
# Note. The EoS tables are loaded internally by the workbench
workbench = planet_workbench.Toolkit()

#######################################################################
# Model specifications
#######################################################################

# define the mass range in earth masses
mass_range = [1., 2.]

# define custom planetary parameters
types = ["telluric", "aqua"]
pps = [{"T_surface_should": 300} for i in range(len(types))]


#######################################################################
# Model creation and execution
#######################################################################

mrd = planet_workbench.MassRadius(tag="example")
mrd.set_up(10, planetary_params=pps, base_types=types)
mrd.create(workbench.iterator)

#######################################################################
# Model inspection
#######################################################################

