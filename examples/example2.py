"""This is an example for the basic usage of the planet_creator and planet_iterator tools.
It initiates an instance of a planetary object of base type TelluricPlanet and construct the
interior model for a set of boundary conditions.
"""

from picse.interiors import planet_iterator
from picse.interiors import planet_creator
from picse.physicalparams import m_earth
import os
import numpy as np

# Load the EoS tables
planet_creator.load_eos_tables()

# Create an instance of the iterator toolkit
iterator = planet_iterator.Toolkit()

#######################################################################
# Model specifications
#######################################################################

# define planetary properties
# parameters that are not specified will be assigned a
# default value for the corresponding base type
planetary_params = {
    "M_surface_should": 1.,  # desired total mass (in earth masses)
    "T_surface_should": 1500,  # desired surface temperature (in kelvin)
    "P_surface_should": 1e5,  # desired surface pressure (in Pa)
    "Mg_number_should": 0.5,  # desired bulk magnesium number
    "Fe_number_mantle": 0.1,  # iron number of the silicates
    "Si_number_mantle": 0.4,  # silicon number of the silicates
    "contents": [[2], [2, 8, 10, 9], [4, 5], [6, 7]],  # composition of each layer
    "fractions":[[1.], [1., 0.1, 0.1, 0.1], [.5, .24564], [.5, .5]], # mole fractions of each layer
    "temperature_jumps":[0, 1500, 0, 0], # temperature discontinuities across layers
   
}

# set up specifications for the iterator
# parameters that are not specified will be assigned a
# default value for the corresponding base type
iterator_specs = {
    "what": ["M_surface", "T_surface"],  # --> target properties
    "how": ["P_center", "T_center"],  # --> adjustable properties
    "val_should": [
        planetary_params["M_surface_should"] * m_earth,
        planetary_params["T_surface_should"],
    ],  # --> target values
    "predictor": ["linear", "linear"],  # --> no effect at this point
    "all_val_should_weights": [
        "log",
        "log",
    ],  # --> log or lin extrapolation for targets
    "all_howval_weights": ["exp", "exp"],  # --> exp or lin prediction for adjustables
    "acc": [1e-4, 1e-3],  # --> desired relative accuracies
    "iterationLimit": 10,  # --> max. number of iterations
    "deltaType": 0,  # --> mode for initial adjustment of adjustables
    "unpredictable": False,  # --> no effect at this point
}

#######################################################################
# Model creation and execution
#######################################################################

# Initialize a telluric planet instance with the specified properties
pl = planet_creator.TelluricPlanet(planetary_params=planetary_params)

# Perform initial structure integration
pl.construct()

# Pass planet instance to iterator to match boundary conditions
# NOTE. planetary objects passed to the iterator must be constructed!
iterator.iterate(planet=pl, iterator_specs=iterator_specs)
print (pl.contents, pl.fractions)
#######################################################################
# Model inspection
#######################################################################

# print fundamental planeatary properties to standard output
pl.print()

# Plot the radial P, T, M, and rho profiles
file_path = os.getcwd()
pl.plot(
    file_name="structure_profiles",
    file_path=file_path,
    write_html=True,
    display=True,
    write_image=True,
    image_extension="pdf",
)
