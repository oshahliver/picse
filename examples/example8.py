"""This is an example for the basic usage of the planet_creator and planet_iterator tools.
It initiates an instance of a planetary object of base type TelluricPlanet and construct the
interior model for a set of boundary conditions.
"""

from picse.interiors import planet_iterator
from picse.interiors import planet_creator
from picse.physicalparams import m_earth, r_earth, G
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
    "M_surface_should": 1.3,  # desired total mass (in earth masses)
    "T_surface_should": 300,  # desired surface temperature (in kelvin)
    "P_surface_should": 1e5,  # desired surface pressure (in Pa)
    "Mg_number_should": 0.5,  # desired bulk magnesium number
    "ener_tot_should": -1.4,  # total planetary energy
    "Fe_number_mantle": 0.0,  # iron number of the silicates
    "Si_number_mantle": 0.4,  # silicon number of the silicates
    "contents": [[2], [2], [4, 5], [6, 7]],  # composition of each layer
}

# set up specifications for the iterator
# parameters that are not specified will be assigned a
# default value for the corresponding base type
iterator_specs = {
    "what": ["M_surface", "ener_tot"],  # --> target properties
    "how": ["P_center", "T_center"],  # --> adjustable properties
    "val_should": [
        planetary_params["M_surface_should"] * m_earth,
        planetary_params["ener_tot_should"],
    ],  # --> target values
    "all_val_should_weights": [
        "log",
        "lin",
    ],  # --> log or lin extrapolation for targets
    "all_howval_weights": ["exp", "lin"],  # --> exp or lin prediction for adjustables
    "acc": [1e-4, 1e-3],  # --> desired relative accuracies
    "iterationLimit": 20,  # --> max. number of iterations
    "deltaType": 0,  # --> mode for initial adjustment of adjustables
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

#######################################################################
# Model inspection
#######################################################################

# print fundamental planeatary properties to standard output
pl.print()

# You can also access individual parameters as attributes. for instance:
# print("total radius (km):", pl.R_surface_is * 1e-3)
# print("mean density (gcc):", pl.mean_density * 1e-3)
# print(
#     "norm. moment of inertia:",
#     pl.moment_of_inertia_is / (pl.R_surface_is**2 * pl.M_surface_is),
# )
# print("total mass (m_earth):", pl.M_surface_is / m_earth)
# print("desired total mass (m_earth):", pl.M_surface_should)
# print("surface pressure (bar):", pl.P_surface_is * 1e-5)
# print("desired surface pressure (bar):", pl.P_surface_should * 1e-5)
# print("surface temperature (K):", pl.T_surface_is)
# print("desired surface temperature (K):", pl.T_surface_should)
# print("core mass fraction", pl.M_core_is / (pl.M_surface_is / m_earth))
# print("desired core mass fraction", pl.M_core_should / pl.M_surface_should)

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
