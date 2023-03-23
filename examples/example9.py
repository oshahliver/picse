"""This is an example for usage of the thermal evolution model. A telluric planet
is created for a set of boundary conditions and its thermal evolution path is
constructed by solving the energy balance equation with a specified incident
stellar flux.
"""

from picse.interiors import planet_iterator
from picse.interiors import planet_creator
from picse.interiors import planet_evolution
from picse.physicalparams import m_earth, r_earth, G
import os
import numpy as np

# Load the EoS tables
planet_creator.load_eos_tables()

# Create an instance of the iterator toolkit
iterator = planet_iterator.Toolkit()
evolver = planet_evolution.Toolkit()

#######################################################################
# Model specifications
#######################################################################

# define planetary properties
# parameters that are not specified will be assigned a
# default value for the corresponding base type
planetary_params = {
    "M_surface_should": 1.,  # desired total mass (in earth masses)
    "T_surface_should": 300,  # desired surface temperature (in kelvin)
    "P_surface_should": 1e5,  # desired surface pressure (in Pa)
    "Mg_number_should": 0.5,  # desired bulk magnesium number
    "ener_tot_should": -1.,  # total planetary energy
    "Fe_number_mantle": 0.0,  # iron number of the silicates
    "Si_number_mantle": 0.4,  # silicon number of the silicates
    "contents": [[2], [2], [4, 5], [6, 7]],  # composition of each layer
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
    "all_val_should_weights": [
        "log",
        "lin",
    ],  # --> log or lin extrapolation for targets
    "all_howval_weights": ["exp", "lin"],  # --> exp or lin prediction for adjustables
    "acc": [1e-4, 1e-4],  # --> desired relative accuracies
    "iterationLimit": 20,  # --> max. number of iterations
    "deltaType": 0,  # --> mode for initial adjustment of adjustables
}

#######################################################################
# Model creation and execution
#######################################################################

# Initialize a telluric planet instance with the specified properties
pl = planet_creator.TelluricPlanet(planetary_params=planetary_params)

# Create initial model
pl.construct()
iterator.iterate(planet=pl, iterator_specs=iterator_specs)

# Perform time evolution
evolver.evolve(planet=pl, iterator=iterator, iterator_specs=iterator_specs)

#######################################################################
# Model inspection
#######################################################################

# print fundamental planeatary properties to standard output
pl.print()
