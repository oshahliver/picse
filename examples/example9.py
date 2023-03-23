"""This is an example for usage of the thermal evolution model. A telluric planet
is created for a set of boundary conditions and its thermal evolution path is
constructed by solving the energy balance equation with a specified incident
stellar flux.
"""

from picse.interiors import planet_iterator
from picse.interiors import planet_creator
from picse.interiors import planet_evolution
from picse.physicalparams import m_earth, r_earth, G, year, day
import os
from picse.utils.plot_tools import plot_tools
import numpy as np
import pandas as pd

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
    "M_surface_should": 1.0,  # desired total mass (in earth masses)
    "T_surface_should": 400,  # desired surface temperature (in kelvin)
    "P_surface_should": 1e5,  # desired surface pressure (in Pa)
    "Mg_number_should": 0.5,  # desired bulk magnesium number
    "ener_tot_should": -1.0,  # total planetary energy
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
        "log",
    ],  # --> log or lin extrapolation for targets
    "all_howval_weights": ["exp", "exp"],  # --> exp or lin prediction for adjustables
    "acc": [1e-5, 1e-3],  # --> desired relative accuracies
    "iterationLimit": 15,  # --> max. number of iterations
    "deltaType": 0,  # --> mode for initial adjustment of adjustables
    "noisy": False,  # --> disable std outputs
}

end_time = 10e6 * year * day  # end time of integration
albedo = 0.3  # planetary albedo
L_in = (
    4 * np.pi * r_earth ** 2 * 1366 * (1 - albedo) / 4 * 0.1
)  # incident radiation power

# set up the specifications for the evolver
# parameters that are not specified will be assigend a
# default value for the corresponding base type
evolver_specs = {
    "start": 0,  # --> start time
    "end": end_time,  # --> end time
    "eps": 0.25,  # --> time step refinement parameter
    "order": 2,  # --> order of integration scheme
    "source": L_in,  # --> external energy source
    "acc": 1e-5,  # --> accuracy of energy iteration
    "write": True,  # --> write timeline data to file
}

#######################################################################
# Model creation and execution
#######################################################################

# Perform time evolution
planetary_params.update({"T_surface_should": 400})
iterator_specs.update(
    {
        "val_should": [
            planetary_params["M_surface_should"] * m_earth,
            planetary_params["T_surface_should"],
        ]
    }
)
# Initialize a telluric planet instance with the specified properties
pl = planet_creator.TelluricPlanet(planetary_params=planetary_params)

# Create initial model
pl.construct()
iterator.iterate(planet=pl, iterator_specs=iterator_specs)

data, instances = evolver.evolve(
    planet=pl,
    iterator=iterator,
    iterator_specs=iterator_specs,
    evolver_specs=evolver_specs,
)
del pl

df = pd.DataFrame(data)

print(df.head())
#######################################################################
# Model inspection
#######################################################################

tl = planet_evolution.TimeLine()
tl.set_up(
    iterator_specs=iterator_specs,
    planetary_params=planetary_params,
    evolver_specs=evolver_specs,
)
tl.create()

df = pd.DataFrame(tl.data)
print (df.head())
# # Plot the time line
# file_path = os.getcwd()
# plot_tools.plot_timeline(df,
#     file_name="timeline",
#     file_path=file_path,
#     write_html=True,
#     display=True,
#     write_image=True,
#     image_extension="pdf",
# )
