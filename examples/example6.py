"""This is an example for the basic usage of the planet_workbench tools. A blind set of planets
is generated from uniform sampling within specified ranges for the central pressure and temperature,
surface pressure, and bulk composition. This data set is then used to train a simple neural
network on to predict the initial conditions for given boundary conditions.
"""

from picse.interiors import planet_workbench
from picse.physicalparams import m_earth
from picse.utils.file_tools import internal_data as id
from picse.utils.plot_tools import plot_tools
import numpy as np

# Create an instance of the workbench toolkit
# Note. The EoS tables are loaded internally by the workbench
workbench = planet_workbench.Toolkit()

#######################################################################
# Model specifications
#######################################################################

# define custom ranges for the planetary parameters
ppr = {
    "Fe_number_mantle": [0.0, 0.5],
    "T_center": [300, 3000],
    "P_center": [1e10, 3e12],
    "P_surface_should": [1e4, 1e9],
}

# define sampling strategy for each parameter
sampling_scales = {
    "Fe_number_mantle": "lin",
    "T_center": "lin",
    "P_center": "log",
    "P_surface_should": "log",
}

# define custom planetary parameters
pps = {}

#######################################################################
# Model creation and execution
#######################################################################

blind = planet_workbench.BlindSet(tag="example")

meta = {
    "base_type": "aqua",
    "planetary_params": pps,
    "planetary_params_ranges": ppr,
    "sampling_scales": sampling_scales,
    "run_params": {},
}

blind.set_up(
    400000, meta=meta, sampling="uni",
)

blind.create()

#######################################################################
# Model inspection
#######################################################################

# Write data to file
filepath = "/home/os18o068/Documents/PHD/Projects/pics_external_data/training_sets/training_1.csv"

blind.export_file(filepath, specs={"conflict": "add"})

blind2 = planet_workbench.BlindSet(tag="reading")
blind2.import_file(filepath)
print(blind2.data.head())

# # Scale data for plotting
# data["mass"] = np.log10(data["mass"] / m_earth)  # in earth masses
# data["pres_center"] = np.log10(data["pres_center"] * 1e-9)  # in GPa
# data["pres_surface"] = np.log10(data["pres_surface"] * 1e-5)  # in bar

# Create 2d scatter plot of selected parameters
# plot_tools.plot_sample(
#     data,
#     specs={"x": "mass", "y": "mg_number", "z": "ocean_mass_fraction"})

# create 3d scatter plot of selected parameters
# plot_tools.plot_sample_3d(
#     data,
#     specs={
#         "x": "temp_surface",
#         "y": "mg_number",
#         "z": "ocean_mass_fraction",
#         "c": "fe_number_mantle",
#     },
# )
