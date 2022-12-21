"""This is a test for the basic usage of the PlanetFactory Toolkits.

In this test we will create a planet using the new input parameter handler.
"""

from pics.interiors import planet_iterator
from pics.interiors import planet_creator
from pics.physicalparams import m_earth

iterator = planet_iterator.Toolkit()

#######################################################################
# Model specifications
#######################################################################

# define planetary properties
# parameters that are not specified will be assigned a
# default value for the corresponding base type
planetary_params = {
    "M_surface_should": m_earth,  # desired total mass (in kg)
    "T_surface_should": 300.0,  # desired surface temperature (in kelvin)
    "P_surface_should": 1.0,  # desired surface pressure (in Pa)
    "Mg_number_should": 0.5,  # desired bulk magnesium number
    "Fe_number_mantle": 0.1,  # iron number of the silicates
    "Si_number_mantle": 0.4,  # silicon number of the silicates
    "contents": [[2], [2, 9, 10, 8, 9], [4, 5], [6, 7]],  # composition of each layer
}

# set up specifications for the iterator
iterator_specs = {
    "what": ["M_surface", "T_surface"],  # --> target properties
    "how": ["P_center", "T_center"],  # --> adjustable properties
    "val_should": [
        planetary_params["M_surface_should"],
        planetary_params["T_surface_should"],
    ],  # --> target values
    "predictor": ["linear", "linear"],  # --> no effect at this point
    "all_val_should_weights": [
        "log",
        "log",
    ],  # --> log or lin extrapolation for targets
    "all_howval_weights": ["exp", "exp"],  # --> exp or lin prediction for adjustables
    "acc": [1e-3, 1e-2],  # --> desired relative accuracies
    "iterationLimit": 20,  # --> max. number of iterations
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
iterator.iterate(planet=pl, **iterator_specs)

#######################################################################
# Model inspection
#######################################################################

# Check if convergence was reached (useful for debugging)
pl.check_convergence()
print("converged:", pl.converged)

# print fundamental planeatary properties to standard output
pl.print()
print("contents =", pl.contents)

# You can also access individual parameters as attributes
print("total radius (km):", pl.R_surface_is * 1e-3)
print("mean density (gcc):", pl.mean_density * 1e-3)
print(
    "norm. moment of inertia:",
    pl.moment_of_inertia_is / (pl.R_surface_is**2 * pl.M_surface_is),
)
print("total mass (m_eart):", pl.M_surface_is / m_earth)
print("desired total mass (m_earth):", pl.M_surface_should)
print("surface pressure (bar):", pl.P_surface_is*1e-5)
print("desired surface pressure (bar):", pl.P_surface_should*1e-5)
print("surface temperature (K):", pl.T_surface_is)
print("desired surface temperature (K):", pl.T_surface_should)
print("core mass fraction", pl.M_core_is / (pl.M_surface_is / m_earth))
print("desired core mass fraction", pl.M_core_should / pl.M_surface_should)
