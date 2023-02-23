"""This is an example for the basic usage of the planet_creator and planet_iterator tools.
It initiates an instance of a planetary object of base type TelluricPlanet and construct the
interior model for a set of boundary conditions.
"""

from pics.interiors import planet_iterator
from pics.interiors import planet_creator
from pics.physicalparams import m_earth
import os
import numpy as np

# Load the EoS tables
planet_creator.load_eos_tables()

# Create an instance of the iterator toolkit
iterator = planet_iterator.Toolkit()

#######################################################################
# Model specifications
#######################################################################

# set up specifications for the iterator
# parameters that are not specified will be assigned a
# default value for the corresponding base type
iterator_specs = {
    "acc": [1e-4, 1e-3],  # --> desired relative accuracies
}

planetary_params = {"ocean_fraction_should": np.log10(0.3)}
#######################################################################
# Model creation and execution
#######################################################################

# Initialize a telluric planet instance with the specified properties

pl = planet_creator.AquaPlanet(planetary_params=planetary_params)

# Perform initial structure integration
pl.construct()

# Pass planet instance to iterator to match boundary conditions
# NOTE. planetary objects passed to the iterator must be constructed!
iterator.iterate(planet=pl, iterator_specs=iterator_specs)

# Check convergence
pl.check_convergence()
print(pl.converged)
#######################################################################
# Model inspection
#######################################################################

# print fundamental planeatary properties to standard output
pl.print()

# Plot the radial P, T, M, and rho profiles
# file_path = os.getcwd()
# pl.plot(
#     file_name="structure_profiles",
#     file_path=file_path,
#     write_html=True,
#     display=True,
#     write_image=True,
#     image_extension="pdf",
# )
