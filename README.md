# The PICS Code

## General Information

The PICS code (Planetary Interior Composition and Structure code) was developed during my joint PhD program at the university of Bern and Zurich. It generates static, multi-layered structure and composition models of the interiors of terrestrial planets and water-worlds under the assumption of hydrostatic equilibrium and adiabatic thermal gradients. 

## Installation

### Manual Installation

For manual installation type:

``` 
make install
```

To permanently add the package path to an anaconda environment activate the corresponding environment and enter the command:

```
conda develop <path_to_so_file>
```

Alternatively:

```
python3 -m pip install .
```

## Basic usage


```python
from pics.interiors import planet_iterator
from pics.interiors import planet_creator
from pics.physicalparams import m_earth

iterator = planet_iterator.Toolkit()
##########################################################

# Initialize a telluric planet instance with the specified properties
pl = planet_creator.TelluricPlanet(planetary_params=planetary_params)

# Perform initial structure integration
pl.construct()

# Pass planet instance to iterator to match boundary conditions
# NOTE. planetary objects passed to the iterator must be constructed!
iterator.iterate(planet=pl, **iterator_specs)

# print fundamental planeatary properties to standard output
pl.print()

# Plot the radial P, T, M, and rho profiles
pl.plot(save=True, file_name="structure_profiles", format="png")

```



