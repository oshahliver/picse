# The PICS Code

## General Information

The PICS code (Planetary Interior Composition and Structure code) was developed during my joint PhD program at the university of Bern and Zurich. It generates one dimensional, static structure and composition models of terrestrial planets and water-worlds.

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

### Base-type planets

#### Example 1: Telluric planets

```python

from pics.interiors import planet_iterator, planet_creator

# initialize iterator instance
iterator = planet_iterator.Toolkit()

# initialize a telluric planet object
planet = planet_creator.TelluricPlanet()

# define planetary properties
# parameters that are not specified will be assigned a
# default value for the corresponding base type
planetary_specs = {
    "M_surface_should":1., # desired total mass (in earth masses)
    "T_surface_should":300., # desired surface temperature (in kelvin)
    "P_surface_should":1., # desired surface pressure (in bar)
    "Mg_number_should":0.5, # desired bulk magnesium number
    "Fe_number_mantle":0.1, # iron number of the silicates
    "Si_number_mantle":0.4, # silicon number of the silicates
    "contents":[[2], [2, 9, 9, 9, 9], [4, 5], [6, 7]], # composition of each layer
    }

# set up specifications for the iterator
iterator_specs = {
    "whats": ["M_surface", "T_surface"], # --> target properties
    "hows": ["P_center", "T_center"], # --> adjustable initials
    "vals_should": [1.0, 300.0], # --> target values
    "predictors": ["linear", "linear"],  # --> no effect at this point
    "should_weights": ["log", "log"], # --> log or lin extrapolation for targets
    "how_weights": ["exp", "exp"], # --> exp or lin prediction for initials
    "accs": [1e-3, 1e-2], # --> desired relative accuracies
    "iterationLimit": 20, # --> max. number of iterations
    "deltaType": 0, # --> mode for initial adjustment of initials
    "unpredictable": False, # --> no effect at this point
}

# perform initial construction of the planet
planet.construct()

# perform iteration on planetary object until desired accuracy for
# all targets is reached
iterator.iterate(planet=pl, **iterator_specs)

# display overview of planetary parameters
planet.show()
```

#### Example 1: Water-Worlds

The creation of an ocean planet is equivalent to the creation of telluric planets with the total water mass fraction being an additional required boundary condition. Using the OcenWorld base-type invokes a different method for the structure iteration under the hood.

```python

from pics.interiors import planet_iterator, planet_creator
import numpy as np

# initialize iterator instance
iterator = planet_iterator.Toolkit()

# initialize a telluric planet object
planet = planet_creator.TelluricPlanet()

# define planetary properties
# parameters that are not specified will be assigned a
# default value for the corresponding base type
planetary_specs = {
    "M_surface_should":1., # desired total mass (in earth masses)
    "T_surface_should":300., # desired surface temperature (in kelvin)
    "P_surface_should":1., # desired surface pressure (in bar)
    "Mg_number_should":0.5, # desired bulk magnesium number
    "Fe_number_mantle":0.1, # iron number of the silicates
    "Si_number_mantle":0.4, # silicon number of the silicates
    "ocean_fraction_should":np.log10(0.25) # log of desired ocean mass fraction
    "contents":[[2], [2, 9, 9, 9, 9], [4, 5], [6, 7], [1]], # composition of each layer. NOTE. additional ocean layer
    }

# set up specifications for the iterator
iterator_specs = {
    "whats": ["M_surface", "T_surface"], # --> target properties
    "hows": ["P_center", "T_center"], # --> adjustable initials
    "vals_should": [1.0, 300.0], # --> target values
    "predictors": ["linear", "linear"],  # --> no effect at this point
    "should_weights": ["log", "log"], # --> log or lin extrapolation for targets
    "how_weights": ["exp", "exp"], # --> exp or lin prediction for initials
    "accs": [1e-3, 1e-2], # --> desired relative accuracies
    "iterationLimit": 20, # --> max. number of iterations
    "deltaType": 0, # --> mode for initial adjustment of initials
    "unpredictable": False, # --> no effect at this point
}

# perform initial construction of the planet
planet.construct()

# perform iteration on planetary object until desired accuracy for
# all targets is reached
iterator.iterate(planet=pl, **iterator_specs)

# display overview of planetary parameters
planet.show()
```

### Customized planets

## Advanced usage

Coming soon (I hope!)