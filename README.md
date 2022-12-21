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

## Usage

```python

from pics.physicalparams import m_earth
from pics.interiors import planet_iterator
from pics.interiors import creator

# initialize iterator instance
iterator = planet_iterator.Toolkit()

# initialize a telluric planet object
planet = creator.TelluricPlanet()

# perform initial construction of the planet
planet.construct()

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

iterator.iterate(planet=pl, **iterator_specs)
```
