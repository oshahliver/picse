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

whats = ["M_surface", "T_surface"]
hows = ["P_center", "T_center"]
vals_should = [m_earth, 300.0]
predictors = [
    "linear",
    "linear",

]  # Only linear predictors possible anyways, this parameter has no effect at this point
should_weights = ["log", "lin"]
how_weights = ["exp", "exp"]
accs = [1e-3, 1e-2]

# pass the planetary object to the iterator to match boundary conditions

iterator.iterate(
    planet=pl,
    what=whats,
    how=hows,
    val_should=vals_should,
    acc=accs,
    all_val_should_weights=should_weights,
    all_howval_weights=how_weights,
    unpredictable=True,
    deltaType=0,
    iterationLimit=20,
    update_val_should=False,
    write_learning_set=False,
)
```
