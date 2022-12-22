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
from pics.interiors import planet_iterator, planet_creator

iterator = planet_iterator.Toolkit()

```

Initialize a telluric planet instance. If no planetary properties are passed default values will be used.

```python
pl = planet_creator.TelluricPlanet()
```
Perform initial structure integration

```python
pl.construct()
```
Pass the planet instance to the iterator to match the boundary conditions. Planetary objects that are passed to the iterator must be constructed!

```python
iterator.iterate(planet=pl, **iterator_specs)
```

If the iterator reached convergence you can inspect the planets properties:

```python
pl.print()
pl.plot(save=True, file_name="structure_profiles", format="png")
```



