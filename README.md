# The PICS Code

## General Information

The PICS code (Planetary Interior Composition and Structure code) was developed during my joint PhD program at the university of Bern and Zurich. It generates static, multi-layered structure and composition models of the interiors of terrestrial planets and water-worlds under the assumption of hydrostatic equilibrium and adiabatic thermal gradients. 

## Installation

Clone the git repository to a local directory via:

```
git clone https://github.com/oshahliver/PICS.git
```
Navigate to the project root ```pics``` and from there install the package from the ```setup.py``` file by typing:

```
python3 -m pip install .
```


## Basic usage


To set up a simple planetary model import the ```planet_creator``` (handles the planetary objects and their properties) and the ```planet_iterator``` (handles the matching of boundary conditions).

```python
from pics.interiors import planet_iterator, planet_creator
```

Next, load the EoS tables into memory for subsequent use during the structure integration. If you forget this step, you'll most likely encounter a segmentation fault.

```python
planet_creator.load_eos_tables() 
```

Create an instance of the iterator toolkit that takes care of matching the boundary conditions of your planets.

```python
iterator = planet_iterator.Toolkit()
```

Initialize a planetary object of a certain base type (only ```TelluricPlanet``` supported so far). If no planetary properties are passed via the ```planetary_params``` argument, default values will be used.

```python
planet_specs = dict(M_surface_should = 1.0, T_surface_should = 300.0)
pl = planet_creator.TelluricPlanet(planetary_params = planet_specs)
```
Perform initial structure integration.

```python
pl.construct()
```
Pass the planet instance to the iterator to match the boundary conditions with the desired precision. Planetary objects that are passed to the iterator must be constructed. If no iterator specifications are passed via the ```iterator_specs``` argument, a default strategy for matching the boundary conditions will be employed for the corresponding base type. The following will iteratively adjust the central values of the pressure and temperature to match the boundary conditions with relative accuracies of 1%:

```python
iterator_specs = dict(acc=[0.01, 0.01])
iterator.iterate(planet=pl, iterator_specs = iterator_specs)
```

If the iterator reached convergence you can inspect the planets properties:

```python
pl.print()
pl.plot()
```



