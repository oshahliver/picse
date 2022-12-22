# The PICS Project

## General Information

The PICS project (Planetary Interior Compositions and Structures) is based on the planetary structure code I developed during my joint PhD program at the University of Bern and the University of Zürich. The general objective of the project is to provide ready-to-use tools for generating interior models of planets or moons. These tools can help us to bridge persisting gaps between observations and theory in the context of classification and characterization of exoplanets and solar system bodies.

PICS generates static, multi-layered structure and composition models of the interiors of terrestrial planets and water-worlds under the assumption of hydrostatic equilibrium and adiabatic thermal gradients. Each layer is described as a multi-component chemical system that is modelled using adequate thermal equations of state and simple mixing laws. For more information about the physics of the model see [[1],[2]](#1).

## Installation

Clone the git repository to a local directory via:

```
git clone https://github.com/oshahliver/PICS.git
```
Navigate to the project root. From there you first have to build the static libraries for the extension modules (requires gfortran to be installed on the system):

```
make static
```

Then the package can be installed from the ```setup.py``` file as:
```
python3 -m pip install .
```

NOTE: The installation process was only tested on a Linux Ubuntu system.

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

Pass the planet instance to the iterator to match the boundary conditions with the desired precision. Planetary objects that are passed to the iterator must be constructed. If no iterator specifications are passed via the ```iterator_specs``` argument, a default strategy for matching the boundary conditions will be employed for the corresponding base type. The following will iteratively adjust the central values of the pressure and temperature using a simple multi-linear predictor to match the boundary conditions with relative accuracies of 1%:

```python
iterator_specs = dict(acc=[0.01, 0.01])
iterator.iterate(planet=pl, iterator_specs = iterator_specs)
```

If the iterator reached convergence you can inspect the planets properties:

```python
pl.print()
pl.plot()
```

## Known issues

...

## Features under development

...

## References
<a id="1">[1]</a> 
O. Shah, et al. (2021).
Internal water storage capacity of terrestrial planets and the effect of hydration on the M-R relation.
A&A 646, A162

<a id="2">[2]</a> 
Oliver Shah, et al. (2022).
Possible Chemical Composition And Interior Structure Models Of Venus Inferred From Numerical Modelling.
ApJ 926 217
