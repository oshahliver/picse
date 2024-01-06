![](assets/chapter_hydro.png)

# The PICSE Project

## General Information

The PICSE ("pic-see") project (Planetary Interior Compositions, Structures, and Evolution) is based on the planetary structure code I developed during my joint PhD program at the University of Bern and the University of Zürich. The general objective of the project is to provide ready-to-use tools for generating interior models of terrestrial and water-rich planets and moons. These tools can help us to bridge persisting gaps between observations and theory in the context of classification and characterization of exoplanets and solar system bodies.

PICSE generates static, multi-layered structure and composition models of the interiors of terrestrial planets and water-worlds under the assumption of hydrostatic equilibrium and adiabatic thermal gradients. Each layer is described as a multi-component chemical system that is modelled using adequate thermal equations of state and simple mixing laws. For more information about the physics of the model see [[1],[2]](#1).

## Installation

Clone the git repository to a local directory via:

```
git clone https://github.com/oshahliver/picse.git
```

Navigate to the project root. The package can be installed from the `setup.py` file as:

```
python3 -m pip install .
```

NOTE: The installation process was only tested on Ubuntu.

## Basic usage

For more detailed documentation, see [main documentation](./docs/README.md).

### Understanding the model

The current version of the PICS code uses a machine-learning enhanced algorithm for matching the boundary conditions. The model was calibrated and optimized for the following ranges for different parameters:

| Parameter                             | Calibration Range     |
| ------------------------------------- | --------------------- |
| Total mass                            | 0.01 - 6 Earth masses |
| Water-mass-fraction                   | 0 - 0.5               |
| Surface pressure                      | 0.1 - 10'000 bar      |
| Surface temperature                   | 100 - 3000 K          |
| Bulk $\rm[Mg]/[Mg + Fe]$              | 0.1 - 0.9             |
| Bulk $\rm[Si]/[Si + Mg]$              | 0.33 - 0.67           |
| Mantle $\rm[FeO]/[FeO + MgO]$         | 0 - 0.5               |
| Mantle $\rm[SiO_2]/[SiO_2 + MgO]$     | 0.33 - 0.67           |
| Mantle $\rm[CaO]/[CaO + MgO]$         | 0 - 0                 |
| Mantle $\rm[Al_2O_3]/[Al_2O_3 + MgO]$ | 0 - 0                 |
| Core $\rm[S]/[S + O + Si + H + Fe]$   | 0 - 0                 |
| Core $\rm[O]/[S + O + Si + H + Fe]$   | 0 - 0                 |
| Core $\rm[Si]/[S + O + Si + H + Fe]$  | 0 - 0                 |
| Core $\rm[H]/[S + O + Si + H + Fe]$   | 0 - 0                 |

The given ranges are only approximate. You may use the code outside of these ranges but I cannot guarantee that your models will converge properly and efficiently if you do so.

### Creating your first planet

To set up a simple planetary model import the `planet_creator` (handles the planetary objects and their properties) and the `planet_iterator` (handles the matching of boundary conditions).

```python
from pics.interiors import planet_iterator, planet_creator
planet_creator.load_eos_tables()
```

Create an instance of the iterator toolkit that takes care of matching the boundary conditions of your planets.

```python
iterator = planet_iterator.Toolkit()
```

Initialize a telluric planet for given boundary conditions:

```python
planet_specs = dict(M_surface_should = 1.0, T_surface_should = 300.0)
pl = planet_creator.TelluricPlanet(planetary_params = planet_specs)
```

Perform initial structure integration and pass the planet to the iterator to match the boundary conditions with desired accuracies.

```python
pl.construct()
iterator_specs = dict(acc=[0.01, 0.01])
iterator.iterate(planet=pl, iterator_specs = iterator_specs)
```

If the iterator reached convergence you can inspect the planets properties:

```python
pl.print()
pl.plot()
```

## Known issues

Currently the structure integration for aqua planets does not always converge properly. Convergence rate should, however, be above 90%.

...

## Coming soon (I hope!)

1. Integration of the single-stage core segregation model from [[2]](#1) into the user interface.

2. A interactive data structure with some basic analysis and visualization capabilities.

3. `SpecificObject` class to sample parameter variability within boundary conditions for specific objects s.a. known exoplanets, the terrestrial planets in the Solar System or the Moons of Jupiter and Saturn.

4. A machine-learning enhanced algorithm for handling initial conditions in order to boost performance of the code.

5. A thermal evolution model that solves the energy balance equation for planetary interiors with customizable energy sources and sinks s.a. radiogenic heating or stellar irradiation.

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
