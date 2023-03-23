# The `pics.interiors` Module: The higher-level tools for creating planetary interior models

The true power of the PICSE package lies in the availability of a sphisticated library of higher-level tools that was specifically created for setting up, creating, and analysing planetary interior models. This library is continuously improved and updated. Here are a few examples of what the `picse.interiors` module achieves:

1. Offers the possibility to use base-type planets (e.g. telluric and aqua planets) for quick and robust creation of planets while maintaining the option to create more user-specific models using more generic classes.

2. Provides Ready-to-use classes to create populations of planets and mass-radius diagrams with ease.

3. Guarantees great flexibility for specifying planetary- and run parameters to explore a wide range of models and modelling strategies.

4. Uses Machine-learning enhanced algorithms for predicting initial conditions and matching user specified boundary conditions over a wide parameter space that boost performance and stability.

## Basic usage

### Create a base-type model

```python
from pics.interiors import planet_iterator, planet_creator

# Load the equations of state
planet_creator.load_eos_tables()

# Create an instance of the iterator
iterator = planet_iterator.Toolkit()

# Define planetary parameters
planet_specs = dict(M_surface_should = 1.0, T_surface_should = 300.0)

# Create and construct planet instance
pl = planet_creator.TelluricPlanet(planetary_params = planet_specs)
pl.construct()

# Pass planet to the iterator to match boundary conditions
iterator_specs = dict(acc=[0.01, 0.01])
iterator.iterate(planet=pl, iterator_specs = iterator_specs)

# Inspect the model
pl.print()
pl.plot()
```

### Create user defined model

Set up all planetary parameters and run parameters for the subsequent structure integration. The example below shows how to manually set up a planet. This is handled internally by PICSE if the `BaseTypePlanet` class is used (recommended).

Example:

```python
from pics.interiors import planet_creator

# Create Planet instance
pl = planet_creator.Planet()

# Define planetary and run parameters if desired
planetary_params = {"M_surface_should":1, "T_surface_should":300.0}
run_params = {"accs":[1e-3, 1e-2]}

# Create instances of parameter handlers
pp = planet_creator.PlanetaryInputParams()
rp = planet_creator.RunInputParams()

pp.set_default_values()
rp.set_default_values()
ls
# update planetary parameters if passed by user
for key, val in planetary_params.items():
    pp.default_values.update({key: val})

# update run parameters if passed by user
for key, val in run_params.items():
    rp.default_values.update({key: val})

# Adopt planetary and run parameters
pl.set_up(planetary_params = pp, run_params = rp)
```

### Mass radius relations

The `MassRadius` class allows the user to create mass-radius diagrams for user-defined compositions and surface conditions over a wide mass range.

```python
from pics.interiors import planet_workbench
workbench = planet_workbench.Toolkit()

mass_range = [0.5, 5.0]
pps = [
    {"T_surface_should": 300, "Mg_number_should": 0.3},
    {"T_surface_should": 300, "Mg_number_should": 0.4},
    {"T_surface_should": 300, "Mg_number_should": 0.5}
]
types = ["telluric" for i in range(len(pps))]

# Create an instance of the MassRadius class
mrd = planet_workbench.MassRadius(tag="my_mass_radius_diagram")

# Set up mass radius curves with 5 data points per curve
mrd.set_up(
    5, planetary_params=pps, base_types=types, mass_range=mass_range, sampling="log"
)

mrd.create(workbench.iterator)

# Extract the data as and plot the mass-radius diagram
mrd.extract_data()
mrd.plot()
```

### Create thermal evolution time line

First you need to import all relvant libraries and define the planetary parameters as well as the specifications for the iterator as we have done before.

```python
from picse.interiors import planet_evolution, planet_creator
from picse.physicalparams import m_earth, r_earth, year, day
import pandas as pd
import numpy as np

planet_creator.load_eos_tables()

planetary_params = {
    "M_surface_should": 1.0,
    "T_surface_should": 400,
    "P_surface_should": 1e5,
    "Mg_number_should": 0.5,
    "ener_tot_should": -1.0,
    "Fe_number_mantle": 0.0,
    "Si_number_mantle": 0.4,
    "contents": [[2], [2], [4, 5], [6, 7]],
}

iterator_specs = {
    "what": ["M_surface", "T_surface"],
    "how": ["P_center", "T_center"],
    "val_should": [
        planetary_params["M_surface_should"] * m_earth,
        planetary_params["T_surface_should"],
    ],
    "all_val_should_weights": [
        "log",
        "log",
    ],
    "all_howval_weights": ["exp", "exp"],
    "acc": [1e-5, 1e-3],
    "iterationLimit": 15,
    "deltaType": 0,
    "noisy": False,
}

```

Next you have to define the start end end points of the time evolution, the energy source term and the specifications for the evolver instance.
Here we use a 2nd order Runge-Kutta scheme with a step size adaption parameter of 1/4 and an accuracy for the total energy of 0.001% to solve the energy balance equation.

```python
start_time = 0.
end_time = 10e6 * year * day

# Define external energy source term
def source(t):
    albedo = 0.3
    F0 = 1366
    return 4 * np.pi * r_earth ** 2 * F0 * (1 - albedo) / 4 * 0.1

evolver_specs = {
    "start": start_time,
    "end": end_time,
    "eps": 0.25,
    "order": 2,
    "source": source,
    "acc": 1e-5,
    "write": True,
}
```

Now you can create a `TimeLine` instance and create the thermal evolution model using the chosen specifications for the `iterator` and the `evolver`.

```python
# Create timeline
tl = planet_evolution.TimeLine()
tl.set_up(
    iterator_specs=iterator_specs,
    planetary_params=planetary_params,
    evolver_specs=evolver_specs,
)
tl.create()

df = pd.DataFrame(tl.data)
print (df.head())
```
