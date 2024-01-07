# The `pics.interiors` Module: The higher-level tools for creating planetary interior models

The true power of the PICSE package lies in the availability of a sphisticated library of higher-level tools that was specifically created for setting up, creating, and analysing planetary interior models. This library is continuously improved and updated. Here are a few examples of what the `picse.interiors` module achieves:

1. Offers the possibility to use base-type planets (e.g. telluric and aqua planets) for quick and robust creation of planets while maintaining the option to create more user-specific models using more generic classes.

2. Provides Ready-to-use classes to create populations of planets and mass-radius diagrams with ease.

3. Guarantees great flexibility for specifying planetary- and run parameters to explore a wide range of models and modelling strategies.

## Basic usage

### Create a base-type model

```python
from picse.interiors import planet_iterator, planet_creator

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
iterator_specs = dict(acc = [0.01, 0.01])
iterator.iterate(planet = pl, iterator_specs = iterator_specs)

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
mrd = planet_workbench.MassRadius(tag = "my_mass_radius_diagram")

# Set up mass radius curves with 5 data points per curve
mrd.set_up(
    5, planetary_params = pps, base_types = types, mass_range = mass_range, sampling = "log"
)

mrd.create(workbench.iterator)

# Extract the mass-radius data and plot the mass-radius diagram
mrd.extract_data()
mrd.plot()
```
