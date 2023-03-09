# The `pics.interiors` Module: The higher-level tools for creating planetary interior models

The true power of the PICSE package lies in the availability of a sphisticated library of higher-level tools that was specifically created for setting up, creating, and analysing planetary interior models. This library is continuously improved and updated. Here are a few examples of what the `picse.interiors` module achieves:

1. Offers the possibility to use base-type planets (e.g. telluric and aqua planets) for quick and robust creation of planets while maintaining the option to create more user-specific models using more generic classes.

2. Provides Ready-to-use classes to create populations of planets and mass-radius diagrams with ease.

3. Guarantees great flexibility for specifying planetary- and run parameters to explore a wide range of models and modelling strategies.

4. Uses Machine-learning enhanced algorithms for predicting initial conditions and matching user specified boundary conditions over a wide parameter space that boost performance and stability.

## Basic usage

Example:

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

## The `pics.interiors.planet_creator` module

### The `planet_creator.PlanetaryInputParams` and `planet_creator.PlanetaryOutputParams` classes

### The `planet_creator.RunInputParams` and `planet_creator.RunOutputParams` classes

### The `planet_creator.Planet` class

The `planet_creator.Planet` class provides the scaffold for planetary interior models in PICSE. It handles the planetary parameters and run parameters specified by the user and provides basic functionalities to set up, create and inspect interior structure models.

#### `Planet.set_values()`:

This method sets up all planetary parameters and run parameters for the subsequent structure integration. The example below shows how to manually set up a planet. This is handled internally by PICSE if the `BaseTypePlanet` class is used (recommended).

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

#### `Planet.construct()`:

This method calls the structure integrator for the specified planetary parameters and run parameters and performs the integration from the center to the surface.

Example:

```python
from pics.interiors import planet_creator

# Create instance of telluric planet using default specifications
pl = planet_creator.TelluricPlanet()

# Perform initial structure integration
pl.construct()
```

#### `Planet.update_initials()`:

#### `Planet.update_values()`:

#### `Planet.compute_ocean_depth()`:

Computes the thickness of the hydrosphere in meters and assigns it to the `Planet.ocean_depth` attribute.

#### `Planet.compute_oxide_fractions()`:

Computes the mole fractions of different oxides in the mantle from the atomic mantle composition. The planet needs to be constracted first.

Example:

```python
pl.compute_oxide_fractions()
print (pl.xi_SiO2_mantle, xi_MgO_mantle, xi_FeO_mantle)
```

#### `Planet.get_hydrosphere_props()`:

Extracts the structure and basic properties of the hydrosphere. These include: the total thickness, the layering into different phase regions, the liquid water mass fraction.

#### Planet.compute_core_mass()

#### Planet.update()

#### Planet.reset()

#### Planet.trim_profiles()

#### Planet.print()

#### Planet.plot()

#### Planet.check_convergence()

### The `BaseTypePlanet` class

### The `TelluricPlanet` class

### The `AquaPlanet` class

### The `CustomPlanet` class

## The `pics.interiors.core_creator` module

## The `pics.interiors.planet_workbench` module

### The `planet_workbench.DataSet` class

### The `planet_workbench.BlindSet` class

The `BlindSet` class generates samples of planetary models within specified parameter ranges without invoking the iterator. This allows for the creation of large data sets for which no pre-defined sets of boundary conditions must be matched with great accuracy. Such samples can, for instance, be used for machine-learning applications and creating new model calibrations for PICSE.

### The `planet_workbench.Population` class

### The `planet_workbench.SpecificObject` class

The `SpecificObject` will provide tools to create samples of planetary models that match specified properties of real objects s.a. observed the terrestrial planets and natural satellites in the Solar System or observed exoplanets.

### The `planet_workbench.MassRadius` class

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

## The `pics.interiors.planet_iterator` module

## The `pics.interiors.planet_tracker` module
