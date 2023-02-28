# The ```pics.interiors``` Module: The higher-level tools for creating planetary interior models

The true power of the PICSE package lies in the availability of a sphisticated library of higher-level tools that was specifically created for setting up, creating, and analysing planetary interior models. This library is continuously improved and updated. Here are a few examples of what the ```picse.interiors``` module has to offer:

1. Offers the possibility to use base-type planets (e.g. telluric and aqua planets) for quick and robust creation of planets while maintaining the option to create more user-specific models using more generic classes.

2. Ready-to-use classes to create populations of planets and mass-radius diagrams with ease.

3. Provides great flexibility for specifying planetary- and run parameters to explore a wide range of models and modelling strategies.

4. Machine-learning enhanced algorithms for predicting initial conditions and matching user specified boundary conditions over a wide parameter space that boost performance and stability. 



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


## The individual components

## The ```pics.interiors.planet_creator``` module

### Input- and output parameters

### The ```planet_creator.Planet``` class

The ```planet_creator.Planet``` class provides the scaffold for planetary interior models in PICSE. It handles the planetary parameters and run parameters specified by the user and provides basic functionalities to set up, create and inspect interior structure models.

#### ```Planet.set_values()```:

This method sets up all planetary parameters and run parameters for the subsequent structure integration. The example below shows how to manually set up a planet. This is handled internally by PICSE if the ```BaseTypePlanet``` class is used (recommended).

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

# update planetary parameters if passed by user
for key, val in planetary_params.items():
    pp.default_values.update({key: val})

# update run parameters if passed by user
for key, val in run_params.items():
    rp.default_values.update({key: val})

# Adopt planetary and run parameters 
pl.set_up(planetary_params = pp, run_params = rp)
```

#### ```Planet.construct()```:

This method calls the structure integrator for the specified planetary parameters and run parameters and performs the integration from the center to the surface.

Example:

```python
from pics.interiors import planet_creator

# Create instance of telluric planet using default specifications
pl = planet_creator.TelluricPlanet()

# Perform initial structure integration
pl.construct()
``` 

#### ```Planet.update_initials()```:

#### ```Planet.update_values()```:

#### ```Planet.compute_ocean_depth()```:

Computes the thickness of the hydrosphere in meters and assigns it to the ```Planet.ocean_depth``` attribute.

#### ```Planet.compute_oxide_fractions()```:

Computes the mole fractions of different oxides in the mantle from the atomic mantle composition. The planet needs to be constracted first.

Example:

```python
pl.compute_oxide_fractions()
print (pl.xi_SiO2_mantle, xi_MgO_mantle, xi_FeO_mantle)
```

#### ```Planet.get_hydrosphere_props()```:

Extracts the structure and basic properties of the hydrosphere. These include: the total thickness, the layering into different phase regions, the liquid water mass fraction.

#### Planet.compute_core_mass()

#### Planet.update()

#### Planet.reset()

#### Planet.trim_profiles()

#### Planet.print()

#### Planet.plot()

#### Planet.check_convergence()

### The ```BaseTypePlanet``` class

### The ```TelluricPlanet``` class

### The ```AquaPlanet``` class

### The ```CustomPlanet``` class

## The ```pics.interiors.core_creator``` module

## The ```pics.interiors.planet_workbench``` module

## The ```pics.interiors.planet_iterator``` module

## The ```pics.interiors.planet_tracker``` module

