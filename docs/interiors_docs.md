# The ```pics.interiors``` Module: The higher-level tools for creating planetary interior models

## Basic Usage

## The ```pics.interiors.planet_creator``` module

### Input- and output parameters

### The ```planet_creator.Planet``` class

The ```planet_creator.Planet``` class provides the scaffold for planetary interior models in PICSE. It handles the planetary parameters and run parameters specified by the user and provides basic functionalities to set up, create and inspect interior structure models.

#### ```Planet.set_values()```:

This method sets up all planetary parameters and run parameters for the subsequent structure integration.

Example:

```python
planetary_params = {"M_surface_should":1, "T_surface_should":300.0}
run_params = {"accs":[1e-3, 1e-2]}

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

pl.set_up(planetary_params = pp, run_params = rp)
```

#### ```Planet.construct()```:

This method calls the structure integrator for the specified planetary parameters and run parameters and performs the integration from the center to the surface.

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

