# The ```pics.interiors``` Module: The higher-level tools for creating planetary interior models

## Basic Usage

## The ```pics.interiors.planet_creator``` module

### Input- and output parameters

### The ```planet_creator.Planet``` class

The ```planet_creator.Planet``` class provides the scaffold for planetary interior models in PICSE. It handles the planetary parameters and run parameters specified by the user and provides basic functionalities to set up, create and inspect interior structure models.

#### ```Planet.set_values()```:

This method sets up all planetary parameters and run parameters for the subsequent structure integration.

#### ```Planet.construct()```:

This method calls the structure integrator for the specified planetary parameters and run parameters and performs the integration from the center to the surface.

#### ```Planet.update_initials()```:

#### ```Planet.update_values()```:

#### ```Planet.compute_ocean_depth()```:

#### ```Planet.compute_oxide_fractions()```:

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

