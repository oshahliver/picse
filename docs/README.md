# Documentation for the ```pics``` package

## Basic Usage

Bla bla

# Interiors

For more detailed documentation, see [interiors](./interiors_docs.md).

## Single Planets

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

Initialize a planetary object of a certain base type (only ```TelluricPlanet``` and ```AquaPlanet``` supported so far). If no planetary properties are passed via the ```planetary_params``` argument, default values will be used.

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


### Blind sets

Blind sets are sets of planetary models for that are integrated once using a given set of initial conditions to which no iterator is applied. This allows for a very time efficient creation of large numbers of structure models over specified parameter spaces without the option to specify all boundary conditions. Blind sets are, for instance, used to create the training data for machine learning algorithms that can replace more primitive methods for predicting initial conditions for given sets of boundary conditions.

```python
from pics.interiors import planet_workbench
workbench = planet_workbench.Toolkit()
```

Create the meta data for the blind set:

```python
planetary_params_ranges = {
    "Fe_number_mantle": [0.0, 0.5],
    "T_center": [300, 3000],
    "P_center": [1e10, 3e12],
    "P_surface_should": [1e4, 1e9],
}

sampling_scales = {
    "Fe_number_mantle": "lin",
    "T_center": "lin",
    "P_center": "log",
    "P_surface_should": "log",
}

meta = {"base_type":"aqua",
    "planetary_params_ranges": planetary_params_ranges,
    "sampling_scales": sampling_scales,
    "run_params":{}
}
```

Initialize, set up, and create a blind set with 1000 planets uniformely sampled from the given parameter ranges:

```python
blindset = planet_workbench.BlindSet(tag="my_blindset")
blindset.set_up(1000, meta = meta, sampling="uni")
blindset.create()
```

Export the main planetary data as a csv file:

```python
filepath = "your/file/path"
blindset.export_file(filepath)
```

Import planetary data from an existing file:

```python
new_blindset = planet_workbench.BlindSet(tag="another_blindset")
new_blindset.import_file(filepath)
print(new_blindset.data.head())
```


### Populations

You can create and manipulate planetary objects on a higher level using the ```planet_workbench``` module. A useful tool is the ```Population``` class that allows you to create populations of planets based on some overall rules and parameter ranges.

```python
from pics.interiors import planet_workbench
workbench = planet_workbench.Toolkit()
```

Create a ```Population``` instance.

```python
pop = planet_workbench.Population()
```

Set up a population of 20 telluric planets uniformely sampled between 1 and 2 Earth masses.

```python
ppr = dict(M_surface_should = [1.0, 2.0])
pop.set_up(20, planetary_params_ranges = ppr)
```

Create the members of the population with the ```create``` method and using the ```iterator``` instance automatically created by the ```planet_workbench.Toolkit``` instance (Note. you may also pass another iterator instance. This allows you to create the same population with different iterators and can be useful for debugging or comparison of different iterator strategies).

```python
pop.create(workbench.iterator)
```

You can inspect members of the population by accessing the individual ```Planet``` instances:

```python
which_planet = 2
pop.planets[which_planet].print()
pop.planets[which_planet].plot()
```

### Mass-Radius diagrams

The ```MassRadius``` class is part of the ```planet_workbench``` module and provides a wrapper to the ```Population``` class to create simple mass-radius diagrams over a specified mass range and for adjustable bulk properties.

```python
from pics.interiors import planet_workbench
workbench = planet_workbench.Toolkit()
```

Define mass range and specify some planetary parameters for each mass-radius curve:

```python
mass_range = [0.5, 5.0]
pps = [
    {"T_surface_should": 300, "Mg_number_should": 0.3},
    {"T_surface_should": 300, "Mg_number_should": 0.4},
    {"T_surface_should": 300, "Mg_number_should": 0.5}
]
types = ["telluric" for i in range(len(pps))]
```

Initialize the mass-radius instance and set up the models using logarithmic spacing between five data points for each curve within the specified mass range:

```python
mrd = planet_workbench.MassRadius(tag="example")
mrd.set_up(
    5, planetary_params=pps, base_types=types, mass_range=mass_range, sampling="log"
)
```

Create the curves as individual populations using the specified iterator:

```python
mrd.create(workbench.iterator)
```

Extract the mass-radius data from the populations for subsequent investigations and plot the mass-radius diagram:

```python
mrd.extract_data()
mrd.plot()
```


# Materials

For more detailed documentation, see [materials](./materials_docs.md).

# Utils

For more detailed documentation, see [utils](./utils_docs.md).
