# The ```pics.materials``` Module: A collection of tools to handle the properties and equations of state of planetary building blocks

The ```pics.materials``` lets you play around with different planetary building blocks and their mixtures using adequate equations of state and simple mixing laws. 

## Basic Usage

| Compound | ID  | EoS      | Pressure range (bar) | Temperature range (K) | Literature |
|----------|-----|----------|----------------------|-----------------------|------------|
| Methane  | 001 | SRK      | 1-100                | 50-200                | Ref 1      |
| Ethane   | 002 | PR       | 10-500               | 100-300               | Ref 2      |
| Propane  | 003 | RK       | 50-1000              | 150-400               | Ref 3      |
| Butane   | 004 | Redlich | 100-2000             | 200-500               | Ref 4      |
| Pentane  | 005 | Peng-Robinson | 500-5000        | 250-600               | Ref 5      |


Example:

```python
from pics.materials import material, eos

# Define the mixture specifications
mix_specs = {"contents": [3, 5], "fractions": [0.25, 0.75], "temp": 300.0, "pres": 2e10}

# Create instance of Mixture class
mix = material.Mixture(specs=mix_specs)

# Compute the mixture properties
mix.compute()

# Update the mixture properties at new conditions
mix.update(temp=350.0, pres=3e10)

# Print the mixture properties
mix.print()

# Compute all EoS parameters and only density for compound 3
all_props = eos.compute(ll = 3, pres = 2e10, temp = 300.)
dens = eos.compute(ll = 3, what = "dens", pres = 2e10, temp = 300.)

# Print the material properties
print(all_props)
```

## The```pics.materials.Material``` module

### The ```Material.Unit``` class

### The ```Material.Mixture``` class

```python
from pics.meterials import material

mix_specs = {"contents": [3, 5], "fractions": [0.25, 0.75], "temp": 300.0, "pres": 2e10}
mix = material.Mixture(specs=mix_specs)
```

## The```pics.materials.EquationsOfState``` module

### ```EquationsOfState.compute()```

## The ```pics.materials.EquationsOfStateTables``` module

## The ```pics.materials.Hydration``` module

Not ready!

## The ```pics.materials.MetalSilicatePartitioning``` module

Not ready!