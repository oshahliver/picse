# The `pics.materials` Module: A collection of tools to handle the properties and equations of state of planetary building blocks

The `pics.materials` lets you play around with different planetary building blocks and their mixtures using adequate equations of state and simple mixing laws. Note that the actual implementation of the equations of state in the structure code uses pre-generated equation of state tables and table interpolation algorithms to boost performance.

## Basic Usage

### Availbale thermal equations of state

1. 3rd order Birch-Murnaghan EoS (BM or BM3)
2. Mie-Grüneisen-Debye EoS (MGD)
3. Vinet EoS (Vinet)
4. Belonoshko EoS (Bel)
5. Vinet-Rydberg EoS (VR)
6. Free energy formulated EoS (FEF)

### Available chemical compounds

The following list contains all compounds for which an equation of state implementation is available.

| Compound                   | ID  | EoS | Pressure range (GPa) | Temperature range (K) | Reference |
| -------------------------- | --- | --- | -------------------- | --------------------- | --------- |
| $\rm H_2O$                 | 000 | FEF |                      |                       | [[1]](#1) |
| $\rm Mg_2SiO_4$            | 001 | BM3 |                      |                       | [[2]](#2) |
| $\rm Mg_2Si_2O_6$          | 002 | BM3 |                      |                       | [[2]](#2) |
| $\rm Fe_2SiO_4$            | 003 | BM3 |                      |                       | [[2]](#2) |
| $\rm Fe_2Si_2O_6$          | 004 | BM3 |                      |                       | [[2]](#2) |
| $\rm MgO$                  | 005 | MGD |                      |                       | [[2]](#2) |
| $\rm MgSiO_3$              | 006 | MGD |                      |                       | [[2]](#2) |
| $\rm FeO$                  | 007 | MGD |                      |                       | [[2]](#2) |
| $\rm FeSiO_3$              | 008 | MGD |                      |                       | [[2]](#2) |
| $\rm Fe(s)$                | 009 | Bel |                      |                       | [[2]](#2) |
| $\rm FeS$                  | 010 | MGD |                      |                       | [[2]](#2) |
| $\rm Mg(OH)_2$             | 011 | BM3 |                      |                       | [[1]](#1) |
| $\rm \alpha-(Mg,Fe)_2Si_4$ | 012 | BM3 |                      |                       | [[1]](#1) |
| $\rm \beta-(Mg,Fe)_2Si_4$  | 013 | BM3 |                      |                       | [[1]](#1) |
| $\rm \gamma-(Mg,Fe)_2Si_4$ | 014 | BM3 |                      |                       | [[1]](#1) |
| $\rm post-(Mg,Fe)SiO_3$    | 015 | MGD |                      |                       | [[1]](#1) |
| $\rm SiO_2$                | 016 | MGD |                      |                       |           |
| $\rm CaCl_2-type \ SiO_2$  | 017 | MGD |                      |                       | -         |
| $\rm FeSi$                 | 018 | MGD |                      |                       | [[2]](#2) |
| $\rm Fe(l)$                | 019 | Bel |                      |                       | -         |
| $\rm Fe_3C$                | 020 | MGD |                      |                       | -         |

The stated temperature and pressure ranges are only rough guidlines. $\rm Al$ and $\rm Ca$ are modelled as substitutes to $\rm Mg$ and $\rm Si$ in the mantle minerals according to [[3]](#3).

### Example 1:

```python
from pics.materials import material, eos

# Define the mixture specifications
mix_specs = {"contents": [3, 5], "fractions": [0.25, 0.75], "temp": 300.0, "pres": 2e10}

# Create instance of Mixture class
mix = material.Mixture(**mix_specs)

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

### Example 2:

The `Mixture` class is a useful tool for investigating and visualizing the properties of mixtures different planetary building blocks and to create equation of state tables.

Note. The `Mixture` class is a separate implementation from the corresponding routines of the structure integrator implemented in fortran.

```python
from pics.meterials import material

mix_specs = {"contents": [3, 5], "fractions": [0.25, 0.75], "temp": 300.0, "pres": 2e10}
mix = material.Mixture(**mix_specs)
mix.compute()
mix.print()
```

### Example 3:

Employ the EOS for the given material and computes the specified quantity at a given conditions.

```python
from pics.materials import eos

# Compute all properties of water at 10^8 Pa and 250 K
params = eos.compute(ll=0, P=1e8, T=250)
print(params)
```

Output:

```pyhton
# density, temperature, pressure, dTdP_S, dPdrho, phase, x_H2O, alpha_th, x_Al
(930.1989519862619, 250, 100000000.0, 1.835892480225174e-08, 10031472.1821997, 0, 0.0, 0.00013061068220559486, 0.0)

```

# References

<a id="1">[1]</a>
O. Shah, et al. (2021).
Internal water storage capacity of terrestrial planets and the effect of hydration on the M-R relation.
A&A 646, A162

<a id="2">[2]</a>
Oliver Shah, et al. (2022).
Possible Chemical Composition And Interior Structure Models Of Venus Inferred From Numerical Modelling.
ApJ 926 217

<a id="3">[3]</a>
C. Sotin, et al. (2007).
