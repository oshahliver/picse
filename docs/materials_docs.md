# The ```pics.materials``` Module: A collection of tools to handle the properties and equations of state of planetary building blocks

The ```pics.materials``` lets you play around with different planetary building blocks and their mixtures using adequate equations of state and simple mixing laws. Note that the actual implementation of the equations of state in the structure code uses pre-generated equation of state tables and table interpolation algorithms to boost performance.

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


| Compound             | ID  | EoS            | Pressure range (GPa) | Temperature range (K) | Literature |
|----------------------|-----|----------------|----------------------|-----------------------|------------|
| $\rm H_2O$           | 000 | FEF            |                |               | Ref 1      |
| $\rm Mg_2SiO_4$      | 001 | BM3             |                |               | Ref 2      |
| $\rm Mg_2Si_2O_6$    | 002 | BM3             |              |               | Ref 3      |
| $\rm Fe_2SiO_4$      | 003 | BM3        |              |                | Ref 4      |
| $\rm Fe_2Si_2O_6$    | 004 | BM3  |            |                | Ref 5      |
| $\rm MgO$            | 005 | MGD               |                      |                       |            |
| $\rm MgSiO_3$        | 006 | MGD               |                      |                       |            |
| $\rm FeO$            | 007 | MGD               |                      |                       |            |
| $\rm FeSiO_3$        | 008 | MGD               |                      |                       |            |
| $\rm Fe(s)$             | 009 |   Bel             |                      |                       |            |
| $\rm FeS$            | 010 |      MGD          |                      |                       |            |
| $\rm Mg(OH)_2$            | 011 |    BM3            |                      |                       |            |
| $\rm \alpha-(Mg,Fe)_2Si_4$            | 012 |   BM3             |                      |                       |            |
| $\rm \beta-(Mg,Fe)_2Si_4$            | 013 |       BM3         |                      |                       |            |
| $\rm \gamma-(Mg,Fe)_2Si_4$            | 014 |       BM3         |                      |                       |            |
| $\rm post-(Mg,Fe)SiO_3$            | 015 |         MGD       |                      |                       |            |
| $\rm SiO_2$            | 016 |        MGD        |                     |                       |            |
| $\rm CaCl_2-type \ SiO_2$            | 017 |         MGD       |                      |                       |            |
| $\rm FeSi$            | 018 |        MGD        |                      |                       |            |
| $\rm Fe(l)$            | 019 |      Bel          |                      |                       |            |
| $\rm Fe_3C$            | 020 |       MGD         |                      |                       |            |


The stated temperature and pressure ranges are only rough guidlines. $\rm Al$ and $\rm Ca$ are modelled as substitutes to $\rm Mg$ and $\rm Si$ in the mantle minerals according to [[1]](#1).

Example:

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

## The ```pics.materials.Material``` module

### The ```Material.Unit``` class

The ```Mixture``` class is a useful tool for investigating and visualizing the properties of isolated planetary building blocks and to create equation of state tables.

Note. The ```Unit``` class is a separate implementation from the corresponding routines of the structure integrator implemented in fortran.

### The ```Material.Mixture``` class

The ```Mixture``` class is a useful tool for investigating and visualizing the properties of mixtures different planetary building blocks and to create equation of state tables.

Note. The ```Mixture``` class is a separate implementation from the corresponding routines of the structure integrator implemented in fortran.

Example:

```python
from pics.meterials import material

mix_specs = {"contents": [3, 5], "fractions": [0.25, 0.75], "temp": 300.0, "pres": 2e10}
mix = material.Mixture(**mix_specs)
mix.compute()
mix.print()
```

### ```Mixture.compute()```:

Initializes a new material instance for each of the
compontents in the mixture and computes volumetric and thermodynamic properties of the individual components at the pressure and temperature conditions specified for the mixture. Uses simple mixing laws to compute the bulk properties of the mixture from the individual components. If values for the pressure or temperature are passed as function arguments the conditions will be updated and the EoS evaluated accordingly.

Parameters:

```float: T=None```: Temperature in kelvin.\
```float: P=None```: Pressure in pascal.\
```int: phase=None```: Phase region signature.


### ```Mixture.update_fractions(new_fractions)```:

Changes the relative abundances of the individual components at fixed P and T and computes the new bulk properties of the mixture.

Parameters:

```list: new_fractions```: list containing the new molar abundances of each component of the mixture.

### ```Mixture.update_weight_fractions()```:

Computes the mass fraction of each component of the mixture from their mole fractions and the molar masses.

### ```Mixture.update()```:

This method updates the material instance of each component of the
mixture individually and then computes the new mean density in the cell
without re-initializing any component objects. If no new pressure or
temperature is passed, nothing will be done. If the density and or density derivative of the pressure are passed the individual contributions of the densities and density derivatives of the pressure of each material will be reconstructed without calling the equation of statee which is more efficient.

Parameters:

```float: P=None```: pressure in pascal.\
```float: T=None```: temperature in kelvin.\
```float: d=None```: density in kilogram per cubic meter.\
```float: dPdrho=None```: density derivative of pressure.

### ```Mixture.plot()```:

Simple plotting routine to visualize basic volumetric and thermodynamic properties of the mixture.

### ```Mixture.print()```:

Prints basic volumetric and thermodynamic properties of the mixture to standard output.

## The ```pics.materials.eos``` module

The ```eos``` module contains all implementations of the different equations of state available in PICSE and links to the corresponding euqation of state parameters for the different planetary building blocks that are stored in a separate file.

### ```eos.compute()```

Employs the EOS for the given material and computes the specified quantity at a given conditions.

Example:

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

### ```eos.P_BM()```:

Computes the pressure using the 3rd order Birch-Murnaghan equation of state.

Parameters:

```float: eta=0```: Compression ratio.\
```float: K=0```: Bulk modulus.\
```float: K0_prime=0```: Pressure derivative of bulk modulus.

### ```eos.rho_BM()```:

Computes the density using the 3rd order Birch-Murnaghan equation of state.

Parameters:

```int: ll=0```: Compound signature.\
```float: T=0```: Temperature in kelvin.\
```float: P=0```: Pressure in Pa.\
```float: dmin=None```: Lower bound for density bisection.\
```float: dmax=None```: Upper bound for density bisection.\
```float: aT=None```: First fitting parameter to thermal expansion coefficient.\
```float: bT=None```: Second fitting parameter to thermal expansion coefficient.\
```float: cT=None```: Third fitting parameter to thermal expansion coefficient.

### ```eos.dPdrho_BM()```:

Computes the density derivative of the pressure using the 3rd order Birch-Murnaghan equation of state.

### ```eos.P_MGD()```:

Computes the pressure using the Mie-Grüneisen-Debye equation of state.

### ```eos.rho_MGD()```:

Computes the density using the Mie-Grüneisen-Debye equation of state.

### ```eos.dPdrho_MGD()```:

Computes the density derivative of the pressure using the Mie-Grüneisen-Debye equation of state.

### ```eos.P_Vinet()```:

Computes the pressure using the Vinet equation of state.

### ```eos.rho_Vinet()```:

Computes the density using the Vinet equation of state.

### ```eos.dPdrho_Vinet()```:

Computes the density derivative of the pressure using the Vinet equation of state.

### ```eos.P_Bel()```:

Computes the pressure using the Belonoshko equation of state.

### ```eos.rho_Bel()```:

Computes the density using the Belonoshko equation of state.

### ```eos.dPdrho_Bel()```:

Computes the density derivative of the pressure using the Belonoshko equation of state.

### ```eos.P_VR()```:

Computes the pressure using the Vinet-Rydberg equation of state.

### ```eos.rho_VR()```:

Computes the density using the Vinet-Rydberg equation of state.

### ```eos.dPdrho_VR()```:

Computes the density derivative of the pressure using the Vinet-Rydberg equation of state.

## The ```pics.materials.eosTables``` module

Not ready!

## The ```pics.materials.Hydration``` module

Not ready!

## The ```pics.materials.MetalSilicatePartitioning``` module

Not ready!

# References

<a id="1">[1]</a> 
C. Sotin, et al. (2007).