"""This is a test for the planet_iterator method. This method takes an interactive planet object as input and interatively reconstructs it until a set of predifined boundary conditions are met.

In this example we create a planet for which we whish to fix the total mass to 1 M_E and the surface temperature to 300 K.
"""
from pics.utils import fortplanet
import numpy as np
from pics.interiors import PlanetFort
from pics.physicalparams import m_earth
from pics.interiors import planet_iterator


# def run():

iterator = planet_iterator.Toolkit()

print("Creating interactive planet instance.")

pl = PlanetFort.Planet(
    contents=[[2], [2, 9, 9, 9, 9], [4, 5], [6, 7]],
    fractions=[[1.0], [1.0, 0.0, 0.0, 0.0, 0.0], [0.5, 0.5], [0.5, 0.5]],
    layermasses=[0.1, 0.3, 1.0, 10.0],
    # layerpres = [1e10, 1e10, 1e10, 1e10],
    # layertemps = [1000, 1000, 1000, 1000],
    T_center=5000,
    P_center=3e11,
    P_surface_should=1e5,
    T_surface_should=300,
    Mg_number_should=0.5,
    Si_number_should=0.5,
    inner_core_frac=0.25,
    adiabatType=1,
    layerConstraint=[1, 1, 3, 1, 1],
    Si_number_layers=[0.0, 0.0, 0.4, 0.4],
)

pl.Construct()

print("Starting iteration to match boundary conditions.")

whats = ["M_surface", "T_surface"]
hows = ["P_center", "T_center"]
vals_should = [m_earth, 300.0]
predictors = [
    "linear",
    "linear",

]  # Only linear predictors possible anyways, this parameter has no effect at this point
should_weights = ["log", "lin"]
how_weights = ["exp", "exp"]
accs = [1e-3, 1e-2]

iterator.iterate(
    planet=pl,
    what=whats,
    how=hows,
    val_should=vals_should,
    acc=accs,
    all_val_should_weights=should_weights,
    all_howval_weights=how_weights,
    unpredictable=True,
    deltaType=0,
    iterationLimit=20,
    update_val_should=False,
    write_learning_set=False,
)

pl.prt()
print ("layers =", pl.layers)
# return pl
print("Exit Code 0")
