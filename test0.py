import pdb
from pics.utils import fortplanet
import numpy as np
from pics.interiors import PlanetFort
from pics.interiors import PlanetFactory as plfac


tlk = plfac.Toolkit()


print("Creating interactive planet instance")

pl = PlanetFort.Planet(
    contents=[[2], [2, 9, 9, 9, 9], [4, 5], [6, 7]],
    fractions=[[1.0], [1.0, 0.0, 0.0, 0.0, 0.0], [0.5, 0.5], [0.5, 0.5]],
    layermasses=[0.1, 0.3, 1.0, 1.0],
    layerpres=[1e10, 1e10, 1e10, 1e10],
    layertemps=[1000, 1000, 1000, 1000],
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

pl.prt()


print("Exit Code 0")
