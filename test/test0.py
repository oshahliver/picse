import pdb
import fortplanet
import numpy as np
import PlanetFort
import PlanetFactory as plfac


tlk = plfac.Toolkit()


print ("Creating interactive planet instance")

pl = PlanetFort.Planet(
    contents = [[2], [2,9,9,9,9], [4,5], [6,7]],
    fractions = [[1.], [1., 0., 0., 0., 0.], [.5, .5], [.5, .5]],
    layermasses = [.1, .3, 1., 1.],
    layerpres = [1e10, 1e10, 1e10, 1e10],
    layertemps = [1000, 1000, 1000, 1000],
    T_center = 5000,
    P_center = 3e11,
    P_surface_should = 1e5,
    T_surface_should = 300,
    Mg_number_should = .5,
    Si_number_should = .5,
    inner_core_frac = .25,
    adiabatType = 1,
    layerConstraint = [1,1,3,1,1],
    Si_number_layers = [0., 0., .4, .4]
)

pl.Construct()

pl.prt()

print ("N Mg =", pl.Mg_count)
print ("N Si =", pl.Si_count)
print ("N Fe =", pl.Fe_count)
print ("Exit Code 0")
