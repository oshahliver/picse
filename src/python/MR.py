import sys
from matplotlib import pyplot as plt
import PlanetFactory as plfac
from PIMPphysicalparams import m_earth
import time

N = int(sys.argv[1])
eps_r = float(sys.argv[2])
T_surf = float(sys.argv[3])
P_surf = float(sys.argv[4])
Mg_numbers = [0.9]

population = []

toolkit=plfac.Toolkit()

t0 = time.time()

for Mg in Mg_numbers:
	print ("\nProcessing Mg# =", Mg)
	mc, devs = toolkit.model_curve(N=N, T_surface_should=T_surf, P_surface_should=P_surf, eps_r=eps_r, Mg_number_should=Mg, M_start=0.1, M_end=4.0, sweeps=25, resweeps=0, modelType=1)


	mc = toolkit.sort_planets(planets=mc, what='M_surface_is')
	population.append(mc)

toolkit.write_population(population=population, loc="./MR_curve_T"+str(int(T_surf))+'test/')
#toolkit.plot_MR(planets=population, save=True)

t=time.time()
print ("elapsed time for model_curve: ", t-t0, "s")

