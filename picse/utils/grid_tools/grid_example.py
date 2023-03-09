#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 15:21:31 2021

@author: os18o068

Parameter keys:
    
    0: T............(surface temperature)
    1: log(P).......(surface pressure)
    2: M_ocean/M....(ocean mass fraction)
    3: R............(total radius)
    4: L............(luminosity)
    5: E_grav.......(gravitational energy)
    6: E_int........(internal energy)

"""

import picse.utils.grid_tools.grid_interpolation as grid_interpolation

# Initiate grid instance
grid = grid_interpolation.Grid()

# Import grid data from file. By default file location is cwd and file name is
#'planets_grid.tab'
grid.read_data_file()
print("\n")

# Example 1
# Compute R and L at T = 254 K, log(P/Pa) = 7.4, M_ocean/M = 0.28
result = grid.interpolate(vals=[254.0, 7.4, 0.28], params=[3, 4])
print("R =", result[0], "m")
print("L =", result[1], "W")

# Example 2
# Compute E_grav and E_int at T = 254 K, log(P/Pa) = 7.4, M_ocean/M = 0.28
result = grid.interpolate(vals=[254.0, 7.4, 0.28], params=[5, 6])
print("E_grav =", result[0], "J")
print("E_int =", result[1], "J")

# Example 3
# Plot R as function of T and P for M_ocean/M = 0.1
grid.plot_slice(axes=[0, 1], slice=0, param=3)

# Example 4
# Plot L as function of T and M_ocean for log(P/Pa) = 5
grid.plot_slice(axes=[0, 2], slice=0, param=4)
