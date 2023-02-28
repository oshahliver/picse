#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 1 2019

@author: oshah

"""

import numpy as np
import os
import sys
import eos
from astropy.io import ascii
from matplotlib import pyplot as plt
from PIMPrunparams import color_list

cwd = os.getcwd()
print(cwd)

path = "./brucite_1.tab"

dataFile = open(path)

# count number of lines
count = 0
for line in dataFile:
    count += 1

# read out length of P and T axis
dataFile = open(path)
c = 0
for line in dataFile:
    if c == 6:
        string = line.split(" ")
        for s in string:
            try:
                dim_x = int(s)

            except ValueError:
                pass

    elif c == 10:
        string = line.split(" ")
        for s in string:
            try:
                dim_y = int(s)

            except ValueError:
                pass
        break

    c += 1

print("dim_x =", dim_x)
print("dim_y =", dim_y)

data_dummy = np.empty([dim_x * dim_y, 3])
data = np.empty([3, dim_x, dim_y])

dataFile = open(path)
print("len =", len(data))

c = 0
for line in dataFile:
    if c >= 13:
        cc = 0
        string = line.split(" ")
        string_list = []
        for s in string:
            if not s == "" and not s == "\n":
                data_dummy[c - 13][cc] = float(s)
                cc += 1

    c += 1

for j in range(dim_y):
    for i in range(dim_x):
        # extract temperature
        data[0][i][j] = data_dummy[i + j * dim_x][0]

        # extract pressure
        data[1][i][j] = data_dummy[i + j * dim_x][1]

        # extract density
        data[2][i][j] = data_dummy[i + j * dim_x][2]

print("min =", np.nanmin(data[2]))
print("max =", np.nanmax(data[2]))

temps = [300, 600, 800, 1000, 2000]
indices = []

# find closest isotherms available in the table
for t in range(len(temps)):
    temp = temps[t]
    diff = np.empty([len(data[0])])
    for i in range(len(data[0])):
        diff[i] = abs(data[0][i][0] - temp)

    ind = np.where(diff == np.amin(diff))[0][0]
    indices.append(ind)
    temps[t] = data[0][ind][0]

# gather pressure and density values
dens = []
pres = []

for ind in indices:
    pres.append(data[1][ind])
    dens.append(data[2][ind])

plot_list = []

# plot isotherms
fig, ax = plt.subplots()

for t in range(len(temps)):
    (pl,) = ax.plot(
        pres[t] * 1.0e-4, dens[t], label=str(int(temps[t])) + " K", color=color_list[t]
    )

    if t == 0:
        plot_list.append(pl)

dens = []
pres = np.linspace(10000, 1000000, 50)

for t in range(len(temps)):
    dens.append([])
    for p in range(len(pres)):
        d = eos.Compute(what="dens", ll=11, T=temps[t], P=pres[p] * 1.0e5)[0]
        dens[t].append(d)

    (pl,) = ax.plot(pres * 1.0e-4, dens[t], linestyle="--", color=color_list[t])

    if t == 0:
        plot_list.append(pl)

legend = ax.legend(plot_list, ["perplex", "MGD"], loc=4)
ax.add_artist(legend)
ax.set_xlabel(r"$P \ [GPa]$")
ax.set_ylabel(r"$\rho \ [kg/m^3]$")
ax.set_title(r"$Mg(OH)_2$")
ax.legend()
nx = 6
ny = 10

x = np.linspace(0, 100, nx)
y = np.linspace(200, 2000, ny)

xx, yy = np.meshgrid(x, y)


# plot density map
fig, ax = plt.subplots()

xmin = np.min(x)
xmax = np.max(x)

ymin = np.min(y)
ymax = np.max(y)

major_tick_labels_x = x

major_tick_labels_y = y

major_ticks_x = np.linspace(0, 1, nx)
major_ticks_y = np.linspace(0, 1, ny)

ax.set_xticks(major_ticks_x)
ax.set_yticks(major_ticks_y)

ax.set_xticklabels(major_tick_labels_x)
ax.set_yticklabels(major_tick_labels_y)

im = ax.imshow(
    data[2], cmap="twilight_shifted", extent=[0.0, 1.0, 0.0, 1.0], origin="lower"
)

cbar = fig.colorbar(im)

cbar.ax.tick_params(size=0.0)

ax.set_ylabel(r"$T \ [K]$")
ax.set_xlabel(r"$P \ [GPa]$")
ax.set_title(r"$Mg(OH)_2  \ (Perplex)$")
cbar.set_label(r"$\rho \ [kg / m^{3}]$")
