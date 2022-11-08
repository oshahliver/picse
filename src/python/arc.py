#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 11:29:54 2022

@author: os18o068
"""

from matplotlib import pyplot as plt
import matplotlib.patches as patches
import numpy as np

fig, ax = plt.subplots()
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
angles = np.linspace(np.pi / 4, np.pi - np.pi / 4)
rad = 1.
x = rad * np.cos(angles)
y = rad * np.sin(angles)

ax.plot(x, y)