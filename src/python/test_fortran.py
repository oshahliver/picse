#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 17:41:12 2019

@author: oshah
"""

import test
import time
import numpy as np

def run(x, y):
    z = x**2 + y**2
    print (z)


def performance(nx, ny):
    t0_p = time.time()
    run(nx, ny)
    t_p=time.time()

    t0_f = time.time()
    test.run(nx, ny)
    t_f = time.time()
    print ("python runtime = ", t_p-t0_p)
    print ("fortran runtime = ", t_f-t0_f)
