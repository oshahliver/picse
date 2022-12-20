#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:11:24 2019

@author: oshah

Routine for EoS data table handling.

After initiating a Table() object use the method(load) to read in a table from
a file. This creates a data array containing the available parameter values
for all the grid points. Use the interpolate() method to use table interpolation
of desired order to interpolate between grid points at a given P and T.
"""

import numpy as np
import os
import time
from astropy.io import ascii


class Table:
    def __init__(self):
        self.material = None
        self.alpha_x = None
        self.alpha_y = None
        self.x_range = None
        self.y_range = None
        self.x_start = None
        self.y_start = None
        self.dim_decade_x = None
        self.dim_decade_y = None
        self.values = None
        self.derivatives_x = None
        self.derivatives_y = None
        self.derivatives_xy = None
        self.b = None
        self.matrix = None
        self.pairs = []
        self.x_axis = []
        self.y_axis = []
        self.dimension = [len(self.x_axis), len(self.y_axis)]
        self.grid_size = 0
        self.data_table = None  # table containing data strings
        self.data = None

    def load(self, file_name, loc="cwd"):
        """Loads existing table and makes data accessable for the current
        python session. This allows quicker data handling as the table content
        need not be read out from the file for every interpolation cycle.
        """
        t0 = time.time()
        self.x_axis = []
        self.y_axis = []
        self.values = []

        if loc == "cwd":
            dir_in = os.getcwd() + "/"

        else:
            dir_in = loc + "/"

        # read ascii file for eos table
        self.data_table = ascii.read(dir_in + file_name, format="ecsv")

        # extract metadata
        self.x_range = self.data_table.meta["x_range"]
        self.y_range = self.data_table.meta["y_range"]
        self.alpha_x = self.data_table.meta["alpha_x"]
        self.alpha_y = self.data_table.meta["alpha_y"]
        self.material = self.data_table.meta["material"]
        self.dimension = self.data_table.meta["dimension"]
        self.grid_size = self.data_table.meta["grid_size"]
        self.x_start = self.data_table.meta["x_start"]
        self.y_start = self.data_table.meta["y_start"]

        self.data = []
        # iterate over rows (x-axis)
        for i in range(len(self.data_table)):
            self.data.append([])
            # iterate over columns (y-axis)
            for j in range(len(self.data_table[i])):
                # gather values
                dat = self.data_table[i][j]
                dat = dat.split(",")
                for d in range(len(dat)):
                    try:
                        dat[d] = float(dat[d])

                    except ValueError:
                        pass

                self.data[i].append(dat)

        self.x_axis = [dat[0][0] for dat in self.data]
        self.y_axis = [dat[1] for dat in self.data[0]]
        print("elapsed time for loading:", time.time() - t0)

    def get_indices(self, x, alpha, start, order=1):
        """This function takes a x or y value as input and gives the indices of
        the neighbouring grid points (left and right) along the corresponding
        grid axis in the eos table. The number of neighbours that are extracted
        depends on the order of the interpolatino polynomial.
        """
        relative_coordinates = [
            [0, 1],  # 1st order
            [-1, 0, 1],  # 2nd order
            [-1, 0, 1, 2],  # 3rd order
            [-2, -1, 0, 1, 2],  # 4th order
            [-2, -1, 0, 1, 2, 3],  # 5th order
            [-3, -2, -1, 0, 1, 2, 3],  # 6th order
            [-3, -2, -1, 0, 1, 2, 3, 4],  # 7th order
            [-4, -3, -2, -1, 0, 1, 2, 3, 4],
        ]  # 8th order

        # compute order of magnitude of input temperature
        exponent = int(np.log10(x))

        # compute scaling factor
        a = 10 ** (exponent - alpha)

        # compute the closest grid point values
        left_value = round((x / a - 0.5), 0) * a
        right_value = round((x / a + 0.5), 0) * a

        # compute distance from left and right grid points to given value
        left_delta = abs(x - left_value)
        right_delta = abs(x - right_value)

        # compute index of left and right values
        left_index = int(left_value / a - 10 ** alpha)
        right_index = left_index + 1  # int(right_value/a+10**alpha)

        # NOTE: if the given value coincides with a existing grid point, the
        # above algorithm might not take the corresponding index as the core
        # point. But since the grid point which coincides with the given value
        # will be taken into account anyways, the interpolation should still yield
        # the precise result for this point

        # take closest grid point as core point
        offset = int(left_value / a - 10 ** alpha)

        # gather indices of all required neighbouring grid points around the
        # core point for n'th order interpolation scheme
        ind = [
            9 * 10 ** alpha * (exponent - start) + offset + r
            for r in relative_coordinates[order - 1]
        ]

        # print ('indices=',ind)
        # print ('values=',[self.pres_axis[i] for i in ind])
        return ind

    def gather_pairs(self, x, y, order=2):
        res = []
        for i in range(len(x)):
            for j in range(len(y)):
                res.append([x[i][j], y[i][j]])
        return res

    def row(self, pnt, order=2):
        """Compute row of coefficient matrix 'M' for 2nd order 2d interpolation"""
        x, y = pnt

        r = []
        for i in range(order + 1):
            for j in range(order + 1):
                r.append(x ** i * y ** j)

        return r

    def construct_matrix(self, x, y, order=1):
        """Construct coefficient matrix for the points (x[i], y[i])"""
        xx, yy = np.meshgrid(x, y)

        self.pairs = self.gather_pairs(xx, yy, order=order)
        self.matrix = [self.row(p, order=order) for p in self.pairs]

    def interpolate(
        self,
        x=None,
        y=None,
        order=1,
        phase_check=True,
        param=2,
        order_check=False,
        warnings=True,
        **kwargs
    ):
        """Takes x=temp and y=pres values as input and gives f(x, y) as output
        using a polynomial nth order interpolation scheme. The scheme solves the
        matrix equation M*a=b where 'M' is the coefficient matrix and 'b'
        the solution vector.
        """
        terminate = False

        if x < self.x_axis[order - 1] or x > self.x_axis[-order]:
            print(
                "WARNING: x value must be in the range",
                self.x_axis[order - 1] + ",",
                self.x_axis[-order] + ". Given:",
                x,
            )
            terminate = True

        if y < self.y_axis[order - 1] or y > self.y_axis[-order]:
            print(
                "WARNING: y value must be in the range",
                self.y_axis[order - 1] + ",",
                self.y_axis[-order] + ". Given:",
                y,
            )

            terminate = True

        if not terminate:
            # Compute the relevant grid indices for the xy coordinates
            ind = [
                self.get_indices(x, self.alpha_x, self.x_start, order=order),
                self.get_indices(y, self.alpha_y, self.y_start, order=order),
            ]

            phase = self.data[ind[0][0]][ind[1][0]][3]

            # gather pT grid points for the coefficient matrix
            try:
                x_list = [self.x_axis[i] for i in ind[0]]
                y_list = [self.y_axis[j] for j in ind[1]]

            except IndexError:
                print(x, y, ind)

            # compute the solutions at the given grid points
            self.b = []
            for i in range(order + 1):
                for j in range(order + 1):
                    # note that here j and i must be interchanged in comparison
                    # to the interpolate_deriv methode, otherwise it doesn't work
                    val = self.data[ind[0][j]][ind[1][i]][param]
                    if not val == "None":
                        self.b.append(val)

                    else:
                        self.b.append(0.0)

            # if grid shift failed and the grid points still spread over
            # multiple phases, just switch to linear interpolation to avoid
            # any kind unphysical results. Of course this is at the cost of
            # accuracy but it only happens for a few grid points along
            # phase transitions
            if order_check:
                if min(self.b) < 0.25 * max(self.b):
                    # print ('Reduziere zu linearer Interpolation')
                    # print ('b-Vektor=',self.b)
                    order = 1

                    # Compute the relevant grid indices for the xy coordinates
                    ind = [
                        self.get_indices(x, self.alpha_x, self.x_start, order=order),
                        self.get_indices(y, self.alpha_y, self.y_start, order=order),
                    ]

                    # gather pT grid points for the coefficient matrix
                    x_list = [self.x_axis[i] for i in ind[0]]
                    y_list = [self.y_axis[j] for j in ind[1]]

                    # compute the solutions at the given grid points
                    self.b = []
                    for i in range(order + 1):
                        for j in range(order + 1):
                            # note that here j and i must be interchanged in comparison
                            # to the interpolate_deriv methode, otherwise it doesnt work
                            val = self.data[ind[0][j]][ind[1][i]][param]
                            self.b.append(val)

            # compute coefficient matrix using the xy-grid points
            self.construct_matrix(x=x_list, y=y_list, order=order)

            # compute the coefficients for the polynom fit
            try:
                self.a = np.linalg.solve(self.matrix, self.b)

            except TypeError:
                print(self.b, self.matrix)
                print(x_list, y_list)

            # compute polynom solution at the given (P, T) point
            result = sum(
                [
                    self.row([x, y], order=order)[i] * self.a[i]
                    for i in range(len(self.a))
                ]
            )

            # print ('Elapsed time:', round(1.0e3*(time.time()-t0), 2), 'ms')
            return result
