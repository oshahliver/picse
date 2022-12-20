#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 13:39:07 2022

@author: os18o068
"""

import numpy as np
import math
import functionTools as ftool
from functionTools import for_recursive


class LinearSystem:
    def __init__(
        self, params, data, target=None, param_weights=None, data_weights=None
    ):
        print("params =", params)
        if param_weights == None:
            self.param_weights = ["lin" for i in range(len(params[0]))]
        else:
            self.param_weights = param_weights

        if data_weights == None:
            self.data_weights = ["lin" for i in range(len(data[0]))]
        else:
            self.data_weights = data_weights

        self.params = params  # current values for desired output parameters
        self.data = data  # current values for input parameters
        self.res = None  # predicted values for input parameters
        self.target = target
        self.count = 0
        self.construct()

    def partial_row(self, iter_list=[], count=0, ind=0):
        # terms = [self.params[ind][i] for i in range(len(iter_list))]
        val = 1e0
        # print ('terms =', terms)
        for i in range(len(iter_list)):
            if self.param_weights[i] == "lin":
                val *= self.params[ind][i] ** iter_list[i]
            elif self.param_weights[i] == "log":
                val *= np.log10(self.params[ind][i]) ** iter_list[i]

        self.body[ind][self.count] = val
        self.count += 1

    def row(self, iter_list=[]):

        iter_ranges = [range(0, 2) for i in range(len(self.params[self.count]))]

        count = self.count
        self.count = 0

        # Compute matrix elements for the current row
        for_recursive(f=self.partial_row, ranges=iter_ranges, ind=count)

        self.count = count

        # Update the current element of the output parameter vector
        # self.b[self.count] = val

        self.count += 1

    def construct_vector(self):
        self.vec = np.empty([2 ** len(self.params[0]), len(self.params[0])])
        for i in range(len(self.data)):
            for j in range(len(self.data[i])):
                if self.data_weights[j] == "lin":
                    self.vec[i][j] = self.data[i][j]

                elif self.data_weights[j] == "exp":
                    self.vec[i][j] = np.log10(self.data[i][j])

    def target_row(self, iter_list=[], count=0, ind=0):
        # terms = [self.params[ind][i] for i in range(len(iter_list))]

        val = 1e0
        # print ('terms =', terms)
        for i in range(len(iter_list)):
            if self.param_weights[i] == "lin":
                val *= self.target[i] ** iter_list[i]
            elif self.param_weights[i] == "log":
                val *= np.log10(self.target[i]) ** iter_list[i]

        for tv in self.target_vector:
            tv[self.count] = val

        self.count += 1

    def construct_target_vector(self):
        self.target_vector = np.zeros([len(self.params[0]), 2 ** (len(self.params[0]))])

        self.count = 0
        iter_ranges = [range(0, 2) for i in range(len(self.params[0]))]
        for_recursive(f=self.target_row, ranges=iter_ranges)
        self.count += 1

    def construct(self):
        self.construct_vector()
        self.body = np.ones([2 ** len(self.params[0]), 2 ** len(self.params[0])])

        self.count = 0

        iter_ranges = [range(0, 2) for i in range(len(self.params[0]))]

        # self.result = 0.
        for_recursive(f=self.row, ranges=iter_ranges)
        self.count = 0

    def create_model(self):
        self.construct_target_vector()
        self.res = np.linalg.solve(self.body, self.vec)

    def predict(self):
        self.pred = [sum(pr) for pr in self.target_vector * self.res.T]

        for i in range(len(self.pred)):
            if self.data_weights[i] == "lin":
                pass
            elif self.data_weights[i] == "exp":
                self.pred[i] = 10 ** self.pred[i]
