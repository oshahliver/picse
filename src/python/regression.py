# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 11:21:42 2021

@author: shlah
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import numpy as np
import random
from scipy.optimize import least_squares, minimize


coefs1 = [ 2.82857837e+00,  4.93387120e-04,  1.05984767e+00, -5.03327217e-04,
  1.03554862e+00, -5.77481881e-04, -1.75655922e+00,  8.88152714e-04,
  1.51865221e+00, -6.39988006e-04, -1.20853459e+00,  6.74354373e-04,
 -1.56195185e+00,  8.46420145e-04,  1.69347895e+00, -9.55141986e-04]
coefs2 = [ 2.04790613, -0.90492109,  0.80672504,  1.45400064,  0.12513535,  1.89508518,
 -0.52088641, -3.15542808, 1.03541819, -1.59623995, -0.17103607,  2.29155182,
 -0.23384214,  2.06838175,  0.29755114, -2.45667181]
coefs3 = [ 1.24446208e+00,  1.08199947e-03, -2.74548003e-01, -2.56267161e-03,
 -3.95392204e-01, -2.20871518e-02, -6.05462506e-02,  5.49056131e-02,
 -8.88769928e-02, -2.37187556e-03, -4.04257214e-01,  5.85667741e-03,
  4.29119899e-01,  3.89283967e-02,  1.53155138e-02, -9.08152620e-02,
 -1.01089550e+00, -1.15070522e-03, -1.10812036e-01,  1.97747982e-03,
 -9.28827664e-01,  2.70495150e-02, -5.81399482e-01, -7.00886805e-02,
 -1.44942032e-01,  3.13011605e-03,  6.37782502e-01, -7.31332439e-03,
  2.10753157e-01, -5.61766236e-02,  4.28251310e-02,  1.26317367e-01,
  1.67226450e-01, -7.00755416e-04, -3.38753955e-01,  1.64383372e-03,
 -2.22626218e-01,  8.97383778e-03, -9.95286377e-02, -1.96545198e-02,
 -1.01776309e-02,  9.99709589e-04, -9.24566196e-02, -2.35511522e-03,
 -1.11178921e-01, -1.38097285e-02, -4.12423482e-02,  3.03780924e-02,
 -1.82729304e-01,  9.01886665e-04,  4.17574947e-01, -2.06982540e-03,
  5.57128973e-01, -1.36589844e-02,  3.29809943e-01,  2.82476046e-02,
 -8.49323131e-02, -1.29489658e-03,  4.18578458e-01,  2.98966457e-03,
  4.53700406e-01,  2.13869177e-02,  2.82532058e-01, -4.45463017e-02]

coefs = [coefs1, coefs2, coefs3]

def for_recursive(c = 0, iter_list = [], f = None, ranges = None,
                  count = 0, n = None, **kwargs):

    if iter_list == []:
        iter_list = [0]*n

    if c == n - 1:
        for iter_list[c] in ranges[c]:
            f(iter_list = iter_list, **kwargs)
            
    else:
        for iter_list[c] in ranges[c]:
            for_recursive(c = c + 1,
                          iter_list = iter_list, 
                          ranges = ranges,
                          f = f,
                          count = count,
                          n = n,
                          **kwargs)


def fct(theta, x):
    #print ('x =', x[:,0])
    return theta[0] + theta[1] * x[1] + theta[2] * x[0] + theta[3] * x[0] * x[1]


def driver(theta, x):
    diff = np.sum( (np.array(fct(theta, x)) - y_data) ** 2 )
    return diff


def function(theta):
    return driver(theta, x_data)


class multi_linear_model():
    def  __init__(self, n_params = 1):
        self.result = []
        self.result_dummy = 0.
        self.count = 0
        self.coefs = []
        self.n_params = n_params
    
    
    def set_data(self, x, y):
        self.x_data = x
        self.y_data = y
    
    
    def set_params(self, coefs = []):
        self.coefs = coefs
        
        if len(coefs) == 0:
            self.coefs = []
            for i in range(self.n_params + 1):
                self.coefs.append(1.)
                
            for i in range(2**self.n_params - self.n_params - 1):
                self.coefs.append(0.)
                
        else:
            self.n_params = int(np.log2(len(self.coefs)))
                

    def driver_function(self, a, x):
        self.set_params(coefs = a)
        self.evaluate(x)
        
        diff = np.sum( (np.array(self.result) - self.y_data) ** 2 )
        
        return diff


    def function(self, theta):
        return self.driver_function(theta, self.x_data)


    def fit(self, theta0):
        """
        Use scipy minimize routine to find the best fit coefficients to the data
        """
        ress = minimize(self.function, theta0)
        return ress

        
    def evaluate(self, x):
        self.result = []
        #Loop over number of params
        for i in range(len(x[0])):
            self.result_dummy = 0.
            self.count = 0
            ranges = [range(0, 2) for xx in range(self.n_params)]
            for_recursive(f = self.test_function, 
                          ranges = ranges, 
                          n = self.n_params, 
                          x = x.T[i])
            
            self.result.append(self.result_dummy)
    
    
    def plot(self):
        if self.n_params > 2:
            print ('WARNING: cannot plot for n_params > 2!')
            
            
        elif self.n_params == 1:
            pass
        
        else:
            
            fig, ax = plt.subplots()
            
            N = 2**6
            
            x1 = np.linspace(min(self.x_data[0]), max(self.x_data[0]), N)
            x2 = np.linspace(min(self.x_data[1]), max(self.x_data[1]), N)
            
            xx, yy = np.meshgrid(x1, x2)
            zz = np.concatenate((xx.reshape(-1,1), yy.reshape(-1,1)), axis = 1)
            self.evaluate(zz.T)
            
            dx1 = max(self.x_data[0]) - min(self.x_data[0])
            dx2 = max(self.x_data[1]) - min(self.x_data[1])
            
            zzz = np.array(self.result).reshape(len(x1), len(x2))
            ax.imshow(zzz, origin = 'lower',
                      aspect = 'equal',
                      extent = [min(self.x_data[0]), max(self.x_data[0]), 
                                min(self.x_data[1]), max(self.x_data[1])])
            
            ax.scatter(self.x_data[0], self.x_data[1], c = self.y_data,
                       transform = ax.transData)
           
            ax.set_aspect(dx1/dx2)
            
            
    def test_function(self, iter_list = None, x = []):
        
        res = 1.
        for i in range(len(iter_list)):
            res *= x[i]**iter_list[i]
        res *= self.coefs[self.count]
        self.count += 1
        self.result_dummy += res
    