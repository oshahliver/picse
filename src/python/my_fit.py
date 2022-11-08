# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 16:20:17 2021

@author: shlah
"""

import regression as regr
from scipy.optimize import least_squares, minimize
from matplotlib.transforms import blended_transform_factory
import matplotlib.pyplot as plt
import numpy as np
import random
import pickle

#Central temperature
coefs1 = [ 2.82857837e+00,  4.93387120e-04,  1.05984767e+00, -5.03327217e-04,
  1.03554862e+00, -5.77481881e-04, -1.75655922e+00,  8.88152714e-04,
  1.51865221e+00, -6.39988006e-04, -1.20853459e+00,  6.74354373e-04,
 -1.56195185e+00,  8.46420145e-04,  1.69347895e+00, -9.55141986e-04]
#Central pressure
coefs2 = [ 2.04790613, -0.90492109,  0.80672504,  1.45400064,  0.12513535,  1.89508518,
 -0.52088641, -3.15542808, 1.03541819, -1.59623995, -0.17103607,  2.29155182,
 -0.23384214,  2.06838175,  0.29755114, -2.45667181]
#Core mass fraction
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

#Central temperature, central pressure, core mass fraction
coefs = [coefs1, coefs2, coefs3]

def fct1(theta, x):
    #print ('x =', x[:,0])
    return theta[0] + theta[1] * x[1] + theta[2] * x[0] + theta[3] * x[0] * x[1]


def fct2(theta, x):
    #print ('x =', x[:,0])
    return theta[0] + theta[1] * x[1]**.9 + theta[2] * x[0]**1.05 + theta[3] * x[0] * x[1]**1.01


def driver(theta, x):
    diff = np.sum( (np.array(fct1(theta, x)) - y_data) ** 2 )
    return diff


def function(theta):
    return driver(theta, x_data)


N = 100
x_data = np.array([[i, i * .5] for i in range(N)]).T
y_data = fct1([0., 1., 1., 1.], x_data) 

'''
for i in range(len(y_data)):
    y_data[i] *= 1. + random.random()*.00001
'''

path = '/home/os18o068/Documents/PHD/Projects/Planets/Data/planet_grid_mars/'
data_all = np.load(path+'planet_data.npy')

#Split data set into different mass regimes
indices = []
for i in range(len(data_all)):
    if data_all[i][45] < 10e1:
        indices.append(i)
        
data = np.empty([len(indices), len(data_all[0])])
for i in range(len(indices)):
    data[i][:] = data_all[indices[i]][:]

#Parameters that should be predicted
predict = [20, 21, 27]
scale = ['log', 'log', 'lin']
labels = ['Central temperature',
          'Central pressure',
          'Core mass fraction']
models = []
cols = ['b', 'g', 'r', 'k']

fig, ax = plt.subplots()

for i in range(len(predict)):
    #Gather relevant data for the fits
    #Central temperature in K
    if i == 0:
        x_data = np.array([#np.log10(data[:,45]), #total mass
                           data[:,0], #Mg number
                           data[:,26], #Fe content in core
                           data[:,17], #T surface (TBL)
                           #data[:,47], #ocean mass frac
                           ])
    #Central pressure in GPa
    elif i == 1:
        x_data = np.array([#np.log10(data[:,45]), #Total mass
                           data[:,0], #Mg number
                           data[:,26], #Fe content in core
                           data[:,5], #Fe content in mantle
                           #data[:,47], #ocean mass frac
                           ])
    #Core mass fraction
    elif i == 2:
        x_data = np.array([#data[:,45], #Total mass
                           data[:,0], #Mg number
                           data[:,26], #Fe content in core
                           data[:,5],  #Fe content in mantle
                           #data[:,4], #Si number mantle
                           #np.log10(data[:,21]), #Central pressure
                           #data[:,47], #ocean mass frac
                           ])

    if scale[i] == 'log':
        y_data = np.log10(data[:,predict[i]])
    
    if scale[i] == 'lin':
        y_data = data[:,predict[i]]
    
    print ()
    print ('Fitting', labels[i])  
    #Initiate model instance
    model = regr.multi_linear_model(n_params = len(x_data))
    #Update data for the model fit
    model.set_data(x_data, y_data)
    #Perform multilinear regression
    model.fit(theta0 = [0. for j in range(2**len(x_data))])
    #model.set_params(coefs[i])
    #print ('coefs =', model.coefs)
    #model.evaluate(np.array([np.array([.5, -.4, 0.])]).T)
    #model.plot()
    
    #N = 2**6
    
    model.evaluate(x_data)
    dev = (model.result - model.y_data) / model.y_data
    mean = np.mean(abs(dev))
    q75 = np.quantile(abs(dev), .75)

    print ('q75 =', q75)
    print ('max =', np.max(abs(dev)))
    print ('mean =', mean)
    print ('std =', np.std(abs(dev)))
    #print ('devs =', dev)
    
    models.append(model)
    
    trafo = blended_transform_factory(ax.transAxes, ax.transData)
    
    ax.plot([0, 1], [np.mean(abs(dev))*100, np.mean(abs(dev))*100], 
            transform = trafo, 
            color = cols[i], 
            label = 'mean = '+str(round(100*mean, 2))+'%')
    ax.plot([0, 1], [np.quantile(abs(dev), .75)*100, 
                     np.quantile(abs(dev), .75)*100], 
            linestyle = '--', 
            color = cols[i], 
            transform = trafo,
            label = r'$q75\% = $'+ str(round(q75*100, 2))+'%')
    ax.scatter(data[:,45], np.abs(dev)*100, label = labels[i],
               s = 2, color = cols[i], zorder = 0)

    #plt.scatter(x_data[0], x_data[1], c = np.log10(dev))
    #plt.colorbar()

ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$M/M_\oplus$')
ax.set_ylabel(r'$\rm Error \ of \ prediction \ [\%]$')

fig.savefig(path+'accuracy_prediction.png', 
            format = 'png', 
            bbox_inches = 'tight',
            dpi = 240)
plt.close(fig)

#Save models to file for later use
filename = path + 'models_calibration_mars.pkl'
pickle.dump(models, open(filename, 'wb'))

test_x = np.array([[0], [.52], [.99], [.1], [.4], [np.log10(450e9)], [40.]])
models[2].evaluate(test_x)
print(models[2].result[0])
'''
'''
x1 = np.linspace(-1, 1.)
x2 = np.linspace(.25, .75, 3)
fig, ax = plt.subplots(1, 3, figsize = (12, 3))
plt.subplots_adjust(wspace = .5)

x3 = [[1., .8], [1., .8], [0., .5]]
linestyles = ['-', '--', ':']
colors = ['r', 'g', 'b']
y_labels = [r'$\rm Central \ temperature \ [K]$',
            r'$\rm Central \ pressure \ [GPa]$',
            r'$\rm Core \ mass \ fraction$']

for axx in range(len(ax)):
    ax[axx].set_xlabel(r'$M/M_\oplus$')
    ax[axx].set_ylabel(y_labels[axx])
    
    for k in range(len(x3[axx])):
        x = np.empty([len(x1), len(x2)])
        
        for i in range(len(x1)):
            for j in range(len(x2)):
                dat = np.array([[x1[i]], [x2[j]], [x3[axx][k]]])
    
                models[axx].evaluate(dat)
                x[i][j] = models[axx].result[0]
                
        for j in range(len(x2)):
            if scale[axx] == 'log':
                xx = 10**x.T[j]
            else:
                xx = x.T[j]
                
            ax[axx].semilogx(10**x1, xx, color = colors[j], linestyle = linestyles[k])
            
            if axx == 1:
                ax[axx].set_yscale('log')

fig.savefig('fig.png', format = 'png', bbox_inches = 'tight')
plt.close(fig)
