#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:03:10 2019

@author: oshah
"""
from matplotlib.colors import ListedColormap, LinearSegmentedColormap 

import Material
import Planet
import PlanetTest
import PlanetFort
import functionTools as ftool
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from sklearn.linear_model import LinearRegression

from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
import os
import shutil
import copy
import sys
import numpy as np
import pandas as pd
import astropy.table
from astropy.io import ascii
import time
import readPREM
import logTrans
import random
import plotTools
from PIMPphysicalparams import m_earth, r_earth, Mg_number_solar, M_solar,\
        R_solar, M_trappist1, R_trappist1, M_trappist1_error, R_trappist1_error,\
        names_trappist1, names_solar, M_others, R_others, names_others,\
        M_others_error, R_others_error, mH, mFe, mMg, mO, mSi, NA, m_jupiter, \
        r_jupiter, mH2O, mOl

import phase_transitions_water_Wagner2002 as waterPhase
from PIMPrunparams import suffix, plot_params
from PIMPrunparams import grid_color, background_color, color_list
import fortplanet

image_loc = '/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/'

plt.rcParams["axes.axisbelow"] = False

list_masses = [.25, .5, 1., 2., 3., 4.]

list_Mg_number = [[0.1, 0.5, 0.9],
                  [0.1, 0.5, 0.9],
                  [0.1, 0.5, 0.9]]

list_temps = [400., 300., 600.]

list_P_center = [[[171, 313, 591, 1135, 1683, 2251],
                 [134, 249, 470, 908, 1344, 1784],
                 [66, 132, 249, 475, 701, 931]],
                    [[172, 315, 593, 1148, 1694, 2263],
                     [135, 252, 475, 914, 1351, 1793],
                     [67, 135, 252, 478, 705, 935]],
                    [[170, 294, 586, 1096,1628, 2185],
                     [130, 238, 450, 877, 1305, 1736],
                     [73, 130, 239, 459, 682, 909]]]

list_T_center = [[[777, 978, 1279, 1661, 3537, 4161],
                 [812, 1326, 2190, 2922, 3371, 3705],
                 [702, 1161, 1826, 2317, 2629, 2845]],
                 [[582, 731, 955, 1598, 2424, 3094],
                  [611, 869, 1504, 2229, 2596, 2863],
                  [526, 766, 1342, 1773, 2032, 2202]],
                 [[3788, 4275, 1932, 7458, 8779, 9952],
                  [3390, 4211, 5317, 6746, 7801, 8694],
                  [2908, 597, 4281, 5280, 6028, 6663]]]


def fit_data():              
    data = [[], [], [], [], []]
    for t in range(len(list_temps)):
        for mg in range(len(list_Mg_number)):
            for m in range(len(list_masses)):
                data[0].append(list_temps[t])
                data[1].append(list_Mg_number[t][mg])
                data[2].append(list_P_center[t][mg][m])
                data[3].append(list_T_center[t][mg][m])
                data[4].append(list_masses[m])
                     
    return data


dat = np.array(fit_data())

mod_temp = ftool.fit_data(dat[4], dat[1], dat[0], dat[3])                 
mod_pres = ftool.fit_data(dat[4], dat[1], dat[0], dat[2])

class Workbench():
    def __init__(self, **kwargs):
        self.contents=contents
        self.fractions=fractions
        self.stuff=[]
        
    def create_planet(self, methode='classic', contents=[], fractions=[]):
        if methode == 'classic':
            planet=Material.Planet(contents=contents, fractions=fractions,
                                   tempType = 'adiabatic', P_center=1.0e12)
            
        if methode == 'bisection':
            pass
            

class Toolkit():
    def __init__(self, a=0, b=0):
        self.blabla = 'blabla'
        self.a = a
        self.b = b
        self.bisection = False
        self.delta=0.
        self.oldvals = []
        self.iteration = False
        self.number = 0


    def all_params(self, planet):
        try:
            all_what = {'Mg_number_is': planet.Mg_number_is,
                    'T_surface_is': planet.T_surface_is,
                    'P_surface_is': planet.P_surface_is,
                    'M_surface_is': planet.M_surface_is,
                    'M_ocean_is': planet.layers[4].indigenous_mass/m_earth,
                    'ocean_frac_is':planet.ocean_frac_is,
                    'xi_H_core':planet.xi_H_core,
                    'xi_FeO_mantle':planet.xi_FeS_core
                    }
        
        #If layer 4 does not exist, no ocean is there
        except IndexError:
            all_what = {'Mg_number_is': planet.Mg_number_is,
                    'T_surface_is': planet.T_surface_is,
                    'P_surface_is': planet.P_surface_is,
                    'M_surface_is': planet.M_surface_is,
                    'M_ocean_is': 0.,
                    'ocean_frac_is':-10,
                    'xi_H_core':planet.xi_H_core,
                    'xi_FeO_mantle':planet.xi_S_core
                    }
            
        except ValueError:
            all_what = {'Mg_number_is': planet.Mg_number_is,
                    'T_surface_is': planet.T_surface_is,
                    'P_surface_is': planet.P_surface_is,
                    'M_surface_is': planet.M_surface_is,
                    'M_ocean_is': planet.layers[4].indigenous_mass/m_earth,
                    'ocean_frac_is':-10,
                    'xi_H_core':planet.xi_H_core,
                    'xi_FeO_mantle':planet.xi_S_core
                    }            
        try:
            all_how = {'M_core': planet.layermasses[1]+planet.layermasses[0],
               'T_center': planet.T_center,
               'P_center': planet.P_center,
               'M_outer_mantle': planet.layers[3].indigenous_mass/m_earth
               }
        
        except IndexError:
            try:
                all_how = {'M_core': planet.layermasses[1]+planet.layermasses[0],
                   'T_center': planet.T_center,
                   'P_center': planet.P_center,
                   'M_outer_mantle': 0.
                   }
                
            except IndexError:
                all_how = {'M_core': planet.layermasses[0],
                   'T_center': planet.T_center,
                   'P_center': planet.P_center,
                   'M_outer_mantle': 0.
                   }
                
        return all_what, all_how


    def bisect(self, val_is, val_should, acc, direction, predictor=True, log=False):
        reldev = (val_should - val_is)/val_should
        #if the original criterion is no longer met, that means that the
        #probed parameter has been overshoot. Then parameter delta must be
        #reduced and the previous step undone
 
        if direction[0]*direction[1]*reldev < -acc:
            print ('overshoot')
            self.iteration = True
            if predictor == 'linear':
                x1 = self.oldhowvals[-2]
                x2 = self.oldhowvals[-1]
                
                y1 = self.oldwhatvals[-2]
                y2 = self.oldwhatvals[-1]

                #Extract values for passive parameter
                p1 = [self.oldpassivevals[i][-2] for i in range(len(self.passives))]
                p2 = [self.oldpassivevals[i][-1] for i in range(len(self.passives))]

                #compute predictor slope for passive parameter
                self.passive_slope = [(p2[i]-p1[i])/(x2-x1) for i in range(len(self.passives))]
                

                if log:
                    x1 = np.log10(x1)
                    x2 = np.log10(x2)
                    y1 = np.log10(y1)
                    y2 = np.log10(y2)
                    val_should = np.log10(val_should)

                #compute predictor slope and intercept
                slope = (y2-y1)/(x2-x1)
                self.delta = (val_should - y2)/slope  

            else:
                self.oldhowvals.pop(-1)
                self.oldwhatvals.pop(-1)

                #perform bisection step
                self.delta = self.delta*.5
                self.bisection = True
                print ('starting bisection')
                
            #print ('overshoot detected')
            
            
        elif abs(reldev) <= acc:
            self.iteration = False
            #print ('desired precission reached')
            
        else:
            self.iteration = True
    
    def fancy_iterate(self, planet=None, what='M_surface', how='P_center',
                val_should = 1.0*m_earth, acc = 1.0e-3, predictor='linear',
                iterationLimit=25, update_val_should=False, echo=False,
                updateType=0, unpredictable = True, deltaType=0,
                write_learning_set=False, start_number=0,
                test_data_dir='test_set_00'):
        pass
        
    
    def iterate(self, planet=None, what=['M_surface'], how=['P_center'],
                val_should = [1.0*m_earth], acc = [1.0e-3], predictor=['linear'],
                iterationLimit=25, update_val_should=False, echo=False,
                updateType=0, unpredictable = True, deltaType=0,
                write_learning_set=False, start_number=0, log=False, 
                test_data_dir='test_set_00', passives=['xi_H_core'], test=False,
                passive_predictors = [2], all_val_should_weights = ['lin'],
                all_howval_weights = ['exp'], deltas = None):
        """Iteratively re-constructs a given planetary object changing a
        specified parameter (how) in each iteration to match a given other 
        parameter (what). By default the total mass is probed by iteratively
        adjusting the central pressure of the planet. The desired value
        for the probed parameter is set by the val_should argument and the 
        desired precission at which the val_should and val_is match is set by 
        the argument acc. This method can be iteratively used to fix several
        surface parameters at once. For instance surface temperature and 
        surface pressure can be fixed by first running the iterate methode for
        the temperature and then for the pressure and repeat this until the
        precission for both variables is reached. This does, however, sometimes
        require large numbers of iteration steps and of course not all composition
        specifications an total masses can be expected to give a solution for
        the desired surface pressure and temperature.
        Note, if what = P_surface_is or what = T_surface_is, then the 
        respective parameter should not be used as majorconstraint in the 
        planet.construct() methode to avoid conflicts between the two.
        
        Note: in order to constrain total mass and surface pressure at the same
        time it is recommended to proceed as follows:
            
        Use P_surface_is as majorconstraint for the planet integration and 
        what = M_surface_is, how = P_center. This methode is very stable and
        efficient.
        
        In principle it is possible, to use M_surface_is as majorconstraint and 
        what = P_surface_is, how = p_center, but as the surface pressure is 
        very sensitive to small changes in P_center, the iteration is very 
        inefficient and can only  find approximate solutions (i.e. reldev >~ 
        10%)
        """
        
        self.iteration = True
        self.bisection = False
        self.passives = passives
        n_points = 2**len(how)
        
        for i in range(len(what)):
            w = what[i]
            if w == 'M_surface' or w == 'M_ocean':
                planet.initials[w+'_should'] = val_should[i]/m_earth
            
            else:
                planet.initials[w+'_should'] = val_should[i]
        
        #perform first construction to determine iteration direction for
        #bisection algorithm
        self.oldhowvals = []
        self.oldwhatvals = []
        self.oldshouldvals = []
        self.oldpassivevals = [[] for i in passives]
        self.oldpassivehowvals = [[] for i in passives]
        self.delta = [None for i in what]
        newval = [None for i in what]
        
        #generate dict containing all possible what and how parameters
        all_what, all_how = self.all_params(planet)
                
        sanity_borders = {'M_core':[1.0e-10, 10.], #in earth masses
                          'T_center':[300., 25000.], #in K
                          'P_center':[1.0e7, 5.0e13], #in Pa
                          'M_outer_mantle':[1.0e-6, 10.] #in earth masses
                          }
        
        #these values are fixed based on runtime experience and have not been
        #optimized by any means
        if deltaType == 0:
            initial_deltas = {'M_core': .1,
                          'T_center': .25,
                          'P_center': .25,
                          'M_outer_mantle':0.1}
            
        elif deltaType == 1:
            initial_deltas = {'M_core': .01,
                          'T_center': .1,
                          'P_center': .1,
                          'M_outer_mantle':0.01}

        self.oldhowvals.append([all_how[h] for h in how])
        self.oldwhatvals.append([all_what[w+'_is'] for w in what])
        self.oldshouldvals.append([v for v in val_should])
        
        for i in range(len(passives)):
            if passives[i] == 'xi_H_core' or passives[i] == 'xi_FeO_mantle':
                self.oldpassivevals[i].append(all_what[passives[i]])
                self.oldpassivehowvals[i].append(None)
            
            else:
                self.oldpassivevals[i].append(all_what[passives[i]+'_is'])
                self.oldpassivehowvals[i].append(all_how['P_center'])
                
        self.passive_slope = [0. for i in passives]
        
        val_is = [all_what[w+'_is'] for w in what]
        
        reldev = [(val_should[i] - val_is[i])/val_should[i] for i in range(len(what))]
        direction = [[0,0] for w in what]
        
        self.already_met = [0 for i in range(len(what))]        
        print ('val_should =', val_should)
        print ('initial howval =', [all_how[h] for h in how])
        print ('initial whatval =', [all_what[w+'_is'] for w in what])
        
        print ('initial passiveval is=', self.oldpassivevals)
        
        #Construct initial grid points for prediction
        for p in range(len(what) - 1):
            for i in range(len(what)):
                print (what[i])
                if what[i] == 'Mg_number':
                    direction[i] = [-1, None]       
                    if how[i] == 'M_core':
                        
                        #initial core mass guess too large -> iteration direction down
                        if direction[i][0]*reldev[i] < -acc[i]:
                            direction[i][1] = -1
                            self.delta[i] = -initial_deltas[how[i]]*(planet.layermasses[1]+\
                                                        planet.layermasses[0])
                        
                        #initial core mass guess too small -> iteration direction up
                        elif direction[i][0]*reldev[i] > acc[i]:
                            direction[i][1] = 1
                            self.delta[i] = initial_deltas[how[i]]*(planet.layermasses[1]+\
                                                        planet.layermasses[0])
                        
                        #accuracy already met
                        else:
                            print ('Condition already met')
                            print ('reldev =', reldev[i])
                            self.delta[i] = initial_deltas[how[i]] * all_how[how[i]]/2.
                            self.already_met[i] = 1
                                                            
                    elif how[i] == 'M_outer_mantle':
                        #initial outer mantle mass guess too small -> iteration direction up
                        if direction[i][0]*reldev[i] < -acc[i]:
                            direction[i][1] = -1
                            self.delta[i] = initial_deltas[how[i]]*planet.layermasses[1]
                        
                        #initial outer mantle mass guess too large -> iteration direction down
                        elif direction[i][0]*reldev[i] > acc[i]:
                            direction[i][1] = 1
                            self.delta[i] = -initial_deltas[how[i]]*planet.layermasses[1]
                        
                        #accuracy already met
                        else:
                            #print ('condition already met')
                            self.delta[i] = 0.0
                            self.already_met[i] = 1
                            
                    elif how[i] == 'P_center':
                        #Mg# too small
                        #initial central pressure too low -> iteration direction up
                        if direction[i][0]*reldev[i] < -acc[i]:
                            direction[i][1] = -1
                            self.delta[i] = initial_deltas[how[i]]*planet.P_center
                            
                        elif direction[i][0]*reldev[i] > acc[i]:
                            direction[i][1] = 1
                            self.delta[i] = -initial_deltas[how[i]]*planet.P_center
                        
                        else:
                            self.delta[i] = 0.0
                            self.already_met[i] = 1
                                            
                elif what[i] == 'T_surface' or what[i] == 'P_surface' or \
                what[i] == 'M_surface' or what[i] == 'M_ocean' or what[i] == \
                    'ocean_frac':
                    #State dependency of target on probed parameter
                    #negative means target increases with increasing parameter
                    if what[i] == 'M_ocean' or what[i] == 'ocean_frac':
                        direction[i] = [-1, None]
                        
                    else:
                        direction[i] = [-1, None]
                        
                    #initial guess for central value too low        
                    if direction[i][0]*reldev[i] < -acc[i]:
                        direction[i][1] = -1
                        self.delta[i] = initial_deltas[how[i]]*all_how[how[i]]
                    
                    #initial guess for central value too high
                    elif direction[i][0]*reldev[i] > acc[i]:
                        direction[i][1] = 1
                        self.delta[i] = -initial_deltas[how[i]]*all_how[how[i]]
                    
                    #accuracy already met
                    else:
                        print ('Condition already met')
                        print ('reldev =', reldev[i])
                        self.delta[i] = initial_deltas[how[i]] * all_how[how[i]]/2.
                        self.already_met[i] = 1
                        
                    print ('initial delta =', self.delta[i])
                        
                all_how[how[i]] += self.delta[i]
                
                #force newval to stay within the defined value ranges for the 
                #currently iterated parameter
                newval[i] = min(all_how[how[i]], sanity_borders[how[i]][1])
                newval[i] = max(newval[i], sanity_borders[how[i]][0])
    
                if abs(reldev[i]) <= acc[i]:
                    print ('Desired precission for ', what[i], ' reached.')
                    
                #M_core has to be treated seperately here as it is not an attribute
                #of the Planet class but rather the first entry of layermasses of a
                #Planet.Planet object. Also the planets properties only need to be
                #updated if the condition is NOT already met
                if how[i] == 'M_core':
                    planet.initials['layermasses'][0] = \
                        newval[i]*planet.initials['inner_core_frac']
                    planet.initials['layermasses'][1] = newval[i] - \
                        planet.initials['layermasses'][0]
                    
                elif how[i] == 'M_outer_mantle':
                    planet.initials['layermasses'][3] = newval[i]
                
                else:
                    planet.initials[how[i]] = newval[i]
            
            print ('vals should =', val_should)
            print ('initial howval =', [all_how[h] for h in how])
            print ('initial whatval =', [all_what[w+'_is'] for w in what])
            print ('initial reldev =', reldev)
            print ('initial deltas =', self.delta)
            
            if sum(self.delta) == 0.:
                print ('All parameters already satisfied')
                self.iteration = False
        
            #Gather initial grid points
        
            planet.Reset()
            planet.Update_initials()
            planet.Construct(echo=echo)
            
            all_what, all_how = self.all_params(planet)
            val_is = [all_what[w+'_is'] for w in what]
            reldev = [(val_should[i] - val_is[i])/val_should[i] 
                      for i in range(len(what))]
            
            self.oldhowvals.append([n for n in newval])
            self.oldwhatvals.append([v for v in val_is])
                        
            if how[p+1] == 'M_core':
                planet.initials['layermasses'][0] = \
                    self.oldhowvals[0][1]*planet.initials['inner_core_frac']
                planet.initials['layermasses'][1] = self.oldhowvals[0][1] - \
                    planet.initials['layermasses'][0]
                
            elif how[p+1] == 'M_outer_mantle':
                planet.initials['layermasses'][3] = self.oldhowvals[0][1]
            
            else:
                planet.initials[how[p+1]] = self.oldhowvals[p][p+1]
            
            print ('howvals =', self.oldhowvals)
            print ('whatvals =', self.oldwhatvals)

            newval[p+1] = self.oldhowvals[p][p+1]
            print ('newvals =', newval)
        
        count = 0
        exitcode = 0
        while self.iteration:
            count += 1
            print ('\n iteration', count)
            planet.Reset()
            planet.Update_initials()
            planet.Construct(echo=echo)
            all_what, all_how = self.all_params(planet)
            val_is = [all_what[w+'_is'] for w in what]
            reldev = [(val_should[i] - val_is[i])/val_should[i] for i in 
                      range(len(how))]
            
            if len(self.passives) > 0:
                self.passive_slope = []
                
                for i in range(len(passives)):
                    if passives[i] == 'xi_H_core' or passives[i] == 'xi_FeO_mantle':
                        self.oldpassivevals[i].append(all_what[passives[i]])
                        self.oldpassivehowvals[i].append(None)
                    
                    else:
                        self.oldpassivevals[i].append(all_what[passives[i]+'_is'])
                        self.oldpassivehowvals[i].append(all_how['P_center'])

                    #Extract dependant parameter for passive parameter
                    x1 = self.oldhowvals[-2][passive_predictors[i]]
                    x2 = self.oldhowvals[-1][passive_predictors[i]]
                    
                    #Extract values for passive parameter
                    p1 = self.oldpassivevals[i][-2]
                    p2 = self.oldpassivevals[i][-1]
                    
                    if np.isnan(p2):
                        p2 = p1
                    
                    #compute predictor slope for passive parameter
                    try:
                        self.passive_slope.append((p2-p1)/(x2-x1))
                    
                    except ZeroDivisionError:
                        self.passive_slope.append(0.)

                '''
                #Extract dependant parameter for passive parameter
                x1 = self.oldhowvals[-2][passive_predictors]
                x2 = self.oldhowvals[-1][passive_predictors]
                
                #Extract values for passive parameter
                p1 = [self.oldpassivevals[i][-2] for i in range(len(passives))]
                p2 = [self.oldpassivevals[i][-1] for i in range(len(passives))]
                
                if np.isnan(p2):
                    p2 = p1
                
                #compute predictor slope for passive parameter
                try:
                    self.passive_slope = [(p2[i]-p1[i])/(x2-x1) for i in 
                                      range(len(passives))]
                
                except ZeroDivisionError:
                    self.passive_slope = [0. for i in range(len(self.passives))]
                '''
            #print ('x1, x2 =', x1, x2)
            #print ('p1, p2 =', p1, p2)
            #print ('passiveslope =', self.passive_slope)
            
            self.oldhowvals.append([n for n in newval])
            self.oldwhatvals.append([v for v in val_is])
    
            #print ('howvals =', self.oldhowvals)
            #print ('whatvals =', self.oldwhatvals)
            
            x = np.log10(np.array(self.oldhowvals[-len(how)-1:]))
            y = np.array(self.oldwhatvals[-len(how)-1:])
            x = np.empty([len(how)+1, len(how)])
            
            for k in range(len(how)+1):
                for l in range(len(how)):
                    if all_howval_weights[l] == 'exp':
                        x[k][l] = np.log10(self.oldhowvals[-len(how)+k-1][l])
                    elif all_howval_weights[l] == 'lin':
                        x[k][l] = self.oldhowvals[-len(how)+k-1][l]
                    '''
                    if l < 2:
                        x[k][l] = np.log10(self.oldhowvals[-len(how)+k-1][l])
                    else:
                        x[k][l] = self.oldhowvals[-len(how)+k-1][l]
                    '''
                    
            for k in range(len(how)+1):
                y[k][1] = np.log10(y[k][1])
            
            matrix = np.empty([len(how)+1, len(how)+1])
            
            for k in range(len(how)+1):
                matrix[k][0] = 1.
                for l in range(len(how)):
                    matrix[k][l+1] = y[k][l]

            #print ('x =', x)
            #print ('y =', y)                
            #print ('matrix =', matrix)
            
            try:
                a = np.linalg.solve(matrix, x)
                
                for k in range(len(how)):
                    the_thing = a[0][k]
                    for l in range(len(how)):
                        if all_val_should_weights[l] == 'lin':
                            the_thing += a[l+1][k] * val_should[l]
    
                        elif all_val_should_weights[l] == 'log':
                            the_thing += a[l+1][k] * np.log10(val_should[l])
                    
                    if all_howval_weights[k] == 'lin':
                        newval[k] = the_thing
                        
                    elif all_howval_weights[k] == 'exp':
                        newval[k] = 10**the_thing
                
                '''
                newval[0] = 10**(a[0][0] + a[1][0] * val_should[0] + \
                                 a[2][0] * np.log10(val_should[1]) + \
                                 a[3][0] * val_should[2])
                    
                newval[1] = 10**(a[0][1] + a[1][1] * val_should[0] + \
                                 a[2][1] * np.log10(val_should[1]) + \
                                 a[3][1] * val_should[2])
                    
                newval[2] = a[0][2] + a[1][2] * val_should[0] + \
                                 a[2][2] * np.log10(val_should[1]) + \
                                 a[3][2] * val_should[2]
                '''
                                 
            except np.linalg.LinAlgError as err:
                if 'Singular matrix' in str(err):
                    print ('Singular matrix')
                    '''The reason for a singular matrix can be the fact that
                    one of the parameters is already met in which case the
                    delta is zero and hence y1 = y2. In this case adopt 
                    separate linear extrapolation for the other parameters
                    keeping the ones that are already met constant.
                    '''
                    ind = []
                    
                    #gather indices for parameters that are not already met
                    for i in range(len(self.delta)):
                        if self.delta[i] != 0.:
                            ind.append(i)
                            print ('i =', i, 'not met')
                            #perform sperate linear extrapolation for those parameters
                            #while keeping the others constant
                            x1 = self.oldhowvals[-2][i]
                            x2 = self.oldhowvals[-1][i]                    
        
                            y1 = self.oldwhatvals[-2][i]
                            y2 = self.oldwhatvals[-1][i]
                            
                            try:
                                slope = (y2-y1)/(x2-x1)
                                self.delta[i] = (val_should[i] - y2)/slope
        
                            #In case the iteration has reached a bottleneck where
                            #it does not proceed anymore, the delta will be zero zero and
                            #hence the new slope aswell. In this case just finish the
                            #iteration and try to match the current parameter in the
                            #next sweep
                            except ZeroDivisionError:
                                pass                              

                            newval[i] = max(self.oldhowvals[-1][i] + self.delta[i], 
                                     sanity_borders[how[i]][0])  
                            
                        else:
                            newval[i] = x2
                    
                else:
                    pass

            for i in range(len(how)):
                #newval[i] = min(newval[i], sanity_borders[how[i]][1])
                #newval[i] = max(newval[i], sanity_borders[how[i]][0])
                if newval[i] >= sanity_borders[how[i]][1] or \
                    newval[i] <= sanity_borders[how[i]][0]:
                        exitcode = 3
                        self.iteration = False
                        'WARNING: Ceiling reached in parameter iteration.'
                        break
                
            #print ('a =', a)
            if len(self.passives) > 0:
                self.passive_delta = [max((newval[passive_predictors[i]] - \
                                       self.oldhowvals[-1][passive_predictors[i]])*\
                                      self.passive_slope[i], 1.0e-10) 
                                      for i in range(len(passives))]
                
                for p in range(len(passives)):
                    if not self.passive_slope[p] == 0. and not np.isnan(self.passive_slope[p]):
                        newpassiveval = [self.oldpassivevals[i][-1] + self.passive_delta[i]
                                         for i in range(len(passives))]
                        
                    else:
                        newpassiveval = [0. for i in range(len(passives))]
            
            print ('newval computed =', newval)
            print ('sanity border =', sanity_borders)
            #Sometimes iteration enters a bottle neck. In this case abort!
            for i in range(len(how)):
                if newval[i] == self.oldhowvals[-1][i] and self.oldhowvals[-1][i] == \
                    self.oldhowvals[-2][i]:
                        self.iteration = False
                        print ('WARNING: Bottleneck for '+how[i]+\
                               ' reached. Aborting iteration...')
                        exitcode = 1
                            
                if np.isnan(newval[i]):
                        self.iteration = False
                        print ('WARNING: '+ how[i]+' is NaN. Aborting iteration...')
                        exitcode = 2
                        
            try:
                print ('newpassiveval =', newpassiveval)
                print ('passivedelta =', self.passive_delta)
                
            except UnboundLocalError:
                pass
            
            print ('new val_is =', val_is)
            print ('new reldev =', reldev)
            for i in range(len(how)):
                if how[i] == 'M_core':
                    planet.initials['layermasses'][0] = \
                        newval[i]*planet.initials['inner_core_frac']
                    planet.initials['layermasses'][1] = newval[i] - \
                        planet.initials['layermasses'][0]
                    
                elif how[i] == 'M_outer_mantle':
                    planet.initials['layermasses'][3] = newval[i]
                
                else:
                    planet.initials[how[i]] = newval[i]

            for i in range(len(passives)):
                if passives[i] == 'xi_H_core' or passives[i] == 'xi_FeO_mantle':
                    planet.initials[passives[i]+'_predicted'] = newpassiveval[i]
                    
                else:
                    pass
                    #planet.initials[passives[i]] = newpassiveval[i]
            
            self.already_met = [0 for i in range(len(how))]
            acc_reached_count = 0
            for h in range(len(how)):
                if abs(reldev[h]) < acc[h]:
                    acc_reached_count += 1
                    self.already_met[h] = 1

                else:
                    pass
                
            if acc_reached_count == len(how):
                print ('Desired precission for all parameters reached!')
                print ('reldev =', reldev)
                print ('acc =', acc)
                self.iteration = False
                
            if count >= iterationLimit:
                self.iteration = False
                print ('WARNING: iteration limit reached after', count, 'iterations !')
                
        return exitcode
