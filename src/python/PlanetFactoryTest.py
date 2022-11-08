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
                    'xi_H_core':planet.xi_H_core
                    }
        
        #If layer 4 does not exist, no ocean is there
        except IndexError:
            all_what = {'Mg_number_is': planet.Mg_number_is,
                    'T_surface_is': planet.T_surface_is,
                    'P_surface_is': planet.P_surface_is,
                    'M_surface_is': planet.M_surface_is,
                    'M_ocean_is': 0.,
                    'ocean_frac_is':-10,
                    'xi_H_core':planet.xi_H_core
                    }
            
        except ValueError:
            all_what = {'Mg_number_is': planet.Mg_number_is,
                    'T_surface_is': planet.T_surface_is,
                    'P_surface_is': planet.P_surface_is,
                    'M_surface_is': planet.M_surface_is,
                    'M_ocean_is': planet.layers[4].indigenous_mass/m_earth,
                    'ocean_frac_is':-10,
                    'xi_H_core':planet.xi_H_core
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
                test_data_dir='test_set_00', passives=['xi_H_core'], test=False):
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
                          'T_center':[100., 50000.], #in K
                          'P_center':[1.0e6, 5.0e13], #in Pa
                          'M_outer_mantle':[1.0e-6, 10.] #in earth masses
                          }
        
        #these values are fixed based on runtime experience and have not been
        #optimized by any means
        if deltaType == 0:
            initial_deltas = {'M_core': .1,
                          'T_center': .5,
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
            if passives[i] == 'xi_H_core':
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
        print ('initial val is =', val_is)
        print ('initial passiveval is=', self.oldpassivevals)
        
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
                    
                print ('inital deltas, allhow =', initial_deltas[how[i]], 
                       all_how[how[i]])
            
            all_how[how[i]] += self.delta[i]
            
            #force newval to stay within the defined value ranges for the 
            #currently iterated parameter
            newval[i] = min(all_how[how[i]], sanity_borders[how[i]][1])
            newval[i] = max(newval[i], sanity_borders[how[i]][0])

            if abs(reldev[i]) <= acc[i]:
                print ('Desired precission for ', what[i], ' reached.')
                
            else:
                print ('initial delta =', self.delta[i])
                
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
                    
        print ('initial howval =', [all_how[h] for h in how])
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
        
        for i in range(len(how)):
            newval[i] = min(all_how[how[i]], sanity_borders[how[i]][1])
            newval[i] = max(newval[i], sanity_borders[how[i]][0])
        
        self.oldhowvals.append([n for n in newval])
        self.oldwhatvals.append([v for v in val_is])
        
        if how[1] == 'M_core':
            planet.initials['layermasses'][0] = \
                self.oldhowvals[0][1]*planet.initials['inner_core_frac']
            planet.initials['layermasses'][1] = self.oldhowvals[0][1] - \
                planet.initials['layermasses'][0]
            
        elif how[1] == 'M_outer_mantle':
            planet.initials['layermasses'][3] = self.oldhowvals[0][1]
        
        else:
            planet.initials[how[1]] = self.oldhowvals[0][1]
            
        print ('howvals =', self.oldhowvals)
        print ('whatvals =', self.oldwhatvals)

        newval[1] = self.oldhowvals[0][1]
        
        count = 0
        while self.iteration:
            print ('')
            count += 1
            planet.Reset()
            planet.Update_initials()
            planet.Construct(echo=echo)
            all_what, all_how = self.all_params(planet)
            val_is = [all_what[w+'_is'] for w in what]
            reldev = [(val_should[i] - val_is[i])/val_should[i] for i in 
                      range(len(how))]

            for i in range(len(passives)):
                if passives[i] == 'xi_H_core':
                    self.oldpassivevals[i].append(all_what[passives[i]])
                    self.oldpassivehowvals[i].append(None)
                
                else:
                    self.oldpassivevals[i].append(all_what[passives[i]+'_is'])
                    self.oldpassivehowvals[i].append(all_how['P_center'])
            
            x1 = self.oldhowvals[-2][1]
            x2 = self.oldhowvals[-1][1]
            
            #Extract values for passive parameter
            p1 = [self.oldpassivevals[i][-2] for i in range(len(passives))]
            p2 = [self.oldpassivevals[i][-1] for i in range(len(passives))]

            #compute predictor slope for passive parameter
            self.passive_slope = [(p2[i]-p1[i])/(x2-x1) for i in 
                                  range(len(passives))]

            self.oldhowvals.append([n for n in newval])
            self.oldwhatvals.append([v for v in val_is])
    
            #print ('howvals =', self.oldhowvals)
            #print ('whatvals =', self.oldwhatvals)
            
            x = np.log10(np.array(self.oldhowvals[-3:]))
            y = np.array(self.oldwhatvals[-3:])
            
            for k in range(len(how)+1):
                y[k][1] = np.log10(y[k][1])
            
            #print ('x =', x)
            #print ('y =', y)
            
            matrix = np.array([[1., y[0][0], y[0][1]], 
                               [1., y[1][0], y[1][1]],
                               [1., y[2][0], y[2][1]]])
                
            #print ('matrix =', matrix)
            try:
                a = np.linalg.solve(matrix, x)
                newval[0] = 10**(a[0][0] + a[1][0] * val_should[0] + \
                                 a[2][0] * np.log10(val_should[1]))
                newval[1] = 10**(a[0][1] + a[1][1] * val_should[0] + \
                                 a[2][1] * np.log10(val_should[1]))
            
            except np.linalg.LinAlgError as err:
                if 'Singular matrix' in str(err):
                    
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
                            
                            #perform sperate linear extrapolation for those parameters
                            #while keeping the others constant
                            x1 = self.oldhowvals[-2][i]
                            x2 = self.oldhowvals[-1][i]                    
        
                            y1 = self.oldwhatvals[-2][i]
                            y2 = self.oldwhatvals[-1][i]
                            
                            slope = (y2-y1)/(x2-x1)
                            try:
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
                
            #print ('a =', a)

            self.passive_delta = [(newval[1] - self.oldhowvals[-1][1])*\
                                  self.passive_slope[i] 
                                  for i in range(len(passives))]
                
            '''
            print ('x1, x2 =', x1, x2)
            print ('p1, p2 =', p1, p2)
            print ('passivedelta =', self.passive_delta)
            print ('passiveslope =', self.passive_slope)'''
            for p in range(len(passives)):
                if not self.passive_slope[p] == 0. and not np.isnan(self.passive_slope[p]):
                    newpassiveval = [self.oldpassivevals[i][-1] + self.passive_delta[i]
                                     for i in range(len(passives))]
                    
                else:
                    newpassiveval = [0. for i in range(len(passives))]
            
            print ('newval computed =', newval)
            print ('newpassiveval =', newpassiveval)
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
                if passives[i] == 'xi_H_core':
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
        
        
    def create_planet(self, acc= 1.0e-3, P_surface_should = 1.0e5,
                      T_surface_should=300., M_surface_should = 1.,
                      contents=[], fractions=[], P_center=1.0e12,
                      layermasses=[], tempType=1,
                      T_center = 2000., Mg_number_should=Mg_number_solar,
                      iterationLimit=10, majorConstraints = [],
                      integrateUpTo='P_surface',
                      Fe_number=[], **kwargs):
        """Employ the iterate method to create a planet that matches three
        external properties at once, e.g. total mass, surface pressure and
        surface temperature. The strategy to achieve this is predefined at
        this point that might change in the future to allow the user to
        specify the methode by which the iteration is performed.
        """
        #convert M_surface_should into kg
        M_surface_should *= m_earth
        
        if len(contents) == 0:
            contents = [[9], [1]]
            fractions = [[1], [1]]
        
        t0 = time.time()
        pl = PlanetTest.Planet(contents=contents, fractions=fractions,
                           majorConstraint=integrateUpTo,
                           layermasses = layermasses, P_center=P_center,
                           P_surface_should=P_surface_should,
                           tempType=tempType, T_center=T_center, 
                           vapor_stop=False, Fe_number=Fe_number, 
                           **kwargs)
                
        pl.Construct()
        
        T_surface_dev = abs(T_surface_should - pl.T_surface_is)/T_surface_should
        M_surface_dev = abs(M_surface_should - pl.M_surface_is)/M_surface_should
        P_surface_dev = abs(P_surface_should - pl.P_surface_is)/P_surface_should
        Mg_number_dev = abs(Mg_number_should - pl.Mg_number_is)/Mg_number_should
        
        devs = {'T_surface_is':T_surface_dev,
                'P_surface_is':P_surface_dev,
                'M_surface_is':M_surface_dev,
                'Mg_number_is':Mg_number_dev}
        
        counter = 0
        all_match = False
        while True:
            
            if counter == iterationLimit:
                break            
            
            #self.iterate(planet=pl, what = 'M_surface_is', how = 'P_center',
             #            val_should=M_surface_should, acc = acc)
            
            if 'T_surface_is' in majorConstraints:
                self.iterate(planet=pl, what='T_surface', how='T_center',
                            val_should=T_surface_should, acc = acc)

            if 'Mg_number_is' in majorConstraints:
                self.iterate(planet=pl, what='Mg_number', how='M_core',
                         val_should=Mg_number_should, acc = acc)

            if 'M_surface_is' in majorConstraints:
                self.iterate(planet=pl, what='M_surface', how='P_center',
                         val_should=M_surface_should*m_earth, acc = acc)

            T_surface_dev = abs(T_surface_should - pl.T_surface_is)/\
                                T_surface_should
            M_surface_dev = abs(M_surface_should - pl.M_surface_is)/\
                                (M_surface_should)
            P_surface_dev = abs(P_surface_should - pl.P_surface_is)/\
                                P_surface_should
            Mg_number_dev = abs(Mg_number_should - pl.Mg_number_is)/\
            Mg_number_should

            devs = {'T_surface_is':T_surface_dev,
                    'P_surface_is':P_surface_dev,
                    'M_surface_is':M_surface_dev,
                    'Mg_number_is':Mg_number_dev}

            breakCheck = 0
            for param in majorConstraints:
                if devs[param] < acc:
                    breakCheck += 1
                                        
            if breakCheck == len(majorConstraints):
                all_match = True
                break
            
            counter += 1
            print ('#########################################################')   
                   
        t = time.time()
        print ()
        print ('> Planet construction finished!')
        if all_match:
            print ('> All specified parameters could be matched!')
            
        else:
            print ('> NOT all specified parameters could be matched!')
            
        print ()
        pl.prt()
        print ()
        print ('T_surface_dev:', round(100*T_surface_dev, 3), '%')
        print ('M_surface_dev:', round(100*M_surface_dev, 3), '%')
        print ('P_surface_dev:', round(100*P_surface_dev, 3), '%')
        print ('Mg_number_dev:', round(100*Mg_number_dev, 3), '%')
        print ('total elapsed time:', round(t-t0, 3), 's')
        return pl

    
    def sweep_again(self, pl, 
                    sweeps=1, 
                    acc_Mg_number = 1.0e-2, 
                    acc_T_surface=1.0e-2, 
                    acc_M_surface=1.0e-5, 
                    echo=False):
        
        t0 = time.time()
        Mg_number = pl.Mg_number_should
        T_surface_should = pl.T_surface_should
        
        for i in range(sweeps):
            
            #print ('\nprocessing Mg_number_is...')
            
            if abs(pl.Mg_number_is - pl.Mg_number_should)/pl.Mg_number_should\
                < acc_Mg_number:
                    pass
                    #print ('Mg_number already satisfied')
                    
            else:
                
                self.iterate(planet=pl, what='Mg_number', how='M_core', 
                     val_should = Mg_number, predictor='linear', 
                     acc=acc_Mg_number)
                    
                    #update layermasses in case they don't match the
                    #previous values anymore
                    #for l in range(len(pl.layers) -1):
                     #  pl.layermasses[l] = pl.layers[l].indigenous_mass/m_earth
                
            #print('\nprocessing M_surface_is...')
            
            if abs(pl.M_surface_is - pl.M_surface_should)/pl.M_surface_should\
                < acc_M_surface:
                    pass
                    #print ('M_surface alredy satisfied')
                
            else:
                try:
                    self.iterate(planet=pl, what='M_surface', how='P_center', 
                         val_should = (pl.layermasses[0] + pl.layermasses[1] + \
                         pl.M_H2O_should)*m_earth,
                         predictor='linear', acc=acc_M_surface, update_type=0,
                         update_val_should=True)
                    
                    #update layermasses in case they don't match the
                    #previous values anymore
                    #for l in range(len(pl.layers) -1):
                     #   pl.layermasses[l] = pl.layers[l].indigenous_mass/m_earth
                
                except TypeError:
                    print ('WARNING: Type Error occured in M_surface iteration')
            
            #print ('\nprocessing T_surface_is...')
            if abs(pl.T_surface_is - pl.T_surface_should)/pl.T_surface_should\
                < acc_T_surface:
                    pass
                    #print ('T_surface already satisifed')
            else:
                try:
                    self.iterate(planet=pl, what='T_surface', how='T_center', 
                         val_should = T_surface_should, predictor='linear',
                         acc=acc_T_surface)
                    
                    #update layermasses in case they don't match the
                    #previous values anymore
                    #for l in range(len(pl.layers) -1):
                     #   pl.layermasses[l] = pl.layers[l].indigenous_mass/m_earth
                
                except TypeError:
                    print ('WARNING: Type Error occured in T_surface iteration')
                    
        t=time.time()
        print ('elapsed time in sweep_again:', t-t0)
        
    
    def model(self, M=0., Mg_number_should=Mg_number_solar, sweeps=3,
              T_surface_should=300., P_surface_should=1.0e5,
              acc_Mg_number = 1.0e-2, acc_T_surface=1.0e-2, 
              acc_M_surface=1.0e-5, eps_r=.25, echo=False):

        t0 = time.time()
        
        M_mantle = M
        
        #estimate M_core in earth masses
        M_core = Material.estimate_Mcore(M_mantle=M_mantle, 
                                         Mg_number=Mg_number_should)
        
        #estimate core pressure to ensure that surface pressure is reached
        P_center = 7.5e11 * (1. + M)
        
        #initiate planet that will be constructed and then iterated to match
        #the given constraints
        pl=Planet.Planet(contents=[[9], [11], [0]], fractions=[[1.], [1.], [1.]], 
                         tempType=1, T_center=2000, P_center=P_center, 
                         brucite_dissoc=True, layermasses=[M_core, M, 20.],
                         majorConstraint='P_surface', eps_r=eps_r,
                         P_surface_should=P_surface_should,
                         Mg_number_should=Mg_number_should,
                         vapor_stop=False)
                
        pl.Construct(echo=echo)
        
        for i in range(sweeps):
            
            #print ('\nprocessing Mg_number_is...')
            
            if abs(pl.Mg_number_is - pl.Mg_number_should)/pl.Mg_number_should\
                < acc_Mg_number:
                    pass
                    #print ('Mg_number already satisfied')
                    
            else:
                
                self.iterate(planet=pl, what='Mg_number', how='M_core', 
                     val_should = Mg_number_should, predictor='linear', 
                     acc=acc_Mg_number, echo=echo)
                    
                    #update layermasses in case they don't match the
                    #previous values anymore
                    #for l in range(len(pl.layers) -1):
                     #  pl.layermasses[l] = pl.layers[l].indigenous_mass/m_earth
                
            #print('\nprocessing M_surface_is...')
            
            if abs(pl.M_surface_is - pl.M_surface_should)/pl.M_surface_should\
                < acc_M_surface:
                    pass
                    #print ('M_surface alredy satisfied')
                
            else:
                try:
                    self.iterate(planet=pl, what='M_surface', how='P_center', 
                         val_should = (pl.layermasses[0] + pl.layermasses[1] +\
                         pl.M_H2O_should)*m_earth,
                         predictor='linear', acc=acc_M_surface, update_type=0,
                         update_val_should=True, echo=echo)
                    
                    #update layermasses in case they don't match the
                    #previous values anymore
                    #for l in range(len(pl.layers) -1):
                     #   pl.layermasses[l] = pl.layers[l].indigenous_mass/m_earth
                
                except TypeError:
                    print ('WARNING: Type Error occured in M_surface iteration')
            
            #print ('\nprocessing T_surface_is...')
            if abs(pl.T_surface_is - pl.T_surface_should)/pl.T_surface_should\
                < acc_T_surface:
                    pass
                    #print ('T_surface already satisifed')
            else:
                try:
                    self.iterate(planet=pl, what='T_surface', how='T_center', 
                         val_should = T_surface_should, predictor='linear',
                         acc=acc_T_surface, echo=echo)
                    
                    #update layermasses in case they don't match the
                    #previous values anymore
                    #for l in range(len(pl.layers) -1):
                     #   pl.layermasses[l] = pl.layers[l].indigenous_mass/m_earth
                
                except TypeError:
                    print ('WARNING: Type Error occured in T_surface iteration')
        
        print ('\nTotal elapsed time:', round(time.time()-t0,2), 's')
            
        return pl
    
    
    def model_hydro(self, 
                    M_core=.32, 
                    Mg_number_should=Mg_number_solar,
                    sweeps=50, 
                    iteration_limit=20,
                    T_surface_should = 300., 
                    P_surface_should=1.0e5,
                    M_ocean_should = 0., 
                    M_surface_should=1., 
                    ocean_frac_should=0.,
                    inner_core_frac=.2,
                    Si_number_should=0.25, 
                    acc_T_surface=5.0e-3,
                    acc_M_surface=1.0e-4, 
                    acc_ocean_frac=1.0e-2,
                    acc_Mg_number=1.0e-4, 
                    eps_r=.25, 
                    eps_H2O=0.,
                    eps_Al=0.,
                    echo=False, 
                    ocean=False,
                    vapor_stop=False, 
                    match_mass=False, 
                    temp_jumps = [0., 800., 0, 0., 0.], 
                    q = [.489, .489, 2.5, 2.9, 1.0],
                    P_center=1.0e11, 
                    T_center=3000.,
                    subphase_res=16,
                    Fe_number_mantle=0.,
                    layerConstraint=[0, 0, 3, 1, 0],
                    contents = [[9], [9], [11, 5, 7, 6], [12, 2, 4], [0]],
                    fractions = [[1.], [1.], [1., 1., 1.], [1., 1., 1.], [1.]],
                    layerradii = [.1915, .0, .0, .0, 0.],
                    layermasses=[0, 0, 1., 10., 10.],
                    layerpres=[0, 0, 3.0e10, 0, 0],
                    adiabatType=1, 
                    predictor_T='none',
                    predictor_P='linear',
                    predictor_M_outer_mantle='linear',
                    predictor_M_core='linear',
                    unpredictable=True,
                    T_zero = [100., 100., 50., 50., 50.],
                    deltaType=0,
                    write_learning_set=False,
                    start_number=None,
                    test_data_dir = 'test_set_00',
                    iterationType=0,
                    log=False, 
                    **kwargs):    
        """Creates planet satisfying the specified constraints using the mantle
        hydration model.
        
        The modelType defines the iteration strategy to match all the boundary
        conditions. 1 means that all parameters are probed. 0 means that only
        the total mass and the surface temperature are probed an the other
        parameters are set initially by analytically computing them from the
        composition.
        """
        
        #If test_data_dir already contains data, in order to not overwrite it,
        #count all the data and start counting above the number of already
        #existing files
        if write_learning_set and start_number == None:
            start_number = 0
            for root, dirs, files in os.walk('./'+test_data_dir):
                for file in files:
                    if file.endswith('.pl'):
                        start_number += 1       
        
        #Reset count number for planet test data set output
        self.number = start_number

        t0 = time.time()
        if ocean:
            contents = [[2], [2], [4,5], [6,7], [1]]
            fractions = [[1.], [1.], [1., 1.], [1., 1.], [1.]]
            Si_number_layers = [.0, .0, Si_number_should, Si_number_should, 0.]
            Fe_number_layers = [1., 1., Fe_number_mantle, Fe_number_mantle, 0.]
            layermasses=[0., 0., 1.0e-6, 10., 100.]
            M_ocean_should = M_surface_should * 10**ocean_frac_should
            
        else:
            contents = [[2], [2], [4,5], [6,7]]
            fractions = [[1.], [1.], [1., 1.], [1., 1.]]
            Si_number_layers = [.0, .0, Si_number_should, Si_number_should]
            Fe_number_layers = [1., 1., Fe_number_mantle, Fe_number_mantle]
            layermasses=[0., 0., 1.0e-1, 100.]
            M_ocean_should = 0.
            ocean_frac_should = -10
            
        if iterationType == 0 or iterationType == 2:
            SiMg = Si_number_should/(1.-Si_number_should)
            M_core = PlanetFort.ComputeCoreMass(Mg_number=Mg_number_should,
                                                M_surface=M_surface_should,
                                                Mg_number_mantle=1.-Fe_number_mantle,
                                                M_ocean=M_ocean_should,
                                                SiMg = SiMg,
                                                contents = contents,
                                                xi_H_core=0.)
            
            M_inner_core = M_core*inner_core_frac
            M_outer_core = M_core - M_inner_core

            if ocean:
                M_outer_mantle_enclosed = M_surface_should - M_ocean_should
                
            else:
                M_outer_mantle_enclosed = 100.
            
        elif iterationType == 1:
            M_inner_core = M_core*inner_core_frac
            M_outer_core = M_core - M_inner_core
            M_outer_mantle_enclosed = layermasses[3]

        layermasses[0] = M_inner_core
        layermasses[1] = M_outer_core
        layermasses[3] = M_outer_mantle_enclosed
        
        pc, mc = self.predict(M_surface_should, Mg_number_should)
        print ('M_core =', M_core)
        print ('M_outer_mantle =', M_outer_mantle_enclosed)
        print ('M_ocean =', M_ocean_should)
        #sys.exit()
        pl=PlanetFort.Planet(contents = contents, 
                         fractions = fractions, 
                         layermasses = layermasses, 
                         layerpres = layerpres, 
                         T_center=T_center, 
                         P_center=pc, 
                         P_surface_should=P_surface_should,
                         T_surface_should=T_surface_should,
                         Mg_number_should=Mg_number_should,
                         Si_number_should=Si_number_should,
                         inner_core_frac=inner_core_frac,
                         layerConstraint = layerConstraint, 
                         ocean_frac_should = ocean_frac_should,
                         M_ocean_should=M_ocean_should,
                         M_surface_should=M_surface_should,
                         Si_number_layers=Si_number_layers,
                         Fe_number_layers=Fe_number_layers,
                         eps_r=eps_r,
                         eps_H2O=eps_H2O,
                         eps_Al=eps_Al,
                         temp_jumps=temp_jumps,
                         xi_H_core_predicted = .0, 
                         subphase_res = subphase_res,
                         **kwargs)

        pl.Construct()
        
        iterationCount = 0
        
        for i in range(sweeps):
            iterationCount += 1
            
            if iterationType == 0:
                Mg_number_match = True
                M_ocean_match = True
            
            elif iterationType == 1 or iterationType == 2:
                Mg_number_match = False
                M_ocean_match = False
            
            if not ocean:
                M_ocean_match = True
            
            T_surface_match = False
            M_surface_match = False
            print ('\n sweep', i)
            print (' M_core =', M_core)
            print ('M =', pl.M_surface_should/m_earth)
            print ('ocean frac =', ocean_frac_should)
            print (' T_surf =', T_surface_should)
            print ('Mg# =', Mg_number_should)
            print (' T_center =', pl.T_center)
            print (' P_center =', pl.P_center*1.0e-9)
            
            print ('\nprocessing T_surface_is...')
            if abs(pl.T_surface_is - pl.T_surface_should)/pl.T_surface_should\
                < acc_T_surface:
                    
                    T_surface_match = True
                    print ('T_surface already satisifed')
                    
            if T_surface_match and Mg_number_match and M_ocean_match:
                if match_mass:
                    if M_surface_match:
                        print ('-> match')
                        break
                    
                else:
                    print ('-> match')
                    break
                
            if not T_surface_match:
                self.iterate(planet=pl, what=['Mg_number', 'M_surface'], how=['M_core', 'P_center'], 
                     val_should = [Mg_number_should, M_surface_should*m_earth],
                     predictor=[predictor_M_core, predictor_P],
                     acc=[acc_Mg_number, acc_M_surface], 
                     echo=echo,
                     update_val_should = False,
                     unpredictable=unpredictable,
                     iterationLimit=iteration_limit,
                     deltaType=deltaType,
                     write_learning_set=write_learning_set,
                     start_number=self.number,
                     test_data_dir = test_data_dir)
     
                
            if match_mass:
                print ('\nprocessing M_surface_is...')
                if abs(pl.M_surface_is - pl.M_surface_should)/pl.M_surface_should\
                    < acc_M_surface:
                        
                        M_surface_match = True
                        print ('M_surface already satisifed')
                
                if M_surface_match:
                    if T_surface_match and Mg_number_match and M_ocean_match:
                        print ('-> match')
                        break
                    
                elif not M_surface_match:
                    self.iterate(planet=pl, what=['M_surface'], how=['P_center'], 
                         val_should = [M_surface_should*m_earth], 
                         predictor=[predictor_P],
                         acc=[acc_M_surface], 
                         echo=echo,
                         update_val_should = False,
                         unpredictable=unpredictable,
                         iterationLimit=iteration_limit,
                         deltaType=deltaType,
                         write_learning_set=write_learning_set,
                         start_number=self.number,
                         test_data_dir = test_data_dir,
                         log=log)
                
            if iterationType > 0:             
                print ('\nprocessing Mg# ...')
                if abs(pl.Mg_number_is - pl.Mg_number_should)/pl.Mg_number_should\
                    < acc_Mg_number:
                        
                        Mg_number_match = True
                        print ('Mg_number already satisifed')
                    
                if T_surface_match and Mg_number_match and M_ocean_match:
                    if match_mass:
                        if M_surface_match:
                            print ('-> match')
                            break
                        
                    else:
                        print ('-> match')
                        break                    
                        
                if not Mg_number_match:
                 
                    try:
                        self.iterate(planet=pl, what=['Mg_number'], how=['M_core'], 
                             val_should = [Mg_number_should], predictor=[predictor_M_core],
                             acc=[acc_Mg_number], 
                             echo=echo,
                             update_val_should = False,
                             unpredictable=unpredictable,
                             iterationLimit=iteration_limit,
                             deltaType=deltaType,
                             write_learning_set=write_learning_set,
                             start_number=self.number,
                             test_data_dir = test_data_dir,
                             passives=['xi_H_core', 'M_surface'])
                        
                        print('Mg# =', pl.Mg_number_is)
                        #update layermasses in case they don't match the
                        #previous values anymore
                        #for l in range(len(pl.layers) -1):
                         #   pl.layermasses[l] = pl.layers[l].indigenous_mass/m_earth
                    
                    except TypeError:
                        print ('WARNING: Type Error occured in T_surface iteration')
                
                '''
                if ocean:
                    print ('\nprocessing ocean_frac_is...')
                    reldev = (pl.ocean_frac_is - pl.ocean_frac_should)/\
                    pl.ocean_frac_should
                    if abs(reldev) < acc_ocean_frac:
                            M_ocean_match = True
                            print ('ocean_frac already satisfied')
                            print ('reldev =', reldev)
    
                    if T_surface_match and Mg_number_match and M_ocean_match:
                        if match_mass:
                            if M_surface_match:
                                print ('-> match')
                                break
                            
                        else:
                            print ('-> match')
                            break  
                     
                    if not M_ocean_match:
                        try:
                            pass
                            
                            self.iterate(planet=pl, what='ocean_frac', 
                                         how='M_outer_mantle', 
                                         val_should=ocean_frac_should,
                                         predictor=predictor_M_outer_mantle,
                                         acc=acc_ocean_frac,
                                         echo=echo,
                                         update_val_should=False,
                                         unpredictable=unpredictable,
                                         iterationLimit=iteration_limit,
                                         deltaType=deltaType,
                                         write_learning_set=write_learning_set,
                                         start_number=self.number,
                                         test_data_dir = test_data_dir)
                            
                        except TypeError:
                            print ('WARNING: Type Error occured in M_ocean iteration')
                '''
            if T_surface_match and Mg_number_match and M_ocean_match:
                if match_mass:
                    if M_surface_match:
                        print ('-> match')
                        break
            
        print ('\nNumber of sweeps:', iterationCount)
        print ('Total elapsed time:', round(time.time()-t0,2), 's')
        print ()
        print ('pc, mc predicted =', pc*1.0e-9, 'GPa', mc)
        return pl        

    
    def read_model_hydro_curve(self, loc='./', planet_dir = ['hydro_curve_01'],
                               discard=False):
        """Scanns given location for model hydro curve outputs and gathers 
        it into dataframe
        """
        for pd in planet_dir:
            try:
                os.mkdir('./'+pd+'/test')
                os.rmdir('./'+pd+'/test')          
            
            except FileNotFoundError:
                print ("WARNING: Output directory '"+'./'+pd+"' not found")
        
        planets = []
        temps = []
        Mg_numbers = []
        temp_count = 0
        Mg_count = 0
        
        #Loop over planet dirs
        for i in range(len(planet_dir)):
            planets.append([])
            #Loop over surface temperatures
            for root, dirs, files in os.walk('./'+planet_dir[i]):
                for dir in dirs:
                    if dir.startswith('Tsurf'):
                        planets[i].append([])
                                                
                        string = dir.split('_')
                        
                        temp = float(string[-1])

                        if i == 0:
                            Mg_numbers.append([])
                            Mg_count = 0
                            temp_count += 1
                            temps.append(temp)
                        
                        for root, subdirs, subfiles in os.walk('./'+\
                                                planet_dir[i]+'/'+dir):
                            for subdir in subdirs:
                                if subdir.startswith('Mg'):
                                    planets[i][-1].append([])
                                    
                                    string = subdir.split('_')
                                    Mg = float(string[-1])/100
                                    
                                    if i == 0:
                                        Mg_count+=1                                    
                                        Mg_numbers[-1].append(Mg)
        
                                    planet_loc = './'+planet_dir[i]+\
                                                        '/'+dir+'/'+subdir                                
                                    
                                    #Count planet outputs in subdirectory
                                    planet_count = 0
                                    for root, subsubdirs, subsubfiles \
                                    in os.walk(planet_loc):
                                        for subsubfile in subsubfiles:
                                            if subsubfile.endswith('.pl'):
                                                planet_count += 1
                                    
                                    numerate = np.arange(0, planet_count)
                                    
                                    #Generate planet target list
                                    planet_names = []
                                    for n in numerate:
                                        if n < 10:
                                            num = '00'+str(n)
                                            
                                        elif n < 100 and n >= 10:
                                            num = '0'+str(n)
                                            
                                        else:
                                            num = str(n)
                                            
                                        planet_names.append('planet_'+num+'.pl')
                                                                        
                                    #Collect planet targets and initiate them
                                    #for outputting
                                    for name in planet_names:
                                        planet_name = name.split('.')
                                        planet_name = planet_name[0]
                                        pl = PlanetTest.Planet(silence=True)
                                        pl.load(loc = planet_loc+'/',
                                                file_name = planet_name)
                                        planets[i][-1][-1].append(pl)
                                
        print (temp_count,'temperatures found')
        print (Mg_count,'Mg_numbers found')
        print (len(planet_dir), ' curves found')
        return planets, temps, Mg_numbers


    def print_model_hydro_curve(self, loc ='./', 
                                planet_dir = ['hydro_curve_XH2O_0']):
    
        planets, temps, Mg_numbers = self.read_model_hydro_curve(loc=loc, 
                                            planet_dir = planet_dir)
        
        for i in range(len(planets[0])):
            for j in range(len(planets[0][i])):
                for pl in planets[0][i][j]:
                    M = pl.finals['M_surface_is']
                    T_S = pl.finals['T_surface_is']
                    P_S = pl.finals['P_surface_is']*1.0e-5
                    R = pl.finals['R_surface_is']
                    T_C = pl.initials['T_center']
                    P_C = pl.initials['P_center']*1.0e-9
                    MH2O_frac = pl.finals['M_H2O_is']/M
                    ocean_frac = 10**pl.finals['ocean_frac_is']
                    M_core = (pl.finals['layer_properties'][1]['indigenous_mass']+\
                    pl.finals['layer_properties'][0]['indigenous_mass'])/m_earth
                    Mg_number = pl.finals['Mg_number_is']
                    
                    digits = 5
                    print (round(M, digits), 
                           round(R, digits), 
                           round(T_S, digits), 
                           round(P_S, digits), 
                           round(T_C, digits), 
                           round(P_C, digits), 
                           '', 
                           round(Mg_number, digits), 
                           round(MH2O_frac, digits), 
                           round(ocean_frac, digits), 
                           round (M_core, digits))
        
    def create_data(self, loc='./', 
                    planet_dir = [],
                    MDWR=True):
        planets, temps, Mg_numbers = self.read_model_hydro_curve(loc=loc, 
                                                    planet_dir = planet_dir)
        
        N_temps = len(temps)
        N_Mg = len(Mg_numbers[0])
        N_curves = len(planet_dir)
        N_planets = len(planets[0][0][0])

        #Data contains: mass, radius, water content, delta R dry vs hyd
        data = np.zeros([N_curves, N_temps, N_Mg, 14, N_planets])
        
        profiles = []
        
        for c in range(N_curves):
            profiles.append([])
            for i in range(N_temps):
                profiles[c].append([])
                print ('\nCurve', c)
                for j in range(N_Mg):
                    profiles[c][i].append([])
                    Mg = round(Mg_numbers[i][j], 2)
                    print ('Mg_number:', Mg)
                    
                    for k in range(len(planets[c][i][j])):
                        pl = planets[c][i][j][k] #hyd planet
                        pl0 = planets[0][i][j][k] #dry planet
                        
                        #Check both planets for convergence
                        Mg_reldev = abs((pl.finals['Mg_number_is']/\
                                         Mg_numbers[i][j]-1.))
                        T_reldev = abs((pl.finals['T_surface_is']/\
                                        pl.initials['T_surface_should']-1.))
                        M_reldev = abs(pl.finals['M_surface_is']/\
                                           pl.initials['M_surface_should'] -1.)
                        
                        skip = False
                        
                        #Extract internal profiles of the planet
                        profiles[c][i][j].append(pl.data_table)
                        
                        if Mg_reldev > 1.0e-2 or \
                        T_reldev > 1.0e-2 or \
                        M_reldev > 1.0e-3:
                            print ('Planet does not match')
                            print ('relM =', M_reldev)
                            print ('relMg =', Mg_reldev)
                            print ('relT =', T_reldev)
                            print ('M_surface_should =', 
                                   pl.initials['M_surface_should'])
                            print ('M_surface_is =', pl.finals['M_surface_is'])
                            skip = True
                            
                        if skip:
                            data[c][i][j][0][k] = None
                            data[c][i][j][1][k] = None
                            data[c][i][j][3][k] = None
                        
                        else:
                            data[c][i][j][0][k] = pl.finals['M_surface_is']
                            data[c][i][j][1][k] = pl.finals['R_surface_is']
                            data[c][i][j][4][k] = pl.initials['P_center']
                            data[c][i][j][5][k] = pl.initials['T_center']
                            data[c][i][j][6][k] = \
                            pl.finals['layer_properties'][-2]['P_out']
                            data[c][i][j][7][k] = \
                            pl.finals['layer_properties'][-2]['T_out']
                            
                            data[c][i][j][8][k] = pl.finals['Mg_number_is']
                            data[c][i][j][9][k] = pl.finals['M_DWR_is']
                            data[c][i][j][10][k] = pl.finals['M_ocean_is']
                            data[c][i][j][11][k] = \
                            pl.finals['layer_properties'][1]['indigenous_mass']
                            data[c][i][j][12][k] = pl.finals['P_surface_is']
                            data[c][i][j][13][k] = pl.finals['T_surface_is']
                            
                            if MDWR:
                                #Extract water content in weight fraction
                                data[c][i][j][2][k] = pl.finals['M_H2O_is']/\
                                                    data[c][i][j][0][k]
                                                    
                            else:
                                #Extract ocean thickness in km
                                try:
                                    R_tot = \
                                    pl.finals['layer_properties'][4]['R_out']
                                    
                                except IndexError:
                                    R_tot = \
                                    pl.finals['layer_properties'][3]['R_out']
                                    
                                R_rocky = \
                                pl.finals['layer_properties'][3]['R_out']
                                R_ocean = (R_tot-R_rocky)/1000.
                                data[c][i][j][2][k] = R_ocean
                                #print ('R_ocean =', R_ocean)
                            deltaR = (pl.finals['R_surface_is']-\
                                        pl0.finals['R_surface_is'])/\
                                        pl0.finals['R_surface_is']
                                        
                            data[c][i][j][3][k] = max(deltaR, .0)
                            
        return data, profiles


    def plot_profiles(self, loc = './',
                               planet_dir = ['hydro_curve_XH2O_01'],
                               fnts1=10, fnts2=10, fnts3=12, 
                               filter=False, 
                               indices = [16, 18, 20],
                               write=False,
                               **kwargs):
        
        data, profiles = self.create_data(loc=loc, planet_dir = planet_dir)
        
        profile_length = np.zeros([len(profiles), len(profiles[0]), 
                                   len(profiles[0][0]), len(profiles[0][0][0])])
        
        #Gather number of shells for each planet
        for c in range(len(profiles)):
            for i in range(len(profiles[0])):
                for j in range(len(profiles[0][0])):
                    for k in range(len(profiles[0][0][0])):
                        profile_length[c][i][j][k] = len(profiles[c][i][j][k])
        
        
        for ind in indices:
        
            linestyles = ['-', '--', ':']
            fig, ax = plt.subplots(3, len(profiles[0]), sharey='row', sharex='col',
                                   figsize=(20, 10))
            plt.subplots_adjust(hspace=.1, wspace=.1)
            
            y_labels = [r'$\rm Temperature \ [K]$', 
                        r'$\rm Pressure \ [GPa]$', 
                        r'$\rm Density \ [kg/m^3]$']
            
            titles = [r'$T_S=200 \ \rm K$', 
                      r'$T_S = 300 \ \rm K$', 
                      r'$T_S = 400 \ \rm K$', 
                      r'$T_S = 500 \ \rm K$']
            
            axis_lims = [[0, 1750], [0, 30], [0, 10000]]
            plot_list1 = []
            plot_list2 = []
            Mg_numbers = [.2, .3, .4, .5, .6, .7]
            x_max = .35
            
            for i in range(3):
                for j in range(len(profiles[0])):
    
                    ax[i][j].set_xlim(0, x_max)
                    ax[i][j].set_ylim(axis_lims[i])  
                    ax[i][j].tick_params(which='both', direction='in',
                      right=True, top=True)
                    ax[i][j].tick_params(which='major', length=8)
                    ax[i][j].set_facecolor(plot_params['backgroundcol'])
                    ax[i][j].grid(which='major', color=plot_params['gridcolor'], 
                      alpha=plot_params['gridalpha'], zorder=0, linewidth=2)
                    
                    ax[i][j].grid(which='minor', color=plot_params['gridcolor'], 
                      alpha=plot_params['gridalpha'], zorder=0, linewidth=1)
                    
                    ax[i][j].xaxis.set_minor_locator(AutoMinorLocator())
                    ax[i][j].yaxis.set_minor_locator(AutoMinorLocator())
    
                    if j == 0:
                        ax[i][j].set_ylabel(y_labels[i], fontsize=16)
                        
                    if i == 2:
                        ax[i][j].set_xlabel(r'$R/R_{\oplus}$', fontsize=16)
                    
                    if i == 0:
                        ax[i][j].set_title(titles[j], fontsize=18, pad=30)
                        ax[i][j].tick_params(labeltop=True)
                        
                    if j == 3:
                        ax[i][j].tick_params(labelright=True)
            
            fig.align_ylabels(ax[:, :])       
            
            #Gather number of shells for each planet
            for c in range(len(profiles)):
                for i in range(len(profiles[0])):
                    for j in range(len(profiles[0][0])):
                        for k in range(len(profiles[0][0][0])):
                            if k == ind:
                                mass=data[c][i][j][0][k]
                                length = len(profiles[c][i][j][k])
                                x = np.zeros([length])
                                y1 = np.zeros([length])
                                y2 = np.zeros([length])
                                y3 = np.zeros([length])
                                for l in range(length):
                                    
                                    x[l] = profiles[c][i][j][k][l][0]
                                    y1[l] = profiles[c][i][j][k][l][2]
                                    y2[l] = profiles[c][i][j][k][l][3]
                                    y3[l] = profiles[c][i][j][k][l][4]
                                    if y3[l] > y3[l-2]:
                                        
                                        y3[l-2] = y3[l]
                                        
                                pl1, = ax[0][i].plot(x, y1, linestyle = linestyles[c],
                                  color=color_list[j])
                                pl2, = ax[1][i].plot(x, y2, linestyle = linestyles[c],
                                  color=color_list[j])
                                pl3, = ax[2][i].plot(x, y3, linestyle = linestyles[c],
                                  color=color_list[j])
                                
                                if i == 0 and j == 0:
                                    plot_list1.append(pl1)
                                    
                                if c == 0 and i == 0:
                                    plot_list2.append(pl2)
        
            legend1 = ax[0][0].legend(plot_list1, ['Type-1', 'Type-2', 'Type-3'],
                                                frameon=True, 
                                                loc=1,
                                                fontsize=12,
                                                framealpha=.5,
                                                edgecolor='None')
    
            legend2 = ax[0][1].legend(plot_list2, Mg_numbers,
                                                frameon=True, 
                                                loc=1,
                                                fontsize=12,
                                                framealpha=.5,
                                                edgecolor='None',
                                                title = r'$\rm Mg\#_{tot}:$',
                                                title_fontsize=12)
    
            ax[0][0].add_artist(legend1)
            ax[0][1].add_artist(legend2)        
            
            mass_scinot = ftool.scinot(num=mass, digits=4).replace('.', '_')
    
            ax[0][-1].text(x_max*.75, axis_lims[0][1]*.75, 
              r'$M/M_{\oplus} = $'+str(round(ftool.fancyround(mass, digits=4),4)), 
              fontsize = 18, bbox=dict(facecolor = 'w', edgecolor='k'))
            
            if write:
                fig.savefig('/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/profiles_M_'+\
                        mass_scinot+'.pdf', format = 'pdf', bbox_inches='tight')
        
    
    def plot_delta_delta_R(self, data = []):
        N_temps = len(data[0][0])
        N_Mg = len(data[0][0][0])
        
        plots = plotTools.Plot(col=N_temps,
                               row=2,
                               figsize=(15, 7.5),
                               logx = True,
                               sharex = True,
                               sharey = True,
                               hspace=0.1, 
                               wspace=0.1, 
                               title_pad=50.,
                               axis_labels = [[[r'$M_{\rm tot}/M_{\oplus}$', 
                                                r'$\delta R_{\rm Ocean} \ [\%]$'], 
                                                [r'$M_{\rm tot}/M_{\oplus}$', '']],
                                              [[r'$M_{\rm tot}/M_{\oplus}$', 
                                                r'$\delta (\delta R_2) \ [\%]$'], 
                                               [r'$M_{\rm tot}/M_{\oplus}$', '']]],
                                plot_titles = [[r'$T_{\rm S} = 200 \ \rm K$', 
                                                r'$T_{\rm S} = 400 \ \rm K$'],
                                                ['', '']],
                               majorlocatory = [[50, 50], [50., 50.]],
                               minorlocatory = [[10, 10], [10, 10]],
                               axis_limits=[[[[-3, np.log10(3.)], [-100., 100.0]], 
                                             [[-3, np.log10(3.)], [-100., 100.0]]],
                                            [[[-3, np.log10(3.)], [-100., 100.0]], 
                                             [[-3, np.log10(3.)], [-100., 100.0]]]]
                               )
        plot_list1 = []
        plot_list2 = []
        linestyles = ['-', '--']
        for c in range(len(data)):
            for i in range(N_temps):
                for j in range(N_Mg):
                    y1 = (data[c][0][i][j][2]-data[c][1][i][j][2])/\
                                data[c][1][i][j][2]
                                
                    y2 = (data[c][0][i][j][3]-data[c][1][i][j][3])/\
                                data[c][1][i][j][3]
                    
                    x = data[c][0][i][j][0]           
                
                    p1, = plots.ax[0][i].plot(np.log10(x), y1*100, marker='s',
                                  color = color_list[j*5],
                                  linestyle = linestyles[c])
                    p2, = plots.ax[1][i].plot(np.log10(x), y2*100, marker='s',
                                  color = color_list[j*5],
                                  linestyle = linestyles[c])
                    
                    if c == 0 and i==0:
                        plot_list2.append(p1)
                    
            plot_list1.append(p1)
        
        
        legend1 = plots.ax[0][0].legend(plot_list1, 
                          [r'$\Delta T_{\rm CMB} = 400 \ \rm K$', 
                           r'$\Delta T_{\rm CMB} = 1200 \ \rm K$'], loc=1)
        legend2 = plots.ax[0][0].legend(plot_list2, 
                          ['Mg# = 0.2', 
                           'Mg# = 0.7'], loc=4)
                
        plots.ax[0][0].add_artist(legend1)
        plots.ax[0][0].add_artist(legend2)
         
        plots.fig.savefig('/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/delta_deltaR2.pdf',
                          bbox_inches='tight', format='pdf')
    
    
    def plot_model_hydro_curve(self, loc = './', 
                               planet_dir = ['hydro_curve_XH2O_01'],
                               fnts1=12, fnts2=14, fnts3=16, 
                               MDWR=True, filter=False, 
                               plot_interior=False, discard=False,
                               planet_label_alpha=.75,
                               write=False,
                               **kwargs):
        
        planets, temps, Mg_numbers = self.read_model_hydro_curve(loc=loc, 
                                                    planet_dir = planet_dir)
        
        N_temps = len(temps)
        N_Mg = len(Mg_numbers[0])
        N_curves = len(planet_dir)
        N_planets = len(planets[0][0][0])
        
        #Data contains: mass, radius, water content, delta R dry vs hyd
        data = np.zeros([N_curves, N_temps, N_Mg, 15, N_planets])
        
        for c in range(N_curves):
            for i in range(N_temps):
                temp = temps[i]
                print ('\nCurve', c)
                for j in range(N_Mg):
                    Mg = round(Mg_numbers[i][j], 2)
                    print ('\nMg_number:', Mg)
                    
                    for k in range(len(planets[c][i][j])):
                        pl = planets[c][i][j][k] #hyd or ocean planet
                        pl0 = planets[0][i][j][k] #dry or hyd planet
                                             
                        #Check both planets for convergence
                        Mg_reldev = abs((pl.finals['Mg_number_is']/Mg_numbers[i][j]-1.))
                        T_reldev = abs((pl.finals['T_surface_is']/pl.initials['T_surface_should']-1.))
                        M_reldev = abs(pl.finals['M_surface_is']/\
                                           pl.initials['M_surface_should'] -1.)
                        
                        skip = False
                        
                        #Check if lower mantle exists
                        if discard:
                            if pl.finals['layer_properties'][1]['P_out']> 25.0e9:
                                skip = True
                        
                        if Mg_reldev > 1.0e-2 or \
                        T_reldev > 1.0e-2 or \
                        M_reldev > 2.0e-3:
                            print ('\n Planet', k, 'does not match')
                            print ('relM =', M_reldev)
                            print ('relMg =', Mg_reldev)
                            print ('relT =', T_reldev)
                            print ('M_surface_should =', pl.initials['M_surface_should'])
                            print ('M_surface_is =', pl.finals['M_surface_is'])
                            print ('M_H2O_is =', pl.finals['M_H2O_is'])
                            print ('M_H2O_should =', pl0.finals['M_DWR_is'])
                            print ('T_surface_should =', pl.initials['T_surface_should'])
                            print ('Mg_number_should =', pl.initials['Mg_number_should'])
                            
                            skip = True
                            
                        if skip:
                            data[c][i][j][0][k] = None
                            data[c][i][j][1][k] = None
                            data[c][i][j][3][k] = None
                            data[c][i][j][8][k] = None
                        
                        else:
                            data[c][i][j][0][k] = pl.finals['M_surface_is']
                            data[c][i][j][1][k] = pl.finals['R_surface_is']
                            data[c][i][j][4][k] = pl.initials['P_center']
                            data[c][i][j][5][k] = pl.initials['T_center']
                            data[c][i][j][9][k] = pl.finals['Mg_number_is']                            
                            data[c][i][j][10][k] = pl.finals['M_DWR_is']
                            data[c][i][j][11][k] = pl.finals['M_ocean_is']
                            data[c][i][j][12][k] = pl.finals['layer_properties'][1]['indigenous_mass']+\
                                                    pl.finals['layer_properties'][0]['indigenous_mass']
                            data[c][i][j][13][k] = pl.finals['P_surface_is']
                            data[c][i][j][14][k] = pl.finals['T_surface_is']                            
                            #If no ocean layer is present the pressure at the
                            #bottom of the ocean is 0 and the temperature difference
                            #between bottom of the ocean and surface is zero
                            if len (pl.finals['layer_properties']) > 4:
                                data[c][i][j][6][k] = \
                                    pl.finals['layer_properties'][3]['P_out']
                                data[c][i][j][7][k] = \
                                    pl.finals['layer_properties'][3]['T_out']
                                data[c][i][j][8][k] = \
                                    pl.finals['layer_properties'][4]\
                                    ['indigenous_mass']/m_earth
                                
                            else:
                                data[c][i][j][6][k] = 0.0
                                data[c][i][j][7][k] = 0.0
                                                            
                            if MDWR:
                                #Extract water content in weight fraction
                                data[c][i][j][2][k] = pl.finals['M_H2O_is']/\
                                                    data[c][i][j][0][k]
                                                    
                            else:
                                #Extract ocean thickness in km
                                try:
                                    R_tot = pl.finals['layer_properties'][4]['R_out']
                                    
                                except IndexError:
                                    R_tot = pl.finals['layer_properties'][3]['R_out']
                                    
                                R_rocky = pl.finals['layer_properties'][3]['R_out']
                                R_ocean = (R_tot-R_rocky)/1000.
                                data[c][i][j][2][k] = R_ocean
                                #print ('R_ocean =', R_ocean)
                                
                            deltaR = (pl.finals['R_surface_is']-\
                                        pl0.finals['R_surface_is'])/\
                                        pl0.finals['R_surface_is']
                                        
                            data[c][i][j][3][k] = max(deltaR, .0)
                            #print ('delta R =', deltaR)
                        
        hspaces = [.1]
        
        plot_params = [[], [], [], [], [], [], [], [], []] 
        plot_params2 = [[], [], [], [], [], [], [], [], []] 
        
        markersize = 5
        planet_label_pos_y = [1.85, 1.85, 1.85, 1.85, 1.85, 10., 1.85, 1.85, 1.85, 10., 
                              1.85, 10., 10., .0]
        
        planet_others_label_pos_y = [2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6,
                                     2.6, 2.6, 2.6, 2.6]
        
        planet_sets_colors = [(.25, .25, .25), (.5, .5, .1), (1., .1, .1)]

        water_col=(.75, .75, .75)
        iron_col = (.1, .1, .1)
        rock_col = (.5, .5, .5)
        
        g = ['Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Europa',
               'Io', 'Ganymede', 'Saturn', 'Titan', 'Uranus', 'Neptune', 'Pluto']

        if MDWR:
            middle_y_lims = [0., 4.5]
            middle_y_majorlocator = 1.
            middle_y_minorlocator = .5
            middle_y_label = r'$M_{\rm H_2 O}/M \ \rm [wt \%]$'
            
            bottom_y_lims = [0., 3.]
            
            bottom_y_label = r'$\delta R_1 \  [\%]$'
            
            line_lims_up = [[2., 6., 3.], [100., 1., 3000., 200.]]
            line_lims_low = [[0., -2., 0.], [0., -1., 1000., 0.]]
            
            planet_trappist1_label_pos_y = [4., 4., 1., 4., 4., 4., 1.]
            
            real_objects_legend_pos = [[-1.8, 3.5], [-1.8, 3.], [-1.8, 2.5]]
            real_objects_legend_box = [(-1.85, 2.4), .74, 1.5]
            
        else:
            middle_y_lims = [0., 250.]
            middle_y_majorlocator = 50.
            middle_y_minorlocator = 25.
            middle_y_label = r'$\rm Ocean \ depth \ [km]$'
            
            bottom_y_lims = [0., 3.]

            bottom_y_label = r'$\delta R_2 \  [\%]$'

            line_lims_up = [[2., 300, 4.], [np.log10(2000), 2., 4000., 200.]]
            line_lims_low = [[0., -50., 0.], [0., -1., 0., 0.]]
            
            planet_trappist1_label_pos_y = [50, 50, 50., 50, 50, 50, 50.]
            
            real_objects_legend_pos = [[-1.8, 175], [-1.8, 150], [-1.8, 125]]
            real_objects_legend_box = [(-1.85, 115), .74, 80.]
    
        real_objects_legend_labels = ['Solar Sys.', 'Trappist-1', 'Kepler']
        x_range = [-2, np.log10(3)]
        major_x = .5
        minor_x = .25
        
        plot_params_vals = [[[x_range, [0., 2.]],
                             [x_range, middle_y_lims],
                             [x_range, bottom_y_lims]], #axis lims
                            [major_x, major_x, major_x], #majorlocator x
                            [.5, middle_y_majorlocator, 1.], #majorlocator y
                            [minor_x, minor_x, minor_x], #minorlocator x
                            [.1, middle_y_minorlocator, .5], #minorlocator y
                            [[r'$M/M_{\oplus}$', 
                              r'$R/R_{\oplus}$'], 
                             [r'$M/M_{\oplus}$',
                              middle_y_label],
                             [r'$M/M_{\oplus}$', 
                              bottom_y_label]], #axis labels
                            ['', '', ''], # Plot titles
                            [True, True, True], #logx
                            [False, False, False], #logy
                            ]
        
        plot_params_vals2 = [[[x_range, [1, np.log10(2000)]],
                             [x_range, [0, 2.0]],
                             [x_range, [1000, 4000]],
                             [x_range, [0, 200]]], #axis lims
                            [.05, .05, .05, .05], #majorlocator x
                            [500, .5, 1000., 50], #majorlocator y
                            [.01, .01, .01, .01], #minorlocator x
                            [100, .1, 500., 10], #minorlocator y
                            [[r'$M/M_{\oplus}$', 
                              r'$P_{Center} \ [\rm GPa]$'], 
                             [r'$M/M_{\oplus}$',
                               r'$P_{Ocean} \ [\rm GPa]$'],
                             [r'$M/M_{\oplus}$', 
                              r'$T_{Center} \ [\rm K]$'],
                              [r'$M/M_{\oplus}$', 
                              r'$\Delta T_{Ocean} \ [\rm K]$']], #axis labels
                            ['', '', '', ''], # Plot titles
                            [True, True, True, True], #logx
                            [True, False, False, False], #logy
                            ]
        
        #Create Plot params
        for i in range(3):
            
            for k in range(len(plot_params)):
                plot_params[k].append([])
            
            for j in range(N_temps):
                temp = temps[j]
                
                #axis limits
                for k in range(len(plot_params)):
                    
                    if k == 6:
                        plot_params[k][i].append(r'$T_{S} \ = \ $'+\
                                   str(int(temp))+' K')
                        
                    else:
                        plot_params[k][i].append(plot_params_vals[k][i])
                        
        #Create Plot params
        for i in range(4):
            
            for k in range(len(plot_params2)):
                plot_params2[k].append([])
            
            for j in range(N_temps):
                temp = temps[j]
                
                #axis limits
                for k in range(len(plot_params2)):
                    
                    if k == 6:
                        plot_params2[k][i].append(r'$T_{S} \ = \ $'+\
                                   str(int(temp))+' K')
                        
                    else:
                        plot_params2[k][i].append(plot_params_vals2[k][i])
                
        #Initiate Plot
        for i in range(1):
            plots = plotTools.Plot(col=N_temps, row=3, 
                        axis_limits=plot_params[0], 
                        majorlocatorx=plot_params[1],
                        majorlocatory=plot_params[2], 
                        minorlocatorx=plot_params[3], 
                        minorlocatory=plot_params[4], 
                        axis_labels = plot_params[5],
                        plot_titles=plot_params[6],
                        logx=plot_params[7],
                        logy=plot_params[8],
                        sharex=True, 
                        sharey=True,
                        hspace=hspaces[i], 
                        wspace=.05,
                        majorlabelsize=fnts2, 
                        axislabelsize=fnts3,
                        titlefontsize=fnts3,
                        title_pad=50,
                        **kwargs)
            
            plots2 = plotTools.Plot(col=N_temps, row=4, 
                        axis_limits=plot_params2[0], 
                        majorlocatorx=plot_params2[1],
                        majorlocatory=plot_params2[2], 
                        minorlocatorx=plot_params2[3], 
                        minorlocatory=plot_params2[4], 
                        axis_labels = plot_params2[5],
                        plot_titles=plot_params2[6],
                        logx=plot_params2[7],
                        logy=plot_params2[8],
                        sharex=True, 
                        sharey=True,
                        hspace=.15, 
                        wspace=.05,
                        majorlabelsize=fnts2, 
                        axislabelsize=fnts3,
                        titlefontsize=fnts3,
                        title_pad=50,
                        **kwargs)
            
            plots3 = plotTools.Plot(col=N_temps, row=2,
                                   axis_limits=[[[[0., .02], [0, 2.0]],
                                                 [[0., .02], [0, 2.0]],
                                                 [[0., .02], [0, 2.0]],
                                                 [[0., .02], [0, 2.0]]],
          
                                  [[[-2, np.log10(3.)], [0, .03]],
                                                 [[-2, np.log10(3.)], [0, .03]],
                                                 [[-2, np.log10(3.)], [0, .03]],
                                                 [[-2, np.log10(3.)], [0, .03]]]],
                                    majorlocatorx = [[.005, .005, .005, .005],
                                                     [1., 1., 1., 1.]],
                                    majorlocatory = [[.5, .5, .5, .5],
                                                     [.01, .01, .01, .01]],
                                    minorlocatorx = [[1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3],
                                                     [.5, .5, .5, .5]],
                                    minorlocatory = [[.1, .1, .1, .1],
                                                     [.005, .005,.005,.005]],
                                    logx = [[False, False, False, False],
                                            [True, True, True, True]],
                                    logy = [[False, False, False, False],
                                            [False, False, False, False]],
                                    axis_labels=[[[ r'$M_{\rm Ocean}/M_{\oplus}$',
                                                   r'$P_{\rm Ocean} \ [\rm GPa]$'],
                                                 [ r'$M_{\rm Ocean}/M_{\oplus}$',
                                                   r'$P_{\rm Ocean} \ [\rm GPa]$'],
                                                  [ r'$M_{\rm Ocean}/M_{\oplus}$',
                                                   r'$P_{\rm Ocean} \ [\rm GPa]$'],
                                                   [ r'$M_{\rm Ocean}/M_{\oplus}$',
                                                   r'$P_{\rm Ocean} \ [\rm GPa]$']],
                                                   [[r'$M/M_{\oplus}$', r'$M_{\rm Ocean}/M_{\oplus}$'], [r'$M/M_{\oplus}$', r'$M_{\rm Ocean}/M_{\oplus}$'], [r'$M/M_{\oplus}$', r'$M_{\rm Ocean}/M_{\oplus}$'], [r'$M/M_{\oplus}$', r'$M_{\rm Ocean}/M_{\oplus}$']]],
                                    plot_titles=[[r'$T_{\rm S} \ = \ 200 \ \rm K$',
                                                  r'$T_{\rm S} \ = \ 300 \ \rm K$',
                                                  r'$T_{\rm S} \ = \ 400 \ \rm K$',
                                                  r'$T_{\rm S} \ = \ 500 \ \rm K$'], ['', '', '', '']],
                                    sharey = True,
                                    sharex = False,
                                    hspace=0.5,
                                    title_pad=50)
                                 
            plots4 = plotTools.Plot(col=N_temps, row=1,
                               axis_limits=[[[[0., .1], [0, .1]],
                                             [[0., .1], [0, .1]],
                                             [[0., .1], [0, .1]],
                                             [[0., .1], [0, .1]]]],
                                majorlocatorx = [[.01, .01, .01, .01]],
                                majorlocatory = [[50, 50, 50, 50]],
                                minorlocatorx = [[5.0e-3, 5.0e-3, 5.0e-3, 5.0e-3]],
                                minorlocatory = [[25, 25, 25, 25]],
                                logx = [[False, False, False, False]],
                                logy = [[False, False, False, False]],
                                axis_labels=[[[ r'$M_{\rm Ocean}/M_{\oplus}$',
                                               r'$R_{\rm Ocean} \ [km]$'],
                                             [ r'$M_{\rm Ocean}/M_{\oplus}$',
                                               r'$R_{\rm Ocean} \ [km]$'],
                                              [ r'$M_{\rm Ocean}/M_{\oplus}$',
                                               r'$R_{\rm Ocean} \ [km]$'],
                                               [ r'$M_{\rm Ocean}/M_{\oplus}$',
                                               r'$R_{\rm Ocean} \ [km]$']]],
                                sharey = True,
                                sharex = False)                                
            
        #Plot masses of real objects            
        for p in range(3):
            for i in range(N_temps):
                
                #Plot lines for solar system planets
                for m_sol in M_solar:
                    if m_sol > 0.01 and m_sol <= 3.:
                        plots.ax[p][i].plot([np.log10(m_sol), np.log10(m_sol)], 
                                [line_lims_low[0][p], 
                                 line_lims_up[0][p]],
                                color=planet_sets_colors[0], clip_on=False,
                                zorder=-100, alpha=.75)
                    
                #Plot lines for trappist1 planets
                for m in range(len(M_trappist1)):
                    m_trappist1 = M_trappist1[m]
                    r_trappist1 = R_trappist1[m]
                    if m_trappist1 > 0.01 and m_trappist1 <= 3. and r_trappist1 < 2.:
                        plots.ax[p][i].plot([np.log10(m_trappist1), 
                                np.log10(m_trappist1)], 
                                [line_lims_low[0][p], 
                                 line_lims_up[0][p]],
                                color=planet_sets_colors[1], clip_on=False,
                                linestyle = '-',
                                zorder=-100, alpha=.75)
                
                #Plot lines for other planets
                for m in range(len(M_others)):
                    m_others = M_others[m]
                    r_others = R_others[m]
                    m_others = m_others*m_jupiter/m_earth
                    r_others = r_others*r_jupiter/r_earth
                    if m_others > 0.01 and m_others <= 3. and r_others < 2.:
                        plots.ax[p][i].plot([np.log10(m_others), 
                                np.log10(m_others)], 
                                [line_lims_low[0][p], 
                                 line_lims_up[0][p]],
                                color=planet_sets_colors[2], clip_on=False,
                                linestyle = '-',
                                zorder=-100, alpha=.75)
        
        for p in range(4):
            for i in range(N_temps):
                for m_sol in M_solar:
                    if m_sol > .01 and m_sol <= 3.:
                        plots2.ax[p][i].plot([np.log10(m_sol), np.log10(m_sol)], 
                                 [line_lims_low[1][p], 
                                 line_lims_up[1][p]],
                                color=planet_sets_colors[0], clip_on=False,
                                zorder=-100)

                #Plot lines for trappist1 planets
                for m in range(len(M_trappist1)):
                    m_trappist1 = M_trappist1[m]
                    r_trappist1 = R_trappist1[m]
                    if m_trappist1 > 0.01 and m_trappist1 <= 3. and r_trappist1 < 2.:
                        plots2.ax[p][i].plot([np.log10(m_trappist1), 
                                np.log10(m_trappist1)], 
                                [line_lims_low[1][p], 
                                 line_lims_up[1][p]],
                                color=planet_sets_colors[1], clip_on=False,
                                linestyle = '-',
                                zorder=-100, alpha=.75)


                #Plot lines for other planets
                for m in range(len(M_others)):
                    m_others = M_others[m]
                    r_others = R_others[m]
                    m_others = m_others*m_jupiter/m_earth
                    r_others = r_others*r_jupiter/r_earth
                    if m_others > 0.01 and m_others <= 3. and r_others < 2.:
                        plots2.ax[p][i].plot([np.log10(m_others), 
                                np.log10(m_others)], 
                                [line_lims_low[1][p], 
                                 line_lims_up[1][p]],
                                color=planet_sets_colors[2], clip_on=False,
                                linestyle = '-',
                                zorder=-100, alpha=.75)
                    
                    
        for i in range(N_temps):
            #Plot solar system planets in upper row
            for s in range(len(R_solar)):
                r_sol = R_solar[s]
                m_sol = M_solar[s]
                plots.ax[0][i].scatter([np.log10(m_sol)], [r_sol],
                            color=planet_sets_colors[0], zorder=400, facecolor='none' )
                
            #Plot Trappist planets in upper row
            for t in range(len(R_trappist1)):
                r_trappist1 = R_trappist1[t]
                m_trappist1 = M_trappist1[t]

                plots.ax[0][i].scatter([np.log10(m_trappist1)], [r_trappist1],
                            color=planet_sets_colors[1], zorder=400, facecolor='none')
                
            #Plot other planets in upper row
            for t in range(len(R_others)):
                r_others = R_others[t]*r_jupiter/r_earth
                m_others = M_others[t]*m_jupiter/m_earth

                plots.ax[0][i].scatter([np.log10(m_others)], [r_others],
                            color=planet_sets_colors[2], zorder=400, facecolor='none')
        
        for i in range(N_temps):
            water_masses, water_radii = self.ReadPureCurve(ll=0, T=int(temps[i]))
            iron_masses, iron_radii = self.ReadPureCurve(ll=9, T=int(temps[i]))
            rock_masses, rock_radii = self.ReadPureCurve(ll=11, T=int(temps[i]))
            
            #Plot Mass-Radius for pure iron and water
            plots.ax[0][i].plot(np.log10(iron_masses), iron_radii, color=iron_col, 
                    linewidth=3)
            
            plots.ax[0][i].plot(np.log10(water_masses), water_radii, color=water_col, 
                    linewidth=3)
            
            plots.ax[0][i].plot(np.log10(rock_masses), rock_radii, color=rock_col, 
                    linewidth=3)
            
            #Add label for pure water and iron curve
            plots.ax[0][i].text(-.55, 0.55, r'$\rmFe$', rotation=25,
                    bbox=dict(facecolor=background_color, zorder=-200,
                            edgecolor='None', boxstyle='square,pad=0.1', alpha=.75,
                            ), fontsize = fnts2, color = iron_col)
            
            plots.ax[0][i].text(-.85, .95, r'$\rmH_2 O$', rotation=30,
                    bbox=dict(facecolor=background_color, zorder=-200,
                            edgecolor='None', boxstyle='square,pad=0.1', alpha=.75,
                            ), fontsize = fnts2, color=water_col)
            
            plots.ax[0][i].text(-.85, .9, r'$\rm Mg(OH)_2$', rotation=25,
                    bbox=dict(facecolor=background_color, zorder=-200,
                            edgecolor='None', boxstyle='square,pad=0.1', alpha=.75,
                            ), fontsize = fnts2, color=rock_col)
                    
        #Plot MR-curve and water content
        for c in range(N_curves):
            for p in range(2):
                for i in range(N_temps):
                   for j in range(N_Mg):
                       if c > 0:
                           
                           if p == 0:
                               marker = 's'
                           else:
                               marker = 's'
                
                           if p == 1:
                               if MDWR:
                                   scale = 100
                                 
                               else:
                                   scale = 1.0

                           else:
                               scale = 1.0e0
                        
                           plots.ax[p][i].plot(np.log10(data[c][i][j][0]), 
                               data[c][i][j][1+p]*scale,
                               marker=marker, 
                               color=color_list[j], 
                               linewidth=2,
                               markersize=markersize, 
                               zorder=200)
        
        #Plot internal pressure and temperature
        for i in range(N_temps):
            for j in range(N_Mg):
                Mg = round(Mg_numbers[i][j], 2)
                plots2.ax[0][i].plot(np.log10(data[1][i][j][0]),
                         np.log10(data[c][i][j][4]*1.0e-9)         ,                      
                         marker=marker, 
                           color=color_list[j], 
                           linewidth=2,
                           markersize=markersize, 
                           zorder=100)
                
                plots2.ax[1][i].plot(np.log10(data[1][i][j][0]),
                         data[c][i][j][6]*1.0e-9,
                         marker=marker, 
                           color=color_list[j], 
                           linewidth=2,
                           markersize=markersize, 
                           zorder=100)
                
                plots2.ax[2][i].plot(np.log10(data[1][i][j][0]),
                         data[c][i][j][5],
                         marker=marker, 
                           color=color_list[j], 
                           linewidth=2,
                           markersize=markersize, 
                           zorder=100)
                
                plots2.ax[3][i].plot(np.log10(data[1][i][j][0]),
                         data[c][i][j][7]-temps[i],
                         marker=marker, 
                           color=color_list[j], 
                           linewidth=2,
                           markersize=markersize, 
                           zorder=100,
                           label = r'$\rm Mg \#_{\rm tot} \ =$' + str(Mg))

                if i == 0:
                    plots2.ax[3][i].legend(fontsize=fnts2-2, loc=2, labelspacing=0)
        
        #Plot relative deviation in radius between hyd vs dry
        for i in range(N_temps):
            for j in range(N_Mg):
                Mg = round(Mg_numbers[i][j], 2)
                print ('Mg# label =', Mg)
                x = data[1][i][j][0]
                y = data[1][i][j][3]
                
                #Compute in %
                y *= 100
                pl, = plots.ax[2][i].plot(np.log10(x), y, color=color_list[j], 
                              label = r'$\rm Mg \#_{\rm tot} \ =$' + str(Mg),
                              marker=marker, markersize=markersize)
                    
                plots3.ax[0][i].plot(data[1][i][j][8], data[1][i][j][6]*1.0e-9,
                         color = color_list[j],
                         marker=marker, markersize=markersize,
                         label = r'$\rm Mg \#_{\rm tot} \ =$' + str(Mg)
                         )
  
                plots3.ax[1][i].plot(np.log10(data[1][i][j][0]), data[1][i][j][8],
                         color = color_list[j],
                         marker=marker, markersize=markersize
                         )

                plots4.ax[0][i].plot(data[1][i][j][8]/data[1][i][j][0], 
                         data[1][i][j][2]/data[1][i][j][1]*1000/r_earth,
                         color = color_list[j]
                         )

                if i == 0:
                    plots.ax[2][i].legend(fontsize=fnts2, loc=2, labelspacing=0)
                    plots3.ax[0][i].legend(fontsize=fnts2, loc='lower right', 
                             labelspacing=0)
            
            #Plot labels for solar system planets
            for j in range(len(M_solar)):
                m=M_solar[j]
                name = names_solar[j]
                if m >= .01 and m <= 10.:
                    scaling = 1.
                    
                    if name == 'Pluto':
                        name = ''
                        
                    if name == 'Titan':
                        scaling = 1.02
                        
                    if name == 'Ganymede':
                        scaling = .97
                    
                    if name == 'Europa':
                        scaling = .1
                        
                    if name == 'Io':
                        scaling = .99
     
                    if name == 'Mars':
                        scaling = 1.
                        
                    if name == 'Moon':
                        scaling = 1.
                        
                    if name == 'Mercury':
                        scaling = 1.
                        
                    if name == 'Venus':
                        scaling = 1.2
                   
                    name = ' '+name
  
                    #Add labels for solar system objects
                    try:
                        plots.ax[0][i].text(np.log10(m)*scaling, planet_label_pos_y[j], 
                            name, rotation=90, 
                           fontsize=fnts1, zorder=300, 
                           bbox=dict(facecolor=background_color,
                            edgecolor='None', boxstyle='square,pad=0.05', 
                            alpha=planet_label_alpha,
                            ), ha='center', color=planet_sets_colors[0],  
                            weight='light')

                        plots2.ax[0][i].text(np.log10(m)*scaling, 3., 
                            name, rotation=90, 
                           fontsize=fnts1, zorder=300, 
                           bbox=dict(facecolor=background_color,
                            edgecolor='None', boxstyle='square,pad=0.05', 
                            alpha=planet_label_alpha,
                            ), color=planet_sets_colors[0], ha='center',  
                            weight='light')
                        
                    except IndexError:
                        pass
            
            #Plot labels for Trappist-1 planets
            for j in range(len(M_trappist1)):
                m=M_trappist1[j]
                name = names_trappist1[j]
                if m >= .01 and m <= 10.:
                    scaling = 1.
                      
                    if name == 'b':
                        scaling = -4.
                        
                    elif name == 'c':
                        scaling = 2.8
                    
                    elif name == 'd':
                        scaling = 1.1

                    elif name == 'e':
                        scaling = 2.5
                        
                    elif name == 'f':
                        scaling = 5.25
                        
                    elif name == 'g':
                        scaling = 1.15
                        
                    elif name == 'h':
                        scaling = .95
                        
                    name = ' '+name

                    #Add labels for solar system objects
                    try:
                        plots.ax[1][i].text(np.log10(m)*scaling, planet_trappist1_label_pos_y[j], 
                            name, rotation=0, 
                           fontsize=fnts1, zorder=300, 
                           bbox=dict(facecolor=background_color,
                            edgecolor='None', boxstyle='square,pad=0.01', 
                            alpha=planet_label_alpha), ha='center',
                            color=planet_sets_colors[1],  weight='light')

                        plots2.ax[1][i].text(np.log10(m)*scaling, .5, 
                            name, rotation=0, 
                           fontsize=fnts1, zorder=300, 
                           bbox=dict(facecolor=background_color,
                            edgecolor='None', boxstyle='square,pad=0.05', 
                            alpha=planet_label_alpha,
                            ), ha='center', color=planet_sets_colors[1],  
                            weight='light')
                        
                    except IndexError:
                        pass
                    
            #Plot labels for other planets
            for j in range(len(M_others)):
                m=M_others[j]*m_jupiter/m_earth
                r=R_others[j]*r_jupiter/r_earth
                name = names_others[j]
                if m >= .01 and m < 3. and r < 2.:
                    scaling = 1.
                    
                    if name == 'K-138 d':
                        scaling = 1.1
                    
                    elif name == 'K-128 b':
                        scaling = .8
                        
                    elif name == 'K-011 b':
                        scaling = .8
                        
                    elif name == 'K-138 c':
                        scaling = 1.15
                        
                    elif name == 'K-114 c':
                        scaling = 1.05
                    
                    name = ' '+name
          
                    #Add labels for solar system objects
                    try:
                        plots.ax[2][i].text(np.log10(m)*scaling, 
                                planet_others_label_pos_y[j], 
                                name, rotation=90, 
                                fontsize=fnts1, zorder=300, 
                                bbox=dict(facecolor=background_color,
                                          edgecolor='None',
                                          boxstyle='square,pad=0.05', 
                                          alpha=planet_label_alpha), 
                                color=planet_sets_colors[2], ha='center',
                                weight='light')

                        plots2.ax[3][i].text(np.log10(m)*scaling, 100., 
                            name, rotation=90, 
                           fontsize=fnts1, zorder=300, 
                           bbox=dict(facecolor=background_color,
                            edgecolor='None', boxstyle='square,pad=0.05', 
                            alpha=planet_label_alpha), 
                            color=planet_sets_colors[2], ha='center',
                            weight='light')
                        
                    except IndexError:
                        pass
        
        #Add labels for real objects (solar system, trappist-1 and kepler planets)
        for i in range(len(real_objects_legend_pos)):
            plots.ax[1][0].text(real_objects_legend_pos[i][0], 
                    real_objects_legend_pos[i][1], 
                    real_objects_legend_labels[i], 
                    color = planet_sets_colors[i],
                    fontsize=fnts2, zorder=200)    
   
        rect = patches.Rectangle(real_objects_legend_box[0],
                                 real_objects_legend_box[1],
                                 real_objects_legend_box[2], facecolor = background_color, 
                                 alpha = planet_label_alpha, zorder=100)
        
        plots.ax[1][0].add_patch(rect)
        
        plots.fig.align_ylabels(plots.ax[:, :]) 
        plots2.fig.align_ylabels(plots2.ax[:, :])
        
        for i in range(3):
            for j in range(len(temps)):
                for k, spine in plots.ax[i][j].spines.items():  #ax.spines is a dictionary
                    spine.set_zorder(100)

        for i in range(4):
            for j in range(len(temps)):
                for k, spine in plots2.ax[i][j].spines.items():  #ax.spines is a dictionary
                    spine.set_zorder(100)
        
        if write:
            if MDWR:
                plots.fig.savefig('/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/Type1_vs_Type2.pdf',
                              format='pdf', bbox_inches='tight')
    
            else:
                plots.fig.savefig('/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/Type2_vs_Type3.pdf',
                              format='pdf', bbox_inches='tight')
                plots2.fig.savefig('/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/Type3_internal.pdf',
                              format='pdf', bbox_inches='tight')
                plots3.fig.savefig('/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/P_ocean_vs_M_ocean.pdf',
                              format='pdf', bbox_inches='tight')

        return data
 

    def plot_model_sensitivity(self, data = None):
        """Takes hydro curve data for different deltT_CMB as input and computes
        the relative difference for the main results with respect to the nominal
        case of deltaT_CMB = 800 K
        
        data[0] is taken to be the nominal case
        """
        x = np.array([.005, .01, .05, .1, .15])
        plots= plotTools.Plot(col=2, row=1, 
                              axis_labels=[[[r'$M/M_{\oplus}$', 
                                             r'$\delta (M_{\rm H2O}/M) \ [\%]$'], 
                            [r'$M/M_{\oplus}$', r'$\delta (M_{\rm H2O}/M)$']]], 
                            sharey=True, 
                            majorlocatory=[[5, 5]], 
                            majorlocatorx=[[.05, .05]], 
                            minorlocatory=[[5, 5]], 
                            minorlocatorx=[[.05, .05]], 
                            axis_limits=[[[[0, .15], [-20, 10]], 
                                          [[.0, .15], [-20, 10]]]], 
                            plot_titles = [[r'$T_S = 200 \rm \ K$', 
                                            r'$T_S = 500 \rm \ K$']], 
                            title_pad=30, 
                            wspace=.1, 
                            axislabelsize=16)
        
        colors = [(.7, .2, .6), (.2, .8, .4)]
        linestyles = ['-', '--']
        TCMB= [400, 1200]
        
        Mg = [.2, .7]
        for c in range(2):
            for i in range(2):
                for j in range(2):
                    plots.ax[0][i].plot(x, data[i][c][j], marker='s',
                            color = colors[j],
                            linestyle = linestyles[c], 
                            label = r'$\rm Mg \# =$'+str(Mg[j])+\
                            r'$, \ \Delta T_{\rm CMB} =$'+str(TCMB[c])+' K')
                            
                plots.ax[0][i].legend()
        
        
    def model_hydro_curve(self, 
                          N=5, 
                          T_center=3000,
                          P_center=10.0e9,
                          M_surface_start = .1,
                          M_surface_end = 1.,
                          Mg_number_should=[Mg_number_solar],
                          T_surface_should=[300], 
                          P_surface_should=1.0e5,
                          acc_T_surface=1.0e-2,
                          acc_Mg_number=1.0e-4,
                          acc_M_surface=1.0e-4,
                          acc_ocean_frac=1.0e-2,
                          T_zero = [100., 100., 50., 50., 50.],
                          sweeps=20, 
                          iteration_limit=30,
                          eps_r=.25,
                          M_core_start=5.0e-4,
                          M_core_end=3.,
                          SiMg=.884,
                          log=False,
                          file_name='', 
                          loc='cwd',
                          temp_jumps=[0., 800, 0., 0., 0.],
                          layerConstraint=['mass', 'mass', 'pres', 
                                           'enclosed mass', 'mass'],
                          FeMg_mantle=0., 
                          X_H2O = 1.,
                          inner_core_frac=.2,
                          unpredictable=True,
                          save_planets = False,
                          planet_dir = 'hydro_curve_01',
                          predictor_T='none',
                          predictor_P='linear',
                          predictor_M_core='linear',
                          predictor_M_outer_mantle='linear',
                          match_mass=False,
                          overwrite = True,
                          ocean=False,
                          ocean_frac_should=-10, 
                          read_data=False,
                          data_dir = None,
                          data = None,
                          deltaType=0,
                          iterationType=0,
                          xiBr = 0.,
                          **kwargs):
        """
        """
        
        if read_data:
            initials = np.load('planet_data.npy')
            
            planets, temps, Mg_numbers = \
                        self.read_model_hydro_curve(planet_dir = [data_dir])

            ocean_fracs = np.zeros([len(temps), len(Mg_numbers[0]), len(planets[0][0][0])])
     
            #Extract initial parameters
            initials = np.empty([len(temps), len(Mg_numbers[0]), len(planets[0][0][0]), 4])
            print ('Mg numbers =', Mg_numbers[0])
            print ('temps = ', temps)
            print ('planets =', planets)
            
            #Gather initial input parameters for hydro model
            for i in range(len(temps)):
                for j in range(len(Mg_numbers[0])):
                    for k in range(len(planets[0][i][j])):
                        planet = planets[0][i][j][k]
                        mc = planet.finals['layer_properties'][0]['indigenous_mass']+\
                        planet.finals['layer_properties'][1]['indigenous_mass']
                        mom = planet.finals['layer_properties'][3]['indigenous_mass']
                        
                        initials[i][j][k] = planet.initials['T_center'],\
                                            planet.initials['P_center'],\
                                            mc/m_earth,\
                                            mom/m_earth
                        
                        #Extract water mass in earth masses
                        M_H2O = planet.finals['M_H2O_is']
                        
                        #Compute logarithmic water mass fraction
                        if M_H2O  > 0.:
                            ocean_fracs[i][j][k] = np.log10(M_H2O/\
                                            planet.finals['M_surface_is'])
                            
                        else:
                            ocean_fracs[i][j][k] = -10
                                                                            
            np.save('planet_data', initials)
            print ('ocean fracs =', ocean_fracs)
            
        if ocean:
            q = [.489, .489, 2.5, 2.9, 1.0]
            layerConstraint=['mass', 'mass', 'pres', 'enclosed mass', 'mass']
            contents = [[9], [9], [11, 5, 7, 6], [12, 2, 4], [0]]
            fractions = [[1.], [1.], [1., 1., 1., 1.], [1., 1., 1.], [1.]]
            layerradii = [.1915, .0, .0, .0, 0.]
            layerpres = [0., 0., 3.0e10, 0., 0.]
            X_H2O_list = [[0.], [.0], [X_H2O, X_H2O, X_H2O, X_H2O], 
                          [X_H2O, X_H2O, X_H2O], [0.]]
            layermasses=[0, 0, 1., 10., 10.]
            
        else:
            q = [.489, .489, 2.5, 2.9]
            layerConstraint=['mass', 'mass', 'pres', 'mass']
            contents = [[9], [9], [11, 5, 7, 6], [12, 2, 4]]
            fractions = [[1.], [1.], [1., 1., 1., 1.], [1., 1., 1.]]
            layerradii = [.1915, .0, .0, .0]
            layerpres = [0., 0., 3.0e10, 0., 0.]
            X_H2O_list = [[0.], [.0], [X_H2O, X_H2O, X_H2O, X_H2O], 
                          [X_H2O, X_H2O, X_H2O]]
            layermasses=[0, 0, 1., 10.]
                    
        subsubdir_names = []
        
        if save_planets:
            #Check if planet output directory exists and create it if not        
            try:
                os.mkdir('./'+planet_dir)
            
            except FileExistsError:
                if overwrite:
                    shutil.rmtree('./'+planet_dir)
                    os.mkdir('./'+planet_dir)            
                
                else:
                    print ('WARNING: the output directory already exists!')
                    sys.exit()
                
            #Now prepare subdirectories for each temperature
            for t in range(len(T_surface_should)):
                temp = T_surface_should[t]
                subdir_name = 'Tsurf_'+str(temp)
                
                os.mkdir('./'+planet_dir+'/'+subdir_name)
                
                subsubdir_names.append([])
                
                #Prepare susubdirectories for total Mg#
                for j in range(len(Mg_number_should)):
                    Mg_round = round(Mg_number_should[j], 2)
                    if round(Mg_round, 1) ==  Mg_round:
                        Mg_string = str(Mg_round) + '0'
                        
                    else:
                        Mg_string = str(Mg_round)
                    
                    subsubdir_name = 'Mg_number_'+Mg_string
                    subsubdir_name = subsubdir_name.replace('.', '_')                
                    os.mkdir('./'+planet_dir+'/'+subdir_name+'/'+subsubdir_name)
                    
                    subsubdir_names[t].append(subsubdir_name)
                    
        t0 = time.time()
        if log:
            core_masses_dummy = np.logspace(np.log10(M_core_start), 
                                      np.log10(M_core_end), N)
            total_masses_dummy = np.logspace(np.log10(M_surface_start), 
                                       np.log10(M_surface_end), N)

            core_masses = []
            for i in range(len(core_masses_dummy)):
                cm = core_masses_dummy[i]
                core_masses.append(cm)
                
                #refine core mass grid to better resolve the structure in 
                #the hydration curve
                if cm > 0.0075 and cm < .4 and not i+1 == len(core_masses_dummy):
                    cm1 = .5*(cm+core_masses_dummy[i+1])
                    cm2 = .5*(cm+cm1)
                    cm3 = .5*(cm1+core_masses_dummy[i+1])
                    core_masses.append(cm1)
                    core_masses.append(cm2)
                    core_masses.append(cm3)
                    
            core_masses = np.asarray(core_masses)
            core_masses.sort()
            
        else:
            core_masses = np.linspace(M_core_start, M_core_end, N)
            total_masses_dummy = np.linspace(M_surface_start, M_surface_end, N)
            
        #Estimate max core mass for Fe + Olivine planets for the maximum
        #planet mass to be roughly 1 m_earth
        core_masses = []        
        M_core_max = np.zeros([len(Mg_number_should)])
        for i in range(len(Mg_number_should)):
            FeMg_tot = 1/Mg_number_should[i]-1
            K = 2.*(FeMg_tot*(1-.2)-.2)*mFe/(mOl+2*(.2*mFe-(1-.2)*mMg))
            
            M_tot = .2
            M_core_max[i] = M_tot*K/(1+K)        
            #core_masses.append(np.linspace(M_core_start, M_core_max[i], N))
            core_masses.append(total_masses_dummy*K/(1+K))
            
        N_curves = len(T_surface_should)
        N_Mg = len(Mg_number_should)
        
        planets = []
        total_masses = np.zeros([N_curves, N_Mg, N])
        total_radii = np.zeros([N_curves, N_Mg, N])
        water_masses = np.zeros([N_curves, N_Mg, N])
        
        #Loop over surface temperatures
        for i in range(len(T_surface_should)):
            planets.append([])
            temp = T_surface_should[i]
            subdir_name = 'Tsurf_'+str(T_surface_should[i])

            for j in range(N_Mg):
                if save_planets:
                    subsubdir_name = subsubdir_names[i][j]
                
                print ('i =', i)
                for k in range(N):
                    cm = core_masses[j][k]
                    print ('\nk =', k)
                    print ('T =', temp)
                    print ('Mcore =', cm)
                    print ('M_tot =', total_masses_dummy[k])
                    print ('Mg# =', Mg_number_should[j])
                    print ('SiMg =', SiMg)
                    print ('predictor T=', predictor_T)
                    print ('predictor P =', predictor_P)
                    print ('predictor M core =', predictor_M_core)
                    print ('T jumps =', temp_jumps)
                    print ('FeMg_mantle =', FeMg_mantle)
                    print ('layer constraint =', layerConstraint)
                    print ('match_mass =', match_mass)
                    print ('eps r =', eps_r)
                    print ('XH2O =', X_H2O)
                    
                    if ocean and read_data:
                        #Ocean frac is log(M_ocean/M_surface)
                        ocean_frac_should = ocean_fracs[i][j][k]
                        #T_center, P_center, cm, mom = initials[i][j][k]
                        #layermasses[3] = mom
                        #print ('T_center, P_center, M_core, ocean frac:')
                        #print (T_center, P_center, cm, mom)
                    
                    if ocean:
                        P_evap = waterPhase.evapPres(temp)
                        P_surface_should = max(P_surface_should, P_evap*1.01)
                        print ('P surface (bar) =', P_surface_should *1.0e-5)
                    
                    pl = self.model_hydro(M_core=cm,
                                          P_center=P_center,
                                          T_center=T_center,
                                          M_surface_should=total_masses_dummy[k],
                                          P_surface_should=P_surface_should,
                                          T_surface_should=temp,
                                          Mg_number_should=Mg_number_should[j],
                                          acc_T_surface=acc_T_surface,
                                          acc_Mg_number=acc_Mg_number, 
                                          acc_M_surface=acc_M_surface,
                                          acc_ocean_frac = acc_ocean_frac,
                                          eps_r = eps_r,
                                          iteration_limit=iteration_limit,
                                          SiMg=SiMg,
                                          inner_core_frac=inner_core_frac,
                                          FeMg_mantle = FeMg_mantle,
                                          sweeps=sweeps,
                                          match_mass=match_mass,
                                          X_H2O=X_H2O_list,
                                          layerConstraint=layerConstraint,
                                          temp_jumps = temp_jumps,
                                          T_zero=T_zero,
                                          predictor_T=predictor_T,
                                          predictor_P=predictor_P,
                                          predictor_M_core=predictor_M_core,
                                          predictor_M_outer_mantle=predictor_M_outer_mantle,
                                          ocean_frac_should=ocean_frac_should,
                                          ocean=ocean,
                                          fractions=fractions,
                                          contents=contents,
                                          q=q,
                                          unpredictable=unpredictable,
                                          deltaType=deltaType,
                                          layermasses=layermasses,
                                          layerradii=layerradii,
                                          layerpres=layerpres,
                                          iterationType=iterationType,
                                          xiBr=xiBr,
                                          **kwargs)
                    
                    planets[i].append(pl)
                    
                    if save_planets:
                        if k < 10:
                            num = '00'+str(k)
                        
                        elif k >= 10 and k < 100:
                            num = '0'+str(k)
                            
                        else:
                            num = str(k)
                        
                        write_loc = './'+planet_dir+'/'+subdir_name+'/'+\
                                    subsubdir_name+'/'
                                    
                        print ('writing planet to:', write_loc)
                        
                        try:
                            pl.write(out='planet_'+num, loc=write_loc)
                        
                        except FileNotFoundError:
                            sys.exit('ERROR: No output directory')
                        
                    total_masses[i][j][k] = pl.M_surface_is
                    total_radii[i][j][k] = pl.R_surface_is
                    water_masses[i][j][k] = pl.H2O_count*mH2O
                
        t=time.time()
        dt = t-t0
        minutes = int(dt/60-.5)
        seconds = round(dt - 60*minutes, 1)
        print ('total elapsed time for model_hydro_curve:', 
               minutes, 'min', seconds, 'sec')

        #prepare data as table and save it to file
        if loc=='cwd':
            dir_out=os.getcwd()+'/'
            
        else:
            dir_out=loc+'/'
        
        #create default file name
        if file_name == '':
            file_name = 'hydro_curve.tab'
        '''
        data = []
        for j in range(len(core_masses)):
            data.append([])

            for i in range(len(T_surface_should)):
                st = T_surface_should[i]
                tm = total_masses[i][j]
                tr = total_radii[i][j]
                wm = water_masses[i][j]
                string = str(st)+','+str(tm)+','+str(tr)+','+str(wm)
                data[j].append(string)
        
        #convert data to ascii file
        data_table=astropy.table.Table(data, meta={'meta': None})
        
        ascii.write(data_table, dir_out+file_name, overwrite=True,
                   format='ecsv')
        '''
        if save_planets:
            return None
        
        else:
            return planets
    
    
    def probe_earth(self, N_Mg=3, N_Si=3, N_cm = 3, N_T = 1):
        planet_list = []
        
        temp_jumps = []
        
        for i in range(N_T):
            for j in range(N_T):
                for k in range(N_T):
                    temp_jumps.append([0., 600+i*200, 250+j*50, 1000+k*200])
                    
        print (temp_jumps)
        
        Mg_numbers = np.linspace(0.5, 0.6, N_Mg)
        SiMg_mantle = np.linspace(.8, .9, N_Si)
        core_masses = np.linspace(.32, .34, N_cm)
        print ('parameter space for Mg#_tot:', Mg_numbers)
        print ('parameter space for SiMg_mantle:', SiMg_mantle)
        print ('parameter space for core mass:', core_masses)
        
        for i in range(len(Mg_numbers)):
            planet_list.append([])
            for j in range(len(SiMg_mantle)):
                planet_list[i].append([])
                for k in range(len(core_masses)):                    
                    for l in range(len(temp_jumps)):
                        mg = Mg_numbers[i]
                        sm = SiMg_mantle[j]
                        mc = core_masses[k]
                        tj = temp_jumps[l]
                        pl = self.model_hydro(T_surface_should=300, M_core=mc,
                                              Mg_number_should=mg, SiMg=sm,
                                              temp_jumps = tj)
                        
                        planet_list[i][j].append(pl)
                    
        return planet_list
        
    
    def analyze_earth_models(self, planets=[], 
                             acc_mass=2.0e-2, 
                             acc_radius=2.0e-2, 
                             shade_color=['silver','silver',' silver'],
                             best_fit_color =['c', 'c', 'c'],
                             Nbins = 50, 
                             normal_border=5, 
                             fat_border=10,
                             rect_coords = [[[.875, 3.3], [2, 0], [.875, 1.75]],
                                            [[.875, 3.3], [2, 0], [.875, 1.75]]],
                             rect_extent = [[[.05, 1], [1, 1,], [.05, .5]],
                                             [[.05, 1], [1, 1,], [.05, .5]]]):

        plot1 = plotTools.Plot(axis_limits=[[[[0., 1.02], [0., 14.]]]],
                              majorlocatorx = [[.2]],
                              majorlocatory = [[2]],
                              minorlocatorx = [[.05]],
                              minorlocatory = [[.5]],
                              axis_labels=[[[r'$\rm Radius \  [R_{\oplus}]$', 
                                             r'$\rm Density \ [g \ cm^{-3}]$']]])
                     
        plot2 = plotTools.Plot(axis_limits=[[[[0., 1.02], [0., 400.]]]],
                              majorlocatorx = [[.2]],
                              majorlocatory = [[100]],
                              minorlocatorx = [[.05]],
                              minorlocatory = [[10]],
                              axis_labels=[[[r'$\rm Radius \  [R_{\oplus}]$', 
                                             r'$\rm Pressure \ [GPa]$']]])

        plot3 = plotTools.Plot(axis_limits=[[[[0., 1.02], [0., 6.]]]],
                              majorlocatorx = [[.2]],
                              majorlocatory = [[1]],
                              minorlocatorx = [[.05]],
                              minorlocatory = [[.1]],
                              axis_labels=[[[r'$\rm Radius \  [R_{\oplus}]$', 
                                             r'$\rm Temperature \ [kK]$']]])

        ax1 = plot1.ax[0][0]
        ax2 = plot2.ax[0][0]
        ax3 = plot3.ax[0][0]
        
        ax = [ax1, ax2, ax3]

        radius_scale = 1/r_earth
        val_scale = [1.0e-3, 1.0e-9, 1.0e-3]
        
        #read PREM data
        #Convention is: [radius [km], density [kg m-3]]
        prem = readPREM.Do()
        
        #Plot PREM
        ax1.plot(prem[0]*1.0e3*radius_scale, prem[1]*1.0e-3, 
                color='k', linestyle='--', label='PREM', zorder=9)
        
        survivors = []
        dev_radius_list = []
        dev_mass_list = []
        dev_list = []

        density = []
        pressure = []
        temperature = []
        radius = []
        water = []
        
        bins = []
        total_radii = []
        
        bin_centers = []
        binned_data = []
        binned_min = []
        binned_max = []
        
        Mg_number_list = []
        SiMg_list = []
        
        for acc in range(len(acc_radius)):
            survivors.append([])
            dev_radius_list.append([])
            dev_mass_list.append([])
            dev_list.append([])
            Mg_number_list.append([])
            SiMg_list.append([])
            
            density.append([])
            radius.append([])
            pressure.append([])
            temperature.append([])
            water.append([])
            
            bin_centers.append([])
            bins.append([])
            total_radii.append([])
            binned_data.append([[], [], []])
            binned_min.append([[], [], []])
            binned_max.append([[], [], []])
            
            #Kill planets with too large deviations in mass and/or radius
            for i in range(len(planets)):
                #iterate over Mg#
                for j in range(len(planets[i])):
                    #iterate over Si/Mg
                    for k in range(len(planets[i][j])):
                        pl = planets[i][j][k]
                        
                        Mg_number_is = pl.Mg_number_is
                        SiMg = pl.SiMg
                        radius_is = pl.R_surface_is/r_earth
                        mass_is = pl.M_surface_is/m_earth
                        
                        if k == 0:
                            Mg_number_list[acc].append(Mg_number_is)
                            SiMg_list[acc].append(SiMg)
                        
                        dev_radius = abs(radius_is - 1.)
                        dev_mass = abs(mass_is - 1.)
                        
                        if dev_radius <= acc_radius[acc] and \
                        dev_mass <= acc_mass[acc]:
                                survivors[acc].append(pl)
                                dev_mass_list[acc].append(dev_mass)
                                dev_radius_list[acc].append(dev_radius)
                                dev_list[acc].append(
                                        np.sqrt(dev_mass**2 + \
                                                dev_radius**2))
                                
            print (len(survivors[acc]), 'survivors for', acc_mass[acc])
            Mg_number_list[acc] = np.asarray(Mg_number_list[acc])
            SiMg_list[acc] = np.asarray(SiMg_list[acc])
                
            best_fit_index = dev_list[acc].index(min(dev_list[acc]))
            print ('dev_list=', min(dev_list[acc]))
            
            min_vals = [[], [], []]
            max_vals = [[], [], []]
            
            min_radius = [[], [], []]
            max_radius = [[], [], []]
            
            mean_density_layers = [[], [], [], []]
            mean_temperature_layers = [[], [], [], []]
            mean_pressure_layers = [[], [], [], []]
            
            min_indices_layers = [[], [], []]
            max_indices_layers = [[], [], []]
            
            mtz_rect_coords = [[[], []], [[], []], [[], []]]
                        
            #Collect profiles and bins
            for i in range(len(survivors[acc])):
                pl = survivors[acc][i]
                water[acc].append(pl.M_H2O_is/pl.M_surface_is)
                density[acc].append([])
                radius[acc].append([])
                pressure[acc].append([])
                temperature[acc].append([])
                
                total_radii[acc].append(pl.R_surface_is/r_earth)
                
                for j in range(len(pl.layers)):
                    lay = pl.layers[j]
                    
                    sum_dens = 0
                    sum_pres = 0
                    sum_temp = 0
                    
                    for sh in lay.shells:
                        density[acc][i].append(sh.dens)
                        radius[acc][i].append(sh.radius)
                        temperature[acc][i].append(sh.temp)
                        pressure[acc][i].append(sh.pres)
                        
                        sum_dens += sh.dens
                        sum_pres += sh.pres
                        sum_temp += sh.temp
                        
                    mean_dens = sum_dens/len(lay.shells)
                    mean_pres = sum_pres/len(lay.shells)
                    mean_temp = sum_temp/len(lay.shells)
                    
                    mean_density_layers[j].append(mean_dens)
                    mean_pressure_layers[j].append(mean_pres)
                    mean_temperature_layers[j].append((lay.shells[0].temp+\
                                           lay.shells[-1].temp)/2.)
            
            lowest_radius = 0.0
            highest_radius = max(total_radii[acc])
            
            bins[acc] = np.linspace(lowest_radius, highest_radius, Nbins+1)

            print('bins[acc] =', bins[acc])
            
            #Compute bin centers from bin borders
            bin_centers[acc] = [(bins[acc][i+1] + bins[acc][i])/2. 
                       for i in range(len(bins[acc])-1)]
            
            #Collect data for each bin
            #iterate over parameters
            for i in range(len(binned_data[acc])):
                #iterate over bins
                for j in range(len(bins[acc])-1):
                    binned_data[acc][i].append([]) 
            
            for i in range(len(survivors)):
                pl = survivors[acc][i]
                
                for lay in pl.layers:
                
                    for sh in lay.shells:
                        r = sh.radius/r_earth
                        
                        #Check if shell radius is between bin borders
                        for bin in range(len(bins[acc])-1):
                            if r >= bins[acc][bin] and \
                            r < bins[acc][bin+1]:
                                binned_data[acc][0][bin].append(sh.dens)
                                binned_data[acc][1][bin].append(sh.pres)                
                                binned_data[acc][2][bin].append(sh.temp)
            
            
            print ('binned data[acc] =', binned_data[acc])
            #Compute max and min values for each bin
            #iterate over parameters
            for i in range(len(binned_data[acc])):
                
                binned_max[acc][i] = np.zeros([len(bins[acc])-1])
                binned_min[acc][i] = np.zeros([len(bins[acc])-1])

                #iterate over bins
                for j in range(len(bins[acc])-1):
                    if not len(binned_data[acc][i][j]) == 0:
                        binned_max[acc][i][j] = max(binned_data[acc][i][j])
                        binned_min[acc][i][j] = min(binned_data[acc][i][j])  
                  
                    else:
                        pass
                    
                    
            print ('binned max =', binned_max[acc])
            print ('binned min =', binned_min[acc])
            print ('len binned max[acc][0] =', len(binned_max[acc][0]))
            print ('len binned min[acc][0] =', len(binned_min[acc][0]))
            
            #Extract index for minimum and maximum mean density for each layer
            for lay in range(4):
                #density
                min_index = \
                mean_density_layers[lay].index(min(mean_density_layers[lay]))
                min_indices_layers[0].append(min_index)
    
                max_index = \
                mean_density_layers[lay].index(max(mean_density_layers[lay]))
                max_indices_layers[0].append(max_index)
                
                #pressure
                min_index = \
                mean_pressure_layers[lay].index(min(mean_pressure_layers[lay]))
                min_indices_layers[1].append(min_index)
    
                max_index = \
                mean_pressure_layers[lay].index(max(mean_pressure_layers[lay]))
                max_indices_layers[1].append(max_index)
                
                #temperature
                min_index = \
                mean_temperature_layers[lay].index(min(mean_temperature_layers[lay]))
                min_indices_layers[2].append(min_index)
    
                max_index = \
                mean_temperature_layers[lay].index(max(mean_temperature_layers[lay]))
                max_indices_layers[2].append(max_index)
                        
            #Collect min and max density profiles for each layer individually
            for lay in range(4):
                for i in range(3):
                    pl_min = survivors[acc][min_indices_layers[i][lay]]
                    pl_max = survivors[acc][max_indices_layers[i][lay]]
                    
                    len_min = len(pl_min.layers[lay].shells)
                    len_max = len(pl_max.layers[lay].shells)
                                
                    N_points = 16
                    for j in range(N_points):
                        ind_min = int((len_min-1)/(N_points-1)*j)
                        ind_max = int((len_max-1)/(N_points-1)*j)
                        
                        if i == 0:
                            min_vals[i].append(pl_min.layers[lay].shells[ind_min].dens)
                            max_vals[i].append(pl_max.layers[lay].shells[ind_max].dens)
                        
                        elif i == 1:
                            min_vals[i].append(pl_min.layers[lay].shells[ind_min].pres)
                            max_vals[i].append(pl_max.layers[lay].shells[ind_max].pres)
    
                        elif i == 2:
                            min_vals[i].append(pl_min.layers[lay].shells[ind_min].temp)
                            max_vals[i].append(pl_max.layers[lay].shells[ind_max].temp)
                        
                        min_radius[i].append(pl_min.layers[lay].shells[ind_min].radius)
                        max_radius[i].append(pl_max.layers[lay].shells[ind_max].radius)
                        
                        if lay == 2 and j == N_points-1:
                            mtz_rect_coords[i][1].append(
                                    pl_max.layers[lay].shells[ind_max].radius)
                            
                            if i == 0:
                                mtz_rect_coords[i][1].append(
                                    pl_max.layers[lay].shells[ind_max].dens)
                                
                            elif i == 1:
                                mtz_rect_coords[i][1].append(
                                    pl_max.layers[lay].shells[ind_max].pres)
    
                            elif i == 2:
                                mtz_rect_coords[i][1].append(
                                    pl_max.layers[lay].shells[ind_max].temp)
                                                            
    
                        if lay == 3 and j == 0:
                            mtz_rect_coords[i][0].append(
                                    pl_min.layers[lay].shells[ind_min].radius)
                            
                            if i == 0:
                                mtz_rect_coords[i][0].append(
                                    pl_min.layers[lay].shells[ind_min].dens)
                        
                            elif i == 1:
                                mtz_rect_coords[i][0].append(
                                    pl_min.layers[lay].shells[ind_min].pres)
                                
                            elif i == 2:
                                mtz_rect_coords[i][0].append(
                                    pl_min.layers[lay].shells[ind_min].temp)
            
            data = np.array([density, pressure, temperature])
            
            for i in range(len(ax)):
                axx = ax[i]
                scale = val_scale[i]
        
            
                #Plot best fit model
                if acc == 0:
                    axx.plot(np.asarray(radius[acc][best_fit_index])*radius_scale, 
                            np.asarray(data[i][acc][best_fit_index])*scale,
                            label ='best fit', color=best_fit_color[i],
                            zorder = 10)
                
                '''
                for s in range(len(survivors[acc])):
                    axx.plot(np.asarray(radius[acc][s])*radius_scale, 
                            np.asarray(data[i][acc][s])*scale,
                            color=best_fit_color[i],
                            zorder = 10)
                '''
                '''
                axx.plot(np.asarray(min_radius[i])*radius_scale, 
                np.asarray(min_vals[i])*scale,
                label ='min', color='k',
                linestyle = '--', marker='s',
                zorder = 10)

                axx.plot(np.asarray(max_radius[i])*radius_scale, 
                np.asarray(max_vals[i])*scale,
                label ='max', color='k',
                linestyle = '--', marker='s',
                zorder = 10)
                '''
                
                axx.fill(
                        np.append(min_radius[i], max_radius[i][::-1])*radius_scale,
                        np.append(min_vals[i], max_vals[i][::-1])*scale,
                        shade_color[acc][i]
                        )
                
                #Plot min and max for each layer individually        
                axx.plot(np.asarray(min_radius[i])*radius_scale,
                        np.asarray(min_vals[i])*scale,
                        linewidth=normal_border, color=shade_color[acc][i], zorder=acc)
                
                axx.plot(np.asarray(max_radius[i])*radius_scale, 
                        np.asarray(max_vals[i])*scale,
                        linewidth=normal_border, color=shade_color[acc][i], zorder=acc,
                        label = r'$\rm \delta = $'+\
                        str(acc_mass[acc]))
 
    
                if acc > 0:
                    axx.plot(np.asarray(min_radius[i])*radius_scale,
                            np.asarray(min_vals[i])*scale,
                            linewidth=fat_border, color=shade_color[acc-1][i], zorder=acc-1)
                    
                    axx.plot(np.asarray(max_radius[i])*radius_scale, 
                            np.asarray(max_vals[i])*scale,
                            linewidth=fat_border, color=shade_color[acc-1][i], zorder=acc-1)
    
    
                '''
                axx.plot(bin_centers[acc],
                        np.asarray(binned_min[acc][i])*scale,
                        linewidth=3, color=shade_color[acc][i], zorder=acc)

                axx.plot(bin_centers[acc],
                        np.asarray(binned_max[acc][i])*scale,
                        linewidth=3, color=shade_color[acc][i], zorder=acc)
                '''
                '''
                x = mtz_rect_coords[i][0][0]*radius_scale
                y = mtz_rect_coords[i][0][1]*scale
                
                dx = (mtz_rect_coords[i][1][0]-mtz_rect_coords[i][0][0])*radius_scale
                dy = (mtz_rect_coords[i][1][1]-mtz_rect_coords[i][0][1])*scale
                '''     
                mtz_rect = patches.Rectangle(rect_coords[acc][i], 
                                             rect_extent[acc][i][0],
                                             rect_extent[acc][i][1], 
                                             edgecolor = shade_color[acc][i], 
                                             facecolor = shade_color[acc][i], 
                                             zorder=0, linewidth=3)
            
                axx.add_patch(mtz_rect)        
                
                axx.legend(loc=1)
            
            dwr_best_fit = water[acc][best_fit_index]
            dwr_min = min(water[acc])
            dwr_max = max(water[acc])
            
            dplus = dwr_max - dwr_best_fit
            dminus = dwr_best_fit - dwr_min
            
            a = round(dwr_best_fit*100, 3)
            b = round(dplus*100, 3)
            c = round(dminus*100, 3)
            print (a, b, c)
            #ax.text(0.05, 2.5, r'$\mathbf {MDWR = {0.028}^{+0.076}_{-0.013} \ wt\%}$')
            ax1.text(0.025, 12.1, r'$\rm inner \ core$')
            ax1.text(0.25, 10.1, r'$\rm outer \ core$')
            ax1.text(0.55, 4.1, r'$\rm lower \ mantle$')
            ax1.text(0.8, 2.1, r'$\rm upper \ mantle$')
                    
            print ('Water content of best fit:', water[acc][best_fit_index]\
                   *100,'wt%')
            
            print ('Min water content:', min(water[acc])*100,'wt%')
            print ('Max water content:', max(water[acc])*100,'wt%')
    
    
    def check_model(self, planets, params = {'P_surface', 'T_surface',
                        'Mg_number', 'M_H2O', 'M_surface'}, 
                        acc = {'P_surface':1.0e-2, 
                        'T_surface':1.0e-2, 'Mg_number':1.0e-4, 'M_H2O':5.0e-2,
                        'M_surface':1.0e-4}):
    
        """Takes planet object and checks if the major constraints are satisfied.
        This is important if several constraints are iteratively matched as it is
        possible that for given conditions no solution to the constraints could be
        found. The argument "constr" is a dict containing the name of the planet
        property to be checked and it's corresponding value.
        """
        
        reldevs = {}
        
        for param in params:
            reldevs[param] = []
        
        for planet in planets:
            planet.match = True
            
            for param in params:
                
                try:
                    val_should = planet.initials[param+'_should']
                    
                except KeyError:
                    try:
                        val_should = planet.finals[param+'_should']
                    
                    except KeyError:
                        print ('WARNING: key', param, 'could not be found')
                        
 
                try:
                    val_is = planet.initials[param+'_is']
                    
                except KeyError:
                    try:
                        val_is = planet.finals[param+'_is']
                    
                    except KeyError:
                        print ('WARNING: key', param, 'could not be found')
                                                
                if val_should == 0.:      
                    try:
                        reldev =  abs(val_should-val_is)
                        
                    except TypeError:
                        print ('type error for ', param)
                
                else:
                    try:
                        if param == 'M_H2O':
                            reldev = abs(val_should-val_is)/val_should
#                            planet.finals['M_surface_is']
                        
                        else:
                            reldev =  abs(val_should-val_is)/val_should
                        
                    except TypeError:
                        print ('type error for ', param)                
                
                if reldev < acc[param]:
                    #print (param+' matches')
                    reldevs[param].append(reldev)
                    
                else:
                    print ('\nreldev',str(param),'=', round(reldev*100, 4), '%')
                    print ('Mg#/M_surface =', 
                           planet.initials['Mg_number_should'],
                           '/',
                           planet.finals['M_surface_is'])
                    reldevs[param].append(reldev)
                    planet.match = False
                    
            dens = []
            for lay in planet.layers:
                for sh in lay.shells:
                    dens.append(sh.dens)
                
            planet.density_inversion=False
            for d in range(len(dens)-1):
                if (dens[d] - dens[d+1])/dens[d] < -.01:
                    planet.density_inversion=True
                    break
        
        for param in params:
            try:
                reldevs[param] = max(reldevs[param])
            
            except  ValueError:
                pass
        
        return reldevs
        
    
    def GeneratePureCurve(self, N=8, ll=0, M_start=.01, M_end=5.,
                          P_surface_should=1.0e5, T_surface_should=300,
                          sweeps=10, save_planets=True, overwrite=False,
                          log=True, fix=False, P_start=1.0e9, P_end=1.0e12):
        """Generates M-R curve for pure substance at given surface conditions
        within the given mass range
        """
        if log:
            if fix:
                masses= np.logspace(np.log10(M_start), np.log10(M_end), N)
            else:
                pressures = np.logspace(np.log10(P_start), np.log10(P_end), N)
        else:
            if fix:
                masses = np.linspace(M_start, M_end, N)
            else:
                pressures = np.linspace(P_start, P_end, N)
        
        planets = []

        if ll == 0:
            P_evap = waterPhase.evapPres(T_surface_should)
            P_surface_should = max(P_surface_should, P_evap*1.01)        

        for i in range(N):
            if fix:
                p = 1.0e9
            else:
                p = pressures[i]
                
            pl = PlanetTest.Planet(contents=[[ll]], layermasses=[100.],
                                   T_surface_should=T_surface_should, 
                                   P_surface_should=P_surface_should,
                                   T_center=3000,
                                   P_center=p)
            
            pl.Construct()
            
            for j in range(sweeps):
                self.iterate(planet=pl, what='T_surface', how='T_center',
                             val_should=T_surface_should)
                
                if fix:
                    self.iterate(planet=pl, what='M_surface', how='P_center',
                             val_should = masses[i]*m_earth)
            
            planets.append(pl)

        if save_planets:
            planet_dir = 'pure_curve_'+str(ll)+'_T_'+str(T_surface_should)

            #Check if planet output directory exists and create it if not        
            try:
                os.mkdir('./'+planet_dir)
            
            except FileExistsError:
                if overwrite:
                    shutil.rmtree('./'+planet_dir)
                    os.mkdir('./'+planet_dir)            
                
                else:
                    pass

            for k in range(N):
                pl = planets[k]
                
                if k < 10:
                    num = '00'+str(k)
                
                elif k >= 10 and k < 100:
                    num = '0'+str(k)
                    
                else:
                    num = str(k)
                
                write_loc = './'+planet_dir+'/'
                            
                print ('writing planet to:', write_loc)
                
                try:
                    pl.write(out='planet_'+num, loc=write_loc)
                
                except FileNotFoundError:
                    sys.exit('ERROR: No output directory')
            
        return planets


    def ReadPureCurve(self, T=300, ll=0):
        """Reads in planets from pure curve directory and extract the M-R
        values
        """
        planet_dir = 'pure_curve_'+str(ll)+'_T_'+str(T)
        masses = []
        radii = []

        planet_loc = './'+planet_dir                              
            
        #Count planet outputs in subdirectory
        planet_count = 0
        for root, subsubdirs, subsubfiles \
        in os.walk(planet_loc):
            for subsubfile in subsubfiles:
                if subsubfile.endswith('.pl'):
                    planet_count += 1
        
        numerate = np.arange(0, planet_count)
        
        #Generate planet target list
        planet_names = []
        for n in numerate:
            if n < 10:
                num = '00'+str(n)
                
            elif n < 100 and n >= 10:
                num = '0'+str(n)
                
            else:
                num = str(n)
                
            planet_names.append('planet_'+num+'.pl')

        #Collect planet targets and initiate them
        #for outputting
        for name in planet_names:
            planet_name = name.split('.')
            planet_name = planet_name[0]
            pl = PlanetTest.Planet(silence=True)
            pl.load(loc = planet_loc+'/',
                    file_name = planet_name)
            
            masses.append(pl.finals['M_surface_is'])
            radii.append(pl.finals['R_surface_is'])
            
        return np.array(masses), np.array(radii)


    def PlotPureCurve(self, ll=0, T=300):
        masses, radii = self.ReadPureCurve(T=T, ll=ll)
        
        fig, ax = plt.subplots()
        ax.set_xlabel(r'$M/M_{\oplus}$')
        ax.set_ylabel(r'$R/R_{\oplus}$')
        
        ax.plot(np.asarray(masses), np.asarray(radii), marker='s')

    
    def model_curve(self, N=10, Mg_number_should=Mg_number_solar,
              T_surface_should=300., P_surface_should=1.0e5, sweeps=20,
              acc_T_surface=5.0e-3, 
              acc_M_surface=1.0e-5, eps_r=.25, M_start=.1, 
              M_end = 1., bisection=True, resweeps=0, 
              refine_curve=0, modelType=0, strategy = 0, **kwargs):
        """Generates a MR diagram for the Mg(OH)2 dissociation model with
        specified boundary conditions. If desired, a bisection algorithm to
        iteratively add a new planet in the gap between two planets where the
        dissociation occures can be employed. This reduces the total number
        of planets that must be computed to resolve the dissociation 
        threshold.
        """
        min_mass = 1.0e-4
        
        if modelType == 0:
            function = self.model
            
        elif modelType == 1:
            function = self.model_next
        
        planets = []
        planets_dummy = []
        devs = []
        
        #employ bisection algorithm for more efficient probing of the
        #dissociation threshold
        if bisection:    
            #start with three planets
            
            #MCreate lowest mass planet
            planet1 = function(M=M_start, 
                                  Mg_number_should=Mg_number_should, sweeps=sweeps,
                                  T_surface_should=T_surface_should, 
                                  P_surface_should=P_surface_should,
                                  acc_T_surface=acc_T_surface, 
                                  acc_M_surface=acc_M_surface, 
                                  eps_r=eps_r)
            
            #Create highest mass planet
            planet2 = function(M=M_end, 
                                  Mg_number_should=Mg_number_should, sweeps=sweeps,
                                  T_surface_should=T_surface_should, 
                                  P_surface_should=P_surface_should,
                                  acc_T_surface=acc_T_surface, 
                                  acc_M_surface=acc_M_surface, 
                                  eps_r=eps_r)

            planets.append(planet1)
            planets.append(planet2)
            
            #check if planets have converged properly and resweep if
            #necessary and desired
            if resweeps > 0 and modelType == 0:
                self.check_model(planets=planets)
                for pl in planets:
                    if not pl.match:
                        self.sweep_again(pl, sweeps=resweeps)
            
            a = M_start
            b = M_end
    
            unused = []
    
            #for the remaining N-3 planets use bisection
            #the last planet will go in the interval that is not beeing used
            #for bisection giving (N-3) + 2 + 1 = N planets in total
            for i in range(N-2):
                
                c = (b+a)/2.0
                print ('c1 =', c)
                olda = a
                oldb = b
                
                #M_mantle_end
                planet = function(M=c, 
                                      Mg_number_should=Mg_number_should, sweeps=sweeps,
                                      T_surface_should=T_surface_should, 
                                      P_surface_should=P_surface_should,
                                      acc_T_surface=acc_T_surface, 
                                      acc_M_surface=acc_M_surface, 
                                      eps_r=eps_r)
                
                #check if planets have converged properly and resweep if
                #necessary and desired
                if resweeps > 0 and modelType == 0:
                    self.check_model(planets=[planet])
                    if not planet.match:
                        self.sweep_again(planet, sweeps=resweeps)
                
                planets.append(planet)
                
                #check how accurately the dissociation transition is probed
                dev = 0.0
                try:
                    dev = abs(planets[-1].M_surface_is - \
                          planets[-2].M_surface_is)/planets[-1].M_surface_is
                              
                    if dev < 1.0e-3:
                        print ('convergance reached after', i+1+3,'iterations')
                        break
                
                except ZeroDivisionError:
                    pass
                
                devs.append(dev)
                
                #if the new planet has non-zero water layer, the dissociation 
                #threshold lies at M_mantle < c --> take (c-a)/2 as new value
                relamount = planet.M_H2O_is/planet.M_surface_is*m_earth
                if relamount > min_mass:
                    b = c
                
                #if the new planet has zero water layer, the dissociation
                #treshold lies at M_mantle > c --> take (b-c)/2 as new value
                else:
                    a = c
                    
                #for the first iteration, remember which interval has been
                #used so that in the end one additional planet can be added
                #to the other interval to make the curve a bit smoother
                if i==0 and refine_curve > 0:
                    if b == c:
                        unused_mass = (c + oldb)/2.
                        
                        if refine_curve > 2:
                            dummyc1 = (c+unused_mass)/2.
                            dummyc2 = (oldb+unused_mass)/2.
                            
                            unused.append(dummyc1)
                            unused.append(dummyc2)
                    
                    else:
                        unused_mass = (olda + c)/2.

                        if refine_curve > 2:
                            dummyc1 = (c+unused_mass)/2.
                            dummyc2 = (olda+unused_mass)/2.
                            
                            unused.append(dummyc1)
                            unused.append(dummyc2)
                        
                    unused.append(unused_mass)
                        
                if i==1 and refine_curve > 1:
                    if b == c:
                        unused_mass = (c + oldb)/2.
                    
                    else:
                        unused_mass = (olda + c)/2.
                        
                    unused.append(unused_mass)                    
                    
            for mm in unused:          
                #now add last planet to unused interval
                planet = function(M=mm, 
                                      Mg_number_should=Mg_number_should, sweeps=sweeps,
                                      T_surface_should=T_surface_should, 
                                      P_surface_should=P_surface_should,
                                      acc_T_surface=acc_T_surface, 
                                      acc_M_surface=acc_M_surface,
                                      eps_r=eps_r) 

                #check if planets have converged properly and resweep if
                #necessary and desired
                if resweeps > 0 and modelType == 0:
                    self.check_model(planets=[planet])
                    if not planet.match:
                        self.sweep_again(planet, sweeps=resweeps)
                    
                planets.append(planet)
                
            #after the dissociation threshold is probed, probe now
            #the threshold for full dissociation
            planets = self.sort_planets(planets=planets)
            
            for i in range(len(planets)):
                planet = planets[i]
                relamount = planet.M_H2O_hidden_is/planet.M_surface_is*m_earth
                if relamount < min_mass:
                    try:
                        previousrelamount = \
                        planets[i-1].M_H2O_hidden_is/planet.M_surface_is*m_earth
                        if previousrelamount > min_mass:
                            if modelType == 0:
                                a = planets[i-1].layermasses[1]
                                b = planets[i].layermasses[1]
                                
                            elif modelType == 1:
                                a = planets[i-1].M_surface_is/m_earth
                                b = planets[i].M_surface_is/m_earth
                    
                    except IndexError:
                        pass            
            
            for i in range(N-len(planets)):
                c = (b+a)/2.0
                print ('c2 =', c)
                olda = a
                oldb = b
                
                #M_mantle_end
                planet = function(M=c, 
                                      Mg_number_should=Mg_number_should, sweeps=sweeps,
                                      T_surface_should=T_surface_should, 
                                      P_surface_should=P_surface_should,
                                      acc_T_surface=acc_T_surface, 
                                      acc_M_surface=acc_M_surface,
                                      eps_r=eps_r)
                
                #check if planets have converged properly and resweep if
                #necessary and desired
                if resweeps > 0 and modelType == 0:
                    self.check_model(planets=[planet])
                    if not planet.match:
                        self.sweep_again(planet, sweeps=resweeps)
                
                planets.append(planet)
                
                #check how accurately the full dissociation transition is probed
                dev = 0.0
                try:
                    dev = abs(planets[-1].M_surface_is - \
                              planets[-2].M_surface_is)/\
                              planets[-1].M_surface_is
                    
                    if dev < 1.0e-3:
                        print ('convergance reached after', i+1+3,'iterations')
                        break
                
                except ZeroDivisionError:
                    pass
                
                devs.append(dev)
                
                #if the new planet has non-zero water layer, the dissociation 
                #threshold lies at M_mantle < c --> take (c-a)/2 as new value
                relamount = planet.M_H2O_hidden_is/planet.M_surface_is*m_earth
                if relamount < min_mass:
                    b = c
                
                #if the new planet has zero water layer, the dissociation
                #treshold lies at M_mantle > c --> take (b-c)/2 as new value
                else:
                    a = c
            
            #now add some planets in big gaps to make curve smoother
            #check M-distance for two adjacent planets. If distance exceeds
            #certain value, add a planet in between
            for s in range(2):
                
                planets = self.sort_planets(planets, what = 'M_surface_is')
                planets_dummy = []
                
                #compute all distances between two adjacent planets
                for i in range(len(planets)-1):
                    if abs(planets[i].M_surface_is-planets[i+1].M_surface_is)/\
                    planets[i].M_surface_is > 0.25:
                        
                        if modelType == 0:
                            c = (planets[i].layermasses[1] + \
                             planets[i+1].layermasses[1])/2.
                            
                        elif modelType == 1:
                            c = (planets[i].M_surface_is + \
                                 planets[i+1].M_surface_is)/2.
                            c = c/m_earth
                        
                        print ('c3 =', c)
    
                        planet = function(M=c, 
                                              Mg_number_should=Mg_number_should, sweeps=sweeps,
                                              T_surface_should=T_surface_should, 
                                              P_surface_should=P_surface_should,
                                              acc_T_surface=acc_T_surface, 
                                              acc_M_surface=acc_M_surface, 
                                              eps_r=eps_r)
                        
                        #check if planets have converged properly and resweep if
                        #necessary and desired
                        if resweeps > 0 and modelType == 0:
                            self.check_model(planets=[planet])
                            if not planet.match:
                                self.sweep_again(planet, sweeps=resweeps)
                        
                        planets_dummy.append(planet)
                        
                for i in range(len(planets_dummy)):
                    planets.append(planets_dummy[i])
                       
        #use equal spacing for mantle mass between given start and end point
        else:
            dM = (M_end - M_start)/(N-1)
            
            for i in range(N):
                planet = function(M=M_start+i*dM, 
                                  Mg_number_should=Mg_number_should, sweeps=sweeps,
                                  T_surface_should=T_surface_should, 
                                  P_surface_should=P_surface_should,
                                  acc_T_surface=acc_T_surface, 
                                  acc_M_surface=acc_M_surface,
                                  eps_r=eps_r)
                
                planets.append(planet)        
        
        return planets, devs
   
    
    def extract_dissoc_positions(self, pops=[]):
        dissoc_positions = []
        fulldissoc_positions = []
        dissoc_positions_error = []
        fulldissoc_positions_error = []
        
        min_mass = 1.0e-4
        
        #iterate over magnesium numbers
        for i in range(len(pops)):
            pop = self.sort_planets(pops[i], what='M_surface_is')
            
            found_dissoc = False
            found_fulldissoc = False
            print ('processing Mg# =', pop[0].initials['Mg_number_should'])
            #iterate over single planets
            for j in range(len(pop)):
                planet = pop[j]
                
                #extract dissociation transition
                relamount = planet.finals['M_H2O_is']/\
                            planet.finals['M_surface_is']
                
                if pop[0].initials['Mg_number_should'] == 0.3:
                    print (planet.finals['M_H2O_is']/\
                           planet.finals['M_surface_is'])
                
                if relamount < min_mass:
                    try:
                        nextrelamount = pop[j+1].finals['M_H2O_is']/\
                            pop[j+1].finals['M_surface_is']
                        
                        nextnextrelamount = pop[j+2].finals['M_H2O_is']/\
                            pop[j+1].finals['M_surface_is'] 
                       
                        if nextrelamount > min_mass and\
                        not found_dissoc:

                            m_this = planet.finals['M_surface_is']
                            m_left = pop[j-1].finals['M_surface_is']
                            m_right = pop[j+1].finals['M_surface_is']

                            x1 = pop[j+1].finals['M_surface_is']
                            x2 = pop[j+2].finals['M_surface_is']
                            
                            y0 = min_mass
                            y1 = nextrelamount
                            y2 = nextnextrelamount
                            
                            #print ('\n', x1, x2, y1, y2)
                            
                            deltax = x2-x1
                            deltay = y2-y1
                            
                            slope = deltay/deltax
                            #print ('slope =', slope)
                            predicted_pos = (y0-y1) * (x2-x1)/(y2-y1) + x1
                            
                            #if the gap between the planets is too big, perform
                            #linear interpolation to estimate the real transition
                            #position
                            predicted_pos = m_this
                            
                            dissoc_positions.append(predicted_pos)
                            
                            #compute uncertainty of dissociation transition
                            #it is estimated as the difference in total mass
                            #between the current planet with 0 MH2O and the
                            #next planet
                            
                            sigma_mass = (m_right-m_this)/2.
                                    
                            sigma_delta = np.sqrt((m_right*0.001)**2 + \
                                                  (m_this*0.001)**2)

                            dissoc_positions_error.append(\
                                    np.sqrt((sigma_delta)**2+sigma_mass**2))
                            
                            found_dissoc = True
                    
                    except IndexError:
                        pass
    
                #extract transition to fully dissociated
                relamount = planet.finals['M_H2O_hidden_is']/\
                            planet.finals['M_surface_is']
                            
                if relamount > min_mass:
                    try:
                        nextrelamount = pop[j+1].finals['M_H2O_hidden_is']/\
                            pop[j+1].finals['M_surface_is']
                            
                        if nextrelamount < min_mass \
                        and not found_fulldissoc:

                            m_this = planet.finals['M_surface_is']
                            m_left = pop[j-1].finals['M_surface_is']
                            m_right = pop[j+1].finals['M_surface_is']
                            
                            x1 = planet.finals['M_surface_is']
                            x2 = pop[j+1].finals['M_surface_is']
                            
                            y0 = min_mass
                            y1 = relamount
                            y2 = nextrelamount
                            
                            #print ('\n', x1, x2, y1, y2)
                            
                            deltax = x2-x1
                            deltay = y2-y1
                            
                            slope = deltay/deltax
                            #print ('slope =', slope)
                            #if the gap between the planets is too big, perform
                            #linear interpolation to estimate the real transition
                            #position
                            predicted_pos = m_this
                                
                            fulldissoc_positions.append(predicted_pos)
    
                            #compute uncertainty of dissociation transition
                            #it is estimated as the difference in total mass
                            #between the current planet with 0 MH2O and the
                            #next planet
                            
                            sigma_mass = (m_right-m_this)/2.
                                    
                            sigma_delta = np.sqrt((m_right*0.001)**2 + \
                                                  (m_this*0.001)**2)
                            
                            fulldissoc_positions_error.append(\
                                    np.sqrt((sigma_delta)**2+sigma_mass**2))
                            
                            found_fulldissoc = True
                                        
                    except IndexError:
                        pass
                    
            if not found_dissoc:
                print ('WARNING: No dissociation transition found')
                dissoc_positions.append(float('NaN'))
                dissoc_positions_error.append(float('NaN'))
            
            if not found_fulldissoc:
                print ('WARNING: No full dissociation transition found')
                fulldissoc_positions.append(float('NaN'))
                fulldissoc_positions_error.append(float('NaN'))
        
        print ('dissoc_positions_error =', dissoc_positions_error)
        return dissoc_positions, fulldissoc_positions, dissoc_positions_error,\
                fulldissoc_positions_error
    
    
    def analyze(self, all_pops=None, loc='./', plot=False, fnts=10,
                lwdth_data=1, lwdth_pure=2, lwdth_frame=2, 
                ticklen_maj=6, ticklen_min=4,  lwdth_ticks_maj=2, 
                lwdth_ticks_min=1, grid_col='b', fnts_title=14, 
                fnts_axislabel=12, fnts_ticklabels=12, labelpad=20,
                lwdth_dissoc=1, background_col = (.95, .95, 1.), **kwargs):
        #all_pops = self.load_all(loc=loc)
        
        dissoc_positions = []
        fulldissoc_positions = []
        dissoc_positions_error = []
        temps = []
        Mg_numbers = []
        
        for i in range(len(all_pops)):
            pops = all_pops[i]
            sorted_planets = self.preprocess_data(planets=pops)
            diss, fulldiss, diss_error, fulldiss_error = \
            self.extract_dissoc_positions(pops=sorted_planets)
            
            temps.append(pops[0][0].initials['T_surface_should'])
            Mg_numbers.append([])
            dissoc_positions.append(diss)
            fulldissoc_positions.append(fulldiss)
            dissoc_positions_error.append(diss_error)
            
            for j in range(len(pops)):
                Mg_numbers[i].append(pops[j][0].initials['Mg_number_should'])
                #dissoc_positions[i].append(diss[j])      
                #fulldissoc_positions[i].append(fulldiss[j])
                #dissoc_positions_error[i].append(diss_error[j])
                    
        if plot:
            fig, ax = plt.subplots()
            
            ax.set_ylim(-1., 5.)
            ax.xaxis.set_major_locator(MultipleLocator(0.1))
     
            ax.yaxis.set_major_locator(MultipleLocator(1.))
            ax.yaxis.set_minor_locator(MultipleLocator(0.1))

            ax.set_facecolor(background_col)
    
            ax.tick_params(top=True, right=True, direction='in', which='both', 
                           pad=5, labelsize=fnts_ticklabels, labelright=True)
            
            ax.tick_params(top=True, right=True, direction='in', which='major',
                           pad=5, length=ticklen_maj, width=lwdth_ticks_maj,
                           labelright=True)
            
            ax.tick_params(top=True, right=True, direction='in', which='minor',
                           pad=5, length=ticklen_min, width=lwdth_ticks_min,
                           labelsize=8)
    
            ax.grid(True, axis='both', which='major', color=grid_col, alpha=.1,
                    linewidth=2)
            ax.grid(True, axis='both', which='minor', color=grid_col, alpha=.1,
                    linewidth=1)  
            
            ax.set_xlabel(r'$\rm magnesium \ number \ Mg \#$')
            ax.set_ylabel(r'$\rm \delta m_{dissoc}^{280K}/m_{dissoc}^{280K} \ [\%]$')     
            ax.set_title (r'$\rm P_{surface} \ = 100 \ bar$')
            
            for i in range(len(all_pops)):
                
                errors = (np.asarray(dissoc_positions[i])/np.asarray(dissoc_positions[0])*\
                    np.asarray(dissoc_positions_error[i]))**2 + \
                        (np.asarray(dissoc_positions[i])/np.asarray(dissoc_positions[0])**2*\
                    np.asarray(dissoc_positions_error[0]))**2
                         
                errors = np.sqrt(errors)
               
                devs_dissoc = (np.asarray(dissoc_positions[i])-\
                        np.asarray(dissoc_positions[0]))/\
                        np.asarray(dissoc_positions[0])

                ax.errorbar(Mg_numbers[i], 100*devs_dissoc, marker='s', 
                        label=str(int(temps[i]))+r'$ \ K$', 
                        yerr=dissoc_positions_error[i])
                
            ax.legend(loc=2)
                
        return dissoc_positions


    def preprocess_data(self, planets=None, sort_by='M_surface_is'):
        """Check input planets and ommit the ones for which one or several
        of the parameters Mg#, T_surface, M_surface, M_H2O, P_surface are
        not matched with the desired precission. Then return an array containing
        the planets sorted by a specified planetary property.
        """
        sorted_planets = []
        
        #pre-process data
        #iterate over magnesium number
        for i in range(len(planets)):
            sorted_planets.append([])
            
            #sort planets by a specified property
            sorted_planets_dummy=self.sort_planets(planets[i], what=sort_by)
            
            #check if planet has converged to specified properties with
            #desired precission
            self.check_model(planets[i])
            
            #iterate over single planets
            for j in range(len(sorted_planets_dummy)):
                planet = sorted_planets_dummy[j]
                
                #if convergence did succeed, ommit planet, else add it to the
                #output array
                if planet.match:       
                    sorted_planets[i].append(planet)
                    
        return sorted_planets
        
    
    def plot_MR(self, planets, others=[], use = 'object', save = False, 
                filename='MR_diagram', format='png', 
                path='/mnt/c/Users/os18o068/Documents/PHD/Abbildungen/', 
                fnts=10, lwdth_data=1, lwdth_pure=2, lwdth_frame=2, 
                ticklen_maj=6, ticklen_min=4,  lwdth_ticks_maj=2, 
                lwdth_ticks_min=1, grid_col=grid_color, fnts_title=14, 
                fnts_axislabel=12, fnts_ticklabels=12, labelpad=20,
                lwdth_dissoc=1, background_col = background_color, 
                fnts_planetlabels=8, **kwargs):
        """Takes as input a population of planets and plots the correspondinig
        MR-diagram. The plot is saved as a image file if desired. By default
        the input must be a list of Planet.Planet() objects but this might change
        in the future as it would be useful to be able to also use a planet
        output file as input instead of python objects (which requires the
        data to have been acquired in the same interactive session as this
        routine is called, this is clearly a restriction of usage)
        """
        pl_cols = ['grey', 'orange', 'red', 'green', 'k', 'blue', 'magenta', 
                   'purple', (0., 0.4, 0.5)]
        bbox_props = dict(boxstyle = 'square,pad=0.2', fc=background_col, 
                          ec='None')
        
        fig, axis = plt.subplots(2, 1, sharex=True)
        fig2, ax2 = plt.subplots()
        fig3, ax3 = plt.subplots()
        fig4, ax4 = plt.subplots()
        
        additional_axes = [ax2, ax3, ax4]
        
        fig.subplots_adjust(hspace=0)
        ax = axis[0]
        axx = axis[1]

        ax.set_xlim(0., 4.)
        ax.set_ylim(.45, 2.0)

        ax4.set_xlim(0., 4.)
        ax4.set_ylim(.45, 2.25)

        axx.set_xlim(0, 4.)
        axx.set_ylim(.0, 1.1) 
        
        x_lims = [[0., 4.], [0., 4.]]
        y_lims = [[.45, 2.], [0., 1.1]]

        x_minor_labels_dummy = [np.around(np.linspace(-0.1, 4., 42), decimals=2),
                                np.around(np.linspace(-0.1, 4., 42), decimals=2)]
       
        y_minor_labels_dummy = [np.around(np.linspace(.4, 2., 17), decimals=1),
                                np.around(np.linspace(-.1, 1.1, 13), decimals=1)]
        
        x_minor_labels = [[], []]
        y_minor_labels = [[], []]
        
        #create customized major and minor tick labels
        for i in range(len(axis)):
            for n in range(len(x_minor_labels_dummy[i])):
                number = x_minor_labels_dummy[i][n]
                if round(number/.5, 0) == number/.5 or \
                n == len(x_minor_labels_dummy[i])-1:
                    x_minor_labels[i].append('')
                    
                else:
                    x_minor_labels[i].append(str(number))
            
            for n in range(len(y_minor_labels_dummy[i])):
                number = y_minor_labels_dummy[i][n]
                
                if round(number/.5, 0) == number/.5 or \
                n == len(y_minor_labels_dummy[i])-1:
                    y_minor_labels[i].append('')
                    
                else:
                    y_minor_labels[i].append(str(number))                
         
        #set width of frame surrounding the plots
        for pos in ['top', 'bottom', 'left', 'right']:
            ax.spines[pos].set_linewidth(lwdth_frame)
            axx.spines[pos].set_linewidth(lwdth_frame)

        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))

        axx.yaxis.set_major_locator(MultipleLocator(0.5))
        axx.yaxis.set_minor_locator(MultipleLocator(0.1))

        ax2.xaxis.set_major_locator(MultipleLocator(0.5))
        ax2.xaxis.set_minor_locator(MultipleLocator(0.1))
 
        ax2.yaxis.set_major_locator(MultipleLocator(0.1))
        ax2.yaxis.set_minor_locator(MultipleLocator(0.01))

        ax3.xaxis.set_major_locator(MultipleLocator(0.1))
        #ax3.xaxis.set_minor_locator(MultipleLocator(0.01))
 
        ax3.yaxis.set_major_locator(MultipleLocator(0.5))
        ax3.yaxis.set_minor_locator(MultipleLocator(0.1))
       
        letters = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
        
        for i in range(len(axis)):
            
            xx = axis[i]
            
            #set plot background color
            xx.set_facecolor(background_col)
            
            #add panel annotation
            xx.text(0.95, 0.1, letters[i], fontsize = fnts_axislabel,
                        transform = xx.transAxes,
                        fontweight='bold')
            
            xx.tick_params(top=True, right=True, direction='in', which='both', 
                           pad=5, labelsize=fnts_ticklabels, labelright=True)
            
            xx.tick_params(top=True, right=True, direction='in', which='major',
                           pad=5, length=ticklen_maj, width=lwdth_ticks_maj,
                           labelright=True)
            
            xx.tick_params(top=True, right=True, direction='in', which='minor',
                           pad=5, length=ticklen_min, width=lwdth_ticks_min,
                           labelsize=8)

            xx.grid(True, axis='both', which='major', color=grid_col, alpha=.1,
                    linewidth=2)
            xx.grid(True, axis='both', which='minor', color=grid_col, alpha=.1,
                    linewidth=1)
            
            xx.grid(True, axis='both', which='major', color=grid_col, alpha=.1,
                    linewidth=2)
            xx.grid(True, axis='both', which='minor', color=grid_col, alpha=.1,
                    linewidth=1)

            xx.set_yticklabels(y_minor_labels[i], minor=True)
            xx.set_xticklabels(x_minor_labels[i], minor=True)
        
        
        for xx in additional_axes:
            xx.set_facecolor(background_col)
    
            xx.tick_params(top=True, right=True, direction='in', which='both', 
                           pad=5, labelsize=fnts_ticklabels, labelright=True)
            
            xx.tick_params(top=True, right=True, direction='in', which='major',
                           pad=5, length=ticklen_maj, width=lwdth_ticks_maj,
                           labelright=True)
            
            xx.tick_params(top=True, right=True, direction='in', which='minor',
                           pad=5, length=ticklen_min, width=lwdth_ticks_min,
                           labelsize=8)
    
            xx.grid(True, axis='both', which='major', color=grid_col, alpha=.1,
                    linewidth=2)
            xx.grid(True, axis='both', which='minor', color=grid_col, alpha=.1,
                    linewidth=1)        

        legend2_list =[None, None]
        legend1_list = [None, None]
        
        surface_temps = []
        surface_pres = []
        surface_temps_mean = []
        surface_pres_mean = []
        surface_temps_error = []
        surface_pres_error = []
        sigam_T = 0.
        sigma_P = 0.
        mean_T = 0.
        mean_P = 0.
    
        dissoc_positions = []
        fulldissoc_positions = []
        dissoc_positions_error = []
        fulldissoc_positions_error = []
        water_contents_analytical = []
        water_contents = []
        water_contents_mean = []
        water_contents_error = []
        
        sorted_planets = []
        
        #load curves for homogenous planets
        water_planets = self.load_population(loc='./water_planets_T300/',
                                             kind='')[0]
        iron_planets = self.load_population(loc='./iron_planets_T300/',
                                            kind='')[0]
        periclase_planets = self.load_population(loc='./periclase_planets_T300/',
                                                 kind='')[0]
        
        water_planets = self.sort_planets(water_planets, what='M_surface_is')
        iron_planets = self.sort_planets(iron_planets, what='M_surface_is')
        periclase_planets = self.sort_planets(periclase_planets, what='M_surface_is')
        
        water_masses = [pl.finals['M_surface_is'] for pl in water_planets]
        iron_masses = [pl.finals['M_surface_is'] for pl in iron_planets]
        periclase_masses = [pl.finals['M_surface_is'] for pl in periclase_planets]

        water_radii = [pl.finals['R_surface_is'] for pl in water_planets]
        iron_radii = [pl.finals['R_surface_is'] for pl in iron_planets]
        periclase_radii = [pl.finals['R_surface_is'] for pl in periclase_planets]
        
        #pre-process data
        sorted_planets = self.preprocess_data(planets=planets)

        dissoc_positions, fulldissoc_positions, dissoc_positions_error, \
        fulldissoc_positions_error = \
        self.extract_dissoc_positions(pops=sorted_planets)
                
        Mg_numbers = []
        
        for i in range(len(planets)):
            water_contents_dummy = []
            surface_temps_dummy = []
            surface_pres_dummy = []
            
            water_contents_analytical.append([])

            Mg_number = round(sorted_planets[i][-1].initials['Mg_number_should'],
                              3)
            Mg_numbers.append(Mg_number)
            
            for j in range(len(sorted_planets[i])):
                planet = sorted_planets[i][j]

                H2Ofraction =  (planet.finals['M_H2O_is']+\
                                planet.finals['M_H2O_hidden_is'])/\
                                planet.finals['M_surface_is']
                
                H2Ofraction_analytical = (mH*2+mO)/\
                        ((1./Mg_number-1.)*mFe + mMg + 2.*(mO+mH))
                
                water_contents_dummy.append(H2Ofraction)
                water_contents_analytical[i].append(H2Ofraction_analytical)
                surface_temps_dummy.append(planet.finals['T_surface_is'])
                surface_pres_dummy.append(planet.finals['P_surface_is'])

            surface_temps.append(surface_temps_dummy)
            surface_pres.append(surface_pres_dummy)
            surface_temps_mean.append(np.mean(np.array(surface_temps_dummy)))
            surface_pres_mean.append(np.mean(np.array(surface_pres_dummy)))
            surface_temps_error.append(np.std(np.array(surface_temps_dummy))) 
            surface_pres_error.append(np.std(np.array(surface_pres_dummy)))
            
            water_contents.append(water_contents_dummy)
            water_contents_mean.append(np.mean(np.array(water_contents_dummy)))
            water_contents_error.append(np.std(np.array(water_contents_dummy)))
            #print ('standard deviation in water content=', 
             #      water_contents_error[i])
        
        sigma_T = np.std(np.array(surface_temps_mean))
        sigma_P = np.std(np.array(surface_pres_mean))

        mean_T = np.mean(np.array(surface_temps_mean))        
        mean_P = np.mean(np.array(surface_pres_mean))        

        #compute how many digits of surface P and T must be displayed
        #according to the uncertainty of the values
        try:
            digits = abs(int(np.log10(sigma_T)))+1
            
        except OverflowError:
            digits=1
            
        sigma_T = sigma_T * 10**digits
        
        mean_T = round(mean_T, digits)
        sigma_T = int(sigma_T)
        
        #convert pressure in bar
        sigma_P *= 1.0e-5
        mean_P *= 1.0e-5
        
        try:
            digits = int(abs(int(np.log10(sigma_P)))+1)
            
        except OverflowError:
            digits = 1
            
        sigma_P = sigma_P * 10**digits
        
        mean_P = round(mean_P, digits)
        sigma_P = int(sigma_P)
        
        format_str_P = '{0:.'+str(digits)+'f}'
        mean_P_string = format_str_P.format(mean_P)
        
        print ('digits=', digits)
        print ('P_mean =', mean_P)
        print ('sigma_P =', sigma_P)
        
        #loop over subpopulations
        for i in range(len(planets)):
                        
            MH2O_in = []
            MH2O_out = []
            MH2O_tot = []

            #extract Mg#
            Mg_number = round(sorted_planets[i][-1].initials['Mg_number_should'],
                              2)
            
            #compute angle for Mg# label using the slope between the last
            #two data points where the label is to be displayed
            #angle for upper panel---------------------------------------------           
            x, y = np.array(((sorted_planets[i][-3].finals['M_surface_is'],
                             sorted_planets[i][-2].finals['M_surface_is']),
                            (sorted_planets[i][-3].finals['R_surface_is'],
                             sorted_planets[i][-2].finals['R_surface_is'])))

            x_pos = sorted_planets[i][-2].finals['M_surface_is']
            y_pos = sorted_planets[i][-2].finals['R_surface_is']
        
            #label position
            xylabel1 = (x_pos+.15, y_pos+.075 - .01*(1.-np.sqrt(Mg_number)))
            xylabel2 = (x_pos-.2, y_pos+.025)
            
            p1 = ax.transData.transform_point((x[0], y[0]))
            p2 = ax.transData.transform_point((x[1], y[1]))
            
            dx = (p2[0] - p1[0])
            dy = (p2[1] - p1[1])
            
            rotn = np.degrees(np.arctan2(dy, dx))*np.sqrt(1.-Mg_number**2)-\
                                0.02*1./Mg_number**3

            error = water_contents_error[i]

            #compute how many digits of water contents must be displayed
            #according to the uncertainty of the value            
            digits = abs(int(np.log10(error)))+1
            error = error * 10**digits

            #add Mg# label to plot
            ax.annotate(r'$\mathbf{\tildeM_{H_2 O} \ =}$'+\
                        str(round(water_contents_mean[i]*100, digits-2))+'('+
                        str(int(error))+')%',
                        xy=xylabel1,
                        color = pl_cols[i], 
                        fontsize = fnts, 
                        rotation = rotn,
                        fontweight='bold',
                        ha='center',
                        va='center')
                    
            ax.annotate(r'$\mathbf{Mg\# \ =}$'+str(Mg_number),
                        xy=xylabel2,
                        color = pl_cols[i], 
                        fontsize = fnts, 
                        rotation = rotn,
                        fontweight='bold',
                        ha='center',
                        va='center')
            
            #angle for lower panel---------------------------------------------
            x, y = np.array(((sorted_planets[i][-3].finals['M_surface_is'],
                             sorted_planets[i][-2].finals['M_surface_is']),
                            (sorted_planets[i][-3].finals['M_H2O_is'],
                             sorted_planets[i][-2].finals['M_H2O_is'])))
            
            x_pos = sorted_planets[i][-2].finals['M_surface_is']
            y_pos = sorted_planets[i][-2].finals['M_H2O_is']
            
            #label position
            xylabel1 = (x_pos+.15, y_pos+.065*Mg_number+.05)
            xylabel2 = (x_pos-.2, y_pos-.01*Mg_number + 0.025 + 1./Mg_number**2*0.001)
            
            p1 = axx.transData.transform_point((x[0], y[0]))
            p2 = axx.transData.transform_point((x[1], y[1]))
            
            dx = (p2[0] - p1[0])
            dy = (p2[1] - p1[1])
            
            rotn = np.degrees(np.arctan2(dy, dx))*np.sqrt(1.-Mg_number**2)-\
                                0.02*1./Mg_number**3
            
            #add Mg# label to plot
            axx.annotate(r'$\mathbf{\tildeM_{H_2 O} \ =}$'+\
                        str(round(water_contents_mean[i]*100, digits-2))+'('+
                        str(int(error))+')%',
                        xy=xylabel1,
                        color = pl_cols[i], 
                        fontsize = fnts, 
                        rotation = rotn,
                        fontweight='bold',
                        ha='center',
                        va='center')
                    
            axx.annotate(r'$\mathbf{Mg\# \ =}$'+str(Mg_number),
                        xy=xylabel2,
                        color = pl_cols[i], 
                        fontsize = fnts, 
                        rotation = rotn,
                        fontweight='bold',
                        ha='center',
                        va='center')
                    
            #loop over planets in subpopulation--------------------------------
            for j in range(len(sorted_planets[i])):
                
                planet = sorted_planets[i][j]
                
                MH2O_in.append(planet.finals['M_H2O_hidden_is'])#/\
                               #planet.finals['M_surface_is'])
                
                MH2O_out.append(planet.finals['M_H2O_is'])#/\
                               #planet.finals['M_surface_is'])
                
                MH2O_tot.append(MH2O_in[j] + MH2O_out[j])

                edge=pl_cols[i]
                face=pl_cols[i]
                
                if planet.vapor_reached:
                    marker = 's'
                    
                else:
                    marker = 'o'
                    
                if planet.density_inversion:
                    face = 'r'
                    
                if planet.match == False:
                    marker = 's'
                    face = 'None'
                
                #plot the individ planets
                try:
                    ax.scatter(planet.finals['M_surface_is'], 
                               planet.finals['R_surface_is'],
                               marker=marker, edgecolor=edge,
                               facecolor=face, s=25)
                    
                except TypeError:
                    pass
            
            
            try:
                #plot MR curve
                ax.plot([planet.finals['M_surface_is'] 
                         for planet in sorted_planets[i]], 
                        [planet.finals['R_surface_is']
                        for planet in sorted_planets[i]],
                        linewidth=lwdth_data,
                        color=edge)
                
                ax4.plot([planet.finals['M_surface_is'] 
                         for planet in sorted_planets[i]], 
                        [planet.finals['R_surface_is']
                        for planet in sorted_planets[i]],
                        linewidth=lwdth_data,
                        color=edge)
                
                #plot position of dissociation onset
                pos = dissoc_positions[i]
                error = dissoc_positions_error[i]
                
                #compute how many digits of pos must be displayed according
                #to the uncertainty of the value
                try:
                    digits = abs(int(np.log10(error)))+1
                    error = error * 10**digits
                    
                    plot_val_string = str(round(pos,digits))
                    plot_err_string = str(int(error))
                    
                except ValueError:
                    digits = 0
                    plot_val_string = ''
                    plot_err_string = ''
                    pass
                
                line1a, = ax.plot([pos, pos],
                        [0., 5.], color = pl_cols[i], linestyle=(0, (1,1)),
                        linewidth=lwdth_dissoc)
                
                axx.plot([pos, pos],
                        [0., 2.], color = pl_cols[i], linestyle=(0, (1,1)),
                        linewidth=lwdth_dissoc)
                
                #add dissoc label to plot
                label_pos = pos
                try:
                    dev = abs(pos-dissoc_positions[i-1])/pos
                    if dev < .25:
                        label_pos -= .02
                        
                except IndexError:
                    pass
                
                axx.annotate(plot_val_string+'('+plot_err_string+')',
                        xy=(label_pos, .91),
                        color = pl_cols[i], 
                        fontsize = fnts, 
                        rotation = 90.,
                        fontweight='bold',
                        ha='center',
                        va='center',
                        bbox = bbox_props)
                    
                #plot position of full dissociation onset
                pos = fulldissoc_positions[i]
                error = fulldissoc_positions_error[i]
                
                #compute how many digits of pos must be displayed according
                #to the uncertainty of the value
                try:
                    digits = abs(int(np.log10(error)))+1
                    error = error * 10**digits
                    
                    plot_val_string = str(round(pos,digits))
                    plot_err_string = str(int(error))
                    
                except ValueError:
                    digits = 0
                    plot_val_string = ''
                    plot_err_string = ''
                    pass
                
                line1b, = ax.plot([pos, pos],
                        [0., 5.], color = pl_cols[i], linestyle=(0, (1,2,5,2)),
                        linewidth=lwdth_dissoc)
            
                axx.plot([pos, pos],
                        [0., 2.], color = pl_cols[i], linestyle=(0, (1,2,5,2)),
                        linewidth=lwdth_dissoc)
                
         
                axx.annotate(plot_val_string+'('+plot_err_string+')',
                        xy=(pos, .91),
                        color = pl_cols[i], 
                        fontsize = fnts, 
                        rotation = 90.,
                        fontweight='bold',
                        ha='center',
                        va='center',
                        bbox=bbox_props)
                
                legend1_list = [line1a, line1b]
                    
                #plot water contents
                pl2a, = axx.plot([planet.finals['M_surface_is'] 
                          for planet in sorted_planets[i]], 
                        MH2O_out, color=edge,
                        linestyle='--',
                        linewidth = lwdth_data,
                        marker='')
                
                pl2b, = axx.plot([planet.finals['M_surface_is']
                          for planet in sorted_planets[i]], 
                         MH2O_tot,
                         color=edge, marker='', 
                         linestyle='-',
                         linewidth=lwdth_data)
                          
                legend2_list = [pl2a, pl2b]
                
            except TypeError:
                print ('error')
                pass
        
        
        #plot solar system bodies
        ax4.scatter(M_solar, R_solar, color='k', zorder=100)
        
        #plot trappist 1
        ax4.errorbar(M_trappist1, R_trappist1, xerr=M_trappist1_error, 
                    yerr=R_trappist1_error, linestyle = '', marker = 'o',
                    color='k', zorder=100)
        
        x_error_others = np.array(M_others_error).T*m_jupiter/m_earth
        y_error_others = np.array(R_others_error).T*r_jupiter/r_earth
        
        #plot other planets
        ax4.errorbar(np.asarray(M_others)*m_jupiter/m_earth, 
                     np.asarray(R_others)*r_jupiter/r_earth, 
                     xerr=x_error_others, 
                    yerr=y_error_others,
                    linestyle = '', marker = 'o',
                    color='k', zorder=100)        

        for i in range(len(R_trappist1)):
            ax4.text(M_trappist1[i]*1., R_trappist1[i]*1.05, names_trappist1[i],
                    size=fnts_planetlabels, color='k')

        for i in range(len(R_solar)):
            if not M_solar[i] < .1 and not R_solar[i] < .1:
                ax4.text(M_solar[i]*1., R_solar[i]*1.05, names_solar[i],
                    size=fnts_planetlabels, color='k', zorder = 10000)
                
        for i in range(len(R_others)):
            ax4.text(M_others[i]*m_jupiter/m_earth, 
                         R_others[i]*r_jupiter/r_earth*1.05, names_others[i],
                        size=fnts_planetlabels, color='k', zorder = 10000)
        
        #plot pure materials lines
        ax.plot(water_masses, water_radii, color = 'blue', linewidth=lwdth_pure,
                alpha=.5)
        ax.plot(iron_masses, iron_radii, color = (.2, .2, .2), linewidth=lwdth_pure,
                alpha=.5)
        ax.plot(periclase_masses, periclase_radii, color = 'brown',
                linewidth=lwdth_pure, alpha=.5)

        ax4.plot(water_masses, water_radii, color = 'blue', linewidth=lwdth_pure,
                alpha=.5)
        ax4.plot(iron_masses, iron_radii, color = (.2, .2, .2), linewidth=lwdth_pure,
                alpha=.5)
        ax4.plot(periclase_masses, periclase_radii, color = 'brown',
                linewidth=lwdth_pure, alpha=.5)
        
        ax.text(1.5, .75, r'$\mathbf{Fe}$', color = (.2, .2, .2), alpha=.5,
                fontsize = fnts+2)
        ax.text(2.82, 1.4, r'$\mathbf{MgO}$', color = 'brown', alpha=.5,
                fontsize = fnts+2)
        ax.text(1., 1.55, r'$\mathbf{H_2 O}$', color = 'blue', alpha=.5,
                fontsize = fnts+2)

        ax4.text(1.5, .75, r'$\mathbf{Fe}$', color = (.2, .2, .2), alpha=.5,
                fontsize = fnts+2)
        ax4.text(2.82, 1.4, r'$\mathbf{MgO}$', color = 'brown', alpha=.5,
                fontsize = fnts+2)
        ax4.text(1., 1.55, r'$\mathbf{H_2 O}$', color = 'blue', alpha=.5,
                fontsize = fnts+2)

        ax.set_title(r'$\rmT_{surface} \ = $'+str(mean_T)+'('+str(sigma_T)+')'+
                     r'$K, \ \rmP_{surface} \ = $'+mean_P_string+'('+
                     str(sigma_P)+')'+r'$ \ bar$', fontsize=fnts_title)

        ax.set_ylabel(r'$\rmR_{Pl}/R_{\oplus}$', fontsize=fnts_axislabel,
                       labelpad=labelpad)
        axx.set_xlabel(r'$\rmM_{Pl}/M_{\oplus}$', fontsize=fnts_axislabel,
                        labelpad=labelpad)
        axx.set_ylabel(r'$\rmwater \ content \ [m/M_{\oplus}]$',
                       fontsize=fnts_axislabel, labelpad=labelpad)

        ax4.set_ylabel(r'$\rmR_{Pl}/R_{\oplus}$', fontsize=fnts_axislabel,
                       labelpad=labelpad)
        ax4.set_xlabel(r'$\rmM_{Pl}/M_{\oplus}$', fontsize=fnts_axislabel,
                        labelpad=labelpad)
        
        legend1 = ax.legend(legend1_list, [r'$\rmdissociation \ onset$', 
                                           r'$\rmfull \ dissociation \ onset$'],
                                            frameon=True, loc=2,
                                            fontsize=fnts_ticklabels,
                                            framealpha=1.,
                                            edgecolor='None')            
        
        legend2 = axx.legend(legend2_list, [r'$\rmm \ = \ m_{H_2 O, out}$', 
                                            r'$\rmm \ = \ m_{H_2 O, out}+m_{H_2 O, in}$'],
                                            frameon=True, loc=(.025, .3),
                                            fontsize=fnts_ticklabels,
                                            framealpha=1.,
                                            edgecolor='None')
        ax.add_artist(legend1)
        axx.add_artist(legend2)
        ax.set_axisbelow(True)          
        
        #save figure as image file if turned on
        if save:
            fig.savefig(path+filename+'.'+format)
    
        
        for i in range(len(planets)):
            masses = np.array([planet.finals['M_surface_is']
            for planet in sorted_planets[i]])
            
            Mg_number = round(sorted_planets[i][0].finals['Mg_number_is'], 3)
     
            diff=np.array((np.array(water_contents[i])-\
                           np.array(water_contents_analytical[i]))/\
                            np.array(water_contents_analytical[i]))
            '''
            diff_trans = []
            for j in range(len(diff)):
                diff_trans.append([diff[j]])
            
            diff_trans, min, max, sign = logTrans.transform(data=diff_trans,
                                                            param = 0, dim=1,
                                                            def_zero=1.0e-6)
            
            exp_max_plot = round(np.nanmax(diff_trans)+.5, 0)
            exp_min_plot = round(np.nanmin(diff_trans)-.5, 0)
            
            if abs(exp_max_plot) > abs(exp_min_plot):
                exp_min_plot = -exp_max_plot
            '''
           
            ax2.plot(masses, diff*100, pl_cols[i])
            
            #add Mg# label to plot
            ax2.annotate(r'$\mathbf{Mg\# \ =}$'+str(Mg_number),
                        xy=(1., 100*(0.0065-0.00075*i)),
                        color = pl_cols[i], 
                        fontsize = fnts, 
                        rotation = 0,
                        fontweight='bold',
                        ha='left',
                        va='bottom')
                    
            ax4.annotate(r'$\mathbf{Mg\# \ =}$'+str(round(Mg_number,2)),
                        xy=(.2, (2.-0.075*i)),
                        color = pl_cols[i], 
                        fontsize = fnts, 
                        rotation = 0,
                        fontweight='bold',
                        ha='left',
                        va='bottom')
            
        ax2.plot([0., 4.], [0., 0.], color='k')
        
        ax2.set_xlim(0., 4.)
        ax2.set_ylim(-.2, .2)        

        ax2.set_ylabel(r'$\rm numerical \ uncertainty \ \delta \tilde{M}_{H_2 O}/ \tilde{M}_{H_2 O} \ [\%]$',
                       fontsize=fnts_axislabel)
        ax2.set_xlabel(r'$\rm M_{Pl}/M_{\oplus}$',
                       fontsize=fnts_axislabel)

        ax2.set_title(r'$\rmT_{surface} \ = $'+str(mean_T)+'('+str(sigma_T)+')'+
                     r'$K, \ \rmP_{surface} \ = $'+mean_P_string+'('+
                     str(sigma_P)+')'+r'$ \ bar$', fontsize=fnts_title)
        
        #plot dissociation position as function of Mg#
        ax3.errorbar(Mg_numbers, dissoc_positions, marker='s', color='b',
                 zorder=100, yerr=dissoc_positions_error)
        
        ax3.errorbar(Mg_numbers, fulldissoc_positions, marker='s', color='orange',
                 zorder=100, yerr=fulldissoc_positions_error)
        
        ax3.set_xlabel(r'$\rm magnesium \ number \ Mg \#$', 
                       fontsize=fnts_axislabel)
        ax3.set_ylabel(r'$\rm dissociation \ position \ [m/M_{\oplus}]$', 
                       fontsize=fnts_axislabel)
        
        ax3.text(.65, .6, r'$\rm m \ = \ m_{onset}$', color = 'blue',
                 fontsize=fnts_axislabel)
        ax3.text(.65, 2.1, r'$\rm m \ = \ m_{full}$', color = 'orange',
                 fontsize=fnts_axislabel)

        ax3.set_title(r'$\rmT_{surface} \ = $'+str(mean_T)+'('+str(sigma_T)+')'+
                     r'$K, \ \rmP_{surface} \ = $'+mean_P_string+'('+
                     str(sigma_P)+')'+r'$ \ bar$', fontsize=fnts_title)                    
        
       
    def write(self, planet, name='planet.out'):
        """Creates data table of the planetary structure and writes it to ascii 
        file. By default, the data is organized as follows:
        
        R           M           T       P       rho         v_esc       g
        (r_earth)   (m_earth)   (K)     (GPa)   (kg/m3)     (km/s)      (m/s2)
        -----------------------------------------------------------------------        
        R_center    M_center    ...
        ...         ...         ...
        R_tot       M_tot       ...
   
        The meta data contains the following parameters by default:
            
            contents
            fractions
            layermasses (m_earth)
            P_center (GPa)
            T_center (K)
            tempType
            majorConstraint
            eps_r
            rhoType
            R_seed (m)
            
        This method is used to write a planet object to a file without using
        it's native write method. This allows to modify the write method after
        the planet object has been initialized, if required for wathever reason.
        """
        
        rad, mass, temp, pres, dens, vesc, grav = [], [], [], [], [], [], []
        
        #gather structure data
        for lay in planet.layers:
            for sh in lay.shells:
                rad.append(sh.radius/r_earth)
                mass.append(sh.mass/m_earth)
                temp.append(sh.temp)
                pres.append(sh.pres*1.0e-9)
                dens.append(sh.dens)
                vesc.append(sh.v_esc/1000.)
                grav.append(sh.gravity)
        
        data = [rad, mass, temp, pres, dens, vesc, grav]
        
        #gather all planetary parameters as meta data
        layer_properties = [{'P_out':None, 'T_out':None, 'rho_out':None,
                                   'indigenous_mass':0.} 
                                    for i in range(len(planet.layers))]
            
        for i in range(len(planet.layers)):
            layer_properties[i]['P_out'] = planet.layers[i].pres
            layer_properties[i]['T_out'] = planet.layers[i].temp
            layer_properties[i]['rho_out'] = planet.layers[i].dens
            layer_properties[i]['indigenous_mass'] = \
                                        planet.layers[i].indigenous_mass
        
        planet.finals['layer_properties'] = layer_properties
        
        meta = {'initials':planet.initials, 'finals':planet.finals}
    
        data_table=astropy.table.Table(data, meta=meta)
            
        names = ('R (r_earth)',
                 'M (m_earth)',
                 'T (K)',
                 'P (GPa)',
                 'rho (kg/m3)',
                 'v_esc (km/s)',
                 'g (m/s2)')    
        
        ascii.write(data_table, './'+name, overwrite=False,
                   format='ecsv', names = names)
        
        
    def sort_planets(self, planets=None, what='M_surface_is'):
        sorted_planets = []
        unsorted_params = []
        
        planets_dummy = copy.deepcopy(planets)
        
        for planet in planets:
            try:
                param = planet.initials[what]
                dont = False
                
            except KeyError:
                try:
                    param = planet.finals[what]
                    dont = False
                    
                except KeyError:
                    print ('ERROR: invalid planet property given')
                    dont = True
                    
            unsorted_params.append(param)
        
        if not dont:
            i = 0
            while len(planets_dummy) > 0:
                i += 1
                min_index = unsorted_params.index(min(unsorted_params))
                
                sorted_planets.append(planets_dummy[min_index])
                planets_dummy.pop(min_index)
                unsorted_params.pop(min_index)
    
        return sorted_planets
            
    
    def write_population(self, population=[], loc='./'):
        """Takes an array containing a population of Planet.Planet objects
        and invokes the Planet.prt() method for each object to write them
        into files in the given location. By default the output files will be
        labeled as planetX.out in the directory loc/populationY where Y is the
        number of the subpopulation the planet belongs to and X is the number
        of the planet within this subpopulation
        
        Example:
        
            population = [[planet00, planet01], [planet10, planet11]]
            loc = './'
            
            will produce the following output by default:
                
                ./population1/planet1.pl
                ./population1/planet2.pl

                &

                ./population2/planet1.pl
                ./population2/planet2.pl   
        """
        
        for p in range(len(population)):
            pop_str = 'Mg0'+\
            str(int(10*population[p][0].initials['Mg_number_should']))
            pop_path = loc+'population'+pop_str
            print ('pop_path=', pop_path)
            os.mkdir(pop_path)
            
            
            for i in range(len(population[p])):
                planet_str = str(i+1)
                
                population[p][i].write(out='planet'+planet_str,
                          loc=pop_path+'/')
                
                
    def load_population(self, loc='./', kind='Mg'):
        """Takes all planet.out files
        """
        
        pop = []
        i = 0
        for root, dirs, files in os.walk(loc):
            for dir in dirs:
                if dir.startswith('population'+kind):
                    pop.append([])
                    i+=1
                    print ('dir=', dir)
                    print ('root=', root)
                    for dd in os.listdir(root+dir):
                        if dd.endswith('.pl'):
                            file = dd.split('.')[0]
                            pl = Planet.Planet(contents=[[9]],
                                               layermasses=[1.])
                            pl.load(loc=root+dir+'/', file_name=file)
                            pop[i-1].append(pl)
                            
        return pop
    
    
    def pure_curve(self, N=10, ll=1, P_min=1.0e9, P_max=1.0e12, T_surface=300.,
                   P_surface=1.0e5, dir_out ='./', file_name=None):
        """Generates a MR-curve with N planets at given surface conditions
        """
        
        if file_name == None:
            file_name = 'pure_curve_'+str(ll)+'.atab'
        
        else:
            pass
        
        data = np.empty([N, 2])
        
        pres = np.logspace(np.log10(P_min), np.log10(P_max), N)
        
        if ll==2 or ll==8:
            Fe_number_layers = [1.]
        else:
            Fe_number_layers = [1.0e-10]
        
        Si_number_layers = [0.]
        
        print (pres)
        for i in range(N):
            
            pl = PlanetFort.Planet(contents=[[ll]],
                                   fractions=[[1.]],
                                   T_surface_should = T_surface,
                                   P_surface_should = P_surface,
                                   T_center = 2000.,
                                   P_center = pres[i],
                                   layermasses=[1.],
                                   layerpres=[0.],
                                   layerradii=[0.],
                                   layertemps=[0.],
                                   temp_jumps = [0.],
                                   layerConstraint = [0],
                                   Fe_number_layers = Fe_number_layers,
                                   Si_number_layers = Si_number_layers,
                                   gammas_layers = [2.],
                                   q = [1.])
            
            pl.Construct()
            
            data[i][0] = pl.M_surface_is
            data[i][1] = pl.R_surface_is
        
        data_table=astropy.table.Table(data, meta={'material':ll,
                                                    'T_surface':T_surface,
                                                    'P_surface':P_surface,
                                                    'N':N})
        
        ascii.write(data_table, dir_out+file_name, overwrite=True,
                   format='ecsv')       
        return data

            
    
    def load_fortpopulation(self, loc='./', file_name=None):
        if not file_name.endswith('.pop'):
            print ('WARNING: Got unexpected file format for planet population.')
        
        t0 = time.time()
        f = open(loc+file_name, 'r')
        dat = []
        contents = f.readlines()
        
        for line in contents:
            line = line.split(' ')
            line.pop(0)
            line[-1] = line[-1][:len(line[-1])-1]
            
            try:
                dat.append([float(l) for l in line])

            except ValueError:
                print ('ValueError')
            
        t = time.time()
        print ('time:', round(t-t0,3), 'sec')
        
        return np.array(dat)
    
    
    def plot_fortpopulation(self, data, meta=4):
        
        pure_curve_labels=['1', '2', '4']
        pure_curves = []        
        for lab in pure_curve_labels:
            pure_curves.append(ascii.read('pure_curve_'+lab+'.atab', format='ecsv'))
        
        c = 0
        pure_curves_dat = np.empty([len(pure_curves), 2, len(pure_curves[0])])
        for i in range(len(pure_curves)):
            c=0
            while c<len(pure_curves[i]):
                pure_curves_dat[i][0][c] = pure_curves[i][c][0]/m_earth
                pure_curves_dat[i][1][c] = pure_curves[i][c][1]/r_earth
                c+= 1
        
        fig, ax = plt.subplots()
        ax.set_xlabel(r'$M/M_{\oplus}$')
        ax.set_ylabel(r'$R/R_{\oplus}$')
        cm = plt.cm.get_cmap('rainbow_r')
        z = 1.-data.T[meta]/data.T[0]
        
        sc = ax.scatter(data.T[0], data.T[1], c=z, cmap=cm)
    
        locs = np.linspace(0., 1., 11, endpoint=True)
         
        cbar = plt.colorbar(sc, ticks=locs)
        cbar.set_label(r'$M_{\rm Ocean}/M$')
        #cbar.set_ticks(locs)
        
        for pc in pure_curves_dat:
            ax.plot(pc[0], pc[1])
        
    
    
    def load_all(self, loc='./'):
        """Finds all MR_curve_T... directories in loc and loads the population
        files they contain into arrays of subpopulations.
        """
        
        #collect all directories containing population outputs in an array
        folders = []
        for root, dirs, files in os.walk(loc):
            for dir in dirs:
                if dir.startswith('MR_curve_T'):
                    folders.append(dir+'/')
                    
        populations = []
        for dir in folders:
            pop = self.load_population(loc=dir)
            populations.append(pop)
            
        return populations
        
        
    def write_latex_table(self, data=None, adendum = ['no', 'no'],
                          skip=1):
        #Number of populations
        N = len(data)
        
        #Number of temperatures
        M = len(data[0])
        
        #Number of Mg#
        K = len(data[0][0])
        
        #Number of planets per MR relation
        L = len(data[0][0][0][0])
        
        out = np.zeros([N*M*K*int(L/skip), 10])
        
        c = 0
        for i in range(N):
            for j in range(M):
                for k in range(K):
                    dat = data[i][j][k]
                    
                    skip_count = 0
                    
                    for l in range(L):
                        print (c)
                        if c == N*M*K*int(L/skip):
                            break
                        
                        skip_count += 1
                        if skip_count == skip:
                            skip_count = 0
                            #total mass
                            out[c][0] = ftool.fancyround(dat[0][l], digits=4)
                            #total radius
                            out[c][1] = dat[1][l]
                            #P_surf
                            out[c][2] = dat[12][l]*1.0e-5
                            #T_surf
                            out[c][3] = dat[13][l]
                            #P_center
                            out[c][4] = dat[4][l]*1.0e-9
                            #T_center
                            out[c][5] = dat[5][l]
                            
                            #Mg_number
                            out[c][6] = dat[8][l]  
        
                            #M_core
                            out[c][7] = dat[11][l]/m_earth
                            
                            #M_H2O
                            out[c][8] = dat[9][l]
                            
                            #M_ocean
                            out[c][9] = dat[10][l]
                        
                            c += 1    
        
        fmt = ['%1.3e',  '%1.3e', '%1.3e',  '%1.1f', '%1.3e',  '%1.3e', 
               '%1.1f',  '%1.3e', '%1.2e',  '%1.2e']
        np.savetxt('table.csv', out, delimiter=' & ', 
                    fmt=fmt, newline='\\\\\n')


    def predict(self, M, Mg):
        a = [12.05797, .9219442, -.245424, .0621981]
        b = [-0.08695424, 1.00000004, -0.06693844, -1.05711204]
        
        Pc = a[0] + a[1]*np.log10(M) + a[2]*np.exp(Mg) +a[3]*np.log10(M)**2
        Pc = 10**Pc
        
        Mc = b[0] + b[1]*np.log10(M) + b[2]*Mg + b[3]*Mg**2
        Mc = 10**Mc
        
        return Pc, Mc
        
        