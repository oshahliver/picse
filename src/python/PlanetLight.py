# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 14:05:05 2018

@author: os18o068
"""

import matplotlib.ticker
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, \
        AutoMinorLocator, LogLocator, FixedLocator
from PIMPphysicalparams import T0_list, K0_list, K0prime_list, rho0_list, \
                                aP_list, molar_mass_list, Rgas, mH, mO, kB, G,\
                                EOS_type, material_list, EOSstring_list, NA, \
                                P_earth, m_earth, r_earth, material_plot_list,\
                                solar_system, abbrevations_solar, P_zero, \
                                T_zero, mFe, mMg, mSi, mH2O, Mg_number_solar,\
                                temperature_jumps, mS, material_list_fort, material_YMg, \
                                material_YSi, material_YO, material_YH, \
                                material_YS
                                
from PIMPrunparams import eps_Psurf, eps_Mtot, eps_layer, param_colors, \
                    plot_units, suffix, eps_Mg_number, plot_params, eps_Tsurf,\
                    color_list
import numpy as np
import time
import functionTools as ftool
import sys
from tabulate import tabulate
from decimal import Decimal
import astropy.table
from astropy.io import ascii
import warnings
import plotTools



class Planet():
    def __init__(self, P_center=1.0e12, 
                 T_center=3000, 
                 R_seed=.1, 
                 contents=[[2]],
                 tempType = 1, 
                 layerType=1, 
                 fractions=[[1.]], 
                 P_surface_should=P_earth, 
                 majorConstraint='P_surface', 
                 layerradii = [], 
                 layerConstraint=[1], 
                 add_atmosphere=False, 
                 eps_r=.5, 
                 layermasses=[], 
                 layerpres=[],
                 M_surface_should=1., 
                 T_surface_should=300, 
                 R_surface_should=1., 
                 Si_number_layers=[],
                 Fe_number_layers=[],
                 rhoType=1, 
                 echo=False, 
                 brucite_dissoc=False, 
                 vapor_stop=False, 
                 Mg_number_should=Mg_number_solar, 
                 Si_number_should=0.,
                 Mg_number_layer = False, 
                 subphase_res = 32,
                 modelType=0,
                 temp_jumps=[0., 0., 0., 0., 0.], 
                 gammas_layer = [2.43, 2.43, 1.96, 1.26, 1.0], 
                 layertemps = [0, 0, 0, 0, 0],
                 T_zero = [-1000, -1000, -1000., -1000., -1000.],
                 inner_core_frac=None,
                 M_ocean_should = 0.,
                 M_core_should = None, 
                 ocean_frac_should=-10.,
                 adiabatType=1,
                 q = [.489, .489, 2.5, 2.9, 1.0],
                 silence = False, 
                 measureTime=False,
                 eps_T_zero = 0.0,
                 eps_Al=0., 
                 eps_H2O=0.,
                 omega= 0.,
                 xi_H_core_predicted=0.,
                 xi_Stv=0., 
                 **kwargs):
        
        self.measureTime = measureTime
        
        #print ('initializing planet ...')
        self.status='shadow of the future'
        self.M_H2 = 0.
        self.xi_H_core_predicted = xi_H_core_predicted

        #dissociation mechnism for Mg(OH)2 -> MgO + H2O in phase 2
        self.brucite_dissoc = brucite_dissoc
        
        #if set to True construction is stopped if vapor phase is reached
        self.vapor_stop = vapor_stop
        self.M_tot_should = None
        self.R_tot_should = None
        self.M_tot_is = None
        self.R_tot_is = None
        self.R_seed = R_seed
        self.rho_center=None
        self.P_center = P_center
        self.T_center = T_center
        self.P_surface_should = P_surface_should
        self.M_surface_should = M_surface_should*m_earth
        self.P_space_should = P_surface_should
        self.P_space_is = None
        self.T_surface_should=T_surface_should
        self.R_surface_should=R_surface_should*r_earth
        self.contents = contents
        self.temp_jumps = temp_jumps
        self.materials = [[material_list[i] for i in con] for con in contents]
        self.fractions = fractions
        self.v_esc_surface = 0.
        self.gravity = 0.
        self.Mg_number_is = 0.
        self.Mg_number_should = Mg_number_should
        self.Si_number_is = 0.
        self.Si_number_should = Si_number_should
        self.Si_number_layers = Si_number_layers
        self.Fe_number_layers = Fe_number_layers
        self.echo = echo
        self.M_H2O_should = 0.
        self.M_H2O_is = 0.
        self.M_H2O_hidden_is = 0.0
        self.delta_M_H2O = 0.
        self.delta_M_H2O_hidden = 0.
        self.H2O_count = 0.0
        self.Fe_count = 0.0
        self.Si_count = 0.
        self.Al_count = 0.
        self.H_count = 0.
        self.O_count = 0. 
        self.Mg_count = 0.
        self.density_inversion=False
        self.match=True
        self.vapor_reached=False
        self.N_shells=0
        self.test_variable=0.0
        self.Mg_count = 0.0
        self.delta_Mg_count = 0.0
        self.Mg_number_layer = Mg_number_layer
        self.subphase_refined = False
        self.subphase_res = subphase_res
        self.subphase_refine_inhibit = False
        self.gammas_layer = gammas_layer
        self.adiabatType = adiabatType
        self.eps_T_zero = eps_T_zero
        self.eps_H2O = eps_H2O
        self.eps_Al = eps_Al
        self.xi_Stv = xi_Stv
        
        self.M_outer_mantle = 0.
        self.M_ocean_is = 0.
        self.M_DWR_is = 0.
        self.M_ocean_should = M_ocean_should
        self.M_core_should = M_core_should
        self.M_core_is = 0.
        self.ocean_frac_should = ocean_frac_should
        self.ocean_frac_is = -10
        self.T_zero = T_zero
        self.omega=omega
        self.ocean_depth = 0.
        
        #defines the physical model that is invoked for the planet
        #e.g Brucite dissociation model without silicates: modelType = 0
        # " " " with hydrated OlivFalseine: modelType = 1
        self.modelType = modelType 

        self.layermasses = layermasses
        self.layerradii = layerradii
        self.layertemps = layertemps
        self.layerpres = layerpres
        self.d0 = []
        self.q = q
        self.xi_H_core = 0.
        self.P_H2_CMB = 0.
        self.M_H2O_core = 0.
        
        #It is important that the inner core fraction is accounted for properly
        #If a BC are to met via toolkit.iterate(), the inner core fraction
        #must be updated each time the core mass is updated otherwise the
        #iteration will in many cases crash and no solution can be found
        if not inner_core_frac == None:
            self.inner_core_frac = inner_core_frac
            self.layermasses[0] = self.layermasses[1]*\
                                    inner_core_frac/(1.-inner_core_frac)
            #self.layermasses[1] = self.M_core_should*(1.-self.inner_core_frac)
        
        else:
            try:
                M_core =  self.layermasses[0]+self.layermasses[1]
                self.inner_core_frac = self.layermasses[0]/M_core
                
            except IndexError:
                self.inner_core_frac = inner_core_frac
            
        self.eps_r = eps_r
        
        if len(self.contents) > 0:
            self.differentiated = True
            
        else:
            self.differentiated = False
        
        self.shellIterationCount=0
        self.layerIterationCount=0
        self.changebisec=True
        self.shellIteration=True
        self.layerIteration=True

        self.tempType = tempType
        self.layerType = layerType
        self.rhoType = rhoType
        self.majorConstraint = majorConstraint
        self.layerConstraint = layerConstraint
        self.layers = []
        self.lay = 0
        self.test = False
        self.atmosphere=None
        self.add_atmosphere = add_atmosphere
        self.have_atmosphere=False
        self.majorComplete=False
        self.layerComplete=False
        self.layer_properties = [{'P_outer':None,
                                  'T_outer':None, 
                                  'rho_outer':None,
                                  'indigenous_mass':0.,
                                  'R_outer':None} 
                                for i in range(len(self.contents))]
        
        self.silence = silence

        #if no fractions have been specified, distribute the components evenly
        #over the mixture
        if len(fractions)==0:
            if not self.silence:
                print ('setting up uniform fractions')
            for c in range(len(self.contents)):
                frac  = 1./len(self.contents[c])
                self.fractions.append([frac for i in self.contents[c]])
        
        '''
        if len(self.layermasses) != len(self.contents):
            if not self.silence:
                print ('\nWARNING: number of layers and layermasses does not match')
                print ('given layermasses:', len(self.layermasses))
                print ('given layers: ', len(self.contents))
        '''
        
        self.status='seed'
        
        self.M_surface_is = 0.
        self.T_surface_is = T_center
        self.P_surface_is = P_center
        self.R_surface_is=R_seed
        self.exit_code = 0
        self.direction = None
        self.constraintValue_is = None
        self.constraintValue_should = None
        
        #gather all initial properties in a dictionary
        self.initials = {'P_center': self.P_center,
                              'T_center': self.T_center,
                              'R_seed': self.R_seed,
                              'contents': self.contents,
                              'tempType': self.tempType,
                              'layerType': self.layerType,
                              'fractions': self.fractions,
                              'P_surface_should': self.P_surface_should,
                              'majorConstraint': self.majorConstraint,
                              'layerConstraint': self.layerConstraint,
                              'add_atmosphere': self.add_atmosphere,
                              'eps_r': self.eps_r,
                              'layermasses':self.layermasses,
                              'layerradii':self.layerradii,
                              'differentiated': self.differentiated,
                              'M_surface_should': self.M_surface_should/m_earth,
                              'T_surface_should': self.T_surface_should,
                              'R_surface_should': self.R_surface_should/r_earth,
                              'rhoType':self.rhoType,
                              'echo':self.echo,
                              'brucite_dissoc':self.brucite_dissoc,
                              'Mg_number_should':self.Mg_number_should,
                              'Si_number_should':self.Si_number_should,
                              'Si_number_layers':self.Si_number_layers,
                              'Fe_number_layers':self.Fe_number_layers,
                              'vapor_stop': self.vapor_stop,
                              'Mg_number_layer':self.Mg_number_layer,
                              'modelType':self.modelType,
                              'temp_jumps':self.temp_jumps,
                              'gammas_layer':self.gammas_layer,
                              'inner_core_frac':self.inner_core_frac,
                              'M_ocean_should':self.M_ocean_should,
                              'ocean_frac_should':self.ocean_frac_should,
                              'layertemps':self.layertemps,
                              'T_zero':self.T_zero,
                              'M_core_should':self.M_core_should,
                              'subphase_res':self.subphase_res,
                              'layerpres':self.layerpres,
                              'eps_H2O':self.eps_H2O,
                              'eps_Al':self.eps_Al,
                              'xi_H_core_predicted':self.xi_H_core_predicted,
                              'xi_Stv':self.xi_Stv,
                              'omega':self.omega}
        
        self.finals = {'P_surface_is':self.P_surface_is,
                       'M_surface_is':self.M_surface_is,
                       'T_surface_is':self.T_surface_is,
                       'R_surface_is':self.R_surface_is,
                       'Mg_number_is':self.Mg_number_is,
                       'Si_number_is':self.Si_number_is,
                       'M_H2O_is':self.M_H2O_is,
                       'M_H2O_should':self.M_H2O_should,
                       'M_H2O_hidden_is':self.M_H2O_hidden_is,
                       'density_inversion':self.density_inversion,
                       'match':self.match,
                       'vapor_reached':self.vapor_reached,
                       'layer_properties':self.layer_properties,
                       'N_shells':self.N_shells,
                       'M_ocean_is':self.M_ocean_is,
                       'ocean_frac_is':self.ocean_frac_is,
                       'M_core_is':self.M_core_is,
                       'M_DWR_is':self.M_DWR_is,
                       'xi_H_core':self.xi_H_core,
                       'M_H2O_core':self.M_H2O_core,
                       'ocean_depth':self.ocean_depth}
        
    
    def prt(self, digits=4, **kwargs):
        print ('=======================================')
        print ('Planet properties:')
        print ('\n--------------------------------')
        print ('Layer details:')
        print ('--------------------------------')
        for i in range(len(self.layers)):
            print ('\nlayer: ', i)
            first_shell = self.layers[i].shells[0]
            for j in range(len(self.contents[i])):
                print (round(first_shell.fractions[j]*100, digits), '%', 
                       material_list[self.contents[i][j]],
                       ' ph =', self.layers[i].shells[-1].mix.mix[j].phase)

            print ('layer mass [M_earth]: ', round(self.layers[i].indigenous_mass/\
                   m_earth, digits))
            print ('outer radius [R_earth]:',round(self.layers[i].radius/\
                   r_earth, digits))
            print ('outer P [GPa]: ', round(self.layers[i].pres*1.0e-9, digits))
            print ('outer T [K]: ', round(self.layers[i].temp, digits))

        
        try:
            print ('\n--------------------------------')
            print ('Major parameters:')
            print ('--------------------------------')
            print ('R_surface_is [R_earth]:', round(self.R_surface_is/r_earth, 
                   digits),
               '\nM_surface_is [M_earth]:', round(self.M_surface_is/m_earth, 
                               digits),
               '\nT_surface_is [K]:', round(self.T_surface_is, digits),
               '\nT_center [K]:', round(self.T_center, digits),
               '\nP_surface_is [bar]:', round(self.P_surface_is*1.0e-5, digits),
               '\nP_center [GPa]:', round(self.P_center*1.0e-9, digits),
               '\nMOI factor:', round(self.MOI_is, digits),
               '\n\nMg_number_is:', round(self.Mg_number_is, digits),
               '\nSi_number_should:', round(self.Si_number_should, digits),
               '\nSi_number_is:', round(self.Si_number_is, digits),
               '\nxi_H_core:', round(self.xi_H_core, digits),
               '\nxi_H_core_predicted:', round(self.xi_H_core_predicted, digits),
               '\nP_H2_CMB [MPa]:', round(self.P_H2_CMB*1.0e-6, digits)
               )
            
            try:
                print('M_H2O [wt%]:', 
                      round((self.H2O_count+self.H_count/2.)*mH2O/self.M_surface_is*100, digits))
                print ('M_H2O core [wt%]:', round(self.H_count/2.*mH2O/self.M_surface_is*100, digits))
                print ('M_H2O mantle [wt%]:', round(self.H2O_count*mH2O/self.M_surface_is*100, digits))
                print ('Ocean frac is:', round(10**self.ocean_frac_is, digits))
                print ('Ocean frac should:', round(10**self.ocean_frac_should, digits))
            except ZeroDivisionError:
                print ('M_H2O [wt%]: NaN')
                    
           
            
        except TypeError:
            print ('WARNING: Type Error in Planet.prt()')

        print ('\n--------------------------------')
        print ('Layer overview:')
        print ('--------------------------------')
        
        dat = []
        for i in range(len(self.layer_properties)):
            lay = self.layer_properties[i]
            material_str = ''
            for c in range(len(self.contents[i])):
                material_str += material_list_fort[self.contents[i][c]-1]
                
            dat.append([i, material_str, ftool.scinot(lay['R_outer']/r_earth, 
                                                      digits=digits),
            ftool.scinot(lay['indigenous_mass']/m_earth, digits=digits), 
            ftool.scinot(lay['P_outer']*1.0e-9, digits=digits),
            ftool.scinot(lay['T_outer'], digits=digits), 
            ftool.scinot(lay['rho_outer'], digits=digits)])
        
        tabl=tabulate(dat, headers=['Layer', 'Contents', "R [R_e]", 
                                    'm [M_e]', "P [GPa]",\
                                        'T [K]','rho [kg m-3]'])
           
        print()
        print (f"{tabl}")
        print () 
        

    def Update_finals(self, **kwargs):
        #gather final properties in a dictionary
        self.finals = {'P_surface_is':self.P_surface_is,
               'M_surface_is':self.M_surface_is/m_earth,
               'T_surface_is':self.T_surface_is,
               'R_surface_is':self.R_surface_is/r_earth,
               'Mg_number_is':self.Mg_number_is,
               'Si_number_is':self.Si_number_is,
               'M_H2O_is':self.M_H2O_is,
               'M_H2O_should': self.M_H2O_should,
               'M_H2O_hidden_is':self.M_H2O_hidden_is,
               'density_inversion':self.density_inversion,
               'match':self.match,
               'vapor_reached':self.vapor_reached,
               'layer_properties':self.layer_properties,
               'N_shells':self.N_shells,
               'M_ocean_is':self.M_ocean_is,
               'ocean_frac_is':self.ocean_frac_is,
               'M_core_is':self.M_core_is,
               'M_DWR_is':self.M_DWR_is,
               'MOI_is':self.MOI_is,
               'xi_H_core':self.xi_H_core,
               'M_H2O_core':self.M_H2O_core,
               'ocean_depth':self.ocean_depth}
    
    
    def Update_initials(self, **kwargs):
        """
        """
        if self.M_surface_should == float('inf'):
            print ('HERE infinity in Planet.update_initials()')
        
        #gather all initial properties in a dictionary
        self.initials = {'P_center': self.P_center,
                              'T_center': self.T_center,
                              'R_seed': self.R_seed,
                              'contents': self.contents,
                              'tempType': self.tempType,
                              'layerType': self.layerType,
                              'fractions': self.fractions,
                              'P_surface_should': self.P_surface_should,
                              'majorConstraint': self.majorConstraint,
                              'layerConstraint': self.layerConstraint,
                              'add_atmosphere': self.add_atmosphere,
                              'eps_r': self.eps_r,
                              'layermasses':self.layermasses,
                              'layerradii':self.layerradii,
                              'differentiated': self.differentiated,
                              'M_surface_should': self.M_surface_should/m_earth,
                              'T_surface_should': self.T_surface_should,
                              'R_surface_should': self.R_surface_should/r_earth,
                              'rhoType':self.rhoType,
                              'echo':self.echo,
                              'brucite_dissoc':self.brucite_dissoc,
                              'Mg_number_should':self.Mg_number_should,
                              'Si_number_should':self.Si_number_should,
                              'Si_number_layers':self.Si_number_layers,
                              'Fe_number_layers':self.Fe_number_layers,
                              'vapor_stop': self.vapor_stop,
                              'Mg_number_layer':self.Mg_number_layer,
                              'modelType':self.modelType,
                              'temp_jumps':self.temp_jumps,
                              'gammas_layer':self.gammas_layer,
                              'inner_core_frac':self.inner_core_frac,
                              'M_ocean_should':self.M_ocean_should,
                              'ocean_frac_should':self.ocean_frac_should,
                              'layertemps':self.layertemps,
                              'T_zero':self.T_zero,
                              'M_core_should':self.M_core_should,
                              'subphase_res':self.subphase_res,
                              'layerpres':self.layerpres,
                              'eps_H2O':self.eps_H2O,
                              'eps_Al':self.eps_Al,
                              'xi_H_core_predicted':self.xi_H_core_predicted,
                              'xi_Stv':self.xi_Stv,
                              'omega':self.omega}
        
        
    def Reset(self, **kwargs):
        """Resets all planetary properties to the initial state. This allows
        to use the exact same specifications for re-constructing the planet.
        """
        
        #delete old properties to release memory
        for lay in self.layers:
            del lay.shells
            del lay

        self.__init__(**self.initials)

        if self.initials['M_surface_should'] > 1.0e30:
            print ('HERE large after Planet.reset() in Planet')
            print ('large value=', self.initials['M_surface_should'])
        
        if self.M_surface_should == float('inf'):
            print ('HERE infinity in Planet.reset() after')

            
    
    def Plot(self, scatter=False, layerColoring = False, axis=[], x_axis='radius',
             save = False, filename='planet_structure', suffix='png', path='./', 
             **kwargs):
        
        if axis == []:
            fig, axis = plt.subplots(2, 3, sharex=True)
            fig.subplots_adjust(hspace=0, wspace = .3)
            
            
        title_list = [r'$\rmPressure \ [GPa]$', r'$\rmMass \ [M_\oplus]$', 
                      r'$\rmDensity \ [10^3 \ kg/m^3]$', r'$\rmTemperature  \ [K]$',
                      r'$\rmGravity \ [m/s^2]$', r'$\rmv_{esc} \ [km/s]$']
        
        radius_list_dummy = np.array([[shell.radius/r_earth 
                                       for shell in layer.shells] 
                                        for layer in self.layers])
        
        
        #collect data for all layers
        plot_list_dummy = [[[shell.pres*1.0e-9 for shell in layer.shells] 
                                            for layer in self.layers ],
                    [[shell.mass/m_earth for shell in layer.shells] 
                                            for layer in self.layers ],
                    [[shell.dens*1.0e-3 for shell in layer.shells] 
                                            for layer in self.layers ],
                    [[shell.temp for shell in layer.shells] 
                                            for layer in self.layers ],
                    [[shell.gravity for shell in layer.shells]
                                            for layer in self.layers],
                    [[shell.v_esc*1.0e-3 for shell in layer.shells]
                                            for layer in self.layers]
                    ]
                    
        #collect thermal velocity data for water vapor to overplot it on
        #the graph for the escape velocity for comparison
        
        v_th_H2O_dummy = [[shell.v_th_H2O*1.0e-3 for shell in layer.shells]
                                            for layer in self.layers]
        
        #extract data of individual layers to one array for each parameter
        plot_list = []
        radius_list = []
        v_th_H2O_list = []
        
        for p in range(len(plot_list_dummy)):
            plot_list.append([])
            
        for p in range(len(plot_list_dummy)):
            for l in range(len(self.layers)):
                for i in range(len(plot_list_dummy[p][l])):
                    plot_list[p].append(plot_list_dummy[p][l][i])
            
        for l in range(len(self.layers)):
            for i in range(len(radius_list_dummy[l])):
                    radius_list.append(radius_list_dummy[l][i])
                    v_th_H2O_list.append(v_th_H2O_dummy[l][i])
                    
        ax = [axis[0][0], axis[0][1], axis[1][0], axis[1][1], axis[0][2], 
              axis[1][2]]
        
        ax[2].set_xlabel(r'$\rmRadius  \ [R_\oplus]$')
        ax[3].set_xlabel(r'$\rmRadius \ [R_\oplus]$')
        ax[5].set_xlabel(r'$\rmRadius \ [R_\oplus]$')
        #for axx in ax:
         #   axx.grid(which='both', axis='both')
        
        for i in range(len(ax)):
            if x_axis == 'radius':
                ax[i].plot(radius_list, plot_list[i], color=param_colors[i], zorder=3)

            elif x_axis == 'pres':
                ax[i].semilogx(plot_list[0], plot_list[i], color=param_colors[i], zorder=3)

            ax[i].tick_params(right=True, top=True, direction='in', which='both')
            ax[i].tick_params(which='major', axis='both', length=8)
            ax[i].grid(True, which='both', zorder=0, 
              color=plot_params['gridcolor'],
              alpha=plot_params['gridalpha'])
            
            ax[i].set_facecolor(plot_params['backgroundcol'])
            ax[i].grid(which='major', axis='both', linewidth=2, zorder=0)
            ax[i].xaxis.set_minor_locator(AutoMinorLocator())
            ax[i].yaxis.set_minor_locator(AutoMinorLocator())
            '''
            for j in range(len(radius_list)):
                if scatter:
                    ax[i].scatter(radius_list[j], plot_list[i][j])
                else:
                    ax[i].plot(radius_list[j], plot_list[i][j])
             '''     
            ax[i].set_ylabel(title_list[i])
        
        ax[-1].plot(radius_list, v_th_H2O_list, color=param_colors[-1],
          linestyle = '--')
        
        #save figure as image file if turned on
        if save:
            fig.savefig(path+filename+'.'+suffix)
            
        fig2, ax2 = plt.subplots()
        ax2.loglog(np.asarray(plot_list[3]), np.asarray(plot_list[0])*1.0e9)
            
            
    def write(self, out='planet', loc='./'):
        """Creates data table of the planetary structure and writes it to ascii 
        file. By default, the data is organized as follows:
        
        R           M           T       P       rho         v_esc       g
        (r_earth)   (m_earth)   (K)     (GPa)   (kg/m3)     (km/s)      (m/s2)
        -----------------------------------------------------------------------        
        R_center    M_center    ...
        ...         ...         ...
        R_tot       M_tot       ...
   
        The meta data contains the properties that are stored in Planet.initials
        and Planet.finas by default
        """
        self.trim_profiles()
        
        #gather all planetary parameters as meta data
        meta = {'initials':self.initials, 'finals':self.finals}
        
        names = ['R (r_earth)',
                 'T (K)',
                 'P (GPa)',
                 'rho (kg/m3)',
                 'M (m_earth)',
                 'g (m/s2)']
        
        data_table=astropy.table.Table(self.profiles.T, meta=meta)
        ascii.write(data_table, loc+out+suffix['planet'], overwrite=True,
                   format='ecsv', names = names)
        
    
    def correct_density_profiles(self):
        i = 0
        
        df = self.data_table.to_pandas()
        
        N_shells = len(df['rho (kg/m3)'])
        
        while True:
            if i == N_shells-2:
                break
            
            else:
                d = df['rho (kg/m3)'][i]
                dd = df['rho (kg/m3)'][i+1]
                
                if d < dd:
                    self.data_table[i][3] = dd
                
            i += 1
                
    
    def load(self, loc='./', file_name='planet'):
        """
        """
        #read ascii file for eos table
        self.data_table=ascii.read(loc+file_name+suffix['planet'], format='ecsv')
         
        self.initials = self.data_table.meta['initials']           
        self.finals = self.data_table.meta['finals']
   
