#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 09:43:30 2021

@author: os18o068

"""
import PlanetInvestigator
from PlanetInvestigator import params as param_labels
import numpy as np
import os
import MOI_model as analytical_models
from matplotlib import pyplot as plt
import functionTools as ftool
import time
import pandas as pd
from matplotlib.patches import Rectangle
import sys
import hist3d
import copy
import matplotlib as mpl
import  matplotlib.transforms as transforms
from PIMPphysicalparams import Mg_earth, Si_earth, SFe_earth, SFe_solar, m_earth, \
    r_earth, R_solar, MoI_solar, M_solar
from PIMPrunparams import layerColors, parameter_mask
import phase_transitions_water_Wagner2002 as phase
import astropy.io
import read_exoplanets
import CircSegPlot

#This mothafucka speeds up saving graphics with mpl like crazy
mpl.use('Agg')

class DataStream():
    def __init__(self, directory = None):
        self.data_root = '/home/os18o068/Documents/PHD/Projects/Planets/Data'
        self.data = None
        self.sample_data = None
        self.filtered_sample_data = {}
        self.filtered_sample_planets = {}
        self.filtered_sample_ranges = {}
        self.filtered_sample_slices = {}
        self.sorted_sample_indices = {}
        self.sorted_sample_ranges = {}
        self.sorted_sample_layer_ranges = {}
        self.analytical_spirits = {}

        self.filtered_data = {}
        self.filtered_planets = {}
        self.filtered_ranges = {}
        self.filtered_slices = {}
        self.sorted_indices = {}
        self.sorted_ranges = {}
        
        if directory != None:
            self.directory = directory
            self.target_directory = self.data_root + '/' + self.directory
    
    
    def get_bulk_structures(self, extract = True):
        #Extract bulk structures for object data
        self.bulk_structures = []
        try:
            if extract:
                for pl in self.planets:
                    for pll in pl:
                        
                        pll.extract_bulk_structure()
               
            for w in self.data:
                self.bulk_structures.append(list(set(w['hyd_struc'].values)))
            self.bulk_structures = list(set(sum(self.bulk_structures, [])))
            
        except TypeError:
            pass
        
        except AttributeError:
            pass
        
        #Extract bulk structures for sample data
        try:
            self.sample_bulk_structures = []
            #Loop over all planets in the sample
            if extract:
                for pl in self.sample_planets:
                    pl.extract_bulk_structure()
            
            self.sample_bulk_structures =\
                [int(x) for x in list(set(self.sample_data['hyd_struc'].values))]
            #self.sample_hydro_structures = list(set(sum(self.sample_hydrostructures, [])))
        
        except TypeError:
            pass
        except AttributeError:
            pass
        
    
    def get_hydro_structures(self, extract = True):
        
        #Extract hydro structures for object data
        self.hydro_structures = []
        try:
            if extract:
                for pl in self.planets:
                    for pll in pl:
                        pll.extract_hydro_structure()
               
        
            for w in self.data:
                self.hydro_structures.append(list(set(w['hyd_struc'].values)))
            self.hydro_structures = list(set(sum(self.hydro_structures, [])))
        except TypeError:
            pass
        except AttributeError:
            pass
        
        #Extract hydro structures for sample data
        try:
            self.sample_hydro_structures = []
            #Loop over all planets in the sample
            if extract:
                for pl in self.sample_planets:
                    pl.extract_hydro_structure()
            
            self.sample_hydro_structures =\
                [int(x) for x in list(set(self.sample_data['hyd_struc'].values))]
            #self.sample_hydro_structures = list(set(sum(self.sample_hydrostructures, [])))
        
        except TypeError:
            pass
        except AttributeError:
            pass
        
        
    def get_object_names(self, max_mass_error = .15):
        #Create Archive instance
        self.archive = read_exoplanets.Archive()
        #Filter the archive data
        self.archive.filter()
        
        #Extract allowed range for equilibrium temperature
        self.temp_range = [[self.archive.df['pl_eqt'][ind] + \
                            self.archive.df['pl_eqterr2'][ind],
                      self.archive.df['pl_eqt'][ind] + \
                            self.archive.df['pl_eqterr1'][ind]]
                      for ind in self.archive.df.index]
        
        #Extract allowed range for total mass
        self.mass_range = [[self.archive.df['pl_bmasse'][ind] + \
                            self.archive.df['pl_bmasseerr2'][ind],
                      self.archive.df['pl_bmasse'][ind] + \
                            self.archive.df['pl_bmasseerr1'][ind]]
                      for ind in self.archive.df.index]
            
        #Extract allowed range for total radius
        self.rad_range = [[self.archive.df['pl_rade'][ind] + \
                            self.archive.df['pl_radeerr2'][ind],
                      self.archive.df['pl_rade'][ind] + \
                            self.archive.df['pl_radeerr1'][ind]]
                      for ind in self.archive.df.index]
        
            
    def get_cleaned_data(self, directory = None):
        """Read in existing pre-processed data. Data is PlanetFort.Planet()
        instances which have already been filtered. Only models which converged
        and for which the mass and radius match the object are contained.
        """
        if directory != None:
            self.directory = directory
            self.target_directory = self.data_root + '/' + self.directory        
        
        #Loads pre-processed data for each object into a list
        try:
            targets = [self.target_directory + '/' + obj + '/' + obj + '_cleaned.pkl'
                   for obj in self.archive.df['pl_cleanname']]
            
        except AttributeError:
            self.get_object_names()
            targets = [self.target_directory + '/' + obj + '/' + obj + '_cleaned.pkl'
                   for obj in self.archive.df['pl_cleanname']]
        
        print ('targets =', targets)
        self.planets = [ftool.load_objects(t) for t in targets]        
        
        
    def get_raw_data(self, directory = None):
        """Read in existing raw data. Raw data is PlanetFort.Planet() instances.
        """
        if directory != None:
            self.directory = directory
            self.target_directory = self.data_root + '/' + self.directory
        
        #Loads all data for each object into a list
        try:
            targets = [self.target_directory + '/' + obj + '_cali.pkl'
                   for obj in self.archive.df['pl_cleanname']]
            
        except AttributeError:
            self.get_object_names()
            targets = [self.target_directory + '/' + obj + '_cali.pkl'
                   for obj in self.archive.df['pl_cleanname']]
        
        print ('targets =', targets)
        self.all_planets = [ftool.load_objects(t) for t in targets]
        
    
    def process_raw_data(self, write = False):
        """Process the raw data. Raw data is PlanetFort.Planet() instances.
        """
        self.planets = []
        not_converged_count = 0
        
        if write:
            #Create list of targets
            try:
                targets = [self.target_directory + '/' + obj + '/' + obj + '_cleaned.pkl'
                       for obj in self.archive.df['pl_cleanname']]
                
            except AttributeError:
                self.get_object_names()
                targets = [self.target_directory + '/' + obj  + '/' + obj + '_cleaned.pkl'
                       for obj in self.archive.df['pl_cleanname']]
        
        #Loop over all data sets
        for i in range(len(self.all_planets)):
            self.planets.append([])
            #Loop over all planets in the set
            for j in range(len(self.all_planets[i])):
                self.all_planets[i][j].check_convergence()
                #Check if planet was actually constructed
                if self.all_planets[i][j].status == 'constructed':
                    #Check if planet converged to boundary conditions
                    if not self.all_planets[i][j].converged:
                        not_converged_count += 1
                        
                    else:
                        R = self.all_planets[i][j].R_surface_is / r_earth
                        
                        #Check if radius is within allowed range for object
                        if R <= self.rad_range[i][1] and R >= self.rad_range[i][0]:
                            self.planets[i].append(self.all_planets[i][j])
            
            #Store processed raw data for each individual planet            
            if write:
                print ('target =', targets[i])
                ftool.save_objects(self.planets[i], targets[i])
    
    
    def create_working_data(self, instream = False, write = True, 
                            directory = None, stats = False, load = False):
        """Creates working data from structure models.
        """
        #Use pre-existing working data. This option can be used if no changes
        #to the data were made but the output format of e.g. the ranges.txt
        #file changed. Otherwise the raw data has to be read in and processed.
        print ("creating working data...", write)
        if load:
            self.get_raw_data(directory = directory)
            self.load_working_data()
            self.process_raw_data(write = write)
            
        else:
            #Use data of current data stream
            if instream:
                pass
            
            #Load data from the objects directories
            else:
                self.get_raw_data(directory = directory)
                self.process_raw_data(write = write)
            
            #Extract relevant parameters for each planet
            self.data = np.array([PlanetInvestigator.extract_params(pl) 
                                  for pl in self.planets])
            self.convert_to_pandas()
        
        if write:
            print ("writing planets to file.")
            #Create list of targets
            try:
                targets = [self.target_directory + '/' + obj + '/' + obj + '_cleaned'
                       for obj in self.archive.df['pl_cleanname']]
                
            except AttributeError:
                self.get_object_names()
                targets = [self.target_directory + '/' + obj  + '/' + obj + '_cleaned'
                       for obj in self.archive.df['pl_cleanname']]
            
            #Save data to numpy array
            for i in range(len(self.data)):
                print ('exporting ', targets[i] + '.npy')
                np.save(targets[i]+'.npy', self.data[i].to_numpy())
               
        if stats:
            for i in range(len(self.archive.df['pl_cleanname'])):
                obj = self.archive.df['pl_cleanname'].values[i]
                
                try:
                    PlanetInvestigator.check_statistics(self.data[i].to_numpy(),
                                 obj = obj,
                                 spec = self.target_directory + '/' + obj,
                                 write = True,
                                 N_all = len(self.all_planets[i]),
                                 N = 1)
                except IndexError:
                    pass
                
                except ValueError:
                    pass
                
                
    def load_working_data(self):
        """Loads existing data for analysis.
        """
        #Create list of targets
        print ('Looking for data in:', self.target_directory)
        try:
            targets = [self.target_directory + '/' + obj + '/' + obj + '_cleaned'
                   for obj in self.archive.df['pl_cleanname']]
            
        except AttributeError:
            self.get_object_names()
            targets = [self.target_directory + '/' + obj  + '/' + obj + '_cleaned'
                   for obj in self.archive.df['pl_cleanname']]
        
        self.data = []
        
        #Save data to numpy array
        for i in range(len(targets)): 
            dat = np.load(targets[i]+'.npy')        
            self.data.append(dat)
            
        self.data = np.array(self.data)
        self.convert_to_pandas()
     
        
    def convert_to_pandas(self):
       #Convert object data to pandas frame
       try:
           for d in range(len(self.data)):       
                self.data[d] = pd.DataFrame(self.data[d],
                                            columns = parameter_mask)
       except TypeError:
           pass
        
       #Convert sample data to pandas frame
       try: 
           self.sample_data = pd.DataFrame(self.sample_data,
                                                  columns = parameter_mask)
       except TypeError:
           pass

            
    def full_sample_process(self, file):
        #self.create_sample(file)
        #self.filter_sample(file, tols = tols)
        #self.get_bulk_structures()
        self.get_sorted_sample_indices()
        self.get_sorted_sample_ranges(file)
        self.create_sample_slices(file)
    
    
    def create_analytical_models(self, obj, N = 10000):
        """ Takes the sorted sample ranges for the layer densities and constructs
        analytical N-layer models within the same density ranges for each layer.
        The resulting parameter space that matches the boundary conditions of the
        objects can then be compared for the two approaches.
        """
        strucs = list(self.sorted_sample_layer_ranges[obj].keys())
        inner = [max(ssr['rho_inner']) for ssr in self.sorted_sample_layer_ranges[obj][strucs[0]]]
        outer = [min(ssr['rho_outer']) for ssr in self.sorted_sample_layer_ranges[obj][strucs[0]]]
        
        dens_ranges = [[a, b] for a, b in zip(outer, inner)]
        dens_ranges[-1] = [1300, 1600]
#        dens_ranges.append([880, 960])

        # Initiate analytical spirit for the specified object
        ams = analytical_models.analyticalModelSpirit(n_models = N,
                                                      total_mass = M_solar[obj.split('_')[0]],
                                                      n_layers = 3,
                                                      dens_ranges = dens_ranges)
        ams.get_layer_masses()
        ams.create_models()
        ams.construct_models()
        ams.filter(obj.split('_')[0], tolerance = [.01, .001, .01])
        self.analytical_spirits.update({obj:ams})
        self.analytical_spirits[obj].get_ranges()
        
        
    def plot_comparison(self, obj):
        #Load profils from literature
        profile_path = "/home/os18o068/Documents/PHD/Projects/Planets/Data/icy_satellites/revived_project/stupid_profiles_from_plots/"
        files_cammarano = ['black_solid', 'black_dashed', 'red_solid', 'red_dashed']
        
        dfs_cammarano = [pd.read_csv('{a}cammarano_2006_{b}.csv'.format(a=profile_path,
                                                                       b= fc),
                                     names = ['density', 'depth'])
                                     for fc in files_cammarano]
                                     
        fig, ax = plt.subplots(1,3, figsize = (9, 4))
        plt.subplots_adjust(wspace = .25)
        
        rads = dfs_cammarano[1]['depth'].max() - dfs_cammarano[1]['depth']
        ax[1].plot(rads * 1e3 / r_earth,
                   dfs_cammarano[1]['density'] * 1e3)
        
        self.analytical_spirits[obj].plot_profiles(filtered=True,
                                                   ax = ax[0])
        
        self.plot_sample_profiles(ax = ax, obj = obj)
        
        obj_name = obj.split('_')[0]
        file_name = '{}_comparison'.format(obj)
        file_path = '{a}/{b}/{c}'.format(a=self.target_directory,
                                         b = obj_name,
                                         c = file_name)
        fig.savefig('{}.pdf'.format(file_path), format = 'pdf', bbox_inches = 'tight')
        fig.savefig('{}.svg'.format(file_path), format = 'svg', bbox_inches = 'tight')
        fig.savefig('{}.png'.format(file_path), format = 'png', bbox_inches = 'tight')
        plt.close(fig)
        
        
    def create_sample(self, file, write = True, stats = True,
                      check_convergence = True,
                      check_params = ["M_surface", "T_surface",
                                      "ocean_frac"],
                      accs = {"M_surface":1e-4,
                              "T_surface":5e-3,
                              "ocean_frac":1e-3}):
        
        print ("loading file:", file + '_cali.pkl')        
        all_planets = ftool.load_objects(self.target_directory + '/' + file + '_cali.pkl')
        self.sample_planets = []
        self.sample_name = file
        self.filtered_sample_planets = {}
        
        #Check conditions for boundary conditions and discard planets for which
        #convergence was not achieved during the structure integration.
        for pl in all_planets:
            pl.check_convergence(accs = accs, check_params = check_params)
            pl.trim_profiles()
            #print ("converged =", pl.converged)
            if check_convergence:
                if pl.converged:
                    self.sample_planets.append(pl)
            else:
                self.sample_planets.append(pl)
                
        conv = len(self.sample_planets) / len(all_planets)
        print ("Number of planets {}:".format(len(self.sample_planets)))
        print ('Convergence rate: {}%'.format(int(100 * conv)))
        
        #Extract parameters from raw data
        data = np.array([PlanetInvestigator.extract_params(self.sample_planets)])[0]
        self.sample_data = data

        if stats:
            PlanetInvestigator.check_statistics(self.sample_data,
                         obj = file + '_cleaned',
                         spec = self.target_directory + '/' + file,
                         write = True,
                         N_all = len(self.sample_planets),
                         N = 1)

        #Convert to pandas frame for further use
        self.convert_to_pandas()
        
        if write:
            try:
                ftool.save_objects(self.sample_planets,
                                   self.target_directory + '/' + file + '/' + file + '_cleaned.pkl')
                np.save(self.target_directory + '/' + file + '/' + file+'_cleaned.npy', 
                        data)    
                self.sample_data.to_csv(self.target_directory + '/' + self.sample_name + '/' + self.sample_name + '.csv')
            except FileNotFoundError:
                os.mkdir(self.target_directory + '/' + file)
                ftool.save_objects(self.sample_planets,
                                   self.target_directory + '/' + file + '/' + file +'_cleaned.pkl')
                np.save(self.target_directory + '/' + file + '/' + file+'_cleaned.npy', 
                        data)
                self.sample_data.to_csv(self.target_directory + '/' + self.sample_name + '/' + self.sample_name + '.csv')    
        
        
    
    def get_sorted_sample_indices(self):
        """Creates list of indices for all possible bulk structures of each
        individual subsample. This list can then be used to specifically target
        models of a specific bulk strucuture from the data set. This method is
        currently applied to the filtered sample data only.
        """
        #Loop over subsamples
        for obj in self.filtered_sample_planets:
            self.sorted_sample_indices.update({obj:{}})
            #Loop over individual models
            for i in range(len(self.filtered_sample_planets[obj])):
                pl = self.filtered_sample_planets[obj][i]
                struc = ''.join([str(s) for s in pl.bulk_structure])
                #If new structure encountered add it to the dict
                if not struc in list(self.sorted_sample_indices[obj].keys()):
                    self.sorted_sample_indices[obj].update({struc:[i]})
                #Add index to existing structure key
                else:
                    self.sorted_sample_indices[obj][struc].append(i)


    def get_sorted_indices(self):
        """Creates list of indices for all possible bulk structures of each
        individual planet. This list can then be used to specifically target
        models of a specific bulk strucuture from the data set. This method is
        currently applied to the filtered data only.
        """
        self.sorted_indices = {}
        #Loop over planets
        for j in range(len(self.planets)):
            obj = self.archive.df['pl_cleanname'].to_list()[j]
            self.sorted_indices.update({obj:{}})
            #Loop over individual models
            for i in range(len(self.planets[j])):
                pl = self.planets[j][i]
                struc = ''.join([str(s) for s in pl.bulk_structure])
                #If new structure encountered add it to the dict
                if not struc in list(self.sorted_indices[obj].keys()):
                    self.sorted_indices[obj].update({struc:[i]})
                #Add index to existing structure key
                else:
                    self.sorted_indices[obj][struc].append(i)
            

    def filter_sample(self, obj, tols = [1e-3, 1e-2], write = False, stats = True):
        """Filter out all models which do not match with the specified objects
        boundary conditions. Currently the object must be a major solar system
        body as specified in PIMPphysicalparams.
        """
        
        tols_string = '_{a:.2g}_{b:.2g}'.format(a=tols[0] * 100,
                                             b=tols[1] * 100)
        tols_string = tols_string.replace('.','_')
        
        R_should = R_solar[obj]
        MoI_should = MoI_solar[obj]
        self.filtered_sample_planets.update({obj + tols_string:[]})
        self.filtered_sample_data.update({obj + tols_string:None})
        data = []
        for i in range(len(self.sample_planets)):
            pl = self.sample_planets[i]
            match = True
            reldevs = [abs(R_should - pl.R_surface_is / r_earth) / R_should,
                       abs(MoI_should - pl.MOI_is) / MoI_should]
            print ("reldev =", reldevs)
            for rd, tl in zip(reldevs, tols):
                
                if rd > tl:
                    match = False
               
            if match:
                self.filtered_sample_planets[obj + tols_string].append(pl)
                data.append(self.sample_data.values[i])

        name = self.sample_name + '_filtered_' + obj + tols_string        
        if write:

            ftool.save_objects(self.filtered_sample_planets[obj + tols_string], 
                               self.target_directory + '/' + self.sample_name + '/' + name + '.pkl')
            np.save(self.target_directory + '/' + self.sample_name + '/' + name + '.npy', 
                    np.array(data))
        
        if stats:
            ranges = PlanetInvestigator.check_statistics(np.array(data),
                         obj = name,
                         spec = self.target_directory + '/' + self.sample_name,
                         write = True,
                         N_all = len(self.filtered_sample_planets[obj + tols_string]),
                         N = 1)
        
            self.filtered_sample_ranges[obj + tols_string] = ranges
        data = pd.DataFrame(data, columns = parameter_mask)
        if write:
            data.to_csv(self.target_directory + '/' + self.sample_name + '/' + name + '.csv')
        
        self.filtered_sample_data[obj + tols_string] = data
        frac = len(self.filtered_sample_planets[obj + tols_string]) / len(self.sample_planets)
        print ("Number of matches:", len(self.filtered_sample_planets[obj + tols_string]))
        print ("Fraction of matches: {}%".format(round(frac * 100, 1)))


    def write_sample(self):
        """Write sample data to file.
        """
        data = self.sample_data.to_numpy()
        self.sample_dir = self.data_root + '/' + self.directory + '/' + self.sample_name
        ranges, std_low, std_high, rmsd = PlanetInvestigator.check_statistics(data,
                                            N = 1, 
                                            red = 1,
                                            obj = self.sample_name,
                                            spec = self.sample_dir, 
                                            write = True,
                                            N_all = len(data))


    def load_filtered_sample(self,  file, obj):
        """Loads filtered sample data unspecific to target objects if present.
        Filters that can be applied are currently restricted to marjo solar
        system objects.
        """
        self.sample_name = file
        
        if not obj in list(self.filtered_sample_data.keys()):
            file_name = '{a}_filtered_{b}'.format(a = file,
                                                  b = obj)
                
            try:
                filtered_planets = ftool.load_objects(self.target_directory + '/' + file + '/' + file_name + '.pkl')
                filtered_data = np.load(self.target_directory + '/' + file + '/' + file_name +'.npy')
                    
            except FileNotFoundError:
                print ("WARNING: Could not find specified file name {}! Doing nothing.".format(file_name))
            
            
            filtered_data = pd.DataFrame(filtered_data, columns = param_labels)
            filtered_ranges = pd.read_pickle(self.target_directory + '/' + file + '/' +file_name + '_ranges.pkl')
            self.filtered_sample_planets.update({obj:filtered_planets})
            self.filtered_sample_data.update({obj:filtered_data})
            self.filtered_sample_ranges.update({obj:filtered_ranges})
            
        else:
            print ('This filter has already been applied to the data set!')
        
        
    def load_sample(self, file, **kwargs):
        """Loads cleaned sample data unspecific to target objects if present.
        """
        self.sample_name = file
        self.filtered_sample_planets = {}
        
        file_name = '{}_cleaned'.format(file)
        loc = self.target_directory + '/' + file + '/' + file_name
        print ("loc = ", loc)
        try:
            self.sample_planets = ftool.load_objects(loc + '.pkl')
            self.sample_data = np.load(loc +'.npy')
            self.sample_ranges = pd.read_pickle(loc + '_ranges.pkl')
        
        except FileNotFoundError:
            print ("WARNING: Could not find specified file name {}! Doing nothing.".format(file_name))
            
        self.convert_to_pandas()
            
    
    def plot_sample(self, xlims = [1., 10.], ylims = [0.5, 3.5], cmap = None,
                    fnts = 14, **kwargs):
        fig, ax = plt.subplots()
        
        if cmap == None:
            cols = np.array([[0.25, 0.3, 0.5], [0.5, 0.6, 1.], [0.8, 0.9, 1.]])
            cmap = ftool.my_dev_cmap(N_col = 64, border_colors = cols)
        
        ax.scatter(self.sample_data['mass'], self.sample_data['radius'],
                   c = self.sample_data['water_mass_frac'], cmap = cmap,
                   zorder = 0)
        sm = plt.cm.ScalarMappable(cmap = cmap)
        cbar = plt.colorbar(sm, ticks = np.linspace(0, 1, 6))
        cbar.set_label(r'$\rm Water \ mass \ fraction$', fontsize = fnts)
        cbar.ax.tick_params(labelsize = fnts)
        #cbar.ax.set_yticklabels(ticklabels.astype(int))
        
        try:
            self.archive.plot(ax = ax, fnts1 = fnts, **kwargs)
        except AttributeError:
            pass
        
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        ax.tick_params(right=True, top=True)
        
        fig.savefig(self.target_directory + '/' + self.sample_name + '/arch.pdf',
                    format = 'pdf', dpi = 480, bbox_inches = 'tight')
    
    
    def get_ranges(self):
        print ("Note that the old version of get_ranges is now called \
               get_ranges_stupid!")
        
        
    def get_ranges_stupid(self, colors = {'supercritical':(.5,0.25,1),
                                             'solid':(.75,.8,1),
                                             'liquid':(.25,.5,1),
                                             'core':(.75,.5,.25),
                                             'mantle':(.4,.25,0)}):
        """
        Plots interior structure vizualization with the ranges for each 
        layer interface for each planet in the data stream.

        Returns
        -------
        None.

        """
        self.all_planet_structures = []
        self.all_planet_colors = []
        self.all_planet_ranges = []
        self.all_planet_locs = []
        self.all_planet_counts = [{} for i in range(len(self.planets))]
        #Loop over all planets
        for p in range(len(self.planets)):
            #Loop over all models for the planet
            all_locs = {}
            all_temps = {}
            all_pres = {}
            all_masses = {}
            all_strucs = []
            all_colors = []
            all_ranges = {}
            #Loop over all models for the current planet
            print ('num. of models for planet {} ='.format(p), len(self.planets[p]))
            for j in range(len(self.planets[p])):
                #Gather locations of phase changes in the hydrosphere
                struc = str(int(self.data[p]['hyd_struc'].values[j]))
                if struc in all_strucs:
                    #Increase counter by one for current structure
                    self.all_planet_counts[p][struc] += 1
                
                else:
                    #Add first counter to for the current structure
                    self.all_planet_counts[p].update({struc:1})
                    all_strucs.append(struc)
                    all_locs.update({struc:[]})
                    all_masses.update({struc:[]})
                    all_pres.update({struc:[]})
                    all_temps.update({struc:[]})
                    all_colors.append([colors['core'], colors['mantle']])
                    for k in range(len(self.planets[p][j].locs) - 1):
                        all_colors[-1].append(colors[self.planets[p][j].phases[k]])

                all_locs[struc].append([0.,
                                        self.data[p]['R_core (km)'].values[j] * 1e3])
                
                all_pres[struc].append([self.data[p]['P_C'].values[j] * 1e9,
                                        self.data[p]['P_CMB'].values[j]])
                            
                all_temps[struc].append([self.data[p]['T_C'].values[j],
                                        self.data[p]['T_CMB'].values[j]])

                all_masses[struc].append([self.data[p]['core_mass_frac'].values[j], 
                                          self.planets[p][j].dms[0] / self.data[p]['M'].values[j] /m_earth- self.data[p]['core_mass_frac'].values[j]])
                
                #Append values to corresponding strcuture
                for k in range(len(self.planets[p][j].locs)):
                    all_locs[struc][-1].append(self.planets[p][j].locs[k])
                    all_pres[struc][-1].append(self.planets[p][j].pres[k])
                    all_temps[struc][-1].append(self.planets[p][j].temps[k])
                    
                    try:
                        #The dms are the enclosed mass at each phase transition
                        mass = self.planets[p][j].dms[k+1]
                        mass -= self.planets[p][j].dms[k]# * (1. - 10**self.planets[p][j].ocean_frac_is)
                        mass /= self.planets[p][j].M_surface_is
                        all_masses[struc][-1].append(mass)
                    except IndexError:
                        pass

            all_ranges ={'layerSizes': [],
                    'layerPressures': [],
                    'layerTemperatures': [],
                    'layerMasses': []}

            #Loop over different structures
            for i in range(len(all_strucs)):
                
                all_ranges['layerSizes'].append([])
                all_ranges['layerPressures'].append([])
                all_ranges['layerTemperatures'].append([])
                all_ranges['layerMasses'].append([])
                                
                locs = np.array(all_locs[all_strucs[i]])
                temps = np.array(all_temps[all_strucs[i]])
                pres = np.array(all_pres[all_strucs[i]])
                masses = np.array(all_masses[all_strucs[i]])
                                
                #Loop over layers of each structure
                for j in range(len(locs.T)):
                    #Append boundary ranges for each layer
                    all_ranges['layerSizes'][i].append([min(locs.T[j]), 
                                                     max(locs.T[j])])
                    all_ranges['layerPressures'][i].append([min(pres.T[j]), 
                                                         max(pres.T[j])])
                    all_ranges['layerTemperatures'][i].append([min(temps.T[j]), 
                                                            max(temps.T[j])])
                    try:
                        all_ranges['layerMasses'][i].append([min(masses.T[j]), 
                                                            max(masses.T[j])])
                    except IndexError:
                        pass     

            self.all_planet_structures.append(all_strucs)
            self.all_planet_ranges.append(all_ranges)
            self.all_planet_colors.append(all_colors)
            self.all_planet_locs.append(all_locs)


    def get_structure_keys(self, obj):
        """Extract the possible interior structures for the sample data.
        Possible structure keys are:
            
            inner core
            outer core
            inner mantle
            outer mantle
            liquid (water)
            supercritical (water)
            solid (water)
            gas (water)
            vapor (water)
        """
        pass


    def get_layer_ranges(self):
        """
        Compute max and min values for layer parameters considering all models
        for the individual planets and for each structure individually.

        Parameters
        ----------
        obj : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.filtered_data = {}
        self.filtered_planets = {}
        self.filtered_ranges = {}
        self.sorted_ranges = {}
        
        #Loop over planets
        for j in range(len(self.planets)):
            obj = self.archive.df['pl_cleanname'].to_list()[j]
            
            if not obj in self.sorted_ranges.keys():
                self.sorted_ranges.update({obj:{}})
            else:
                self.sorted_ranges[obj] = {}
                
            #Loop over structures
            for struc in list(self.sorted_indices[obj].keys()):
                ind = self.sorted_indices[obj][struc]
                #gather all models for current structure
                planets = [self.planets[j][i] for i in ind]
                
                if not struc in self.sorted_ranges[obj].keys():
                    self.sorted_ranges[obj].update({struc:[{} for i in range(len(struc))]})
                
                #Loop over individual layers of current structure
                for i in range(len(struc)):
                    #Loop over layer properties
                    for param in planets[0].cleaned_layer_properties[i].keys():
                        vals = [pl.cleaned_layer_properties[i][param] for pl in planets]
                        
                        try:
                            self.sorted_ranges[obj][struc][i].update({param:[min(vals), max(vals)]})
                        except TypeError:
                            pass
                    
                self.sorted_ranges[obj]
    
    
    def get_sorted_sample_ranges(self, obj):
        """
        Compute max and min values for layer and main parameters considering all models
        in the subsample and for each structure individually.

        Parameters
        ----------
        obj : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if not obj in self.sorted_sample_ranges.keys():
            self.sorted_sample_ranges.update({obj:{}})
        else:
            self.sorted_sample_ranges[obj] = {}
            
        if not obj in self.sorted_sample_layer_ranges.keys():
            self.sorted_sample_layer_ranges.update({obj:{}})
        else:
            self.sorted_sample_layer_ranges[obj] = {}

        #Loop over structures
        for struc in list(self.sorted_sample_indices[obj].keys()):
            ind = self.sorted_sample_indices[obj][struc]
            #gather all models for current structure
            planets = [self.filtered_sample_planets[obj][i] for i in ind]
            #gather all data for current structure
            data = [self.filtered_sample_data[obj].values[i] for i in ind]
            df = pd.DataFrame(data, columns=self.filtered_sample_data[obj].columns)
        
            if not struc in self.sorted_sample_layer_ranges[obj].keys():
                self.sorted_sample_layer_ranges[obj].update({struc:[{} for i in range(len(struc))]})
            
            if not struc in self.sorted_sample_ranges[obj].keys():
                self.sorted_sample_ranges[obj].update({struc:{}})
            
            #Extract layer ranges
            #Loop over individual layers of current structure
            for i in range(len(struc)):
                #Loop over layer properties
                for param in planets[0].cleaned_layer_properties[i].keys():
                    vals = [pl.cleaned_layer_properties[i][param] for pl in planets]
                    
                    try:
                        self.sorted_sample_layer_ranges[obj][struc][i].update({param:[min(vals), max(vals)]})
                    except TypeError:
                        pass
                    
            #Extract main ranges
            #Loop over main parameters
            for param in self.filtered_sample_data[obj].columns:
                vals = [df[param].min(), df[param].max()]
                self.sorted_sample_ranges[obj][struc].update({param:vals})
                
    
    def create_slices(self, label_font_size = 8, split_mantle = True):
        """Creates CircSegPlot planet_chart instances for each planet and each
        individual structure based on the layer ranges.
        """
        self.filtered_slices = {}
        middleKeys = {0:["mass_fraction"],
                      1:["mass_fraction"],
                      2:["mass_fraction"],
                      3:[],
                      4:[],
                      5:[],
                      6:[]}
        
        for obj in list(self.sorted_ranges.keys()):
            
            ranges = self.sorted_ranges[obj]
            structures = list(ranges.keys())
            structures_dummy = copy.deepcopy(structures)
            combined_structures = {}
            # Check for different mantle structures
            # Loop over structures
            for i in range(len(structures_dummy)):
                structures_dummy[i] = structures_dummy[i].replace('56', '6')
                structures_dummy[i] = structures_dummy[i].replace('5', '6') 
                if not structures_dummy[i] in combined_structures.keys():
                    combined_structures.update({structures_dummy[i]:[]})
                    
                combined_structures[structures_dummy[i]].append(structures[i])
            
            # Collect combined data
            data = {}
            structure_ranges = {}
            for cs in structures_dummy:
                sr = [ranges[key] for key in combined_structures[cs]]
                structure_ranges.update({cs:sr})
            
            #Loop over structures
            for cs in list(combined_structures.keys()):
                dat = structure_ranges[cs]
                print ('----------')
                
                combo = [{}] * len(cs)
                
                #Loop over different mantle structures for structure
                for i in range(len(dat)):
                    
                    param_keys = list(dat[i][0].keys())
                    
                    for param in param_keys:
                        if not param in combo[i]:
                            combo[i].update({param:[]})
                        #Add minimum and maximum
                        #combo[i][param].append(dat[i][0][param])
                #{key:[min(d[key][0]), max(d[key][1])]}
                print ('combo =', combo)
            print ('structures =', structures)
            print ('structures dummy =', structures_dummy)
            print ('\ncombined strucs =', combined_structures)
            #sprint ('\nstructure ranges =', structure_ranges)
            #print ('\ndata =', data)
            if len(structures) > 2: 
                sys.exit()
                
            fig, ax = plt.subplots(1, len(list(ranges.keys())), figsize = (6,6))
            plt.subplots_adjust(wspace = 1.)
            
            for i in range(len(structures)):
                #struc = list(ranges.keys())[i]
                struc = structures[i]
                data = {struc:ranges[struc]}
                print ('data =', data)
                pc = CircSegPlot.PlanetChart(data,
                                             on_scale = False,
                                             phi1 = np.pi / 2.5)
                #pc.get_ghosts()
                pc.get_slices()
                self.filtered_slices.update({obj:pc})
                try:
                    pc.draw(ax=ax[i])
                except TypeError:
                    pc.draw(ax=ax)
                    
                pc.add_labels(middleKeys = [middleKeys[int(s)] for s in struc],
                              outer = False,
                              labelFontSize = label_font_size)
            
            file_name = '{a}/{b}/{c}_strucs'.format(a=self.target_directory,
                                        b=obj,
                                        c=obj)
            
            print ('file name =', file_name)
            fig.savefig('{}.pdf'.format(file_name), 
                        format = 'pdf', 
                        bbox_inches = 'tight')
            fig.savefig('{}.svg'.format(file_name), 
                        format = 'svg', 
                        bbox_inches = 'tight')
            
            plt.close(fig)
            
    
    def create_sample_slices(self, obj):
        middleKeys = {0:["mass_fraction"],
                      1:["mass_fraction"],
                      2:["mass_fraction"],
                      3:["mass_fraction"],
                      4:["mass_fraction"],
                      5:["mass_fraction", "x_SiO2", "x_MgO"],
                      6:["mass_fraction", "x_SiO2", "x_MgO"]}
        
        ranges = self.sorted_sample_layer_ranges[obj]
        n = len(list(ranges.keys())[0])
        fig, ax = plt.subplots(1, len(list(ranges.keys())), figsize = (6,6))
        for i in range(len(list(ranges.keys()))):
            struc = list(ranges.keys())[i]
            data = {struc:ranges[struc]}
        
            pc = CircSegPlot.PlanetChart(data)
            pc.get_ghosts()
            pc.get_slices()
            self.filtered_sample_slices.update({obj:pc})
            try:
                pc.draw(ax=ax[i])
            except TypeError:
                pc.draw(ax=ax)
                
            pc.add_labels(middleKeys = [middleKeys[int(s)] for s in struc],
                          labelFontSize = 8)
        
        file_name = '{a}/{b}/{c}_strucs'.format(a=self.target_directory,
                                    b=obj.split('_')[0],
                                    c=obj)
        
        fig.savefig('{}.pdf'.format(file_name), 
                    format = 'pdf', 
                    bbox_inches = 'tight')
        fig.savefig('{}.svg'.format(file_name), 
                    format = 'svg', 
                    bbox_inches = 'tight')
        fig.savefig('{}.png'.format(file_name), 
                    format = 'png', 
                    bbox_inches = 'tight')
    
    def create_slices_stupid(self, obj, label = None):
        """Creates CircSegPlot pie slice instance for given DataStream instance.
        The slice can then be plotted directly or combined to an overview plot
        for a series of different data streams within a WorkFlow instance.
        """
        if label == None:
            label = obj.capitalize()
        
        layerKeys = [['R_core (km)', 'T_CMB', 'P_CMB', 'xi_MgO m', 'xi_SiO2 m', 'xi_FeO m', 'X_S', 'X_O', 'X_Si', 'fO2'],
                     ['R_HMI (km)', 'T_HMI', 'P_HMI', 'xi_MgO m', 'xi_SiO2 m', 'xi_FeO m', 'X_S', 'X_O', 'X_Si', 'fO2'],
                     ['R', 'T_TBL (K)', 'P_surf', 'xi_MgO m', 'xi_SiO2 m', 'xi_FeO m', 'X_S', 'X_O', 'X_Si', 'fO2']]
        centralKeys = ['R_core (km)', 'T_C', 'P_C', 'xi_MgO m', 'xi_SiO2 m', 'xi_FeO m', 'X_S', 'X_O', 'X_Si', 'fO2']
        keys = ['R_core (km)', 'R_HMI (km)',
        'R']
        layer_vals_max = np.array([[self.filtered_sample_ranges[obj]['{}_max'.format(k)].values[0]
                for k in key]
                          for key in layerKeys])        
        layer_vals_min = np.array([[self.filtered_sample_ranges[obj]['{}_min'.format(k)].values[0]
                for k in key]
                          for key in layerKeys])        
        
        central_vals_max = np.array([self.filtered_sample_ranges[obj]['{}_max'.format(k)].values[0]
                for k in centralKeys])        
        central_vals_min = np.array([self.filtered_sample_ranges[obj]['{}_min'.format(k)].values[0]
                for k in centralKeys])        #layer_vals_min = np.insert(layer_vals_min, [0], [0,0,0, 0, 0, 0, 0, 0, 0], axis=0)
        
        #stupid pressure is in GPa -> convert to Pa
        central_vals_min[2] *= 1e9
        central_vals_max[2] *= 1e9
        #layer_vals_max = np.insert(layer_vals_max, [0], [0,0,0, 0, 0, 0, 0, 0, 0], axis=0)
        #print ('layermax =', layer_vals_max)            
        #print ('layermin =', layer_vals_min)        
        #layer_vals_min[0][1] = self.filtered_sample_ranges[obj]['T_C_min'].values[0]
        #layer_vals_max[0][1] = self.filtered_sample_ranges[obj]['T_C_max'].values[0]
        #layer_vals_min[0][2] = self.filtered_sample_ranges[obj]['P_C_min'].values[0] * 1e9
        #layer_vals_max[0][2] = self.filtered_sample_ranges[obj]['P_C_max'].values[0] * 1e9
        print ('centralmin =', central_vals_min)
        vals_max = [self.filtered_sample_ranges[obj]['{}_max'.format(k)].values[0]
                for k in keys]
        vals_min = [self.filtered_sample_ranges[obj]['{}_min'.format(k)].values[0]
                for k in keys]
        
        #stupid total radius is in earth radii -> convert to km
        vals_min[-1] *= r_earth * 1e-3
        vals_max[-1] *= r_earth * 1e-3
        layer_vals_max[-1][0] *= r_earth * 1e-3
        layer_vals_min[-1][0] *= r_earth * 1e-3


        layerSizes = [v * 1e3 / R_solar[obj] / r_earth for v in vals_max]
        ghostSizes = [v * 1e3 / R_solar[obj] / r_earth for v in vals_min]

        pc = CircSegPlot.PlanetChart(nLayers = 3,
                                     layerSizes = layerSizes,
                                     structure = [4,6, 2])
        
        pc.get_slices()
        pc.get_ghosts(positions = ghostSizes)
        
        self.filtered_sample_slices.update({obj:pc})
        
        fig, ax  = plt.subplots()
        pc.draw(ax = ax, title =  label)
      
        data = {"layerSizes": [[vmin * 1e3, vmax * 1e3]
                               for vmin, vmax in zip(layer_vals_min.T[0], layer_vals_max.T[0])], 
                "layerPressures": [[vmin, vmax]
                               for vmin, vmax in zip(layer_vals_min.T[2], layer_vals_max.T[2])], 
                "layerMasses": [[0, 0] for k in keys], 
                "layerTemperatures": [[vmin, vmax]
                               for vmin, vmax in zip(layer_vals_min.T[1], layer_vals_max.T[1])],
                "S_core": [[vmin, vmax]
                               for vmin, vmax in zip(layer_vals_min.T[6], layer_vals_max.T[6])],
                "Si_core": [[vmin, vmax]
                               for vmin, vmax in zip(layer_vals_min.T[8], layer_vals_max.T[8])],
                "O_core": [[vmin, vmax]
                               for vmin, vmax in zip(layer_vals_min.T[7], layer_vals_max.T[7])],
                "MgO_mantle": [[vmin, vmax]
                               for vmin, vmax in zip(layer_vals_min.T[3], layer_vals_max.T[3])],
                "FeO_mantle": [[vmin, vmax]
                               for vmin, vmax in zip(layer_vals_min.T[5], layer_vals_max.T[5])],
                "SiO2_mantle": [[vmin, vmax]
                               for vmin, vmax in zip(layer_vals_min.T[4], layer_vals_max.T[4])],
                "logfO2_mantle": [[vmin, vmax]
                               for vmin, vmax in zip(layer_vals_min.T[9], layer_vals_max.T[9])]}
        
        central_data = {"Temperature": [central_vals_min.T[1], central_vals_max.T[1]], 
                        "Pressure": [central_vals_min.T[2], central_vals_max.T[2]]}
        
        middleKeys = [ 
                       [], 
                       ['MgO_mantle', 'SiO2_mantle'],
                       ['MgO_mantle']
                       ]
                      
        
        pc.addData(data, centralKeys = centralKeys, middleKeys = middleKeys,
                   central_data = central_data )
        fig_path = self.target_directory + '/' + self.sample_name
        fig.savefig(fig_path + '/{}_slices.pdf'.format(label),
                    format = 'pdf',
                    bbpox_inches = 'tight')


    def plot_sample_profiles(self, obj, ax = [], fnts = 8, save = False):
        if ax == []:
            fig, ax = plt.subplots(1,3, figsize=(16, 4))
            
        else:
            pass
        
        df1 = pd.read_csv("/home/os18o068/Documents/PHD/Projects/Planets/Data/icy_satellites/cammarano_2006_hot_1/hot_1.csv",
                          delimiter=';', decimal=',',
                          names = ["temperature", "depth"])
        
        df2 = pd.read_csv("/home/os18o068/Documents/PHD/Projects/Planets/Data/icy_satellites/cammarano_2006_hot_1/hot_chondritic.csv",
                          delimiter=';', decimal=',',
                          names = ["density", "depth"])
        
        df3 = pd.read_csv("/home/os18o068/Documents/PHD/Projects/Planets/Data/icy_satellites/cammarano_2006_hot_1/hot_pyrolite.csv",
                          delimiter=';', decimal=',',
                          names = ["density", "depth"])
        
        cmap = plt.cm.get_cmap('plasma')
        sm = plt.cm.ScalarMappable(cmap = cmap)
        inset_axis = ax[-1].inset_axes((.6, .6, .05, .3))
        prof = '4602'
        tickvals = [self.sorted_sample_ranges[obj][prof]['moment_of_inertia'][0],
                          sum(self.sorted_sample_ranges[obj][prof]['moment_of_inertia']) / 2,
                          self.sorted_sample_ranges[obj][prof]['moment_of_inertia'][1]]
        
        ticklabels = ['{:.4g}'.format(tv) for tv in tickvals]
        cbar = plt.colorbar(sm, cax = inset_axis, ticks=[0., .5, 1.])
        cbar.ax.set_yticklabels(ticklabels)
        cbar.ax.tick_params(labelsize = fnts)
        cbar.set_label(r'$C/MR^2$', fontsize = fnts)
        ax[0].plot([714, 714], [0, 1500], color = 'k', linestyle = ':')
        for pl in self.filtered_sample_planets[obj]:
            pl.trim_profiles()
            c = (pl.MOI_is - self.sorted_sample_ranges[obj][prof]['moment_of_inertia'][0]) /\
                    (self.sorted_sample_ranges[obj][prof]['moment_of_inertia'][1] -\
                     self.sorted_sample_ranges[obj][prof]['moment_of_inertia'][0])
            
            ax[0].plot(pl.profiles[0] * 1e-3, pl.profiles[1],
                    color = cmap(c))
            ax[1].plot(pl.profiles[0] * 1e-3, pl.profiles[2] * 1E-9,
                       color = cmap(c))
            ax[2].plot(pl.profiles[0] * 1e-3, pl.profiles[3],
                       color = cmap(c))
        
        #Plot literature
        ax[0].plot(df1['depth'], df1['temperature'], color = 'k',
                   linestyle = '-', label ="Cammarano et al. 2006")
        ax[2].plot(1561. - df2["depth"], df2["density"] * 1e3, color = 'k',
                   linestyle = '--', label ='chondritic')
        ax[2].plot(1561. - df3["depth"], df3["density"] * 1e3, color = 'k',
                   linestyle = '-', label="pyrolite")
        ax[0].legend()
        ax[2].legend(loc = 3)
        
        
        for axx in ax:
            axx.set_xlabel(r'${\rm Radius} \ (\rm km)$')
        ax[0].set_ylabel(r'$\rm Temperature \ (K)$')
        ax[1].set_ylabel(r'$\rm Pressure \ (GPa)$')
        ax[2].set_ylabel(r'$\rm Density \ (kg \ m^{-3})$')
        
        if save:
            fig.savefig('{a}/_{b}_profiles.pdf'.format(a=self.target_directory,
                                                       b=obj), 
                        bbox_inches = 'tight',
                        format ='pdf')
            fig.savefig('{a}/_{b}_profiles.png'.format(a=self.target_directory,
                                                       b=obj), 
                        bbox_inches = 'tight',
                        format ='png',
                        dpi=320)        
    
    def plot_sample_structures(self):
        
        self.sample_planet_colors = [[layerColors['outer core'], 
                                      layerColors['outer mantle']] +\
                                     [list(layerColors.values())[int(col)]
                                     for col in str(struc)]
                                     for struc in self.sample_hydro_structures]
            
        self.sample_planet_colors = []
        plot_strucs = []
        for i in range(len(self.sample_hydro_structures)):
            struc = self.sample_hydro_structures[i]
            l = [s for s in str(struc)]
            print ('l =', l)
            if not '3' in l and not '4' in l:
                self.sample_planet_colors.append([layerColors['outer core'],
                                                 layerColors['outer mantle']])
                print ('adding...')
                plot_strucs.append(struc)
                for j in range(len(l)):
                    
                    self.sample_planet_colors[-1].append(list(layerColors.values())[int(l[j])])
        print ('plot strucs =', plot_strucs)
        fig, ax = plt.subplots(len(plot_strucs), figsize = (10, 5))
        plt.subplots_adjust(wspace=0, hspace=0.)
        for k in range(len(plot_strucs)):
            
            struc = plot_strucs[k]
            nLayers = len(str(struc)) + 2
            print ('len =', len(self.sample_planet_colors[k]), nLayers)
            pc = CircSegPlot.PlanetChart(nLayers,
                                         phi1 = np.pi / 2.5,
                                         colors = self.sample_planet_colors[k])
            pc.get_slices()
            try:
                pc.draw(ax = ax[k])
            
            except TypeError:
                pc.draw(ax = ax)
            
        obj = self.sample_name
        objLabel = ''
        
        fig.suptitle(f'{objLabel}', fontsize = 20, ha = 'center')
        fig.savefig(self.target_directory + '/' + obj  + '/' + obj + '_pies.pdf',
                    format = 'pdf',
                    bbox_inches = 'tight')
        
        plt.close(fig)
        

    #This was just a shity and quick fix to plot these stupid structures for
    #the ocean planets. Not usable for general samples.
    def plot_sample_structures_stupid(self, colors = {
                                             'liquid':(.25,.5,1),
                                             'supercritical':(.5,0.25,1),
                                             'solid':(.75,.8,1),
                                             'vapor': (1,0,0),
                                             'gas':(0,1,0),
                                             'core':(.75,.5,.25),
                                             'mantle':(.4,.25,0)}):
        
        fig, ax = plt.subplots(3, 6, 
                               figsize = (10, 5))
        plt.subplots_adjust(wspace=0, hspace=0.)
        
        self.sample_planet_colors = [[colors['core'], colors['mantle']] +\
                                     [list(colors.values())[int(col)]
                                     for col in str(struc)]
                                     for struc in self.sample_hydro_structures]
            
        self.sample_planet_colors = []
        plot_strucs = []
        for i in range(len(self.sample_hydro_structures)):
            struc = self.sample_hydro_structures[i]
            l = [s for s in str(struc)]
            print ('l =', l)
            if not '3' in l and not '4' in l:
                self.sample_planet_colors.append([colors['core'],
                                                 colors['mantle']])
                print ('adding...')
                plot_strucs.append(struc)
                for j in range(len(l)):
                    
                    self.sample_planet_colors[-1].append(list(colors.values())[int(l[j])])
        
        
        for i in range(3):
            for j in range(6):
                k = 6 * i + j
                struc = plot_strucs[k]
                nLayers = len(str(struc)) + 2
                print ('len =', len(self.sample_planet_colors[k]), nLayers)
                pc = CircSegPlot.PlanetChart(nLayers,
                                             phi1 = np.pi / 2.5,
                                             colors = self.sample_planet_colors[k])
                
                try:
                    pc.draw(ax = ax[i][j])
                
                except TypeError:
                    pc.draw(ax = ax)
                
        obj = self.sample_name
        objLabel = ''
        
        fig.suptitle(f'{objLabel}', fontsize = 20, ha = 'center')
        fig.savefig(self.target_directory + '/' + obj  + '/' + obj + '_pies.pdf',
                    format = 'pdf',
                    bbox_inches = 'tight')
        
        plt.close(fig)
    

    def plot_structures(self, colors = {'supercritical':(.5,0.25,1),
                                             'solid':(.75,.8,1),
                                             'liquid':(.25,.5,1),
                                             'core':(.75,.5,.25),
                                             'mantle':(.4,.25,0)}):
                            
        for p in range(len(self.planets)):
            fig, ax = plt.subplots(1, len(self.all_planet_structures[p]), 
                                   figsize = (10, 5))
            for i in range(len(self.all_planet_structures[p])):
                nLayers = len(list(self.all_planet_locs[p].values())[i][0]) - 1
                pc = CircSegPlot.PlanetChart(nLayers,
                                             phi1 = np.pi / 2.5,
                                             colors = self.all_planet_colors[p][i])
                
                try:
                    pc.draw(ax = ax[i])
                
                except TypeError:
                    pc.draw(ax = ax)
                
                title = len(list(self.all_planet_locs[p].values())[i]) / len(self.planets[p])
                pc.addData({k:self.all_planet_ranges[p][k][i] 
                            for k in list(self.all_planet_ranges[p].keys())},
                           title = f'{100*title:.0f}\%')
                
            obj = self.archive.df['pl_cleanname'].values[p]
            objLabel = self.archive.df['pl_name'].values[p]
            fig.suptitle(f'{objLabel}', fontsize = 20, ha = 'center')
            fig.savefig(self.target_directory + '/' + obj  + '/' + obj + '_pies.pdf',
                        format = 'pdf',
                        bbox_inches = 'tight')
            
            plt.close(fig)
    
    
    def plot_data(self):
        """

        Returns
        -------
        None.

        """
        for i in range(len(self.planets)):
            try:
                xmin = np.nanmin(self.data[i]['M'].values)
                xmax = np.nanmax(self.data[i]['M'].values)
                xmid = (xmin + xmax)*.5
                x_range = (xmid - xmin)/xmid
                
                ymin = np.nanmin(self.data[i]['R'].values)
                ymax = np.nanmax(self.data[i]['R'].values)
                ymid = (ymin + ymax)*.5
                y_range = (ymid - ymin)/ymid
                
                indices = np.arange(0, len(self.data[i]))
        
                cols = np.array([[0, 0., 0.25], [0.1, 0.1, .5], [0.5, 0.5, 1.]])
                cmap = ftool.my_dev_cmap(N_col = 32, border_colors = cols)
                
                cols = np.array([[.8, .8, 0.0], [0.8, 0.5, .0], [1., 0., 0.]])
                cmap_arch = ftool.my_dev_cmap(N_col = 32, border_colors = cols)
                
                obj = self.archive.df['pl_cleanname'].values[i]
                spec = self.target_directory + '/' + obj
                PlanetInvestigator.plot_models(self.data[i].values, 
                                               self.planets[i], 
                                               indices, 
                                               y_real= ymid, 
                                               obj = obj,
                                               reldev_x = 1., 
                                               reldev_y = 1., 
                                               plot_label = 'All models',
                                               spec = spec, 
                                               x_real = xmid, 
                                               x_range = x_range*1.1, 
                                               y_range = y_range*1.1,
                                               plot_rect = False, 
                                               plot_cross = False, 
                                               mass = 1., 
                                               xy = [46, 1],
                                               x_lims = [0., 10.], 
                                               y_lims = [.75, 3.], 
                                               plot_PT = False,
                                               plot_PREM = False, 
                                               plot_CMB = False, 
                                               profile_cbar_loc = 1,
                                               xy_scale = ['linear', 'linear'], 
                                               profile_x_lims = [0., 2.5],
                                               profile_y_lims = [[0., 15000.], [0, 5500], [0., 30.]], 
                                               plot_real = True, 
                                               colormap = cmap, 
                                               colormap_arch = cmap_arch,
                                               real_filter = True, 
                                               max_mass_err = .15)
                '''
                PlanetInvestigator.check_statistics(all_all_out,
                                                    obj = '',
                                                    spec = obj,
                                                    write = True,
                                                    N_all = len(all_planets),
                                                    N = 1)
                '''
                plt.close('all')
            except ValueError:
                pass
    
    
    def box_plot(self, bar_params = ['w_liq', 'w_liq/w_H2O', 'h_liq', 'V_liq'], 
                 bar_facs = [1., 1., 1e-3, 1.], whis = 100,
                 bar_cols = ['r', 'g', 'b', 'k'], dx = 0., top_label = 'pl_eqt',
                 **kwargs):
        if 'ax' in kwargs:
            ax = kwargs.get('ax')
            
        else:
            fig, ax = plt.subplots(1, len(bar_params))
        
        y_labels = [r'$w_{\rm liq}$',
                    r'$w_{\rm liq}/w_{\rm H_2O}$',
                    r'$h_{\rm liq} \ {\rm \ [km]}$',
                    r'$V_{\rm liq} / V_{\rm liq \oplus}$']
        
        column_labels = [r'${\rm Total \ mass \ fraction}$',
                         r'${\rm Relative \ mass \ fraction}$',
                         r'${\rm Cummulative \ layer \ depth}$']
        fnts2 = 14
        dy = .05
        delta_y = 1.
        top_label_props = {'pl_eqt':{'format':'', 
                                     'unit':'K'},
                           'pl_bmasse':{'format':'', 
                                        'unit':'M_\oplus'}}
        
        for i in range(len(ax)):
            for k in range(len(ax[0])):
                #ax[i].set_title(column_labels[i], fontsize = fnts2)
                ind = 2 * i + k
                ax[i][k].set_xticklabels(self.archive.df['pl_name'], 
                                         rotation = 45, 
                                         minor = False,
                                         ha = 'right')
                ax[i][k].set_xticks(np.arange(0, 
                                              len(self.archive.df['pl_cleanname'])))
                ax[i][k].set_xlim(-dx * 1.5 * 2, 
                                  len(self.archive.df['pl_cleanname']) - 1 + dx * 1.5 * 2)
                ax[i][k].set_ylabel(y_labels[ind])
                ax[i][k].set_yscale('log')
                trans = transforms.blended_transform_factory(ax[i][k].transData, 
                                                         ax[i][k].transAxes) 

                
                for j in range(len(self.data)):
                    xmin = j - dx
                    
                    bplot = ax[i][k].boxplot(self.data[j][bar_params[ind]].values \
                                             * bar_facs[ind],
                                     positions = [xmin], showfliers = False,
                                     whiskerprops = {'color' : bar_cols[i]},
                                     patch_artist = True, whis = whis)
    
                    for patch in bplot['boxes']:
                        patch.set_facecolor(bar_cols[i])
                        
                    #Plot guide bars
                    ax[i][k].add_patch(Rectangle((j-dx*(2), 0), dx * 4, delta_y, 
                                          facecolor = 'k', fill = True,
                                          transform = trans, zorder = 0,
                                          clip_on = False, alpha = .1))
                    
                    #Generate strings for top labels (value, +error, -error)
                    top_label_str0 =\
                        str((self.archive.df[top_label][self.archive.df.index[j]]))
                    top_label_str1 =\
                        str((self.archive.df['{}err1'.format(top_label)][self.archive.df.index[j]]))
                    top_label_str2 =\
                        str((self.archive.df['{}err2'.format(top_label)][self.archive.df.index[j]]))
                    
                    top_label_str =\
                        r"${{{a}}}^{{+{b}}}_{{{c}}} \ \rm {u}$".format(a=top_label_str0, 
                                                                     b = top_label_str1, 
                                                                     c = top_label_str2,
                                                                     u = top_label_props[top_label]['unit'])
                    if i == 0:
                        ax[i][k].text(j, 1.05, top_label_str,
                                   transform = trans, rotation = 45, ha = 'left')
                ax[i][k].set_xticks(range(0, len(self.archive.df['pl_cleanname'])))
            
            
################################################################################ 
class WorkFlow():
    def  __init__(self, directories = [], create = False, stats = False, 
                  sortBy = 'pl_eqt', load = True, write = False):
        self.data_streams = []
        self.ranges = None
        self.combined_ranges = None
        #Initiate data streams for each target directory
        for dir in directories:
            ds = DataStream(directory = dir)
            #Create working data from raw data
            if create:
                ds.create_working_data(stats = stats, load = load, write = write)
            #Use pre-processed data as input
            else:
                ds.load_working_data()
                ds.get_cleaned_data()
                
            ds.get_hydro_structures()
            self.data_streams.append(ds)
            
            
    def sortData(self, sortBy = 'pl_eqt'):
        for ds in self.data_streams:
            #Sort data according to equilibrium temperature
            ds.data = [x for _, x in sorted(zip(ds.archive.df[sortBy], ds.data))]
            ds.archive.df = ds.archive.df.sort_values(by = [sortBy])


    def create_table(self):
        """Compiles the data into csv tables and saves them to files for 
        subsequent use in Latex tables.
        """
        
        selectedParams = {'ocean_frac':[1e2, 2], 
                          'w_liq':[1e2,2], 
                          'w_liq/w_H2O':[1e2,2], 
                          'h_liq':[1e-3,2], 
                          'V_liq':[1,2],
                          'R_core (km)':[1,2],
                          'R_HMI (km)':[1,2]}
        
        hydroTypes = {'0':'a', '10':'b', '210':'c', '1210':'d'}
        '''
        hydroTypes = {'21020':'a',
                      '2120': 'b',
                      '1210': 'c',
                      '1202': 'd',
                      '121020': 'e',
                      '121': 'f',
                      '210': 'g',
                      '1020': 'h',
                      '120': 'i',
                      '202': 'j',
                      '12120': 'k',
                      '21': 'l',
                      '20': 'm',
                      '2': 'n',
                      '0': 'o',
                      '10': 'p',
                      '12': 'q',
                      '1': 'r'}
        '''
        #l = len(self.data_streams)
        #m = len(self.data_streams[0].data)
        #n = len(selectedParams)
        #all_rows = [[[[None] * 2] * l] * n] * m
        all_rows = []
        all_hyd_types = []
        #Create txt file for LaTex tables of the data of each DataStream
        for i in range(len(self.data_streams)):
            ds = self.data_streams[i]
            rows = []
            all_rows.append([])
            all_hyd_types.append([])
            #Loop over planets
            for j in range(len(ds.data)):
                name = ds.archive.df['pl_name'].values[j]
                vals = [name]
                all_rows[i].append([])
                #Loop over different parameters
                for k in range(len(list(selectedParams.keys()))):
                    all_rows[i][j].append([])
                    param = list(selectedParams.keys())[k]
                    try:
                        #extract value range
                        min_ = np.nanmin(ds.data[j][param].values)
                        max_ = np.nanmax(ds.data[j][param].values)
                        min_ *= selectedParams[param][0]
                        max_ *= selectedParams[param][0]
                        all_rows[i][j][k] = [min_, max_]
                        p = selectedParams[param][1]
                        
                        string1 = '{:g}'.format(float('{:.{p}g}'.format(min_, p=p)))
                        string2 = '{:g}'.format(float('{:.{p}g}'.format(max_, p=p)))
                        vals.append(string1 + '-' + string2)
                    
                    #handle empty data
                    except ValueError:
                        all_rows[i][j][k] = [np.nan, np.nan]
                        vals.append('...')
                
                #append hydro structure
                strucs = []
                for s in list(set(ds.data[j]['hyd_struc'].values)):
                    strucs.append(hydroTypes[str(int(s))])
                
                print ('struc =', strucs)
                all_hyd_types[i].append(strucs)
                strucString = ','.join(strucs)
                vals.append(strucString)
                row = ' & '.join(vals) + ' \\\\' 
                rows.append(row)
                
            #write to file
            path = ds.target_directory + '/ranges_table.txt'
            with open(path, 'w+') as file:
                for row in rows:
                    file.write(row + '\n')
        
        
        #Create file for the data of all DataStreams
        #jesus that's ugly as fuck
        all_rows = np.array(all_rows)
        rows = []
        vals = np.empty((len(self.data_streams[0].data), 
                         len(selectedParams) + 2), 
                        dtype='object')
    
        #loop over params
        for k in range(len(all_rows.T[0])):
            param = list(selectedParams.keys())[k]
            
            #loop over planets
            for j in range(len(all_rows.T[0][k])):
                name = self.data_streams[0].archive.df['pl_name'].values[j]
                vals[j][0] = name
    
                try:
                    #extract value range
                    
                    min_ = np.nanmin(np.asarray(all_rows.T[0][k][j]))
                    max_ = np.nanmax(np.asarray(all_rows.T[1][k][j]))
                    p = selectedParams[param][1]
                    string1 = '{:g}'.format(float('{:.{p}g}'.format(min_, p=p)))
                    string2 = '{:g}'.format(float('{:.{p}g}'.format(max_, p=p)))
                    string = string1 + '-' + string2
                    
                    vals[j][k + 1] = string
                    
                #handle empty data
                except ValueError:
                    vals[j][k + 1] = '...'     
        
        # Add hydro structure
        hyd_types = []
        #Loop over planets
        for i in range(len(all_rows.T[0][0])):
            hyd_types.append([])
            #Loop over over data streams
            for j in range(len(all_hyd_types)):
                #Loop over structures
                for k in range(len(all_hyd_types[j][i])):
                    
                    hyd_types[i].append(all_hyd_types[j][i][k])
                    
            hyd_types[i] = list(set(hyd_types[i]))
            
            #append hydro structure
            strucs = []
            vals[i][-1] = ','.join(hyd_types[i])
            row = ' & '.join(vals[i]) + ' \\\\' 
            rows.append(row)
        
        #write to file
        path = ds.data_root + '/ranges_table.txt'
        with open(path, 'w+') as file:
            for row in rows:
                file.write(row + '\n')
                 
        return vals

    
    def get_ranges(self):
        for ds in self.data_streams:
            ds.get_ranges_stupid()
            
        self.params = list(self.data_streams[0].all_planet_ranges[0].keys())
        self.ranges = [{} 
                       for i in range(len(self.data_streams[0].planets))]
        self.combined_ranges = [{} 
                                for i in range(len(self.data_streams[0].planets))]
        
        #Collect ranges for each data stream        
        for ds in self.data_streams:
            ds.get_ranges()
            #Loop over all planets in the data stream
            for i in range(len(ds.all_planet_structures)):
                  
                #Loop over all possible structures for the planet
                for j in range(len(ds.all_planet_structures[i])):
                    struc = ds.all_planet_structures[i][j]

                    #If specific structure was not already identified for the
                    #current planet, add it.
                    if not struc in self.ranges[i].keys():
                        self.ranges[i].update({struc:{p:[] for p in self.params}})
                        self.combined_ranges[i].update({struc:{p:[] for p in self.params}})
                    
                    #Get ranges for all layer parameters
                    for param in self.params:
                        d = ds.all_planet_ranges[i][param][j]
                        self.ranges[i][struc][param].append(d)            
        
        #Combine ranges of the different data streams to overall ranges
        #Loop over planets
        for i in range(len(self.ranges)):
            strucs = self.ranges[i].keys()
            #Loop over possible structures for each planet
            for struc in strucs:
                #Extract all relevant layer parameters
                for param in self.params:
                    #Note that the array shape needed for the structure chart
                    #plot needs to be inverted.
                    arr = np.array(self.ranges[i][struc][param]).T
                    dummy = np.zeros([len(arr), len(arr[0])])
                    
                    #Compute min and max for each parameter
                    for k in range(len(dummy[0])):
                        dummy[0][k] = min(arr[0][k])
                        dummy[1][k] = max(arr[1][k])
                    
                    #Add overall min and max to array
                    for d in dummy.tolist():
                        self.combined_ranges[i][struc][param].append(d)
    

    def plot_structures(self, ):
        if self.combined_ranges == None:
            self.get_ranges()
        
        #Loop over planets
        for p in range(len(self.data_streams[0].planets)):
            nStrucs = len(list(self.combined_ranges[p].keys()))
            fig, ax = plt.subplots(1, nStrucs, 
                                   figsize = (10, 5))
            
            #Loop over structures for the planet
            for i in range(nStrucs):
                struc = list(self.combined_ranges[p].keys())[i]
                nLayers = len(self.combined_ranges[p][struc]['layerSizes'][0]) -1
                data = {k:np.array(self.combined_ranges[p][struc][k]).T for k in self.params}
 
                pc = CircSegPlot.PlanetChart(data,
                                             phi1 = np.pi / 2.5)
                
                try:
                    pc.draw(ax = ax[i])
                
                except TypeError:
                    pc.draw(ax = ax)
                    
                data = {k:np.array(self.combined_ranges[p][struc][k]).T for k in self.params}
                #Get number of models for the current planet and structure for
                #each data stream individually
                numStruc = 0
                numTot = 0
                for ds in self.data_streams:
                    try:
                        numStruc += ds.all_planet_counts[p][struc]
                    #current structure does not occur in the current data stresam
                    except KeyError:
                        pass
                    
                    numTot += sum(ds.all_planet_counts[p].values())
                
                percentage = numStruc / numTot
                if nStrucs > 1:
                    if percentage > .99:
                        title = r'>99\%'
                    elif percentage < .01:
                        title = r'<1\%'
                    else:
                        title = f'{100*percentage:.0f}\%'
                else:
                    title = f'{100*percentage:.0f}\%'
            
            obj = self.data_streams[0].archive.df['pl_cleanname'].values[p]
            objLabel = self.data_streams[0].archive.df['pl_name'].values[p]
            fig.suptitle(f'{objLabel}', fontsize = 20, ha = 'center')
            fig.savefig(self.data_streams[0].data_root + '/' + obj + '_pies.svg',
                        format = 'svg',
                        bbox_inches = 'tight')
            fig.savefig(self.data_streams[0].data_root + '/' + obj + '_pies.pdf',
                        format = 'pdf',
                        bbox_inches = 'tight')            
            plt.close(fig)    
            

    def box_plot(self, dx = .2, x0 = .55, y0 = .9, dy = .05, d_stripe = .1,
                               sortBy = 'pl_eqt', whis=100, name = 'hydro_bars',
                               pres = False, top_label = 'pl_eqt'):
        fig, ax = plt.subplots(2, 2, figsize = (9, 7), sharex = True)
        '''
        for ds in self.data_streams:
            #Sort data according to specified parameter
            ds.data = [x for _, x in sorted(zip(ds.archive.df[sortBy], ds.data))]
            ds.archive.df = ds.archive.df.sort_values(by = [sortBy])
        '''
        plt.subplots_adjust(hspace = .1)
        bar_cols = [[(.5, .5, .5) for i in range(3)], 
                    [(.5, .0, .5) for i in range(3)], 
                    [(.25, .5, 1.) for i in range(3)]]
        
        if pres:
            bar_labels = [r'$\rm 1 \ bar$',
                          r'$\rm 10 \ bar$',
                          r'$\rm 100 \ bar$']
            name = "{}_pres".format(name)
        
        else:
            bar_labels = [r'$\rm Mg/Fe = 0.5$',
                          r'$\rm Mg/Fe = 0.57$',
                          r'$\rm Mg/Fe = 0.8$']
        
        bar_facs = [1.,
                    1.,
                    1e-3,
                    1.]
        
        for i in range(len(self.data_streams)):
            ds = self.data_streams[i]
        
            ds.box_plot(ax = ax, dx = -dx + dx * i, bar_cols = bar_cols[i],
                        bar_facs = bar_facs, whis = whis,
                        top_label = top_label)
            
        #Legend
        trans = transforms.blended_transform_factory(ax[0][0].transAxes, 
                                                     ax[0][0].transAxes) 
        
        for i in range(len(self.data_streams)):
            print (x0, y0 + dy * i)
            #Add legend bar
            ax[0][0].add_patch(Rectangle((x0, y0 - dy * i), d_stripe, dy / 2, 
                                      facecolor = bar_cols[i][0], 
                                      fill = True, 
                                      zorder = 30,
                                      transform = trans))
            
            #Add legend label
            ax[0][0].text(x0 + d_stripe + .05, y0 - dy * i, bar_labels[i], 
                          zorder = 30,
                          transform = trans)

        # Add legend box background
        ax[0][0].add_patch(Rectangle((x0 - .01, y0 - dy * 3 + dy / 2), .425, 7 * dy / 2, 
                                     zorder = 20,
                color = 'white', alpha = .75, transform = trans))
        
        fig.savefig('/home/os18o068/Documents/PHD/Projects/Planets/Data/'+\
                    name + '.pdf', format = 'pdf', bbox_inches = 'tight')
        plt.close(fig)
