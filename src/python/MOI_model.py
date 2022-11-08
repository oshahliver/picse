# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 20:11:41 2020

@author: shlah
"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from PIMPphysicalparams import m_earth, r_earth, G, R_solar, M_solar, MoI_solar
import random
import functionTools as ftool
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import pandas as pd
import plotly.express as px
import plotly.io as pio
pio.renderers.default = 'browser'

matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

R_earth = 6.371e6
M_earth = 5.972e24
omega = 1.6098e-5


layer_colors = np.array([[0,0, .75], [0, 0, 1], [.5, .25, 0], [1, 0, 0]])

layer_colors_finals = np.array([[.25, .25, .75],[.5, .5, 1], [.75, .5, 0], [1., 1., 1.]])

im_dir = '/home/os18o068/Documents/PHD/Projects/Planets/Data/icy_satellites/revived_project/'

#Values from Schubert et al. 2004 and Iess 2010
#Io, Europa, Ganymede, Callisto, Titan, Enceladus
j2 = [1859.5, 435.5, 127.53, 32.7, 33.46, 5435.2]
c22 = np.array([558.8, 131.5, 38.26, 10.2, 10.022, 1549.8])
dj2 = np.array([2.7, 8.2, 2.9, .8, .632, 34.9])
dc22 = np.array([.8, 2.5, .87, .3, .071, 15.6])

references_names = ['Schubert et al. 2004', 
              'Iess et al. 2010 (SOL2)', 
              'Iess et al. 2014']

references = [0, 0, 0, 0, 1, 2]

colors = ['red', 'green', 'blue']

#Data for comparison between hydrous and anhydrous mantles for core-less 
#planets and T=300

masses = np.array([.001, .005, .01, .02, .03, .04, .05, .06, .07, .08, .09, .1])
radii = np.array([[.1193, .2037, .2563, .3223, .3684, .4049, .4356, .4618,
                   .4851, .5062, .5256, .5436], #hydrous
                  [.1193, .2037, .2563, .3221, .368, .4044, .4347, .4606,
                   .4834, .5041, .5231, .5406]]) #anhydrous
MoI = np.array([[.3978, .3973, .3969, .3963, .3959, .3955, .3948, .3934, .3923, 
                 .3915, .3909, .3904],
                [.3978, .3972, .3967, .396, .3954, .3948, .3938, .3919, .3903,
                 .3892, .3883, .3876]])


#Data for solar composition and total masses of 0.0225 (~Titan) and 0.025 (Ganymede) and 0.008 (Europa)
data = np.array([[[.3025, .2956, .2894, .2843, .2802, .2775, .2763, .2769, .2798, .2855, .295, .31, .3342], 
                  [.422, .4151, .4077, .3999, .3921, .3837, .375, .3658, .356, .3455, .3341, .3213, .3063],
                  [.6, .55, .5, .45, .4, .35, .3, .25, .2, .15, .1, .05, 0.]], 
                 [[.3025, .2955, .2895, .2844, .2804, .2777, .2765, .2771, .2799, .2855, .295, .3099, .334],
                  [.436, .4287, .4213, .4135, .4052, .3968, .3878, .3784, .3683, .3575, .3457, .3326, .3171],
                  [.6, .55, .5, .45, .4, .35, .3, .25, .2, .15, .1, .05, 0.]], 
                 [[.3027, .2952, .2886, .2831, .2788, .276, .2747, .2754, .2785, .2847, .2949, .3109, .3357],
                  [.3044, .299, .2935, .2878, .2818, .2757, .269, .2621, .2547, .2469, .2383, .2288, .218],
                  [.6, .55, .5, .45, .4, .35, .3, .25, .2, .15, .1, .05, 0.]]])



def compute_curve_difference(x, y1, y2, method = "area"):
    """Computes the difference between two curves using the specified method.
    The input curves must have the same shape and the same x values.
    """
    
    #Use linear estimation of area between the curves over discrete intervals
    if method == "area":
        areas = []
        #Loop over discrete invervals
        for i in range(len(x) - 1):
            dx = x[i + 1] - x[i]
            dy1 = y1[i] - y2[i]
            dy2 = y1[i + 1] - y2[i + 1]
            dy = .5 * (dy1 + dy2)
            areas.append(abs(dy * dx))
        
        return sum(areas)
    
    #Use distance between y values
    elif method == "distance":
        dists = []
        #Loop over points
        for i in range(len(x)):
            dy = y1[i] - y2[i]
            dists.append(abs(dy))
            
        return sum(dists)




class analyticalModelPopulation():
    def __init__(self, total_mass = 1., res = [0, 1, 2], n_layers = 3,
                 n_objects = 10, core_mantle_ratio = None,
                 n_models = 10):
        self.n_layers = n_layers
        self.total_mass = total_mass
        self.n_objects = n_objects
        self.res = res
        self.ocean_mass_fractions = np.linspace(0, .6, n_objects)
        self.core_mantle_ratio = core_mantle_ratio
        self.n_models = n_models
        
        
    def create_curves(self):
        self.objects = []
        for r in self.res:
            self.objects.append([])
            for i in range(self.n_objects):
                amo = analyticalModelObjects(n_layers = self.n_layers,
                                         res = [r],
                                         total_mass = self.total_mass,
                                         ocean_mass_fraction = self.ocean_mass_fractions[i],
                                         core_mantle_ratio = self.core_mantle_ratio,
                                         n_models = self.n_models)
                amo.create_spirits()
                self.objects[-1].append(amo)
                
                
    def extract_best_fits(self):
        pass
    
                
    def plot(self):
        x = [[obj.spirits[0].data['ocean'] for obj in objs]
             for objs in self.objects]
        y = [[obj.spirits[0].data['MOI'] for obj in objs]
              for objs in self.objects]
             
        fig, ax = plt.subplots()
        ax.set_xlabel('ocean mass fraction')
        ax.set_ylabel('MOI factor')
        for i, j in zip(range(len(x)), range(len(y))):
            ax.scatter(x[i], y[j], marker = 'o', label = 'res = {}'.format(self.res[i]))
             
        ax.legend()


class analyticalModelObjects():
    def __init__(self, n_layers = 3, n_models = 1000, total_mass = 1., res = [0, 2, 4],
                 ocean_mass_fraction = None,
                 core_mantle_ratio = None):
        self.n_layers = n_layers
        self.total_mass = total_mass
        self.res = res
        self.n_models = n_models
        self.ocean_mass_fraction = ocean_mass_fraction
        self.core_mantle_ratio = core_mantle_ratio
    
    
    def create_spirits(self):
        self.spirits = []
        for r in self.res:
            sp = analyticalModelSpirit(n_models = self.n_models,
                                       n_layers = self.n_layers,
                                       total_mass = self.total_mass,
                                       res = r)
            
            sp.get_layer_masses(core_mantle_ratio = self.core_mantle_ratio,
                                ocean_mass_fraction = self.ocean_mass_fraction)
            sp.create_models()
            sp.construct_models()
            
            self.spirits.append(sp)
            
    
    def plot(self):
        fig, ax = plt.subplots(1, len(self.spirits), figsize = (12, 3))
        plt.subplots_adjust(wspace = .75)
        for axx in ax:
            axx.set_xlabel(r"$\rm Water \ mass \ fraction$")
            axx.set_ylabel(r"$\rm Moment \ of \ inertia \ factor$")
        try:
            for sp, axx in zip(self.spirits, ax):
                sp.plot(ax = axx, surface = False)
        except TypeError:
            for sp in self.spirits:
                sp.plot(ax = ax)
                
        fig.savefig("./analytical_model_n_layer.png", format="png", bbox_inches = "tight")


class analyticalModelSpirit():
    def __init__(self, n_models = 100, n_layers = 3, res = 0, total_mass = 1.,
                 dens_ranges = np.array([[1e4, 1.2e4],
                        [4.0e3, 5.0e3],
                        [1e3, 1.5e3],
                        [9e2, 1e3]])):
        
        self.dens_ranges = dens_ranges
        self.n_models = n_models
        self.n_layers = n_layers
        self.total_mass = total_mass
        self.res = res
        self.layer_mass_fractions = [np.zeros([n_layers]) for i in range(n_models)]
    
    
    def filter(self, obj, tolerance = [.01, .01, .01]):
        """Filter out objects for which the mass, radius and MoI do not match
        with the observed values of the specified object within the specified
        tolerance. If the mass was already set to the objects mass it will
        automatically match and only the radius and MoI matter.
        """
        self.filtered_models = []
        for mod in self.models:
            reldev_R = abs(R_solar[obj] - mod.total_radius) / R_solar[obj]
            reldev_M = abs(M_solar[obj] - mod.total_mass) / M_solar[obj]
            reldev_MoI = abs(MoI_solar[obj] - mod.MOI) / MoI_solar[obj]
            
            if reldev_R < tolerance[0] and reldev_M < tolerance[1] and reldev_MoI < tolerance[2]:
                self.filtered_models.append(mod)
                print (reldev_R, reldev_M, reldev_MoI)    
                
            else:
                pass

        data = np.array([[mod.total_mass, 
                          mod.total_radius, 
                          mod.MOI, 
                          sum(mod.layer_masses[-2**self.res:]) /\
                              self.total_mass,
                            mod.layer_masses[0] / sum(mod.layer_masses)] 
                         for mod in self.filtered_models])
        labels = ["mass", "radius", "MOI", "water_fraction", 
                  "core_mass_fraction"]
        self.filtered_data = pd.DataFrame(data, columns = labels)
        
    
    def fix_ocean_mass(self, mass_fraction):
        for i in range(self.n_models):
            self.layer_mass_fractions[i][-1] = mass_fraction
    
    
    def get_layer_masses(self, core_mantle_ratio = None,
                         ocean_mass_fraction = None):
        """
        Computes layer masses. If no ocean mass fraction and/or core mantle mass
        ratio is given the layer masses will be chosen randomly to add up to the
        total specified mass of the spirit.

        Parameters
        ----------
        core_mantle_ratio : TYPE, optional
            DESCRIPTION. The default is None.
        ocean_mass_fraction : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        None.

        """
        for i in range(self.n_models):
            off = 0
            #Ocean mass fraction already fixed
            if not ocean_mass_fraction == None:
                # Fix outer layer mass
                self.layer_mass_fractions[i][-1] = ocean_mass_fraction
                off += 1
            
            else:
                # Probe outer layer mass
                self.layer_mass_fractions[i][-1] = random.random() 
                off += 1
            
            # Core mantle ratio already fixed
            if not core_mantle_ratio == None:
                m = self.total_mass - sum(self.layer_mass_fractions[i])
                self.layer_mass_fractions[i][1] = m / (1. + core_mantle_ratio)
                off += 1
                
            for j in range(self.n_layers - 1 - off):
                m = random.random() * (1 - sum(self.layer_mass_fractions[i]))
                self.layer_mass_fractions[i][j] = m
            
            self.layer_mass_fractions[i][-1 - off] = 1. - sum(self.layer_mass_fractions[i])
                
    
    def create_models(self):
        self.models = [analyticalModel(n_layers = self.n_layers,
                                       res = self.res,
                                       total_mass = self.total_mass,
                                       layer_mass_fractions = self.layer_mass_fractions[i])
                       for i in range(self.n_models)]    
    
    
    def construct_models(self):
        for mod in self.models:
            mod.draw_random(dens_ranges = self.dens_ranges)
            mod.get_layer_masses()
            mod.get_layer_radii()
            mod.compute_moi()
            mod.compute_ocean()
        
        
        omf = sum(mod.layer_masses[-2**self.res:]) / self.total_mass
        
        data = np.array([[mod.total_mass, 
                          mod.total_radius, 
                          mod.MOI, 
                          mod.ocean_mass_fraction,
                            mod.layer_masses[0] / sum(mod.layer_masses)] 
                         for mod in self.models])
            
        labels = ["mass", "radius", "MOI", "water_fraction", 
                  "core_mass_fraction"]
        self.data = pd.DataFrame(data, columns = labels)
        
    
    def get_ranges(self):
        self.ranges = [self.data['MOI'].min(), self.data['MOI'].max()]
        self.filtered_ranges = [self.filtered_data['MOI'].min(), self.filtered_data['MOI'].max()]
        
        
    def plot(self, ax = None, fig = None, surface = True, obj = ''):
    
        if ax == None:
            if surface:
                fig = plt.figure()
                ax = fig.add_subplot(projection = '3d')
            else:
                fig, ax = plt.subplots(figsize = (6,6))
        
            
        X = self.data['core_mass_fraction']
        Y = self.data['water_fraction']
        Z = self.data['radius']
        C = self.data['MOI']
        axes_labels = [r"$\rm Core \ mass \ fraction$", 
                       r"${\rm Radius} \ (R_\oplus)$",
                       r"$\rm Water \ mass \ fraction$"]
        
        if surface:
            surf = ax.scatter(Z, Y, X, c = C, cmap = 'jet', vmin = .25, vmax = .4)
            ax.set_xlabel(axes_labels[1])
            ax.set_ylabel(axes_labels[2])
            ax.set_zlabel(axes_labels[0])
            plt.colorbar(surf, ax = ax, label = r"$\rm Moment \ of \ inertia \ factor$")
            
            
        else:
            sc = ax.scatter(self.data['radius'], self.data['MOI'], 
                       c = self.data['water_fraction'],
                       cmap = 'ocean')
            
            inset_axis = ax.inset_axes((1.05, 0., .05, 1.))
            cbar = plt.colorbar(sc, cax = inset_axis)
            cbar.set_label(r'$\rm Water \ mass \ fraction$')
            ax.set_xlabel(r'${\rm Radius} \ (R_\oplus)$')
            ax.set_ylabel(r'$\rm MoI$')
        
        fig.savefig("{}_layers.png".format(self.n_layers), format = "png", bbox_inches = "tight")
        fig1 = px.scatter_3d(self.data, x = 'water_fraction',
                            y = 'MOI',
                            z = 'radius',
                            color = 'core_mass_fraction')
        
        fig1.show()
        fig1.write_html('{a}/{b}/{b}.html'.format(a = im_dir,
                                         b = obj))
        
        
    def plot_profiles(self, ax = None, filtered = False, fnts = 8):
        if ax == None:
            fig, ax = plt.subplots()
        
        else:
            pass
        
        cmap = ftool.my_dev_cmap()
        cmap = plt.cm.get_cmap('jet')
        sm = plt.cm.ScalarMappable(cmap = cmap)
        inset_axis = ax.inset_axes((.7, .6, .05, .3))
        
        if filtered:
            tickvals = [self.filtered_ranges[0],
                          sum(self.filtered_ranges) / 2,
                          self.filtered_ranges[1]]
        else:
            tickvals = [self.ranges[0],
                          sum(self.ranges) / 2,
                          self.ranges[1]]
        
        ticklabels = ['{:.4g}'.format(tv) for tv in tickvals]
        cbar = plt.colorbar(sm, cax = inset_axis, ticks=[0., .5, 1.])
        cbar.ax.set_yticklabels(ticklabels)
        cbar.ax.tick_params(labelsize = fnts)
        cbar.set_label(r'$C/MR^2$', fontsize = fnts)

        
        if filtered:
            for mod in self.filtered_models:
                c = (mod.MOI - self.filtered_ranges[0]) /\
                    (self.filtered_ranges[1] - self.filtered_ranges[0])
                mod.plot_profile(ax = ax,
                                 color = cmap(c))
        else:
            for mod in self.models:
                mod.plot_profile(ax = ax)

        ax.set_xlabel(r'${\rm Radius} \ (R_\oplus)$')
        ax.set_ylabel(r'$\rm Density \ (kg \ m^{-3})$')


class analyticalModel():
    def __init__(self, n_layers = 3, res = 0, total_mass = 1., 
                 layer_mass_fractions = [.3, .3, .4]):
        self.n_layers = n_layers
        self.res = res
        self.total_mass = total_mass
        self.surface_pressure = 0e0
        self.layer_mass_fractions = layer_mass_fractions
        
    
    def plot_profile(self, ax = None, color = 'k'):
        if ax == None:
            fig, ax = plt.subplots()
        
        x = [0, self.layer_radii[0]]
        for i in range(len(self.layer_radii) - 1):
            x.append(self.layer_radii[i])
            x.append(self.layer_radii[i + 1])
            
        y = np.array([2*[d] for d in self.layer_densities]).flatten()
        
        ax.plot(x, y, color = color)
            
    
    def get_pressure_profile(self):
        """ Computes pressure profile as function of the radial distance
        from the center if the layer sizes and layer densities are knwon and
        assuming hydrostatic equilibrium.
        """
        
        # Compute pressure increments for each layer
        dRsq = [self.layer_radii[i + 1]**2 - self.layer_radii[i]**2 for i in range(self.n_layers - 1)]
        dRsq.append(0.)
        dRsq = np.array(dRsq)
        dRsq *= r_earth**2
        dP = np.array([d**2 for d in self.layer_densities])
        dP *= dRsq * 2 * np.pi * G / 3.
        #dP += self.surface_pressure
        self.layer_pressures = np.array([sum(dP[i:-1]) for i in range(self.n_layers)])
        
        self.layer_pressures = [0.]
        for i in range(self.n_layers-1,-1,-1):
            print ('i =', i, self.n_layers - i - 1)
            f = 0.
            for j in range(2, i):
                ff = self.layer_radii[j]**3 - self.layer_radii[j-1]**3
                ff *= self.layer_densities[j]
                f += ff
    
            a = f * (1. / self.layer_radii[i] - 1. / self.layer_radii[i-1])
            a += .5 * (self.layer_radii[i-1]**2 - self.layer_radii[i]**2) * self.layer_radii[i]
            a += self.layer_radii[i-1]**3 * (1. / self.layer_radii[i] - 1. / self.layer_radii[i]) * self.layer_radii[i]
            dP = G * self.layer_densities[i-1] * 4. * np.pi / 3. * a * r_earth**2
            print ('dP =', dP, self.layer_pressures[self.n_layers - i - 1])
            self.layer_pressures.append(dP + self.layer_pressures[self.n_layers - i - 1])
        
        self.layer_pressures = np.array(self.layer_pressures)
    
    
    def get_layer_densities(self):
        """Computes the density of each layer for fixed layer massas and
        layer sizes.
        """
        
        self.layer_densities = 3. / (4. * np.pi)
        self.layer_densities *= self.layer_masses
        #(self.layer_radii[i]**3 - self.layer_radii[i - 1]**3)
    
    
    def get_layer_masses(self):       
        self.layer_masses = [[self.total_mass * lmf / (2**self.res) 
                              for i in range(2**self.res)] 
                             for lmf in self.layer_mass_fractions]
        self.layer_masses = np.array(self.layer_masses).flatten()
        

    def draw_random(self, dens_ranges = np.array([[7.50e3, 8.0e3],
                        [3.0e3, 4.0e3],
                        [1e3, 1.2e3],
                        [9e2, 1e3]])):
        """Get density and mass for each layer by randomly sampling within
        predefined ranges.
        """
        self.layer_densities = np.empty([self.n_layers * 2**self.res])

        c = 0
        #Loop over layers
        for i in range(self.n_layers):
            #Loop over sublayers
            for j in range(2**self.res):
                p = random.random()
                
                low = dens_ranges[i][0]
                
                #Make sure that density is lower than in previous layer
                if c > 0:
                    high = min(dens_ranges[i][1], self.layer_densities[c-1])
                
                else:
                    high = dens_ranges[i][1]
                
                #Compute value for density within ranges
                rho = low + p * (high - low)
        
                self.layer_densities[c] = rho
                
                #Make sure density is lower than in previous sublayer
                if c > 0:
                    self.layer_densities[c] = min(rho, self.layer_densities[c-1])
                    
                c += 1
                

    def get_layer_radii(self):
        """Compute radius of each layer for given masses and densities
        """
        N = len(self.layer_masses)
        R1 = (3./(4.*np.pi)*self.layer_masses[0]*m_earth/self.layer_densities[0])**(1./3.)
        
        self.layer_radii = [R1 / r_earth]
        
        for i in range(N-1):
            M, rho = self.layer_masses[i+1]*m_earth, self.layer_densities[i+1]
            
            R = (3./(4.*np.pi)*M/rho + (r_earth*self.layer_radii[i])**3)**(1./3.)
            self.layer_radii.append(R / r_earth)
        
        self.total_radius = self.layer_radii[-1]


    def compute_moi(self):
        
        I = 2./5. * self.layer_masses[0]*m_earth*(self.layer_radii[0]*r_earth)**2
        
        N = len(self.layer_radii)
        for i in range(N - 1):
            M = self.layer_masses[i+1] * m_earth
            R = self.layer_radii[i+1] * r_earth
            R_before = self.layer_radii[i] * r_earth
            
            I += 2./5. * M *(R**5-R_before**5)/(R**3-R_before**3)
        
        self.MOI = I / (self.total_mass * m_earth * self.total_radius**2 * r_earth**2)
    

    def compute_ocean(self):
        if self.n_layers > 2:
            self.ocean_mass_fraction = sum(self.layer_masses[2:]) / sum(self.layer_masses)
            
        else:
            self.ocean_mass_fraction = 0.


class numericalModel():
    def __init__(self):
        pass


def plot_hydro_effect():
    
    omegas = np.array([1.0e-2])
    
    j2 = [J_2(o, radii, masses, MoI) for o in omegas]
    
    
    
    fig, ax = plt.subplots(figsize=(4, 3))
    
    ax.tick_params(right=True, top=True, direction='in', which='both',
                   )
    
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    
    ax.plot(masses, (radii[0]-radii[1])/radii[1]*100, color='k')
    ax.plot(masses, (MoI[0]-MoI[1])/MoI[1]*100, color='k', linestyle = '--')
    
    for j in j2:
        ax.plot(masses, (j[0]-j[1])/j[1]*100, linestyle = ':', color='k')

    ax.text(.05, .01*100 , r'$J_2$', rotation=60)    
    ax.text(.06, .005*100, 'MoI factor', rotation=20)
    ax.text(.06, .001*100 , 'radius', rotation=15)
    
    ax.set_xlabel(r'$M/M_\oplus$')
    ax.set_ylabel(r'$\rm max. \ rel. \ effect \ of \ hydration \ \rm [\%]$')
            
    fig.savefig(im_dir+'hydro_effect.pdf', format='pdf', bbox_inches='tight')
    plt.close(fig)    
    return j2


def plot_j2toc22():
    xticks = np.linspace(1, 9, len(j2))
    xlabels = ['Io', 'Europa', 'Ganymede', 'Callisto', 'Titan', 'Enceladus']
    
    sigma = sigma_j2toc22(j2, c22, dj2, dc22)

    fnts1 = 14
    fnts2 = 18
    
    fig, ax = plt.subplots(figsize=(8,6))
    
    for i in range(len(j2)):
        col = colors[references[i]]

        ax.errorbar(xticks[i], j2[i]/c22[i], yerr=sigma[i],
                 marker='o', linestyle = 'None',
                 color=col)
        
        ax.text(xticks[i]+.1, j2[i]/c22[i], str(references[i]+1), color = col)
        
    for i in range(len(references_names)):
        ax.text(.5, 3.15-i*.03, str(i+1)+': '+references_names[i],
                color = colors[i], fontsize=fnts1)
        
    bbox_props = dict(edgecolor='None', facecolor='white',
              linewidth = .5, pad=.1, alpha=1.)
    
    ax.text(1.5, 10/3, '10/3', va='center', color='grey',
            bbox=bbox_props, fontsize=fnts1)
    
    ax.plot([0, 10], [10/3, 10/3], color='grey')
    
    ax.set_xlim(0, 10)
    ax.tick_params(which='both')
    ax.set_ylabel(r'$J_2/C_{22}$', fontsize = fnts2, labelpad=20)
    
    plt.rcParams['ytick.labelsize'] = fnts1 
    
    plt.xticks(xticks, xlabels, rotation=60, fontsize=fnts1)
    
    fig.savefig(im_dir+'j2toc22.pdf', format='pdf', bbox_inches='tight')
    plt.close(fig)
    
    
def sigma_j2toc22(a, b, sa, sb):
    return np.sqrt(1/b**2*sa**2+(a/b**2)**2*sb**2)


def h_f(MoI):
    """Love number hf from the Radau equation (Fortes 2004, PhD thesis)
    """
    return 5./(((5/2-15/4*MoI))**2+1.)


def k_f(MoI):
    return h_f(MoI)-1.


def J_2(omega, a, M, kf):
    """J2 gravitational coefficient for a satellite in synchronous orbit
    assuming the body to be in hydrostatic equilibrium (Schubert et al. 1994)
    
    omega: rotation or orbital period (equal as synchronous orbit is assumed)
    a: radius of satellite
    M: mass of satellite
    kf: Love number
    
    """
    
    return 5./6.*kf*(omega**2*a**3/(G*M))


def C_22(omega, a, M, kf):
    """C22 gravitational coefficient for a satellite in synchronous orbit
    assuming the body to be in hydrostatic equilibrium (Schubert et al. 1994)
    """
    
    return 1./4.*kf*(omega**2*a**3/(G*M))


def rho(r):
    return 5500.


def P2(N):
    
    res = 0.
    
    dtheta = np.pi/(N-1)
    theta = 0.
    
    for i in range(N):
        res += .5*(3.*np.cos(theta)**2-1.) * dtheta
        
        theta += dtheta
        
    return res
        
'''
def J2(R, M, N=10000, Q=10, eps=.01):
    
    res = 0.
    
    theta = 0.
    
    dtheta = np.pi/2/(N-1)
    dr = R/(Q-1)
    
    for i in range(N):
        r = 0.
        #print ('theta =', theta/np.pi) 
        for j in range(Q):
            
            rbar = r/np.sqrt(1.-eps**2*np.cos(theta)**2)
            P2 = .5*(3.*np.cos(theta)**2 - 1.)

            #print ('r=', r, 'rbar =', rbar)
            #res += rho(r) * .5*(3*np.cos(theta)**2-1.) * np.sin(theta) * \
             #   dtheta*dr/R*(r/R)**4
            
            res += P2 * np.sin(theta) * dtheta * dr/R * (r/R)**4
            
            #print ((r/R)**4)
            
            r += dr
                       
        theta += dtheta
        

        
    res *= -4.*np.pi*R**3/M
    
    return res
'''


def compute_diff_area(x, y1, y2):
    """Compute 
    """
    N = len(x)
    
    A = 0.
    
    for i in range(N-1):
        
        delta_x = abs(x[i] - x[i+1])

        delta_y1 = abs(y1[i] - y2[i])
        delta_y2 = abs(y1[i+1] - y2[i+1])
        
        delta_y = (delta_y1+delta_y2)*.5
        
        delta_A = delta_x*delta_y
        
        A += delta_A
        
    return A
        

def draw_densities(N_layer, N_res):
    dens = np.empty([N_layer*N_res])

    c = 0
    for i in range(N_layer):
        for j in range(N_res):
            p = random.random()
            
            low = dens_ranges[i][0]
            
            if c > 0:
                high = min(dens_ranges[i][1], dens[c-1])
            
            else:
                high = dens_ranges[i][1]
            
            rho = low + p*(high-low)
    
            dens[c] = rho
            
            if c > 0:
                dens[c] = min(rho, dens[c-1])
                
            c += 1

    return dens


def draw_pressure_enhancement(N_layer, N_res, omf):
    a = np.empty([N_layer*N_res])
    
    for i in range(N_layer*N_res):
        sign = (random.random()-.5)/.5
        
        #Ocean
        #The density of water increases monotonously with the ocean mass fraction
        if i >= N_layer*N_res-N_res:
            a[i] = random.random()/5.
        
        #Core or Mantle
        else:
            a[i] = random.random()/100.*sign
            
    return a


def compute_mean_densities(N_layer, N_res, omf):
    
    dens = draw_densities(N_layer, N_res)
    
    a = draw_pressure_enhancement(N_layer, N_res, omf)
    
    result = np.empty([len(omf), N_layer*N_res])
    
    for i in range(len(omf)):
        for j in range(N_layer*N_res):
            result[i][j] = dens[j]#*omf[i]**a[j]#(1.+a[j]*omf[i])
        
    return result, a
        

def Ri(layer_masses, layer_densities):
    print ('layer masses', layer_masses)
    N = len(layer_masses)
    R1 = (3./(4.*np.pi)*layer_masses[0]/layer_densities[0])**(1./3.)
    
    layer_radii = [R1]
    
    for i in range(N-1):
        M, rho = layer_masses[i+1], layer_densities[i+1]
        
        R = (3./(4.*np.pi)*M/rho + layer_radii[i]**3)**(1./3.)
        layer_radii.append(R)
        
    return layer_radii


def MOI_factor(layer_masses, layer_radii):
    
    M_tot = sum(layer_masses)
    R_tot = layer_radii[-1]
    I = 2./5.*layer_masses[0]/M_tot*(layer_radii[0]/R_tot)**2
    
    N = len(layer_masses)
    
    for i in range(N-1):
        M = layer_masses[i+1]
        R = layer_radii[i+1]
        R_before = layer_radii[i]
        
        I += 2./5. * M *(R**5-R_before**5)/(R**3-R_before**3)
    
    return I/(M_tot*R_tot**2)


def initiate_layer_masses(M, fractions):
    
    masses = []
    
    for i in range(len(fractions)):
        masses.append(M*fractions[i])
        #masses.append(.5*M*(1.-omf)*fractions[i])
            
    return masses


def rms_dev(a, b):
    
    pass


def initiate_colors(N_layer, N_res):
    colors = np.empty([N_layer*N_res, 3])
    
    c = 0
    for i in range(N_layer):
        dcolor =  (layer_colors_finals[i] - layer_colors[i])/(N_res-1)
        for j in range(N_res):

            colors[c] = layer_colors[i] + dcolor*j
            
            c += 1
        
    return colors
    

def run():
    N_res = 2
    
    N = 13
    N_sample = 1
    N_layer = 3
    cmf = .33 #m_core/(m_core+m_mantle)
    
    lc = initiate_colors(N_layer, N_res)
    
    ocean_mass_fractions = np.linspace(.0001, .6, N)
    
    MOIs = np.empty([N])
    
    sample_curves = np.empty([N_sample, N])
    sample_densities = np.empty([N_sample, N, N_layer*N_res])
    sample_rms1 = np.empty([N_sample])
    sample_rms2 = np.empty([N_sample])
    sample_radii = np.empty([N_sample, N])
    sample_J2 = np.empty([N_sample, N])
    
    sample_layer_masses = np.empty([N, N_layer*N_res])
    
    
    for j in range(len(ocean_mass_fractions)):
        omf = ocean_mass_fractions[j]
        #Core
        for i in range(N_res):
            sample_layer_masses[j][i] = cmf/N_res*(1.-ocean_mass_fractions[j])
            
        #Mantle
        for i in range(N_res):
            sample_layer_masses[j][N_res+i] = (1.-cmf)/N_res*(1.-ocean_mass_fractions[j])
        
        #Ocean
        for i in range(N_res):
            sample_layer_masses[j][2*N_res+i] = ocean_mass_fractions[j]/N_res
    
    
    
    for i in range(N_sample):
    
        ld, a = compute_mean_densities(N_layer, N_res, ocean_mass_fractions)
        
        
        for o in range(N):
            for j in range(N_layer*N_res):
                sample_densities[i][o][j] = ld[o][j] 
        
        
        for o in range(N):
            omf = ocean_mass_fractions[o]
    
            lm = sample_layer_masses[o]*M_earth*.008
            
            #print ('masses =', lm)
            
            lr = Ri(lm, ld[o])
            #print ('radii =', [r/max(lr) for r in lr])
        
            MOI = MOI_factor(lm, lr)
            MOIs[o] = MOI
            
            kf = k_f(MOI)
            J2 = J_2(omega, lr[-1], .008*M_earth, kf)        
            
            sample_curves[i][o] = MOI
            sample_radii[i][o] = lr[-1]
            sample_J2[i][o] = J2
            
        diff_area_MOI = compute_diff_area(data[0][2], data[0][0], MOIs)
        diff_area_radius = compute_diff_area(data[0][2], data[0][1], sample_radii[i]/R_earth)
        
        rms1 = diff_area_radius#sum(abs(data[1][0]-MOI)/data[1][0])#**2/data[1][0]**2)
        rms2 = diff_area_MOI#np.sqrt(sum(abs(data[1][0]-MOI)**2/data[1][0]**2) + \
                       #sum(abs(data[1][1]-sample_radii[:][-1]/R_earth)**2/data[1][1])**2)
        
        #print ('rms2 =', rms2)
        
        sample_rms1[i] = rms1
        sample_rms2[i] = rms2
        
        
    
    ind1 = np.argmin(sample_rms1)
    ind2 = np.argmin(sample_rms2)
    
    print ('min rms dev =', sample_rms1[ind2])
    print ('dens =', sample_densities[ind2])
    print ('a =', a)
    
    
    fig1, ax1 = plt.subplots(1,3)
    
    
    ax1[0].set_xlim(0., .6)
    ax1[1].set_xlim(0., .6)
    
    ax1[0].text(.1, .378, 'Io')
    ax1[1].text(.1, .286, 'Io')
    ax1[0].plot([0., 1.], [.378, .378])
    ax1[1].plot([0., 1.], [.286, .286])
    
    ax1[0].text(.1, .341, 'Titan')
    ax1[0].plot([0., 1.], [.3414, .3414])
    ax1[1].plot([0., 1.], [.404, .404])
    
    ax1[0].text(.1, .331, 'Enceladus')
    ax1[0].plot([0., 1.], [.3305, .3305])
    ax1[1].plot([0., 1.], [.0395, .0395])
    
    ax1[0].text(.1, .346, 'Europa')
    ax1[0].plot([0., 1.], [.346, .346])
    ax1[1].plot([0., 1.], [.245, .245])
    
    ax1[0].text(.1, .3929, 'Moon')
    ax1[0].plot([0., 1.], [.3929, .3929])
    ax1[1].plot([0., 1.], [.2727, .2727])
    
    ax1[0].text(.1, .3115, 'Ganymede')
    ax1[0].plot([0., 1.], [.3115, .3115])
    ax1[1].plot([0., 1.], [.413, .413])
    
    ax1[0].text(.1, .3549, 'Callisto')
    ax1[0].plot([0., 1.], [.3559, .3549])
    ax1[1].plot([0., 1.], [.378, .378])
    
    
    ax1[0].plot(ocean_mass_fractions, sample_curves[ind2], label='analytical', 
                color=[0., 1., 0.])
    ax1[1].plot(ocean_mass_fractions, sample_radii[ind2]/R_earth, label='analytical', 
                color=[0., 1., 0.])
    ax1[2].plot(ocean_mass_fractions, sample_J2[ind2]*1.0e6, label='analytical', 
                color=[0,1,0])
    
    
    linestyles = ['-', '--', ':']
    
    for i in range(len(data)):
        ax1[0].plot(data[i][2], data[i][0], marker='s', color='k', label='numerical',
           linestyle = linestyles[i])
        ax1[1].plot(data[i][2], data[i][1], marker='s', color='k', label='numerical',
           linestyle = linestyles[i])
        
    
    ax1[0].legend()
    
    for ax in ax1:
        ax.tick_params(which='both', top=True, right=True, direction='in')
        ax.set_xlabel('Ocean mass fraction')
        
    ax1[0].set_ylabel('MOI factor')
    ax1[1].set_ylabel(r'$R/R_\oplus$')
    ax1[2].set_ylabel(r'$J_2 \ [10^{-6}]$')
    
    fig, ax = plt.subplots()
    
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    
    ax.set_aspect('equal')
    ax.set_xlabel('normalized radius')
    ax.set_ylabel('normalized radius')
    
    cols = ('aqua', 'blue', 'mediumblue', 'midnightblue', 
            'brown', 'indianred', 'lightcoral', 'goldenrod', 'gold', 'white')
    
    print ('')
    for i in range(len(lr)):
        circle = plt.Circle((0, 0), lr[-i-1]/max(lr), color=lc[i])
        
        ax.add_artist(circle)
        
