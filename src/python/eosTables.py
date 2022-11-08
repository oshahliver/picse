#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:11:24 2019

@author: oshah
"""

from tqdm import tqdm
from PIMPphysicalparams import material_list
import numpy as np
import os
from matplotlib import pyplot as plt
import Feistel2006EOS
import French2015EOS
import Wagner2002EOS
import Mazevet2018EOS
import eos
import sys
import time
import math
import matplotlib as mpl
import functionTools as ftool
import astropy.table
import warnings
from astropy.io import ascii
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PIMPphysicalparams import T_triple, T_critical, P_triple, P_critical, NA,\
                                mH, mO, kB, material_plot_list, mSi, mMg, mFe, \
                                mH2O, mPer, mPerov, mOl, mWs

from PIMPrunparams import param_labels, param_lims, water_phases

import phaseCheck
import Wagner2002EOS
import eoswater
import logTrans


def s2min(val, r=3):
    m=int(round(val/60-.5, 0))
    s=round(val-m* 60, r)
    return m, s

def print_time(t0, t, **kwargs):
        print ('\nelapsed time: '+str(s2min(t-t0, **kwargs)[0])+ 'min '+\
                           str(s2min(t-t0, **kwargs)[1])+ 's')

def fct(x=None, y=None, param=2, table=False, **kwargs):
    param_scaling = [0, 0, 0, 0, 3, 6]
    #compute density
    t0=time.time()
    result = eos.Compute(what='all', T=x, P=y, table=table, **kwargs)[param_scaling[param]]
    actual_time = time.time()-t0
    
    return result, actual_time

def deriv(x=None, y=None, whicharg=None, **kwargs):
    if whicharg=='x':
        d=ftool.deriv(f=fct, whicharg='x', x0=x, y=y, **kwargs)
    
    elif whicharg=='y':
        d=ftool.deriv(f=fct, whicharg='y', x0=y, x=x, **kwargs)
    return d
    
def derivderiv(x=None, y=None, **kwargs):
    """Compute partial derivative with respect to x and y
    """
    #first differentiate with respect to x
    def fx(x=x, y=y, arg='x'):
        d=deriv(x=x, y=y, whicharg=arg, **kwargs)
        return d
    
    #then differentiate this with respect to y
    #here the **kwargs mustn't be passed in order to avoid argument conflicts
    d=ftool.deriv(f=fx, whicharg='y', arg='x', x=x, x0=y)
    return d

class Table():
    def __init__(self, x_start=None, y_start=None, x_range=None, y_range=None,
                 alpha_x = None, alpha_y = None):
        self.material=None
        self.ll=None
        self.Fe_number = 0
        self.saturation = False
        self.alpha_x = alpha_x
        self.alpha_y = alpha_y
        self.x_range= x_range
        self.y_range= y_range
        self.x_start= x_start
        self.y_start= y_start
        self.dim_decade_x = None
        self.dim_decade_y = None
        self.values=None
        self.derivatives_x=None
        self.derivatives_y=None
        self.derivatives_xy=None
        self.b=None
        self.matrix=None
        self.pairs=[]
        self.x_axis=[]
        self.y_axis=[]
        self.dimension=[len(self.x_axis), len(self.y_axis)]
        self.grid_size=0
        self.data_table=None #table containing data strings
        self.data=None
        self.phase=None
        self.data_contents=['x', 'y', 'function value', 'x-derivative',
                            'y-derivative', 'xy-derivative']
    
    
    def test_plot(self, param=2, def_zero=1.0e-10):
        
        #compute x and y axis range
        xmin=self.x_axis[0]
        xmax=self.x_axis[-3]
        ymin=self.y_axis[0]
        ymax=self.y_axis[-3]
        
        #generate log scale for axis
        nx = int(np.log10(xmax)-np.log10(xmin))+1
        ny = int(np.log10(ymax)-np.log10(ymin))+1
        
        major_ticks = [np.linspace(0, 1, nx), np.linspace(0, 1, ny)]
        major_tick_labels = [np.linspace(int(np.log10(xmin)), 
                                        int(np.log10(xmax)), nx),
                            np.linspace(int(np.log10(ymin)), 
                                        int(np.log10(ymax)), ny)]
        
        a,b,c,d = logTrans.transform(data=self.data, param=param, 
                                     def_zero=def_zero)
        
        logTrans.plot(a,b,c,d, major_ticks=major_ticks,
                      major_tick_labels=major_tick_labels,
                      cbar_label = param_labels[param],
                      material=material_plot_list[self.material])
    
    
    def Plot(self, param=2, dim=2, log=False, **kwargs):
        
        #gather data points in dim dimensional array
        #ignore the last two points as they are extra points
        dat_plot=np.empty([self.dimension[0]-3, self.dimension[1]-3])
        #leave out ghost cells
        for i in range(len(self.x_axis)-3):
            for j in range(len(self.y_axis)-3):
                try:
                    val = float(self.data[i][j][param])
                
                except TypeError:
                    val = None
                    
                except ValueError:
                    val = None
                    
                #plot density
                if param == 2:
                    if log:
                        try:
                            dat_plot[i][j]=np.log10(val)
                        
                        except TypeError:
                            dat_plot[i][j]=None
                            
                    else:
                        dat_plot[i][j] = val
                
                elif param==3 or param==4 or param==5 or param==6 or param==7:
                    if log:
                        try:
                            val = np.log10(val)#max(np.log10(val), -11)
                        
                        except TypeError:
                            pass
                    
                    dat_plot[i][j] = val
                            
        #compute x and y axis range
        xmin=self.x_axis[0]
        xmax=self.x_axis[-3]
        ymin=self.y_axis[0]
        ymax=self.y_axis[-3]
        
        plot_title = material_plot_list[self.material]
        
        if self.material == 0:
            try:
                plot_title = plot_title + ' ' + str(water_phases[self.phase])
            
            except AttributeError:
                pass
            
            except TypeError:
                pass
            
        #generate log scale for axis
        nx = int(np.log10(xmax)-np.log10(xmin))+1
        ny = int(np.log10(ymax)-np.log10(ymin))+1
        print (nx, ny)
        
        major_ticks_x = np.linspace(0, 1, nx)
        major_ticks_y = np.linspace(0, 1, ny)

        major_tick_labels_x = np.linspace(int(np.log10(xmin)), 
                                        int(np.log10(xmax)), nx)
        
        major_tick_labels_y = np.linspace(int(np.log10(ymin)), 
                                        int(np.log10(ymax)), ny) 
        
        min_plot = np.nanmin(dat_plot)
        max_plot = np.nanmax(dat_plot)
        
        #major_tick_labels_x = 10**major_tick_labels_x
        #major_tick_labels_y = 10**major_tick_labels_y
            
        bounds=[0,1,2,3,4,5,6]
        if param == 3:
            cmap ='rainbow'
            norm = norm = mpl.colors.Normalize(vmin=min_plot, vmax=max_plot)
            '''
            cmap=mpl.colors.ListedColormap([(.6, .4, .6),
                                            (.2, .8, .8),
                                            (.2, .8, .2),
                                            (.5, .5, 0.),
                                            (.8, .4, 0.),
                                            (1., .2,.2)])
            norm=mpl.colors.BoundaryNorm(bounds, cmap.N)
            '''
        else:
            cmap = 'rainbow'
            norm = norm = mpl.colors.Normalize(vmin=min_plot, vmax=max_plot)
            
        fig, ax = plt.subplots()
        ax.set_xlim(left=0., right=1.)
        ax.set_ylim(bottom=0., top=1.)
        
        ax.set_xticks(major_ticks_x)
        ax.set_yticks(major_ticks_y)
        
        ax.set_xticklabels(major_tick_labels_x)
        ax.set_yticklabels(major_tick_labels_y)
        
        ax.set_xlabel(r'$log(T) \ [K]$')
        ax.set_ylabel(r'$log(P) \ [Pa]$')
        
        im=ax.imshow(dat_plot.T, cmap=cmap, extent=[0.,1.,0.,1.],
        norm=norm, origin='lower', vmin= param_lims[param][0],
        vmax =param_lims[param][1])
        
        divider=make_axes_locatable(ax)
        cax=divider.append_axes('right', size='5%', pad=.5)
        if param == 33:
            cbar=fig.colorbar(im, cax=cax, cmap=cmap, norm=norm, boundaries=bounds,
                          ticks=[.5, 1.5, 2.5, 3.5, 4.5, 5.5])
            cbar.ax.set_yticklabels(['solid', 'liquid', 'vapor', 'gas', 
                                     'supercritical', 'high pressure'])        
        
        else:
             cbar=fig.colorbar(im, cax=cax, cmap=cmap, norm=norm)
        
        cbar.ax.tick_params(size=0.)
        cbar.ax.set_ylabel(param_labels[param])
        
        ax.set_title(plot_title)
        
        plot_color = (.2,.3,.5)
        phase_color = (.2,.3,.5)
        
        font1 = {'family': 'sans',
                'color':  'grey',
                'weight': 'bold',
                'size': 10,
                }
        
        font2= {'family': 'sans',
                'color':  'k',
                'weight': 'normal',
                'size': 14,
                }
        
        font3= {'family': 'sans',
                'color':  'k',
                'weight': 'normal',
                'size': 10,
                }
        '''
        exp_T_critical = int(np.log10(T_critical))
        exp_P_critical = int(np.log10(P_critical))
        exp_T_triple = int(np.log10(T_triple))
        exp_P_triple = int(np.log10(P_triple))
        
        ax.scatter((exp_T_critical+T_critical/10**(exp_T_critical+1))/(nx-1),
                   (exp_P_critical+P_critical/10**(exp_P_critical+1))/(ny-1), 
                   color='black', s=80)
        
        ax.scatter((exp_T_triple+T_triple/10**(exp_T_triple+1))/(nx-1),
                   (exp_P_triple+P_triple/10**(exp_P_triple+1))/(ny-1),
                   color='white', s=80)
        
        #ax.scatter(np.log10(1)/(nx-1), np.log10(1)/(ny-1))
        
        ax.plot([2/(nx-1), 2/(nx-1)], [0, 1], 
                 color='white', 
                 linestyle='--')
        
        ax.plot([0, 5/(nx-1)], [9/(ny-1), 9/(ny-1)], color='white', 
                 linestyle='--')
        
        '''
        
        
    def add_data(self, param=2, x_start=0, y_start=0, x_range=1, y_range=1, 
                 **kwargs):
        print ('computing grid values...\n')

        try:
            ll=kwargs['ll']
        
        except KeyError:
            print ('WARNING: no material given! Setting ll=0.')
            ll=0
        
        t_start=time.time()
        dt_list=[]
        time_list=[]
        total_length=50
        count=0
        IPS=0 #iterations per second
        
        percentage=int(count/self.grid_size*100)
        space=(total_length-int(percentage/100*total_length))*' '
        bar=int(percentage/100*total_length)*'=O'
        sys.stdout.write('\rProgress: ['+bar+space+']'+str(percentage)+\
                         ' % , ETA=\u221e min')
        sys.stdout.flush()

        for i in range(len(self.data)):
            #iterate over columns (y-axis)
            for j in range(len(self.data[i])):
                t0=time.time()
                count+=1
                #gather values
                dat = self.data[i][j]
                new = self.data[i][j][param]
                x, y, val, ph, dPdrho, dPdT = dat
                
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        if param == 2 and y >= 1.0e8:
                            new = fct(x=x, y=y, ll=0)
 
                        elif param == 3:
                            new = phaseCheck.getPhase(T=x, P=y, double_check=True,
                                                   d=val, ll=ll)
                        
                        elif param == 4:
                            new = eos.dPdrho(P=y, d=val, T=x, ll=ll)
                        
                        elif param == 5:
                            new = eos.dTdP(P=y, d=val, T=x, ll=ll)
                        
                        else:
                            new = dat[param]
                            
                except OverflowError:
                    pass
                    
                self.data[i][j][param] = new
                
                if not len(dt_list)==0:
                    dt_max=max(dt_list)
                    dt_mean=sum(dt_list)/len(dt_list)
                    eta_mean=(self.grid_size-count+1)*dt_mean
                    eta_max=(self.grid_size-count+1)*dt_max
                    
                    percentage=round(count/self.grid_size*100, 1)
                    space=(total_length-int(percentage/100*total_length))*' '
                    bar=(int(percentage/100*total_length)-1)*'='
                    
                    bar+='>'
                    time_max=s2min(eta_max)
                    eta_mean = (eta_mean+eta_max)/2
                    time_mean=s2min(eta_mean)
                    IPS=1/(time.time()-t0)
                    sys.stdout.write('\rProgress: ['+bar+space+'] '+\
                                     str(percentage)+'%, IPS='+str(int(IPS))+\
                                     ', ETA_mean='+ \
                                     str(int(time_mean[0]))+'min '+\
                                     str(int(time_mean[1]))+'s'+\
                                     ' / ETA_max='+str(int(time_max[0]))+\
                                     'min '+str(int(time_max[1]))+'s'+\
                                     100*' ')
                    
                    sys.stdout.flush()
                    time_list.append(time.time()-t_start)

                dt=time.time()-t0
                dt_list.append(dt)
    
    
    def write_fortran_grid(self, file_name, loc='cwd'):
        """Writes text file containing a plain PT gridin the format that is 
        used by the fortran routine
        """
        file = open(file_name, 'w')
        
        for i in range(len(self.x_axis)):
            for j in range(len(self.y_axis)):
                file.write('\n'+str(self.x_axis[i])+','+str(self.y_axis[j]))
        
    
    def write_fortran(self, file_name, loc='cwd'):
        file = open(file_name, 'w')
        file.write('material='+str(self.material))
        file.write('\nnx='+str(self.dimension[0]))
        file.write('\nny='+str(self.dimension[1]))
        file.write('\nalpha_x='+str(self.alpha_x))
        file.write('\nalpha_y='+str(self.alpha_y))
        file.write('\nx_start='+str(self.x_start))
        file.write('\ny_start='+str(self.y_start))
        file.write('\nx_range='+str(self.x_range))
        file.write('\ny_range='+str(self.y_range))
        #file.write('\Fe_number='+str(self.Fe_number))
        
        for i in range(len(self.x_axis)):
            for j in range(len(self.y_axis)):
                file.write('\n'+str(self.data[i][j][0])+','+\
                           str(self.data[i][j][1])+','+\
                           str(self.data[i][j][2])+','+\
                           str(self.data[i][j][3])+','+\
                           str(self.data[i][j][4])+','+\
                           str(self.data[i][j][5])+','+\
                           str(self.data[i][j][6])+','+\
                           str(self.data[i][j][7])+','+\
                           str(self.data[i][j][8])+','+\
                           str(self.data[i][j][9]))
        
        
    def write(self, file_name, loc='cwd'):
        """Saves eos table which has been generated with the <generate>
        methode to ascii file
        """
        t0=time.time()
        data=[]
        
        print ('gathering grid values...')
        for j in range(len(self.y_axis)):
            data.append([])
            #note that rows correspond to x-axis and columns to y-axis
            for i in range(len(self.x_axis)):
                dat=self.data[i][j]
                string=str(dat[0])
                for d in range(len(dat)-1):
                    string+=','+str(dat[d+1])
                
                data[j].append(string)
        
        if loc=='cwd':
            dir_out=os.getcwd()+'/'
            
        else:
            dir_out=loc+'/'
        
        print ('generating table...')
        self.data_table=astropy.table.Table(data, meta={
                                                'material':self.material,
                                                'saturation':self.saturation,
                                                'Fe_number':self.Fe_number,
                                                'phase':self.phase,
                                                'alpha_x': self.alpha_x,
                                                'alpha_y': self.alpha_y,
                                                'x_range': self.x_range,
                                                'y_range': self.y_range,
                                                'dimension':self.dimension,
                                                'grid_size':self.grid_size,
                                                'x_start':self.x_start,
                                                'y_start':self.y_start})
    
        print ('writing eos table to:', dir_out+file_name)
        ascii.write(self.data_table, dir_out+file_name, overwrite=True,
                   format='ecsv')
        print ('table successfully written to file')
        t=time.time()
        print_time(t0, t)
    
    
    def load(self, file_name, loc='cwd'):
        """Loads existing table and makes data accessable for the current
        python session. This allows quicker data handling as the table content
        need not be read out from the file for every interpolation cycle.
        """
        t0=time.time()
        self.x_axis=[]
        self.y_axis=[]
        self.values=[]
        
        if loc=='cwd':
            dir_in=os.getcwd()+'/'
            
        else:
            dir_in=loc+'/'
        
        #read ascii file for eos table
        self.data_table=ascii.read(dir_in+file_name, format='ecsv')
        
        #extract metadata
        self.x_range=self.data_table.meta['x_range']
        self.y_range=self.data_table.meta['y_range']
        self.alpha_x=self.data_table.meta['alpha_x']
        self.alpha_y=self.data_table.meta['alpha_y']
        self.material=self.data_table.meta['material']
        self.dimension=self.data_table.meta['dimension']
        self.grid_size=self.data_table.meta['grid_size']
        self.x_start=self.data_table.meta['x_start']
        self.y_start=self.data_table.meta['y_start']
        self.saturation=self.data_table.meta['saturation']
        self.Fe_number=self.data_table.meta['Fe_number']
        
        print ('extracting grid values...')
        self.data=[]
        #iterate over rows (x-axis)
        for i in range(len(self.data_table)):
            self.data.append([])
            #iterate over columns (y-axis)
            for j in range(len(self.data_table[i])):
                #gather values
                dat = self.data_table[i][j]
                dat=dat.split(',')
                for d in range(len(dat)):
                    try:
                        dat[d]=float(dat[d])
                        
                    except ValueError:
                        pass
                    
                self.data[i].append(dat)
                
        print ('-> grid data successfully extracted')
        print ('generating grid axes...')
        self.x_axis = [dat[0][0] for dat in self.data]
        print ('-> x axis successfully generated')
        self.y_axis = [dat[1] for dat in self.data[0]]
        print ('-> y axis successfully generated')
        t=time.time()
        print_time(t0=t0, t=t)
    
    
    def test(self, T=None, d=None, P=None):
        return Wagner2002EOS.cp_spec(T=T, d=d, P=P)
       
        
    def construct_grid(self):
        """Construct PT grid according to the specification given in the
        Table instance initiation
        """
        print ('initializing x-y-grid...')
        print ('processing x-axis...')
        
        #construct x-axis
        for i in range(self.x_range):
            #for the last decade include endpoint to cover entire decade            
            if i == self.x_range-1:
                arr = np.linspace(10**(i+self.x_start), 10**(i+1+self.x_start), 
                                  9*10**self.alpha_x+1, endpoint=True)
            
            #otherwise don't
            else:
                arr = np.linspace(10**(i+self.x_start), 10**(i+1+self.x_start), 
                                  9*10**self.alpha_x, endpoint=False)     
                
            #append decade to temperature axis
            for a in arr:
                self.x_axis.append(a)

        #Include two additional grid points to allow at least 2nd order
        #interpolation up to desired upper limit of the temperature axis        
        self.x_axis.append(10**(self.x_range+self.x_start)+10**(self.x_range+\
                           self.x_start)*\
                           10**(-self.alpha_x))
        
        self.x_axis.append(10**(self.x_range+self.x_start)+10**(self.x_range+\
                           self.x_start)*\
                           2*10**(-self.alpha_x))
        
        #construct y-axis
        print ('processing y-axis...')
        for j in range(self.y_range):
            
            #for the last decade include endpoint to cover entire decade
            if j == self.y_range-1:
                arr = np.linspace(10**(j+self.y_start), 10**(j+1+self.y_start), 
                                  9*10**self.alpha_y+1, endpoint=True)

            else:
                arr = np.linspace(10**(j+self.y_start), 10**(j+1+self.y_start), 
                                  9*10**self.alpha_y, endpoint=False) 
            
            #append decade to pressure axis
            for a in arr:
                self.y_axis.append(a)

        #Include two additional grid point to allow at least 2nd order
        #interpolation up to desired upper limit of the pressure axis        
        self.y_axis.append(10**(self.y_range+self.y_start)+10**(self.y_range+\
                           self.y_start)*\
                           10**(-self.alpha_y))

        self.y_axis.append(10**(self.y_range+self.y_start)+10**(self.y_range+\
                           self.y_start)*\
                           2*10**(-self.alpha_y))
                
        #update dimension of the grid
        self.dimension=[len(self.x_axis), len(self.y_axis)]
        
        #grid size
        self.grid_size=self.dimension[0]*self.dimension[1]  
    
    
    def Generate(self, alpha_x=0, alpha_y=0, x_range=4, y_range=12, 
                 x_start=0, y_start=0, ll=None, n_params=10, phase=None,
                 saturation=False, SiMg=None, Fe_number=0.0,
                 pres_thres = False, new_data=True):
        """Generates a logarithmic 2d x-y-grid. The spacing between two grid
        points is set by the alpha parameter. Conventions are:
            
        alpha (type=int, range=[0, infinity]): sets the order of grid
            refinement. Each decade on the p and T axis is divided into
            10**(alpha+1) grid points
            
        T_min, T_max (type=int): defines temperature range for the table as
            [10**T_min, 10**T_max]
            
        xstart, ystart (type=int): defines the starting exponent for the x and y
            axis. Axes are then generated betweend 10**start and 10**(start+range)
            for x and y respectively and using the corresponding alpha value
        
        
        Convention for the eos table data is:
            T, P, rho, phase, dP/drho, dP/dT
            
        where derivatives are carried out at constant entropy. The derivative
        drho/dT can be computed as (dP/dT)/(dP/drho)
        """
        t_start=time.time()
        if ll == None:
            print ('WARNING: No material specified!')
        
        #set refinement level for x and y axis
        self.alpha_x=alpha_x
        self.alpha_y=alpha_y
        
        self.x_range=x_range
        self.y_range=y_range
        
        self.x_start=x_start
        self.y_start=y_start
        
        #first set up x-y-grid
        self.x_axis=[]
        self.y_axis=[]
        
        self.material=ll
        self.phase = phase
        self.Fe_number = Fe_number
        self.saturation = saturation
        
        print ('initializing x-y-grid...')
        print ('processing x-axis...')
        
        #construct x-axis
        for i in range(x_range):
            #for the last decade include endpoint to cover entire decade            
            if i == x_range-1:
                arr = np.linspace(10**(i+x_start), 10**(i+1+x_start), 
                                  9*10**self.alpha_x+1, endpoint=True)
            
            #otherwise don't
            else:
                arr = np.linspace(10**(i+x_start), 10**(i+1+x_start), 
                                  9*10**self.alpha_x, endpoint=False)     
                
            #append decade to temperature axis
            for a in arr:
                self.x_axis.append(a)

        #Include two additional grid points to allow at least 2nd order
        #interpolation up to desired upper limit of the temperature axis        
        self.x_axis.append(10**(x_range+x_start)+10**(x_range+x_start)*\
                           10**(-alpha_x))
        
        self.x_axis.append(10**(x_range+x_start)+10**(x_range+x_start)*\
                           2*10**(-alpha_x))
        
        #construct y-axis
        print ('processing y-axis...')
        for j in range(y_range):
            
            #for the last decade include endpoint to cover entire decade
            if j == y_range-1:
                arr = np.linspace(10**(j+y_start), 10**(j+1+y_start), 
                                  9*10**self.alpha_y+1, endpoint=True)

            else:
                arr = np.linspace(10**(j+y_start), 10**(j+1+y_start), 
                                  9*10**self.alpha_y, endpoint=False) 
            
            #append decade to pressure axis
            for a in arr:
                self.y_axis.append(a)

        #Include two additional grid point to allow at least 2nd order
        #interpolation up to desired upper limit of the pressure axis        
        self.y_axis.append(10**(y_range+y_start)+10**(y_range+y_start)*\
                           10**(-self.alpha_y))

        self.y_axis.append(10**(y_range+y_start)+10**(y_range+y_start)*\
                           2*10**(-self.alpha_y))
                
        #update dimension of the grid
        self.dimension=[len(self.x_axis), len(self.y_axis)]
        
        #grid size
        self.grid_size=self.dimension[0]*self.dimension[1]        
        
        print ('computing grid values...\n')

        nx = len(self.x_axis)
        ny = len(self.y_axis)
        
        if new_data:
            self.data = np.zeros([nx, ny, n_params])
       
        dt_list=[]
        time_list=[]
        total_length=50
        count=0
        IPS=0 #iterations per second
        
        percentage=int(count/self.grid_size*100)
        space=(total_length-int(percentage/100*total_length))*' '
        bar=int(percentage/100*total_length)*'=O'
        sys.stdout.write('\rProgress: ['+bar+space+']'+str(percentage)+\
                         ' % , ETA=\u221e min')
        
        sys.stdout.flush()
        
        #if specific phase is given, enforce this phase for all P and T
        if not phase == None:
            update_phase = False
            
        else:
            update_phase = True
        
        for i in range(len(self.x_axis)):
            #print ('Progress {:2.1%}'.format(count/self.grid_size), end='\r')
            for j in range(len(self.y_axis)):
                t0=time.time()
                count+=1
                x=self.x_axis[i]
                y=self.y_axis[j]
             
                '''
                val=Mazevet2018EOS.density(t=x, p=y)#fct(x=x, y=y, param=2, table=False, **kwargs)[0]
                #ph=phaseCheck.getPhase(P=y, T=x, double_check=True, d=val, ll=ll)
                dTdP_S = Mazevet2018EOS.dTdP_S(t=x, p=y, d=val)#eos.dPdrho(P=y, d=val, T=x, ll=ll)
                #print ('dTdP_S=', dTdP_S)
                alpha = Mazevet2018EOS.alpha_th_p(T=x, P=y, d=val)
                cp = Mazevet2018EOS.cp_spec(T=x, P=y, d=val)
                s = 0.0
                u = Mazevet2018EOS.u_spec(t=x, p=y, d=val)
                '''
                
                '''
                val = Feistel2006EOS.density(P=y, T=x)
                dTdP_S = Feistel2006EOS.dTdP_s(T=x, P=y, d=val)
                #sif x==50. and y==1.:
                    #print ('a=', a)
                    #self.test()
                    #sys.exit()
                alpha = Feistel2006EOS.alpha_th_p(T=x, P=y, d=val)
                cp = Feistel2006EOS.cp_spec(T=x, P=y)
                s = Feistel2006EOS.s_spec(T=x, P=y)
                u = Feistel2006EOS.u_spec(T=x, P=y)
                '''
    
                if True:
                #if y < 1.0e4 and x > 800.:
                    if ll == 0:
                        dens, dTdPS, dPdrhoT, alpha, cp, s, u, phase =\
                        None, None, None, None, None, None, None, None
                        xiH2O = 0.0
                        
                        #French at very low temperatures takes forever to run
                        #and is not valid anyways
                        if not y > 3.0e12:
                            if x < 100 and y > 2.0e8:
                                pass
                            
                            else:
                                dens, dTdPS, dPdrhoT, alpha, cp, s, u, phase = \
                                eoswater.together(T=x, P=y, ph = None)
                            
                        #For some weird reasons, calling the water EoS from
                        #outside produces weird results at low P and calling
                        #it manually again to overwrite seems to work, wtf
                        if x>800. and y<1.0e4:
                            dTdPS = Wagner2002EOS.dTdP_S(T=x, P=y, d=dens)
                        
                    else:
                        try:
                            if update_phase:
                                ph = None
                                
                            else:
                                ph = phase
                            
                            P_eval = y
                            
                            if pres_thres:
                                if ll == 12 or ll==13 or ll==14:
                                    P_eval = min(y, 2.5e10)
                            
                            dens, xx, yy, dPdrhoT, ph, XH2O, alpha = \
                                                        eos.Compute(what='all', 
                                                        ll=ll, T=x, P=P_eval, 
                                                        saturation=saturation,
                                                        SiMg=SiMg, phase=ph,
                                                        Fe_number=Fe_number)
               
                            #If no phase is enforced, update the phase according
                            #to the phase diagram
                            if update_phase:
                                phase = ph
                            
                            #Compute water content as mass fraction
                            XH2O *= 0.01
                            
                            #Convert water content into mole fraction
                            if ll == 11:
                                m1 = mPer
                                m2 = mH2O
                                
                            elif ll == 12:
                                m1 = (1-.01*Fe_number)*mOl + 0.01*Fe_number*\
                                (mOl - 2*mMg + 2*mFe)
                                m2 = mH2O
                                
                            else:
                                m1 = 1.
                                m2 = 1.
                            
                            eta = XH2O
                            xiH2O = m1*eta/(eta*m1-eta*m2+m2)
                            
                        except ZeroDivisionError:
                            print ('zero division at T (K) / P (GPa):', x, y*1.0e-0)
                            xiH2O = 0.0
                            
                        except OverflowError:
                            print ('overflow at T (K) / P (GPa):', x, y*1.0e-9)
                            xiH2O = 0.0
                            
                        #here dTdP is dTdP/gammaG
                        dTdPS = x/(dens*dPdrhoT)
                        cp = None
                        s = None
                        u = None
                    
                    #except OverflowError:
                     #   print ('vals=', x, y, val)
                        #sys.exit()
                    
                    dat=[x, y, dens, dTdPS, dPdrhoT, alpha, cp, s, u, phase]
                    
                    #loop over all parameters
                    for k in range(n_params):
                        self.data[i][j][k] = dat[k]
                
                    if not len(dt_list)==0:
                        '''
                        dt_max=max(dt_list)
                        dt_mean=sum(dt_list)/len(dt_list)
                        eta_mean=(self.grid_size-count+1)*dt_mean
                        eta_max=(self.grid_size-count+1)*dt_max
                        '''
                        percentage=round(count/self.grid_size*100, 3)
                        
                        space=(total_length-int(percentage/100*total_length))*' '
                        bar=(int(percentage/100*total_length)-1)*'='
                        
                        bar+='>'
                        #time_max=s2min(eta_max)
                        #eta_mean = (eta_mean+eta_max)/2
                        #time_mean=s2min(eta_mean)
                        IPS=1./(time.time()-t0)
                        sys.stdout.write('\rProgress: ['+bar+space+'] '+\
                                         str(percentage)+'%, IPS='+str(int(IPS))+\
                                         'T (K):'+str(x)+', P (GPa):'+\
                                         str(round(y*1.0e-9,9)))
                        '''+\
                                         ', ETA_mean='+ \
                                         str(int(time_mean[0]))+'min '+\
                                         str(int(time_mean[1]))+'s'+\
                                         ' / ETA_max='+str(int(time_max[0]))+\
                                         'min '+str(int(time_max[1]))+'s'+\
                                         100*' ')
                        '''
                        sys.stdout.flush()
                        time_list.append(time.time()-t_start)
                        
                    dt=time.time()-t0
                    dt_list.append(dt)
                
        print ('\n-> grid has succesfully been generated')       
        print_time(t_start, time.time())
    
    
    def GenerateOlivine(self, Fe_numbers=[0., 5., 10., 15., 20., 25.],
                        alpha_x=1, alpha_y=1, x_range=3, y_range=11,
                        start=0, end=-1, x_start=1, y_start=3):
        
        for i in range(len(Fe_numbers)-start):
            print ('---------')
            
            #Two seeps: dry and hyd
            for j in range(2):
                
                if j == 0:
                    sat = False
                
                else:
                    sat = True
                
                #Three phases: alpha, beta, gamma
                for k in range(3):
                
                    count = (12+6*start+3*(2*i+j)+k)
                    print(str(count), sat, Fe_numbers[i+start])
                    
                    if sat == False:
                        self.Generate(ll=12 + k, x_start=x_start, y_start=y_start,
                                      x_range=x_range, y_range=y_range, 
                                      alpha_x=alpha_x, alpha_y=alpha_y,
                                      Fe_number=Fe_numbers[i+start], 
                                      saturation=sat,
                                      pres_thres=True,
                                      phase = k)
                        
                        self.write(file_name='eos_'+str(count)+'.tab')
                        self.write_fortran(file_name='eos_'+\
                                           str(count)+'_fortran.tab')
                        
                    else:
                        pass
        
        
    def reprocess(self, param=2, x_range=[], y_range=[], **kwargs):
        x1, x2 = x_range
        y1, y2 = y_range
        i00, i10 = self.get_index(x1, self.alpha_x), \
                    self.get_index(x2, self.alpha_x)
        i01, i11 = self.get_index(y1, self.alpha_y), \
                    self.get_index(y2, self.alpha_y)        
        
        print (i00, i10, i01, i11)
        #adjust x and y values to existing grid points
        
        x1, x2 = self.x_axis[i00], self.x_axis[i10]
        y1, y2 = self.y_axis[i01], self.y_axis[i11]

        #gather target x and y values
        x_axis_dummy = [self.x_axis[i+i00] for i in range(i10-i00+1)]        
        y_axis_dummy = [self.y_axis[j+i01] for j in range(i11-i01+1)]
        
        #reprocess all grid points in between these indices
        
        print ('computing grid values...\n')
        
        t_start=time.time()
        dt_list=[]
        time_list=[]
        total_length=50
        count=0
        IPS=0 #iterations per second
        grid_size=len(x_axis_dummy)*len(y_axis_dummy)
        
        percentage=int(count/grid_size*100)
        space=(total_length-int(percentage/100*total_length))*' '
        bar=int(percentage/100*total_length)*'=O'
        sys.stdout.write('\rProgress: ['+bar+space+']'+str(percentage)+\
                         ' % , ETA=\u221e min')
        
        sys.stdout.flush()
        
        for ii in range(len(x_axis_dummy)):
            #print ('Progress {:2.1%}'.format(count/self.grid_size), end='\r')
            for jj in range(len(y_axis_dummy)):
                i = i00+ii
                j = i01+jj
                t0=time.time()
                count+=1
                x=x_axis_dummy[ii]
                y=y_axis_dummy[jj]
                val=fct(x=x, y=y, param=param, **kwargs)[0]
                #print ('new=',val)
                #print ('old=', self.data[i][j][param])
                self.data[i][j][param]=val
                
                if not len(dt_list)==0:
                    dt_max=max(dt_list)
                    dt_mean=sum(dt_list)/len(dt_list)
                    eta_mean=(grid_size-count+1)*dt_mean
                    eta_max=(grid_size-count+1)*dt_max
                    
                    percentage=round(count/grid_size*100, 1)
                    space=(total_length-int(percentage/100*total_length))*' '
                    bar=(int(percentage/100*total_length)-1)*'='
                    
                    bar+='>'
                    time_max=s2min(eta_max)
                    eta_mean = (eta_mean+eta_max)/2
                    time_mean=s2min(eta_mean)
                    IPS=1/(time.time()-t0)
                    sys.stdout.write('\rProgress: ['+bar+space+'] '+\
                                     str(percentage)+'%, IPS='+str(int(IPS))+\
                                     ', ETA_mean='+ \
                                     str(int(time_mean[0]))+'min '+\
                                     str(int(time_mean[1]))+'s'+\
                                     ' / ETA_max='+str(int(time_max[0]))+\
                                     'min '+str(int(time_max[1]))+'s'+\
                                     100*' ')
                    
                    sys.stdout.flush()
                    time_list.append(time.time()-t_start)
                    
                dt=time.time()-t0
                dt_list.append(dt)
                
        
    def get_index(self, val, alpha, start):
        #compute order of magnitude of input temperature
        exponent=int(np.log10(val))
        #print ('exponent=', exponent)
        a=10**(exponent-alpha)
        #print ('a=', a)
       
       #compute the closest grid point values
        left_value=round((val/a-0.5), 0)*a
        right_value=round((val/a+0.5), 0)*a
        
        left_delta = abs(val-left_value)
        right_delta = abs(val-right_value)
        
        left_index = int(left_value/a-10**alpha)
        right_index = left_index+1#int(right_value/a+10**alpha)
        
        if left_delta < right_delta:
            offset=int(left_value/a-10**alpha)
            
        else:
            offset=int(left_value/a-10**alpha)+1
        
        
        print ('left value=', left_value)
        print ('right value=', right_value)  
        print ('left index=', left_index)
        print ('right_index=', right_index)
        
        #compute grid index that is closest to the given value
        ind=9*10**alpha*(exponent-start)+offset
        return ind

        
    def get_indices(self, x, alpha, start, order=2):
        """This function takes a x or y value as input and gives the indices of
        the neighbouring grid points (left and right) along the corresponding 
        grid axis in the eos table. The number of neighbours that are extracted
        depends on the order of the interpolatino polynomial.
        """
        relative_coordinates=[[0, 1], #1st order
                              [-1, 0, 1], #2nd order
                              [-1, 0, 1, 2], #3rd order
                              [-2, -1, 0, 1, 2], #4th order
                              [-2, -1, 0, 1, 2, 3], #5th order
                              [-3, -2, -1, 0, 1, 2, 3], #6th order
                              [-3, -2, -1, 0, 1, 2, 3, 4], #7th order
                              [-4, -3, -2, -1, 0, 1, 2, 3, 4]] #8th order
        
        
        #first check if value coincides with existing grid point
        
        
        #compute order of magnitude of input temperature
        exponent=int(np.log10(x))
        #print ('exponent=', exponent)
        
        #compute scaling factor
        a=10**(exponent-alpha)
        #print ('a=', a)
       
        #compute the closest grid point values
        left_value=round((x/a-0.5), 0)*a
        right_value=round((x/a+0.5), 0)*a
        
        #compute distance from left and right grid points to given value
        left_delta = abs(x-left_value)
        right_delta = abs(x-right_value)
        
        left_index = int(left_value/a-10**alpha)
        right_index = left_index+1#int(right_value/a+10**alpha)

        #NOTE: if the given value coincides with a existing grid point, the
        #above algorithm might not take the corresponding index as the core
        #point. But since the grid point which coincides with the given value
        #will be taken into account anyways, the interpolation will still yield
        #the precise result for this point
        
        
        #take left grid point as core
        offset = int(left_value/a-10**alpha)

        #gather indices of all required neighbouring grid points around the
        #core point for n'th order interpolation
        ind=[9*10**alpha*(exponent-start)+offset+r
             for r in relative_coordinates[order-1]]

        #print ('indices=',ind)
        #print ('values=',[self.pres_axis[i] for i in ind])
        return ind
    
    def gather_pairs(self, x, y, order=2):
        res=[]
        for i in range(len(x)):
            for j in range(len(y)):
                res.append([x[i][j], y[i][j]])
        return res
    
    def gather_values(self, x, y, order=2):
        res=[]
        for i in range(len(x)):
            for j in range(len(y)):
                res.append([self.values[i][j], self.values[i][j]])
        
        return res
    
    def row(self, pnt, order=2):
        """Compute row of coefficient matrix 'M' for 2nd order 2d interpolation
        """
        x, y = pnt
        r=[]
        for i in range(order+1):
            for j in range(order+1):
                r.append(x**i*y**j)
                
        return r
    
    def row1d(self, pnt, order=2):
        """Compute row of coefficient matrix 'M' for 2nd order 2d interpolation
        """
        x, y = pnt
        r=[]
        for i in range(order+1):
            r.append(x**i)
        return r
    
    def deriv_x_row(self, pnt, order=2):
        x, y = pnt
        r=[]
        for i in range(order+1):
            for j in range(order+1):
                r.append(i*x**(i-1)*y**(j))
        
        return r
    
    def deriv_x_row1d(self, pnt, order=2):
        x, y = pnt
        r=[]
        for i in range(order+1):
            r.append(i*x**(i-1))
        
        return r

    def deriv_y_row(self, pnt, order=2):
        x, y = pnt
        r=[]
        for i in range(order+1):
            for j in range(order+1):
                r.append(x**(i)*j*y**(j-1))
                
        return r
    
    def deriv_xy_row(self, pnt, order=2):
        x, y = pnt
        r=[]
        for i in range(order+1):
            for j in range(order+1):
                r.append((i)*x**(i-1)*j*y**(j-1))
                
        return r
        
    
    def construct_matrix(self, x, y, order=2):
        """Construct coefficient matrix for the points (x[i], y[i])
        """
        xx, yy = np.meshgrid(x, y)
        self.pairs = self.gather_pairs(xx, yy, order=order)
#        print ('pairs=', self.pairs)
        self.matrix = [self.row(p, order=order) for p in self.pairs]

    def interpolate_deriv(self, x=None, y=None, order=3, **kwargs):
        """Takes x and y values as input and gives f(x, y) as output using a
        predifined interpolation scheme. A 3rd order polynomial is fitted
        using the function values and derivatives at two bracketing grid points
        on the x-, and y-axis respectively.
        """
        t0=time.time()
        
        #Compute the relevant grid indices for the xy coordinates
        #for this scheme only two neighbouring points and the corresponding
        #derivatives are used which is why the argument order=1 is used although
        #the interpolation scheme here is 3rd order
        ind = [self.get_indices(x, self.alpha_x, order=1), 
               self.get_indices(y, self.alpha_y, order=1)]

        #gather pT grid points for the coefficient matrix
        x_list = [self.x_axis[i] for i in ind[0]]
        y_list = [self.y_axis[j] for j in ind[1]]

        print ('x_list:',x_list)
        print ('y_list:',y_list)
        #construct coefficient matrix
        self.matrix=[]
        self.b=[]
    
        #the first 4 rows correspond to the four grid point values
        for a in range(2):
            for b in range(2):
                self.b.append(self.data[ind[0][a]][ind[1][b]][2])
                self.matrix.append(self.row(pnt=[x_list[a], y_list[b]], 
                                            order=order))
                #print ('fct(a, b)=',fct(x_list[a], y_list[b]))
        
        #the next four rows correspond to the x-derivatives on six four points
        for a in range(2):
            for b in range(order-1):
                self.b.append(self.data[ind[0][a]][ind[1][b]][3])
                self.matrix.append(self.deriv_x_row(pnt=[x_list[a], y_list[b]], 
                                                    order=order))
                #print ('deriv(a,b)=', deriv(x=x_list[a], y=y_list[b], whicharg='x'))
    
        #the next four rows correspond to the y-derivatives at four grid points
        for a in range(2):
            for b in range(order-1):
                self.b.append(self.data[ind[0][a]][ind[1][b]][4])
                self.matrix.append(self.deriv_y_row(pnt=[x_list[a], y_list[b]], 
                                                    order=order))
                #print ('deriv(a,b)=', deriv(x=x_list[a], y=y_list[b], whicharg='y'))
        
        #the last four rows correspond to the mixed xy-derivatives at four 
        #grid points
        for a in range(order-1):
            for b in range(order-1):
                self.b.append(self.data[ind[0][a]][ind[1][b]][5])
                self.matrix.append(self.deriv_xy_row(pnt=[x_list[a], y_list[b]], 
                                                     order=order))
                #print ('deriv(a,b)=', deriv(x=x_list[a], y=y_list[b], whicharg='y'))
        
        
        #print ('matrix=',self.matrix)
        #print ('b vector=', self.b)
        self.a=np.linalg.solve(self.matrix, self.b)
        #print ('a vecotr=', self.a)
        result=sum([self.row([x, y], order=order)[i]*self.a[i] 
                        for i in range(len(self.a))])
        
        t=time.time()
        dt_int=t-t0                
        t0=time.time()
        f=fct(x=x, y=y, **kwargs)[0]
        dt_fct=time.time()-t0
        #result = self.a[0]+self.a[1]*x+self.a[2]*x**2+self.a[3]*x**3
        
        print ('true result=', f)
        print ('relative deviation=', round((f-result)/f*100, 3),'%')
        
        fig, ax = plt.subplots(1,2)

        x_plot_list=np.linspace(x_list[0], x_list[-1], 20)
        y_plot_list=np.linspace(y_list[0], y_list[-1], 20)
        
        ax[0].set_xlabel('x')
        ax[1].set_xlabel('y')
        for o in range(order-2):
            ax[0].scatter(x_list, [fct(x=t, y=y_list[o], **kwargs)[0]
            for t in x_list], color='r')
            ax[1].scatter(y_list, [fct(x=x_list[o], y=p, **kwargs)[0] 
            for p in y_list], color='b')
            
            ax[0].plot(x_plot_list, [fct(x=t, y=y_list[o], **kwargs)[0] 
            for t in x_plot_list], color='r')
            ax[1].plot(y_plot_list, [fct(x=x_list[o], y=p, **kwargs)[0]
            for p in y_plot_list], color='b')
            
            
            temp_results=[]
            for t in x_plot_list:
                temp_results.append(sum([self.row([t, y], order=order)[i]*\
                                         self.a[i] for i in 
                                         range(len(self.a))]))

            pres_results=[]
            for p in y_plot_list:
                pres_results.append(sum([self.row([x, p], order=order)[i]*\
                                         self.a[i] for i in 
                                         range(len(self.a))]))
            
            ax[0].plot(x_plot_list, temp_results, linestyle='--', color='r')
            ax[1].plot(y_plot_list, pres_results, linestyle='--', color='b')

        
            ax[0].scatter(x, fct(x=x, y=y_list[o], **kwargs), color='r', 
              zorder=10, marker='x')
            ax[1].scatter(y, fct(x=x_list[o], y=y, **kwargs), color='b', 
              zorder=10, marker='x')
            
            ax[0].scatter(x, fct(x=x, y=y, **kwargs), color='k', zorder=10, 
              marker='x')
            ax[1].scatter(y, fct(x=x, y=y, **kwargs), color='k', zorder=10, 
              marker='x')
            
            ax[0].scatter(x, result, color='k', zorder=10, marker='o')
            ax[1].scatter(y, result, color='k', zorder=10, marker='o')
            
            print ('elapsed time for interpolation:', 
                   round(dt_int*1000,3),'ms')
            print ('elapsed time for direct evaluation:', 
                   round(dt_fct*1000,3),'ms')
            print ('time gain factor:', 
                   round(dt_fct/dt_int, 3))
            
        return result
        
        
    def interpolate(self, x=None, y=None, order=1, plot=False, compare=False,\
                    phase_check=True, param=2, order_check=True, 
                    warnings=True, **kwargs):
        """Takes x and y values as input and gives f(x, y) as output using a
        polynomial nth order interpolation scheme. The scheme solves the
        matrix equation M*a=b where 'M' is the coefficient matrix and 'b'
        the solution vector.
        """
        t0=time.time()
        terminate=False
        if x<self.x_axis[order-1] or x>self.x_axis[-order]:
            print ('WARNING: x value must be in the range',\
                   self.x_axis[order-1],',', self.x_axis[-order])
            print ('Got', x)
            terminate=True
            
        if y<self.y_axis[order-1] or y>self.y_axis[-order]:
            print ('WARNING: y value must be in the range',\
                   self.y_axis[order-1], ',',self.y_axis[-order])
            print ('Got', y)
            terminate=True

        if not terminate:
            #Compute the relevant grid indices for the xy coordinates
            ind = [self.get_indices(x, self.alpha_x, self.x_start, 
                                    order=order), 
                   self.get_indices(y, self.alpha_y, self.y_start,
                                    order=order)]
                   
            phase = self.data[ind[0][0]][ind[1][0]][3]
            
            #gather pT grid points for the coefficient matrix
            try:
                x_list = [self.x_axis[i] for i in ind[0]]
                y_list = [self.y_axis[j] for j in ind[1]]                
            
            except IndexError:
                print (x,y,ind)

            #compute the solutions at the given grid points
            t1=time.time()
            self.b=[]
            for i in range(order+1):                
                for j in range(order+1):
                    #note that here j and i must be interchanged in comparison
                    #to the interpolate_deriv methode, otherwise it doesn't work
                    val = self.data[ind[0][j]][ind[1][i]][param]
   
                    if not val == 'None':
                        self.b.append(val)
                        
                    else:
                        self.b.append(0.)
            #print ('time for b vector:', round(1000*(time.time()-t1), 3))
            #if grid shift failed and the grid points still spread over
            #multiple phases, just switch to linlear interpolation to avoid
            #any kind unphysical results. Of course this is at the cost of
            #accuracy but it only happens for a few grid points along 
            #phase transitions
            if order_check:
                if min(self.b) < .25*max(self.b):
                    #print ('Reduziere zu linearer Interpolation')
                    #print ('b-Vektor=',self.b)
                    order = 1

                    #Compute the relevant grid indices for the xy coordinates
                    ind = [self.get_indices(x, self.alpha_x, self.x_start, 
                                            order=order), 
                           self.get_indices(y, self.alpha_y, self.y_start,
                                            order=order)]         
                    
                    #gather pT grid points for the coefficient matrix
                    x_list = [self.x_axis[i] for i in ind[0]]
                    y_list = [self.y_axis[j] for j in ind[1]]

                    #compute the solutions at the given grid points
                    self.b=[]
                    for i in range(order+1):                
                        for j in range(order+1):
                            #note that here j and i must be interchanged in comparison
                            #to the interpolate_deriv methode, otherwise it doesnt work
                            val = self.data[ind[0][j]][ind[1][i]][param]
                            self.b.append(val)
            
            #print ('bin_list=', bin_list)
            #compute coefficient matrix using the xy-grid points
            t2=time.time()
            self.construct_matrix(x=x_list, y=y_list, order=order)            
            #print ('time for matrix:', round(1000*(time.time()-t2), 3))
            #print ('b=',self.b)                    
            #print ('dim matrix=', len(self.matrix), len(self.matrix[0]))
    
            #compute the coefficients for the 2nd order polynom fit
            try: 
                t3=time.time()
                self.a=np.linalg.solve(self.matrix, self.b)
                #print ('time for linalg:', round(1000*(time.time()-t3), 3))
            
            except TypeError:
                print (self.b, self.matrix)
                print (x_list, y_list)
            #print ('intermediate time 2:', round((time.time()-t)*1000, 3), 'ms')
            #Evaluate 2nd order polynomial at x, y = T, p
            
            #compute 2nd order polynom solution at the given (P, T) point
            result=sum([self.row([x, y], order=order)[i]*self.a[i]
                        for i in range(len(self.a))])

            t=time.time()
            dt_int=t-t0      
            
            #ftool.printTime(sec=dt_int, digits=3, ms=True, where='old table interpolation')
            
            if plot:
                #print ('x_list=', x_list) 
                #print ('y_list=', y_list)
                #print ('x indices=', ind[0])
                #print ('y indices=', ind[1])
                fig, ax = plt.subplots(1,2)
    
                x_plot_list=np.linspace(x_list[0], x_list[-1], 15)
                y_plot_list=np.linspace(y_list[0], y_list[-1], 15)
                
                for o in range(order+1):
                    #gather grid points
                    valx_list = [self.data[ind[0][o]][j][param] for j in 
                                 range(order+1)]
                    valy_list = [self.data[i][ind[1][o]][param] for i in 
                                 range(order+1)]
                    
                    #print ('valx = ' ,valx_list)
                    #print ('valy = ' ,valy_list)
                    #plot used grid points
                    ax[0].scatter(x_list, valx_list, color='r')
                    ax[1].scatter(y_list, valy_list, color='b')
                    
                    #plot true function values between used grid points
                    ax[0].plot(x_plot_list, [fct(x=t, y=y_list[o], param=param,
                               **kwargs) for t in x_plot_list], color='r')
                    ax[1].plot(y_plot_list, [fct(x=x_list[o], y=p, param=param,
                               **kwargs) for p in y_plot_list], color='b')
                    
                    #compute function value between used grid points using
                    #the solution polynomial at input x, y
                    temp_results=[]
                    for t in x_plot_list:
                        temp_results.append(sum([self.row([t, y], order=order)[i]*\
                                                 self.a[i] for i in 
                                                 range(len(self.a))]))
    
                    pres_results=[]
                    for p in y_plot_list:
                        pres_results.append(sum([self.row([x, p], order=order)[i]*\
                                                 self.a[i] for i in 
                                                 range(len(self.a))]))
                    
                    #plot the solution polynomial between the used gridpoints
                    #at input x, y
                    ax[0].plot(x_plot_list, temp_results, linestyle='--', 
                      color='r')
                    ax[1].plot(y_plot_list, pres_results, linestyle='--', 
                      color='b')
    
                    #plot true values at grid points
                    #ax[0].scatter(x, fct(x=x, y=y_list[o], **kwargs), color='g', zorder=10, marker='s')
                    #ax[1].scatter(y, fct(x=x_list[o], y=y, **kwargs), color='g', zorder=10, marker='s')
                    
                #plot true value at x, y
                ax[0].scatter(x, fct(x=x, y=y, param=param, **kwargs), 
                  color='k', zorder=10, marker='x')
                ax[1].scatter(y, fct(x=x, y=y, param=param, **kwargs), 
                  color='k', zorder=10, marker='x')
                
                ax[0].scatter(x, result, color='k', zorder=10, marker='o')
                ax[1].scatter(y, result, color='k', zorder=10, marker='o')
            
            if compare:
                t0=time.time()
                f, actual_time=fct(x=x, y=y, param=param, **kwargs)
                dt_fct=(time.time()-t0)
                print ('estimated value=', result)
                print ('true value=', f)
                print ('deviation=', round((f-result)/f*100, 3),'%')
                print ('elapsed time for interpolation:', 
                       round(dt_int*1000, 3),'ms')
                print ('elapsed time for direct evaluation:', 
                       round(actual_time*1000, 3),'ms')
                print ('time gain factor:', 
                       round(dt_fct/dt_int, 2))
            
            return result
        