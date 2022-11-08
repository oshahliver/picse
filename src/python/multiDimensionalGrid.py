
from tqdm import tqdm
import pandas as pd
from PIMPphysicalparams import material_list
import numpy as np
import os
import eosTables
import Material
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
from astropy.nddata import NDData, CCDData
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PIMPphysicalparams import T_triple, T_critical, P_triple, P_critical, NA,\
                                mH, mO, kB, material_plot_list, mSi, mMg, mFe, \
                                mH2O, mPer, mPerov, mOl, mWs

from PIMPrunparams import param_labels, param_lims, water_phases

total_length=50

def s2min(val, r=3):
    m=int(round(val/60-.5, 0))
    s=round(val-m* 60, r)
    return m, s

def print_time(t0, t, **kwargs):
        print ('\nelapsed time: '+str(s2min(t-t0, **kwargs)[0])+ 'min '+\
                           str(s2min(t-t0, **kwargs)[1])+ 's')

def fct(x=None, y=None, param=2, table=False, **kwargs):
    #compute density
    if param == 2:
        t0=time.time()
        result = eos.Compute(what='dens', T=x, P=y, table=table, **kwargs)[0]
        actual_time = time.time()-t0
    
    #compute dP/drho
    elif param == 4:
        d = eos.Compute(what='dens', T=x, P=y, table=table, **kwargs)[0]
        t0=time.time()
        result = eos.dPdrho(T=x, P=y, d=d, table=table, **kwargs)
        actual_time=time.time()-t0
        
    #compute dP/dT    
    elif param == 5:
        pass
    
    return result, actual_time

    

class Table():
    def __init__(self,
                 scalings = ['log', 'lin', 'lin', 'lin'],
                 alphas = [0, 0, 0, 0],
                 starts = [100, 5, 1, 0.1],
                 ends = [1000, 9, 10, 1.],
                 n_out = 2,
                 n_phase = 1,
                 ):
        
                                 
        self.param_strings = [r'$T \ [\rm K]$',
                              r'$P \ [\rm GPa]$',
                              r'$\rm Fe\# \ [\rm mol\%]$',
                              r'$\rm Si\# \ [\rm mol\%]$',
                              r'$\epsilon_{\rm H_2 O}$',
                              r'$\epsilon_{\rm Al}$',
                              r'$\xi_{\rm Si O_2 }$']
        self.scalings = scalings
        self.alphas = alphas
        self.starts = starts
        self.ends = ends
        self.axes_lengths = [n_phase]
        self.n_phase = n_phase
        self.ranges = []
        self.single_axes = []
        self.n_out = n_out
        self.count = 0
        self.result = 0
        self.indices = None
        self.Fe_number = 0.
        self.ll = None
        self.tab = None
        self.N_pnts = 0
        self.phase = 0
        
        #Compute axis length for each input parameter
        for i in range(len(self.alphas)):
            
            if self.scalings[i] == 'log':
                length = (int(self.ends[i])-int(self.starts[i]))*9*10**self.alphas[i]+1
                
            elif self.scalings[i] == 'lin':
                length = 2**(self.alphas[i]+1)
    
            self.axes_lengths.append(length)
            self.ranges.append(range(length))
                
        self.axes_lengths.append(len(self.alphas)+self.n_out)
        
        #Define empty axes with the right dimensions
        self.axes = np.zeros([l for l in self.axes_lengths])
        
        self.construct_single_axes()
        self.N_pnts = 1
        
        #The first entry of axes_lengths corresponds to n_phase and the
        #last one to n_out. These two entries need to be ommited here
        for i in range(len(self.axes_lengths)-2):
            self.N_pnts *= self.axes_lengths[i+1]

    
    def construct_single_axes(self):
        self.single_axes = []
        for p in range(self.n_phase):
            self.single_axes.append([])
            for i in range(len(self.alphas)):
                self.single_axes[p].append([])
                if self.scalings[i] == 'log':        
                    
                    for j in range(int(self.ends[i]-self.starts[i])):
                        #for the last decade include endpoint to cover entire decade            
                        if j == int(self.ends[i]-self.starts[i]) -1:
                            arr = np.linspace(10**(j+self.starts[i]), 
                                              10**(j+1+self.starts[i]), 
                                              9*10**self.alphas[i]+1, 
                                              endpoint=True)
                        
                        #otherwise don't
                        else:
                            arr = np.linspace(10**(j+self.starts[i]), 
                                              10**(j+1+self.starts[i]), 
                                              9*10**self.alphas[i], 
                                              endpoint=False)     
                            
                        #append decade to temperature axis
                        for a in arr:
                            self.single_axes[p][i].append(a)
                            
                elif self.scalings[i] == 'lin':
                    arr = np.linspace(self.starts[i], self.ends[i], 
                                      2**(self.alphas[i]+1), endpoint=True)
           
                    #append decade to temperature axis
                    for a in arr:
                        self.single_axes[p][i].append(a)                        
            
        
    def get_index(self, val=None, which = 0):
        if self.scalings[which] == 'log':
            #compute order of magnitude of input temperature
            exponent=int(np.log10(val))
            #print ('exponent=', exponent)
            a=10**(exponent-self.alphas[which])
            #print ('a=', a)
           
            #compute the closest grid point values
            left_value=round((val/a-0.5), 0)*a
            right_value=round((val/a+0.5), 0)*a
            
            left_delta = abs(val-left_value)
            right_delta = abs(val-right_value)
            
            left_index = int(left_value/a-10**self.alphas[which])
            right_index = left_index+1#int(right_value/a+10**alpha)
            
            offset=int(left_value/a-10**self.alphas[which])
                
            
            print ('left value=', left_value)
            print ('right value=', right_value)  
            print ('left index=', left_index)
            print ('right_index=', right_index)
            
            #compute grid index that is closest to the given value
            ind=9*10**self.alphas[which]*(exponent-self.starts[which])+offset
            
            
        elif self.scalings[which] == 'lin':
            
            ind = int(val*(len(self.single_axes[which]) - 1)/\
                      (self.ends[which]-self.starts[which])+.5)
            
        return ind
   
    
    def get_indices(self, val = None, which=0, order=1):
        """This function takes a x or y value as input and gives the indices of
        the neighbouring grid points (left and right) along the corresponding
        grid axis in the eos table. The number of neighbours that are extracted
        depends on the order of the interpolation polynomial.
        """
        
        relative_coordinates=[[0, 1], #1st order
                              [-1, 0, 1], #2nd order
                              [-1, 0, 1, 2], #3rd order
                              [-2, -1, 0, 1, 2], #4th order
                              [-2, -1, 0, 1, 2, 3], #5th order
                              [-3, -2, -1, 0, 1, 2, 3], #6th order
                              [-3, -2, -1, 0, 1, 2, 3, 4], #7th order
                              [-4, -3, -2, -1, 0, 1, 2, 3, 4]] #8th order

        
        if self.scalings[which] == 'log':
            #compute order of magnitude of input temperature
            exponent=int(np.log10(val))
            #print ('exponent=', exponent)
            a=10**(exponent-self.alphas[which])
            #print ('a=', a)
           
            #compute the closest grid point values
            left_value=round((val/a-0.5), 0)*a
            right_value=round((val/a+0.5), 0)*a
                        
            left_index = int(left_value/a-10**self.alphas[which])
            right_index = left_index+1#int(right_value/a+10**alpha)
            
            offset=int(left_value/a-10**self.alphas[which])
            '''
            #compute grid index that is closest to the given value
            ind=[left_index+r
             for r in relative_coordinates[order-1]]            
            '''
            #gather indices of all required neighbouring grid points around the
            #core point for n'th order interpolation scheme
            ind=[9*10**self.alphas[which]*(exponent-self.starts[which])+offset+r
                 for r in relative_coordinates[order-1]]
        
        elif self.scalings[which] == 'lin':
            delta = (self.ends[which]-self.starts[which])/\
            (len(self.single_axes[self.phase][which])-1)
            ind = [int((val-self.starts[which])/delta) + r               
            for r in relative_coordinates[order-1]]
            
        return ind
    
    
    def for_recursive(self, c = 0, iter_list=[], f = None, ranges=None, 
                      count = 0, n=None, **kwargs):
        #print (iter_list, ranges, c, n)
        if n==None:
            n = len(self.alphas)
            
        if iter_list == []:
            iter_list = [0]*n
        
        if c == n-1:
             
            for iter_list[c] in ranges[c]:    
                f(iter_list = iter_list,
                  **kwargs)
                
        else:
            for iter_list[c] in ranges[c]:
                
                self.for_recursive(c=c+1, 
                                   iter_list=iter_list, 
                                   ranges=ranges, 
                                   f=f, 
                                   count=count, 
                                   n=n,
                                   **kwargs)
                
    
    def compute_input_grid(self, iter_list=None, **kwargs):
        for p in range(self.n_phase):
            for i in range(len(self.alphas)):
                self.axes[p][tuple(iter_list)][i] = self.single_axes[p][i][iter_list[i]]
                #print (p, i, iter_list[i])
                
    
    def construct(self):
        self.construct_single_axes()
        self.construct_grid()
        
    
    def construct_grid(self, **kwargs):
        """Construct input parameter grid according to the specification given
        in the Table instance initiation
        """
        print ('\ninitializing input grid...')
        print ('processing axes...')
         
        if not len(kwargs) == 0:
            self.__init__(**kwargs)
        
        print ('ranges =', self.ranges)
        self.for_recursive(f=self.compute_input_grid, 
                           ranges=self.ranges)

        self.N_pnts = 1
        
        #The first entry of axes_lengths corresponds to n_phase and the
        #last one to n_out. These two entries need to be ommited here
        for i in range(len(self.axes_lengths)-2):
            self.N_pnts *= self.axes_lengths[i+1]
            
            
    def compute_output_grid(self, iter_list=None, ll=None, Fe_number=0., 
                            table_data=False, mix=False, ph=0, **kwargs):
        
        x, y, z, a, b = self.axes[ph][tuple(iter_list)][0:5]

        #Only P and T as input
        if len(iter_list) == 2:
            z = Fe_number
            
        #P, T, and Fe# as input
        elif len(iter_list) == 3:
            pass
        
        #print (tuple(in_list))
        #print (self.axes[tuple(in_list)][0:3])
        if table_data:
            rho = self.tab.interpolate(x=x, y=y, param=2)
            if ll==0:
                #For pure water, instead of dPdrho, dTdP_S is stored directly
                dTdP_S = self.tab.interpolate(x=x, y=y, param=3)
            else:
                dTdP_S = 0.
                
            dPdrho = self.tab.interpolate(x=x, y=y, param=4)
            alpha = self.tab.interpolate(x=x, y=y, param=5)
            pnt = [rho, None, None, dTdP_S, dPdrho, None, 0., alpha, 0.]
            
        else:
            if self.eps_H2O > .9999:
                sat = True
            else:
                sat = False

            pnt = eos.Compute(what='all', 
                              ll = 
                              ll, 
                              T=x, 
                              P=y, 
                              Fe_number = z * 100.,
                              saturation = sat, 
                              phase = ph, 
                              **kwargs)

        if mix:
            self.mix.fractions=[1.0 - z, z]
            self.mix.Compute(T=x, P=y, phase=ph)
            
            rho = self.mix.dens
            dPdrho = self.mix.dPdrho
            alpha = self.mix.alpha
            
            pnt = [rho, None, None, .0, dPdrho, None, 0., alpha, 0.]
        
        self.axes[ph][tuple(iter_list)][-6] = pnt[0] #rho #x + 2*y + x*y + z*x*3
        self.axes[ph][tuple(iter_list)][-5] = pnt[3] #dTdP_S only for pure water
        self.axes[ph][tuple(iter_list)][-4] = pnt[4] #dPdrho #3*z*y + x*y*2 + z*x*4
        self.axes[ph][tuple(iter_list)][-3] = pnt[7] #alpha_th  #np.sqrt(x**2+y**2+z**2)
        self.axes[ph][tuple(iter_list)][-2] = pnt[6] #X_H2O
        self.axes[ph][tuple(iter_list)][-1] = pnt[8] #xi_Al
        
        self.it_counter += 1
        
        percentage=round(self.it_counter/self.N_pnts*100, 3)
        space=(total_length-int(percentage/100*total_length))*' '
        bar=(int(percentage/100*total_length)-1)*'='
        
        bar+='>'
        #time_max=s2min(eta_max)
        #eta_mean = (eta_mean+eta_max)/2
        #time_mean=s2min(eta_mean)
        sys.stdout.write('\rProgress: ['+bar+space+'] '+\
                         str(percentage)+'%')
        '''+\
                         ', ETA_mean='+ \
                         str(int(time_mean[0]))+'min '+\
                         str(int(time_mean[1]))+'s'+\
                         ' / ETA_max='+str(int(time_max[0]))+\
                         'min '+str(int(time_max[1]))+'s'+\
                         100*' ')
        '''
        sys.stdout.flush()        

            
    def generate(self, ll = None, Fe_number=0., table_data=False, tab=None,
                 mix=False, ll_plus=None, **kwargs):
        self.ll = ll
        self.Fe_number = Fe_number
        
        for ph in range(self.n_phase):
            print ('\nProcessing phase', ph+1, 'of', self.n_phase,':')
            if table_data:
                if tab == None:
                    print ('WARNING: No table given.')
                    print ('Trying to locate table file...')
                    self.tab = eosTables.Table()
                    self.tab.load('eos_'+str(ll)+'.tab')
                    print ('Table file loaded successfully.')
                    
                else:
                    self.tab = tab
                    
            if mix:
                self.mix = Material.Mixture(contents=[ll, ll_plus], 
                                            eos_table=False)
            
            self.it_counter = 0
            
            dt_list=[]
            time_list=[]
            IPS=0 #iterations per second
            
            percentage=int(self.it_counter/self.N_pnts*100)
            space=(total_length-int(percentage/100*total_length))*' '
            bar=int(percentage/100*total_length)*'=O'
            sys.stdout.write('\rProgress: ['+bar+space+']'+str(percentage))
            
            sys.stdout.flush()
            
            self.for_recursive(f=self.compute_output_grid,
                               ranges=self.ranges,
                               ll=ll,
                               Fe_number=Fe_number,
                               table_data=table_data,
                               mix=mix,
                               ph=ph, 
                               **kwargs)
            self.it_counter = 0
    
    
    def partial_row(self, iter_list, row_index=0,
                    vals=None):
        """Computes the individual matrix elements for the current row. In each
        row there are (order+1)**how_many elements. The elements are permutations
        of the products of the input parameters.
        """
        
        vec = [vals[i]**iter_list[i] for i in range(len(iter_list))]
        res = 1.
        for v in vec:
            res *= v
        
        #print ('iter_list', iter_list,' vals', vals,' res =', res)
        
        self.matrix[row_index][self.count] = res
        self.count += 1
        
    
    def index_row(self, iter_list=None, order=1):
        """Computes a row of the index_matrix. The row corresponds to the
        indices of all input parameters at the given row index. There are a 
        total of (order+1)**how_many rows. Each row has how_many entries.
        """
        #print (iter_list)
        self.index_matrix[self.count] = [self.indices[i][iter_list[i]]
                                            for i in range(len(self.alphas))]
        self.count += 1
    
    
    def construct_index_matrix(self, vals=None, order = 1):
        """Computes all indices of the grid points that are needed for
        the interpolation and collects them in a matrix. The number of rows
        corresponds to (order+1)**n, i.e. the total amount of grid points
        that is needed for the interpolation at given order and with n
        input parameters. The number of columns corresponds to the number of
        input parameters, i.e. n because for each grid point, each
        input parameter has a distinct index.
        """
        self.index_matrix = np.zeros([(order+1)**len(self.alphas), 
                                      len(self.alphas)], dtype='int')
    
        self.indices = [self.get_indices(val=vals[i], which=i) 
        for i in range(len(self.alphas))]
        
        self.count = 0
        iter_ranges = [range(0, order+1) for i in range(len(self.alphas))]
        self.for_recursive(f=self.index_row, 
                           ranges=iter_ranges)
        self.count = 0        
        
        
    def row(self, iter_list=None, which=0, order=1):
        """Compute row of coefficient matrix 'M' for n-d interpolation
        """
        #Extract values for input parameters at all relevant grid points
        vals = self.axes[self.phase][tuple(self.index_matrix[self.count])][0:len(iter_list)]
        
        #Extract value for output parameter(s) at these grid points
        val =self.axes[self.phase][tuple(self.index_matrix[self.count])][which]

        #x, y, z = self.axes[tuple(iter_list)][0:3]
        #val = self.axes[tuple(iter_list)][which]
        '''
        r=[]
        for i in range(order+1):
            for j in range(order+1):
                for k in range(order+1):
                    r.append(x**i*y**j*z**k)
        '''
        count = self.count
        self.count = 0
        
        iter_ranges = [range(0,order+1) for i in range(len(self.alphas))]
        
        #Compute matrix elements for the current row
        self.for_recursive(f=self.partial_row, 
                           ranges=iter_ranges,
                           row_index=count,
                           vals = vals)
        
        self.count = count
        
        #Update the current element of the output parameter vector
        self.b[self.count] = val
        
        self.count += 1
    

    def construct_matrix(self, which=[],  order=1):
        """Construct coefficient matrix for the grid points in n-d according
        to the interpolation order
        """
        self.count = 0
        self.matrix = np.zeros([2**len(self.alphas), 2**len(self.alphas)])
        self.b = np.zeros([2**len(self.alphas), len(which)])
        self.a = np.ones([2**len(self.alphas), len(which)])
        #self.matrix = [self.row(p, order=order) for p in self.pairs]
        
        iter_ranges = [range(0,order+1) for i in range(len(self.alphas))]
        
        #self.result = 0.
        self.for_recursive(f=self.row, ranges=iter_ranges, which=which)
        self.count = 0
    
    
    def compute_result(self, iter_list = None, vals=None):
        """Returns F(vals) where F is the desired parameter(s) which is given
        as a function(s) of all input parameters, which are contained in the
        vals array.
        """
        prod = 1.
        for i in range(len(iter_list)):
            prod *= vals[i]**iter_list[i]
            #print ('prod=', prod)
        
        delta = self.a[self.count]*prod
        self.result += delta
        self.count += 1
    
        
    def interpolate(self, vals = None, order=1, phase_check=True, params=[2], 
                    order_check=False, warnings=True, proc=0, phase=0):
        """Takes x=temp and y=pres values as input and gives f(x, y) as output
        using a polynomial n'th order interpolation scheme. The scheme solves the
        matrix equation M*a=b where 'M' is the coefficient matrix and 'b'
        the solution vector.
        """
        t0 = time.time()
        
        terminate = False
        print ('interpolate: vals, params =', vals, params)
        self.phase = phase
        if not terminate:            
            if proc == 0:
                #Construct matrix containing the indices of all grid points that
                #are needed for interpolation along all relevant axes
                self.construct_index_matrix(vals=vals, order=order)
                
                #Construct n-linear coefficient matrix
                self.construct_matrix(which=params)
                
                #Compute the solution vector which gives the coefficients
                #of the n-linear interpolation
                self.a=np.linalg.solve(self.matrix, self.b)
                self.count = 0
                self.result = 0.
                
                #Use the solution to compute the result(s) at the given values for
                #the input parameters
                iter_ranges = [range(0,order+1) for i in range(len(self.alphas))]
                self.for_recursive(f=self.compute_result, ranges=iter_ranges,
                                   vals = vals)
                
                self.count = 0
                
            elif proc==1:
                self.construct_index_matrix(vals=vals, order=order)
                
                #Compute slopes for first slice
                slopes = []
                pnts = []
                for i in range(4):
                    slopes.append([])
                    s = (self.axes[phase][tuple(self.index_matrix[2*i+1])][3] -
                        self.axes[phase][tuple(self.index_matrix[2*i])][3])/\
                        (self.axes[phase][tuple(self.index_matrix[2*i+1])][2] -
                        self.axes[phase][tuple(self.index_matrix[2*i])][2])
                    
                    print ('s =', s)
                    
                    pnt = self.axes[phase][tuple(self.index_matrix[2*i])][3]+\
                        s*(vals[2]-self.axes[phase][tuple(self.index_matrix[2*i])][2])
                        
                    print ('pnt =', pnt)
                    pnts.append(pnt)
                    
                print ('pnts =', pnts)
                
                sy1 = (pnts[1] - pnts[0])/\
                    (self.axes[phase][tuple(self.index_matrix[2])][1]-\
                     self.axes[phase][tuple(self.index_matrix[0])][1])
                            
    
                sy2 = (pnts[3] - pnts[2])/\
                    (self.axes[phase][tuple(self.index_matrix[3])][1]-\
                     self.axes[phase][tuple(self.index_matrix[1])][1])
                    
                print ('sy1, sy2 =', sy1, sy2)
                
                dy1 = (vals[1]-self.axes[phase][tuple(self.index_matrix[0])][1])
                dy2 = (vals[1]-self.axes[phase][tuple(self.index_matrix[2])][1])
                
                q11 = pnts[0] + sy1*dy1
                q12 = pnts[2] + sy2*dy2
                print ('dy =', dy1, dy2)
                print ('q11, q12 =', q11, q12)
                
                        
        t = time.time()
        ftool.printTime(sec=t-t0, digits=4, ms=True, where='new table interpolation')
        return self.result
  

    def write(self, file_name, loc='cwd'):
        """Saves eos table which has been generated with the <generate>
        methode to ascii file
        """
        t0=time.time()
        '''
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
        '''
        if loc=='cwd':
            dir_out=os.getcwd()+'/'
            
        else:
            dir_out=loc+'/'
        
        print ('generating data table...')
        self.data_table=astropy.table.Table(self.axes, meta={})
        
        fh = open(file_name, 'bw')
        self.axes.tofile(fh)
        
        
    def load(self, file_name, loc='./eos_tables/'):
        fh = open(loc+file_name, 'rb')
        self.axes = np.fromfile(fh).reshape(self.axes_lengths[1:])
    
    
    def write_individual_params(self, iter_list=None, ph=0):
        for n in range(self.axes_lengths[-1]):
            self.data[ph][n][self.count] = self.axes[ph][tuple(iter_list)][n]
            self.data_str[ph][n] += str(self.data[ph][n][self.count])+','
        
        self.count += 1
    

    def read_individual_params(self, iter_list=None, ll=None, cc=None, ph=0):
        """
        """
        #print ('count =', self.count)
        self.axes[ph][tuple(iter_list)][cc] = float(ll[self.count])
        self.count += 1
        
       
    def create_data_file(self, file_name = None, meta = True):
        """Write table data to file which can be read back for interpolation by
        the read_data_file methode or loaded into the fortran module eosfort.f95
        """
        t0 = time.time()
        new_ll = [1, 6, 7, 0, 0, 4, 5, 9, 0, 2, 8, None, None, None, None, None,
                  3, None, 10, 0, 11]
        
        self.material = new_ll[self.ll]
        ll_str = str(self.material)
        H2O_str = str(round(self.eps_H2O, 3)).replace('.', '_')
        Al_str = str(self.eps_Al).replace('.', '_')
        
        #Generate default name
        if file_name == None:
            if len(self.single_axes[0]) == 2:
                Fe_str = str(self.Fe_number/100).replace('.', '_')
            else:
                Fe_str = 'var'
            
            file_name = 'eos_'+ll_str+'_Fe_'+Fe_str+'_H2O_'+H2O_str+'_Al_'+Al_str+'.tab'
        
        file_name = 'eos_tables/'+file_name
        print ('Creating table: ', file_name)
        #Compute number of grid points
        n_grid = 1
        for i in range(len(self.axes_lengths)-2):
            n_grid *= self.axes_lengths[i+1]
        
        #Rearange data
        self.data = np.zeros([self.n_phase, self.axes_lengths[-1], n_grid])
        self.data_str = [['']*self.axes_lengths[-1] for p in range(self.n_phase)]
        
        for ph in range(self.n_phase):
            #c=0
            self.count = 0
            iter_ranges = [range(0, al) for al in [self.axes_lengths[i+1] 
                            for i in range(len(self.axes_lengths)-2)]]

            self.for_recursive(f=self.write_individual_params, ranges=iter_ranges,
                               ph=ph)
        
        '''
        for i in range(self.axes_lengths[0]):
            for j in range(self.axes_lengths[1]):
                for k in range(self.axes_lengths[2]):
                    for n in range(self.axes_lengths[-1]):
                        data[n][c] = self.axes[i][j][k][n]
                        data_str[n] += str(data[n][c])+','
                    
                    c += 1
        '''
        
        self.count = 0
        
        f = open(file_name, 'w')
        
        #Write meta data
        f.write('material='+ll_str+'\n')
        f.write('eps_H2O='+str(self.eps_H2O)+'\n')
        f.write('eps_Al='+str(self.eps_Al)+'\n')
        f.write('n_input_params='+str(len(self.alphas))+'\n')
        f.write('n_phase='+str(self.n_phase)+'\n')
        
        dims_string = 'dims='
        for dd in range(len(self.axes_lengths)-1):
            d = self.axes_lengths[dd+1]
            dims_string += str(d)+','
        
        dims_string = dims_string[:-1]
        f.write(dims_string+'\n')
        
        alphas_string = 'alphas='
        for a in self.alphas:
            alphas_string += str(a)+','
            
        alphas_string = alphas_string[:-1]
        f.write(alphas_string+'\n')
        
        scaling_string = 'scaling='
        for s in self.scalings:
            scaling_string += s+','
        
        scaling_string = scaling_string[:-1]
        f.write(scaling_string+'\n')
        
        start_string ='starts='
        for s in self.starts:
            start_string += str(s)+','
        
        start_string = start_string[:-1]
        f.write(start_string+'\n')

        ends_string ='ends='
        for e in self.ends:
            ends_string += str(e)+','
        
        ends_string = ends_string[:-1]
        f.write(ends_string+'\n')        

        #Write data
        for p in range(self.n_phase):
            print ('writing phase ', p)
            for n in range(self.axes_lengths[-1]):
                f.write('\n'+self.data_str[p][n][:-1])
            
        f.close()
        
        t = time.time()
        print ('Elapsed time for creating table file:', round(t-t0,2), 'sec')


    def read_data_file(self, file_name='table.tab', loc='./eos_tables/'):
        f = open(loc+file_name, 'r')
    
        meta_strings = ['material', 'eps_H2O', 'eps_Al', 
                        'n_input_params', 'n_phase', 'dims', 
                        'alphas', 'scalings', 'starts', 'ends']
        
        meta_data = []
        
        print ('extracting meta data...')
        
        #Extract meta data
        for i in range(len(meta_strings)):
            line = f.readline()
            meta = line.split('=')[-1]
            meta = meta.split(',')

            #Extract string values
            if i == 7:
                meta = [str(m) for m in meta]
                meta[-1] = meta[-1][:-1]
            
            #Extract float values
            elif i == 8 or i == 9 or i==1 or i==2:
                meta = [float(m) for m in meta]
                
                if i==1 or i==2:
                    meta = meta[0]
            
            #Extract integer values
            else:
                meta = [int(m) for m in meta]
                
                if i == 0 or i == 3 or i == 4:
                    meta = meta[0]
                    
            print (i, ':', meta_strings[i]+' =', meta)               
            meta_data.append(meta)
        
        print ('----------')
        print ('extracting table data...')
        #Prepare empty data array according to the axes dims. The axes dims
        #are stored in the first meta data entry
        self.n_phase = meta_data[4]
        self.eps_H2O = meta_data[1]
        self.eps_Al = meta_data[2]
        self.material = meta_data[0]
        self.n_input_params = meta_data[3]
        self.dims = meta_data[5]
        
        axes_dims = [self.n_phase]
        for i in range(len(meta_data[5])):
            axes_dims.append(meta_data[5][i])
        
        self.axes = np.zeros(axes_dims)
        self.axes_lengths = axes_dims
        self.alphas = meta_data[6]
        self.scalings = meta_data[7]
        self.starts = meta_data[8]
        self.ends = meta_data[9]
        
        self.construct_single_axes()
        self.ranges = [range(0, r) for r in self.dims]
        self.construct_grid()       
 
        #Extract data
        #Continue at the next line after the meta data has been extracted
        
        
        #scip one line
        line = f.readline()
        
        #loop over individual parameters
        for p in range(self.n_phase):
            c = 0
            while True:
                print ('Extracting parameter:', c)
                line = f.readline()
                ll = line.split(',')
                ll[-1] = ll[-1].split('\n')[0]
    
                self.count = 0
                iter_ranges = [range(0, al) for al in axes_dims[1:-1]]
                print ('iter_ranges =', iter_ranges)
                self.for_recursive(f=self.read_individual_params, 
                                   ranges=iter_ranges,
                                   ll=ll, cc=c, ph=p)
                self.count = 0
                
                '''
                cc=0
                
                for i in range(dims[0]):
                    for j in range(dims[1]):
                        for k in range(dims[2]):
                            try:
                                data[i][j][k][c] = float(ll[cc])
                            except ValueError:
                                pass
                            
                            except IndexError:
                                pass
                            
                            cc += 1
                '''         
                c += 1
                
                if c == meta_data[5][-1]:
                    break
            
        print ('table data extraction succesful.')
            
            
    def plot_slice(self, axes = [0, 1], params = None, vals = [], ph = 0,
                   extent = None, scalings = ['log', 'log'], cmap = 'jet', 
                   addCbar = True, label = 'A Colorbar', 
                   scale = 1., fnts = 14, **kwargs):
        
        
        if extent is None:
            extent = [self.starts[0], self.ends[0], self.starts[1], self.ends[1]]
        
        if not len(self.alphas) - len(axes) == len(vals):
            print ('WARNING: invalid number of values given!')
            
        else:
            #Gather data for specified parameter along the given axes
            #Note that 0th axis corresponds to phase so 1 + axes value
            plot_data = np.zeros([self.axes_lengths[1 + axes[0]], 
                                self.axes_lengths[1 + axes[1]]])
            
            indices = np.zeros([len(self.alphas)], dtype='int')
            
            '''
            iii = 0
            for ii in range(len(indices)):
                if not ii in axes:
                    indices[ii] = vals[iii]
                    iii += 1
            '''
            c = 0
            #Insert index of slices along correct axes (not the main axes given)
            for k in range(len(indices)):
                if k not in axes:
                    indices[k] = vals[c]
                    c += 1
            
            #Gather table data along specified axes across given slice
            for i in range(self.axes_lengths[1 + axes[0]]):
                for j in range(self.axes_lengths[1 + axes[1]]):
                    indices[axes[0]] = i
                    indices[axes[1]] = j
                    plot_data[i][j] = self.axes[ph][tuple(indices)][params]
        
        try:
            ax = kwargs['ax']
            fig = kwargs['fig']
        except KeyError:
            fig, ax = plt.subplots()
            
        xTicks = np.linspace(self.starts[0], self.ends[0], 3)
        yTicks = np.linspace(self.starts[1], self.ends[1], 5)
        if scalings[0] == 'log':
            xString = f'$log({self.param_strings[axes[0]]})$'
        else:
            xString = self.param_strings[axes[0]]
        
        if scalings[1] == 'log':
            fore = r'\rm log'
            yString = f'${fore}({self.param_strings[axes[1]]})$'
        else:
            yString = self.param_strings[axes[1]]
        
        ax.set_xticklabels([f'${t:.0f}$' for t in xTicks],
                           fontsize = fnts)
        ax.set_yticklabels([f'${t:.0f}$' for t in yTicks * 1e-9],
                           fontsize = fnts)
            
        ax.set_xlabel(xString, fontsize = fnts)
        ax.set_ylabel(yString, fontsize = fnts)
        
        extent = [self.starts[0], self.ends[0], self.starts[1], self.ends[1]]
        aspect = (extent[1] - extent[0]) / (extent[3] - extent[2])

        im = ax.imshow(plot_data.T * scale, 
                       origin='lower', 
                       extent = extent, 
                       aspect = aspect,
                       cmap = cmap)
        
        ax.set_xticks(xTicks)
        ax.set_yticks(yTicks)
        
        if addCbar:
            axin = ax.inset_axes([1.1, .0, 0.05, 1.])
            cbar = plt.colorbar(im, cax = axin)
            cbar.ax.tick_params(labelsize = fnts)
            cbar.set_label(label, size = fnts)
            
        