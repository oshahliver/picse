"""
#Example: 
    
import grid_interpolation

#Initiate grid instance
grid = grid_interpolation.Grid()

#Inport grid data from file
grid.read_data_file(file_name = <file name>, loc = <file location>)

#1
#Compute radius and luminosity at T = 254 K, log(P/Pa) = 7.4, M_ocean/M = 0.28
result = grid.interpolate(vals = [254., 7.4, 0.28], params = [3,4])

#2
#Compute E_grav and E_int at T = 254 K, log(P/Pa) = 7.4, M_ocean/M = 0.28
result = grid.interpolate(vals = [254., 7.4, 0.28], params = [5,6])

#3
#Plot radius as function of T and P for M_ocean/M = 0.1
grid.plot_slice(axes = [0, 1], slice = 0, param = 3)

#3
#Plot luminosity as function of T and M_ocean for log(P/Pa) = 5
grid.plot_slice(axes = [0, 2], slice = 0, param = 4)
"""


import numpy as np
from matplotlib import pyplot as plt
import eos
import sys
import time
#from PlanetInvestigator import run_planet

total_length=50
r_earth = 6371e3

class Grid():
    def __init__(self, 
                 scalings = ['lin', 'lin', 'lin'],
                 alphas = [0, 0, 0],
                 starts = [100, 5, .1],
                 ends = [500, 8, .5],
                 n_out = 4,
                 n_phase = 1
                 ):
        """Convention is: T (K), log(P/Pa), M_ocean/M, R/R_E, Lum. (W), E_grav (J), E_int (J)
        """
        
        self.param_scales = [1., 1., 1., 1. / r_earth, 1., 1., 1.]
        self.param_strings = [r'$T \ [\rm K]$',
                              r'$\log\left(P / \rm Pa\right)$',
                              r'$M_{\rm Ocean}/M$',
                              r'$R/R_\oplus$',
                              r'$L \ [\rm W]$',
                              r'$E_{\rm grav} \ [\rm J]$',
                              r'$E_{\rm int} \ [\rm J]$']
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
                              [-1, 0, 1, 2]] #3rd order

        
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
                
        
    def construct_grid(self, **kwargs):
        """Construct input parameter grid according to the specification given
        in the Table instance initiation
        """
        print ('initializing input grid...')
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
            
            
    def compute_output_grid(self, iter_list=None, ph = 0, **kwargs):
        
        params = self.axes[ph][tuple(iter_list)][0:7]

        planet = run_planet(T_surface = params[0],
                            P_surface = 10**params[1],
                            ocean_mass_frac = params[2])
            
        self.axes[ph][tuple(iter_list)][-4] = planet.R_surface_is
        self.axes[ph][tuple(iter_list)][-3] = planet.luminosity
        self.axes[ph][tuple(iter_list)][-2] = planet.E_grav
        self.axes[ph][tuple(iter_list)][-1] = planet.E_int
        
        self.it_counter += 1
        
        percentage=round(self.it_counter/self.N_pnts*100, 3)
        space=(total_length-int(percentage/100*total_length))*' '
        bar=(int(percentage/100*total_length)-1)*'='
        
        bar+='>'
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

            
    def generate(self, **kwargs):
        
        for ph in range(self.n_phase):
            print ('\nProcessing phase', ph+1, 'of', self.n_phase,':')            
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
    
        
    def interpolate(self, vals = None, order=1, params=[2], warnings='on', 
                    phase=0):
        """Takes x=temp and y=pres values as input and gives f(x, y) as output
        using a polynomial n'th order interpolation scheme. The scheme solves the
        matrix equation M*a=b where 'M' is the coefficient matrix and 'b'
        the solution vector.
        """
        
        terminate = False
        
        self.phase = phase
        if not terminate:            
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
    
        return self.result

    
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
                
        #Generate default name
        if file_name == None:
            
            file_name = 'planets_grid.tab'
            
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
        
        self.count = 0
        
        f = open(file_name, 'w')
        
        #Write meta data
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


    def read_data_file(self, file_name='planets_grid.tab', loc='./'):
        f = open(loc+file_name, 'r')
    
        meta_strings = [
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
            if i == 4:
                meta = [str(m) for m in meta]
                meta[-1] = meta[-1][:-1]
            
            #Extract float values
            elif i == 5 or i == 6:
                meta = [float(m) for m in meta]
            
            #Extract integer values
            else:
                meta = [int(m) for m in meta]
                
                if i == 0 or i == 1:
                    meta = meta[0]
                    
            print (i, ':', meta_strings[i]+' =', meta)               
            meta_data.append(meta)
        
        print ('----------')
        print ('extracting table data...')
        #Prepare empty data array according to the axes dims. The axes dims
        #are stored in the first meta data entry
        self.n_phase = meta_data[1]
        self.n_input_params = meta_data[0]
        self.dims = meta_data[2]
        
        axes_dims = [self.n_phase]
        for i in range(len(meta_data[2])):
            axes_dims.append(meta_data[2][i])
        
        self.axes = np.zeros(axes_dims)
        self.axes_lengths = axes_dims
        self.alphas = meta_data[3]
        self.scalings = meta_data[4]
        self.starts = meta_data[5]
        self.ends = meta_data[6]
        
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
                         
                c += 1
                
                if c == meta_data[2][-1]:
                    break
            
        print ('table data extraction succesful.')
            
            
    def plot_slice(self, axes = [0, 1], slice = 0, param = 3, ph = 0,
                   extent = None, cmap = 'viridis'):
        
        axx = np.arange(self.n_input_params)
        for a in axes:
            axx = np.delete(axx, np.argwhere(axx == a))
        #axx = np.delete(axx, drop)
        if extent is None:
            extent = [self.starts[0], self.ends[0], self.starts[1], self.ends[1]]

        fig, ax = plt.subplots()
        
        ax.set_xlabel(self.param_strings[axes[0]])
        ax.set_ylabel(self.param_strings[axes[1]])
        
        x = np.meshgrid(self.single_axes[ph][axes[0]],self.single_axes[ph][axes[0]])
        y = np.meshgrid(self.single_axes[ph][axes[1]],self.single_axes[ph][axes[1]])

        #Gather data for specified parameter along the given axes
        #Note that 0th axis corresponds to phase so 1 + axes value
        plot_data = np.zeros([self.axes_lengths[1 + axes[0]], 
                            self.axes_lengths[1 + axes[1]]])
        
        indices = np.array([slice, slice, slice])
        
        for i in range(self.axes_lengths[1 + axes[0]]):
            for j in range(self.axes_lengths[1 + axes[1]]):
                indices[axes[0]] = i
                indices[axes[1]] = j
                plot_data[i][j] = self.axes[ph][tuple(indices)][param]
    
        sc = ax.scatter(x[0], y[0].T,
                   c = plot_data.T * self.param_scales[param],
                   marker = 's',
                   s = 100,
                   cmap = cmap)
        
        ax.set_title(self.param_strings[axx[0]] + ' = ' + str(self.single_axes[ph][axx[0]][slice]))
        cbar = fig.colorbar(sc, ax=ax, label = self.param_strings[param])
        
        
        
