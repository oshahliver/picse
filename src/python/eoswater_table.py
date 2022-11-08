#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:11:24 2019

@author: oshah
"""

import numpy as np
from matplotlib import pyplot as plt
import functionTools as ftool
import eoswater as water

def fct(x=0, y=0):
    d=water.density(t=x, p=y)
    return d

class Table():
    def __init__(self, alpha):
        self.alpha = alpha
        self.values=None
        self.b=None
        self.matrix=None
        self.pairs=[]
        self.temp_axis=[]
        self.pres_axis=[]
        self.dimension=[len(self.temp_axis), len(self.pres_axis)]
    
    def generate(self):
        
        #first set up T-P grid
        self.temp_axis=[]
        self.pres_axis=[]
        
        print ('generating PT-grid...')
        print ('processing T-range...')
        
        for i in range(5):
            arr = np.linspace(10**i, 10**(i+1), 9*10**self.alpha, endpoint=False)
            for a in arr:
                self.temp_axis.append(a)
        
        print ('processing P-range...')
        for i in range(15):
            arr = np.linspace(10**i, 10**(i+1), 9*10**self.alpha, endpoint=False)
            for a in arr:
                self.pres_axis.append(a)
                
         
        #update dimension
        self.dimension=[len(self.temp_axis), len(self.pres_axis)]
        
        #compute value for each grid point
        self.values=[[fct(x=t, y=p) for p in self.pres_axis] for t in self.temp_axis]


     
    def get_T_index(self, T):
        """This function takes a temperature value as input and computes the 
        corresponding indices along the temperature axis in the eos table. The 
        alpha parameter specifies the refinement of the table that is used. See
        Table.generate() methode for further information.
        """
        
        #compute order of magnitude of input temperature
        exponent=int(np.log10(T))
        #print ('exponent=', exponent)
        a=10**(exponent-self.alpha)
        #print ('a=', a)
        left_value=round((T/a-0.5), 0)*a
        right_value=round((T/a+0.5), 0)*a
        #print ('left value=', left_value)
        #print ('right value=', right_value)
        left_index=exponent
        right_index=int(left_value/a-10**self.alpha)
        
        #print ('left index=', left_index)
        #print ('right index=', right_index)
        
        ind = [9*10**self.alpha*exponent+right_index-1,
               9*10**self.alpha*exponent+right_index,
               9*10**self.alpha*exponent+right_index+1]     
                                                                                                                                                            
        #print (ind)
        #print ('values=',[self.temp_axis[i] for i in ind])
        return ind
 
    def get_P_index(self, P):
        #compute order of magnitude of input temperature
        exponent=int(np.log10(P))
        #print ('exponent=', exponent)
        a=10**(exponent-self.alpha)
       # print ('a=', a)
        left_value=round((P/a-0.5), 0)*a
        right_value=round((P/a+0.5), 0)*a
       # print ('left value=', left_value)
       # print ('right value=', right_value)
        left_index=exponent
        right_index=int(left_value/a-10**self.alpha)
        
        #print ('left index=', left_index)
        #print ('right index=', right_index)
        
        ind = [9*10**self.alpha*exponent+right_index-1,
               9*10**self.alpha*exponent+right_index,
               9*10**self.alpha*exponent+right_index+1]     
                                                                                                                                                            
        #print ('indices=',ind)
        #print ('values=',[self.pres_axis[i] for i in ind])
        return ind
    
    def get_index(self, T, P):
        return [self.get_T_index(T), self.get_P_index(P)]
    
    def gather_pairs(self, x, y):
        res=[]
        for i in range(len(x)):
            for j in range(len(y)):
                res.append([x[i][j], y[i][j]])
        
        return res
    
    def b_vector(self):
        self.b=[fct(p[0], p[1]) for p in self.pairs]
    
    def row(self, pnt):
        x, y = pnt
        return [1, y, y**2, x, x*y, x*y**2, x**2, x**2*y, x**2*y**2]
    
    def construct_matrix(self, x, y):
        """construct coefficient matrix for the points (x[i], y[i])
        """
        xx, yy = np.meshgrid(x, y)
        self.pairs = self.gather_pairs(xx, yy)
#        print ('pairs=', self.pairs)
        self.matrix = [self.row(p) for p in self.pairs]

    
    def interpol_auto(self, T=None, P=None):
        if T<self.temp_axis[1] or T>self.temp_axis[-1]:
            print ('WARNING: temperature must be in the range',\
                   self.temp_axis[1],',', self.temp_axis[-2])
            
        if P<self.pres_axis[1] or P>self.pres_axis[-1]:
            print ('WARNING: pressure must be in the range',\
                   self.pres_axis[1], ',',self.pres_axis[-2])

        else:
            ind = self.get_index(T, P)
            temp = [self.temp_axis[i] for i in ind[0]]
            pres = [self.pres_axis[j] for j in ind[1]]
            self.construct_matrix(x=temp, y=pres)
            self.b_vector()
            self.a=np.linalg.solve(self.matrix, self.b)
            result=sum([self.row([T, P])[i]*self.a[i] for i in range(len(self.a))])
            f=fct(x=T, y=P)
            print ('true value=', f)
            print ('deviation=', round((f-result)/f*100,6),'%')
            return result
        
    def interpolate(self, T=None, P=None, order=1):
        """Compute the value at T using a polynom interpolation of order 
        <order> using an appropriate amount of grid point values in the 
        vicinity of T. Input format is:
            
            T (type=float, range=[0, 1.0e5]): temperature at which the function
                should be evaluated
                
            order (type=int, range=[1,6]): order of the polynomial 
                interpolation scheme
        """
        #define grid coordinates
        igrid, ii = self.get_T_index(T)
        jgrid, jj = self.get_P_index(P)
        
        #gather neighboring grid points
        grid_indices_T=[ii[0], ii[1], ii[1]+1, ii[0]-1, ii[1]+2, ii[0]-2]
        grid_indices_P=[jj[0], jj[1], jj[1]+1, jj[0]-1, jj[1]+2, jj[0]-2]
        
        grid_indices=[[ [grid_indices_T[i], grid_indices_P[j]] for i in range(2)] for j in range(2)]
        
        #extract x, y and dy/dx input for interpolation
        x_list = [self.temp_axis[igrid][tar] for tar in grid_indices_T]
        y_list = [self.values[igrid][tar] for tar in grid_indices_T]
        deriv_list = [self.derivatives[igrid][tar] for tar in grid_indices_T]
        
        #coefficient matrix. Format is:
        #(x1**order           x1**(order-1)           ....x1  1)
        #(x2**order           x2**(order-1)           ....x2  1)
        #(order*x1**(order-1) (order-1)*x1**(order-2) ....1   0)
        #(order*x2**(order-1) (order-1)*x2**(order-2) ....1   0)
        #(x3**order           x3**(order-1)           ....x3  1)
        #(.                   .                            .  .)
        #(.                   .                            .  .)
        #(.                   .                            .  .)
        matrix_dummy=[]
        
        #solution vecotr
        y_dummy=[]
        
        #parameter to control the dimension of linear system that must be
        #solved to find the polynom coefficients
        check=order
        
        #counts the rows added to the coefficient matrix
        counter=0
        
        #c-parameters to account for the fact that as the coefficient matrix
        #alternatingly contains two rows for the grid values and two for 
        #the corresponding derivatives for the interpolation up to the 
        #specified order
        c1=0
        c2=0
        
        #compute coefficient matrix (order x order)
        #the while loop allows for automatic adaption to desired order of the
        #interpolated polynomial
        while check>=0:
            #until the specified order of the matrix is reached, continue to
            #add rows to it
            if order-counter>=0:
                matrix_dummy.append([x_list[c1]**(order-o) 
                        for o in range(order+1)])
                counter+=1
                y_dummy.append(y_list[c1])
            
            if order-counter>=0:
                matrix_dummy.append([x_list[c1+1]**(order-o) 
                        for o in range(order+1)])
                y_dummy.append(y_list[c1+1])
                counter+=1
                c1+=2
                
            if order-counter>=0:
                matrix_dummy.append([(order-o)*x_list[c2]**(order-o-1)
                        for o in range(order+1)])
                y_dummy.append(deriv_list[c2])
                counter +=1
                
            if order-counter>=0:
                matrix_dummy.append([(order-o)*x_list[c2+1]**(order-o-1) 
                        for o in range(order+1)])
                y_dummy.append(deriv_list[c2+1])
                counter +=1
                c2+=2
            
            
            check-=counter
            #print ('\ncounter=', counter)
            #print ('check=', check)
      
        #compute coefficients of polynom of order "order"
        x_dummy=np.linalg.solve(matrix_dummy, y_dummy)
        
        #compute result at point T with the new polynom coefficients
        result=sum([x_dummy[o]*T**(order-o) for o in range(order+1)])
        #print ('true value=', fct(T))
        #print ('result=', result)
        print ('relative deviation=', round((result-fct(T))/fct(T)*100,4),'%')
        return result
