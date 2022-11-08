# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 10:56:16 2020

@author: shlah
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from samternary.ternary import Ternary

import Material

x = np.linspace(0., 1., 8)
y = np.linspace(0., 1., 8)
z = np.linspace(0., 1., 8)



# OP's data                                                             
A = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0,
              0.1, 0.2, 0.3, 0.4, 0.2, 0.2, 0.05, 0.1])
B = np.array([0.9, 0.7, 0.5, 0.3, 0.1, 0.2, 0.1, 0.15, 0, 0.1,
              0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
C = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.2, 0.2, 0.05, 0.1, 0.9,
              0.7, 0.5, 0.3, 0.1, 0.2, 0.1, 0.15, 0])

D = np.zeros([len(A)])

for i in range(len(D)):
    SiMg = Material.SiMg_UM(xiOl=A[i], xiPyr=B[i], xiStv=C[i])
    D[i] = SiMg/(1.+SiMg)

print (D)
# note that the array C above is not necessary since A+B+C=1            

titles = ['Lower mantle', 'Upper mantle']

# plot the data in two ways, in cartesian coordinates (ax_norm)         
# and in ternary-plot coordinates (ax_trans)                            

# create the figure and the two sets of axes                            
#fig, ax = plt.subplots(1,3, figsize=[5,2.8])
fig = plt.figure(figsize=[5, 3])

gs = gridspec.GridSpec(1, 2, width_ratios=[3, 3]) 
# plot data in normal way first using tricontourf                       

fnts1=10
fnts2=12
fnts3=14

ax = []

for i in range(2):
    axx = plt.subplot(gs[i])
    ax.append(axx)


for i in range(len(ax)):
    ax_trans = ax[i]

    # transform ax_trans to ternary-plot style, which includes              
    # building axes and labeling the axes                                   
    cob = Ternary(ax_trans, bottom_ax = 'bottom', left_ax = 'left',
                  right_ax = 'right',
                  labelpad=20)
    
    # use change of bases method within Ternary() to                        
    points = cob.B1_to_B2(A, B)
    
    # affine transform x,y points to ternary-plot basis                     
    cs = ax_trans.tricontourf(points[0], points[1],D)
    
    
    #ax_norm.set_title("Cartesian "
    #                 "(basis " + r"$\mathcal{B}_1$" + ")")
    ax_trans.set_title(titles[i])
    

cbaxes = fig.add_axes([.85, .3, .1, .5])
cbar = fig.colorbar(cs,ax=cbaxes, shrink=1., pad=0)
cbar.set_label(r'$\rm Si \#$', rotation=270, labelpad=20)

cbaxes.remove()

fig.subplots_adjust(bottom=0.2,hspace=0.01)
plt.show()