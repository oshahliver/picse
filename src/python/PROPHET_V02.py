#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 12:38:17 2020

@author: oshah
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap 
from matplotlib import cm

#Define rotatin matrix along angle alpha in 2d
def R(alpha):
    c, s = np.cos(alpha), np.sin(alpha)
    
    return np.array(((c, -s), (s, c)))


def matrix_product(M=None, X=None):
    return np.array([sum([M[0][i]*X[i] for i in range(len(X))]),
                    sum([M[1][i]*X[i] for i in range(len(X))])])
    
def vector_product(X=None, Y=None):
    return sum(X[i]*Y[i] for i in range(len(X)))
    

def rotation(alpha, X):
    """Rotate the vector x around the angle alpha
    """
    return np.array([sum([R(alpha)[0][i]*X[i] for i in range(len(X))]),
                    sum([R(alpha)[1][i]*X[i] for i in range(len(X))])])
    
def D_tilde(t=None, alpha=None, S=None):
    
    unity = np.array(((1., 0.), (0., 1.)))
    
    matrix = unity - t*R(alpha)
    
    return matrix_product(M=matrix, X=S)


def r(theta=None, alpha=None, S=np.array([1., 0.])):
    """Compute radial distance from planet center of a point defined by the
    azimuth theta and the incident angle alpha
    """
    sx, sy = S
    D = np.sqrt(sx**2 + sy**2)
    r = abs(D/np.cos(theta)/(1.+np.tan(theta)/np.tan(alpha)))
    return r


def GetIndex(r=None, theta=None, rs=None, n=2):
    
    thetas_dummy = np.linspace(0, np.pi, 2**n - 2**(n-1) + 1)
    delta_theta = (thetas_dummy[-1]-thetas_dummy[0])/(len(thetas_dummy)-1)
    delta_r = (rs[-1] - rs[0])/(len(rs)-1)
    i = int((theta-thetas_dummy[0]+delta_theta/2.)/delta_theta\
        -(abs(theta)/theta -1)*2**(n-1))
    
    j = int((r-rs[0]+delta_r)/delta_r)
    return i, j


def alpha_incident(r=None, theta=None, S=np.array([1., 0.])):
    """Compute the incident angle alpha as a function of the polar coordinates
    in the planets frame
    """
    sx, sy = S
    D = np.sqrt(sx**2 + sy**2)
    return np.arctan(np.tan(theta)/(D/(r*np.cos(theta))-1.))


def CartToPolar(x):
    """Convert vector x in cartesian coordinates into polar coordinates where
    the angle theta goes from - pi to pi
    """
    radius = np.sqrt(sum(x[i]**2 for i in range(len(x))))
    
    if not radius == 0. and not x[1] == 0.:
        theta = np.arccos(x[0]/radius)*x[1]/abs(x[1])
        
    else:
        if not x[0] == 0:
            theta = np.pi * x[0]/abs(x[0])
            
        else:
            theta = 0.
        
    return np.array([theta, radius])


def CrossIntersect(Q1 = None, Q2 = None, P1=None, P2=None):
    """Computes intersection point between straight line Q1-Q2 and P1-P2
    """
    
    K1 = P1[1]-Q1[1] + (P2[1]-P1[1])*(Q1[0]-P1[0])/(P2[0]-P1[0])
    K2 = Q2[1]-Q1[1] - (P2[1]-P1[1])*(Q2[0]-Q1[0])/(P2[0]-P1[0])
    tQ = K1/K2
    return Q1 + tQ*(Q2-Q1)


def t12(r=None, alpha=None, S=None):
    
    RS = rotation(alpha, S)
    
    SS = np.sqrt(sum(S**2))
    
    RSS = vector_product(X=RS, Y=S)
    
    RSRS = vector_product(X=RS, Y=RS)
    
    return [(RSS + np.sqrt(RSS**2 - RSRS*(SS**2 - r**2)))/RSRS,
            (RSS - np.sqrt(RSS**2 - RSRS*(SS**2 - r**2)))/RSRS]


def Qs(alpha=None, r=None, S=None):
    """Compute coordinates of intersection points between circle of radius r
    and straight line at anlge alpha originating at S
    """
    t1, t2 = t12(r=r, alpha=alpha, S=S)
    Q1 = D_tilde(t=t1, alpha=alpha, S=S)
    Q2 = D_tilde(t=t2, alpha=alpha, S=S)
    return Q1, Q2


def Propagate(phi2=None, r=None, n1=None, n2=None,  S=None, offset=None, signsign=None):
    Q1, Q2 = Qs(alpha = phi2, r = r, S = S)
    #If phi2 is zero that means that the ray is along the x-axis. To avoid
    #errors the theta, delta_theta and sign must be set manually here as they
    #are simply 0, pi and 1 anyways
    if phi2 == 0:
        delta_theta = np.pi
        sign = 1.
        theta = 0.
        
    else:
        #Account for spurious numerical behaviour for Q/r close to 1. If it is 
        #1 numerical errors can lead it to be slightly greater than 1 which will
        #give nan in arccos. Therefore it has to be forced to remain bounded
        #between -1 and 1
        argument = Q2[0]/r
        
        argument = min(abs(argument), 1.)*argument/abs(argument)
        theta = np.arccos(argument)*Q2[1]/abs(Q2[1])
        delta_theta = np.arccos(Q1[0]/r)*Q1[1]/abs(Q1[1])
        
        sign = abs(phi2)/phi2
        
    phi1 = sign*(abs(theta)*signsign + abs(phi2) - offset)
    phi2 = np.arcsin(n1/n2*np.sin(phi1))
    return Q1, Q2, theta, delta_theta, phi1, phi2


def PropagateFull(phi2=None, r=None, n1=None, n2=None,  S=None, offset=None, signsign=None):
    Q1, Q2 = Qs(alpha = phi2, r = r, S = S)
    #If phi2 is zero that means that the ray is along the x-axis. To avoid
    #errors the theta, delta_theta and sign must be set manually here as they
    #are simply 0, pi and 1 anyways
    if phi2 == 0:
        delta_theta = np.pi
        sign = 1.
        theta = 0.
        
    else:
        #Account for spurious numerical behaviour for Q/r close to 1. If it is 
        #1 numerical errors can lead it to be slightly greater than 1 which will
        #give nan in arccos. Therefore it has to be forced to remain bounded
        #between -1 and 1
        argument = Q2[0]/r
        
        argument = min(abs(argument), 1.)*argument/abs(argument)
        theta = np.arccos(argument)*Q2[1]/abs(Q2[1])
        delta_theta = np.arccos(Q1[0]/r)*Q1[1]/abs(Q1[1])
        
        sign = abs(phi2)/phi2
        
    phi1 = sign*(abs(theta)*signsign + abs(phi2) - offset)
    phi2 = np.arcsin(n1/n2*np.sin(phi1))
    return Q1, Q2, theta, delta_theta, phi1, phi2
    


def Qs_refracted(alphas=None, rs=None, ns=None, S=None, forbidden=False,
                 n_pizza = 3):
    """Compute for a set of circles around the origin the intersection points 
    with a straight line at angle alpha originating at S while changing the
    the angle of S at each intersection point. Returns an array containing
    the coordinates of the nodes (intersections of ray and circles) 
    for all alphas
    """
    
    alpha_min = abs(np.arctan(rs[0]/S[0]))
    alpha_max = abs(np.arctan(rs[-1]/S[0]))
    
    Nx = 2**n_pizza
    Ny = len(rs)
    
    x_axis = np.linspace(-np.pi, np.pi, 2**n_pizza+1)
    y_axis = rs
    delta_x = (x_axis[-1]-x_axis[0])/(len(x_axis)-1)
    
    thetas = np.linspace(0, 2*np.pi, 2**n_pizza+1)
    
    xx, yy = np.meshgrid(x_axis, y_axis)
    angle_correct = 1.0+1.0e-8
    intersect_limit = 2**n_pizza
    
    #paths
    data = [[S] for a in range(len(alphas))]
    nodes = [[] for a in range(len(alphas))]
    for a  in range(len(alphas)):
        print ('\n------------------------------')
        print ('alpha = ', alphas[a]/np.pi)
        alpha = alphas[a]
    
        if abs(alpha) >= alpha_max:
            far = S[0]*5*np.array([np.cos(alpha), np.sin(alpha)])
            data[a].append(S-far)
            
        else:
            signsign = 1
            offset = 0.
            blabla = False
            surface_hit = False
            node_count = 0
    
            Q2 = S
            Q1 = np.zeros([2])
            phi2 = alpha
            phi1 = alpha
            theta = 0.
            theta_tot = 0.
            Q = Q2
            index = len(rs) - 1

            #Gather front points
            for i in range(len(rs)):
                node_count += 1
                iterator = len(rs) - 1 - i
                index = abs(iterator)
                phi2_old = phi2
                
                #Compute next intersection points Q1 and Q2 and corresponding
                #refraction angles phi1 and phi2 aswell as the new polar angle
                Q1, Q2, theta, delta_theta, phi1, phi2 = Propagate(phi2=phi2,
                                                      r=rs[index],
                                                      n1=ns[index+1],
                                                      n2=ns[index],
                                                      S=Q,
                                                      offset=offset,
                                                      signsign=signsign)
                
                #Rotate Q2 onto x-axis
                Q2 = rotation(alpha=-theta, X=Q2)
                
                if i == len(rs)-1:
                    if abs(alpha) < alpha_max:
                        blabla = True
                        
                    stop_index = index
                    offset = np.pi
                
                if np.isnan(Q2[0]) or np.isnan(Q2[1]):
                    if abs(alpha) < alpha_max:
                        blabla = True
    
                    phi2 = phi2_old
                    stop_index = index+1
                    offset = np.pi
                    
                else:             
                    #Update theta to track the total rotation that has been performed
                    theta_tot += theta
                    Q=Q2
                    Q_plot = rs[index]*np.array([np.cos(theta_tot), np.sin(theta_tot)])
                    
                    data[a].append(Q_plot)
                    
                    polar1 = CartToPolar(data[a][-2])
                    polar2 = CartToPolar(data[a][-1])
                    node_iteration = True
                    k=0
                    
                    while node_iteration:
                        k+= 1
                        grid_indices1 = GetIndex(rs = rs,
                                                theta = polar1[0],
                                                r = polar1[1],
                                                n=n_pizza)
                        
                        grid_indices2 = GetIndex(rs = rs,
                                                theta = polar2[0],
                                                r = polar2[1],
                                                n=n_pizza)
                        
                        if grid_indices1[0] == grid_indices2[0] or polar1[1] > \
                        rs[-1] or k >= intersect_limit :
                            node_iteration = False
                            
                        else:
                            
                            pizza_angle = thetas[grid_indices1[0]]
                            
                            if not abs(pizza_angle) == 0.0:
                                pizza_angle += delta_x/2. * \
                                abs(pizza_angle)/pizza_angle*polar1[0]/\
                                abs(polar1[0])
                                
                            else:
                                pizza_angle += delta_x/2. *polar1[0]/abs(polar1[0])
                            node = CrossIntersect(Q1=data[a][-1], Q2=data[a][-2], 
                                                  P1=np.zeros([2]),
                                                  P2 = 10*np.array([np.cos(pizza_angle), 
                                                                    np.sin(pizza_angle)]))
            
                            nodes[a].append(node)
                            
                            if not np.isnan(node[0]):
                                polar1 = CartToPolar(node)
                                polar1[0] *= angle_correct
                                
                    #If surface is hit, do not proceed with the ray propagation
                    if index == 0:
                        surface_hit = True
                    
                if blabla and not surface_hit:
                    #Compute next intersection points Q1 and Q2 and corresponding
                    #refraction angles phi1 and phi2 aswell as the new polar angle
                    Q1, Q2, theta, delta_theta, phi1, phi2 = Propagate(phi2=phi2,
                                                          r=rs[stop_index],
                                                          S=Q,
                                                          n1=ns[stop_index],
                                                          n2=ns[stop_index+1],
                                                          offset=offset,
                                                          signsign=signsign)
                    
                    #Rotate Q1 onto x-axis
                    Q1 = rotation(alpha=-delta_theta, X=Q1)
    
                    Q1_plot = rotation(alpha=theta_tot+delta_theta, X=Q1)
                    Q2_plot = rotation(alpha=theta_tot, X=Q2)
    
                        
                    data[a].append(Q2_plot)
                    data[a].append(Q1_plot)
                    
                    #Update total theta angle
                    theta_tot += delta_theta
                    
                    theta = delta_theta
                    
                    polar1 = CartToPolar(Q2_plot)
                    polar2 = CartToPolar(Q1_plot)
                        
                    node_iteration = True
                    k=0
                    while node_iteration:
                        k+= 1
                        grid_indices1 = GetIndex(rs = rs,
                                                theta = polar1[0],
                                                r = polar1[1],
                                                n=n_pizza)
                        
                        grid_indices2 = GetIndex(rs = rs,
                                                theta = polar2[0],
                                                r = polar2[1],
                                                n=n_pizza)
                        
                        if grid_indices1[0] == grid_indices2[0] or polar1[1] > rs[-1] or k >= intersect_limit :
                            node_iteration = False
                            
                        else:
                            
                            pizza_angle = thetas[grid_indices1[0]]
                            if not abs(pizza_angle) == 0.0:
                                pizza_angle += delta_x/2. * abs(pizza_angle)/pizza_angle*polar1[0]/abs(polar1[0])
                                
                            else:
                                pizza_angle += delta_x/2. *polar1[0]/abs(polar1[0])
                            
                            node = CrossIntersect(Q1=data[a][-1], Q2=data[a][-2], P1=np.zeros([2]),
                                                  P2 = 10*np.array([np.cos(pizza_angle), 
                                                                 np.sin(pizza_angle)]))
            
                            nodes[a].append(node)
                            
                            if not np.isnan(node[0]):
                                polar1 = CartToPolar(node)
                                polar1[0] *= angle_correct
                    break
            
                    
            print ('continuing on the other side...')
            #Rotate Q1 such that initial theta is zero again
            iterator = stop_index + 1
            index = iterator
            Q2 = Q1
            signsign = -1
            offset = 0.

            if not surface_hit:
                #Gather front points
                for i in range(len(rs) - stop_index - 1):
                    node_count += 1
                    index = iterator + i
                    Q = Q2    
                    #Compute next intersection points Q1 and Q2 and corresponding
                    #refraction angles phi1 and phi2 aswell as the new polar angle
                    Q1, Q2, theta, delta_theta, phi1, phi2 = Propagate(phi2=phi2,
                                                          r=rs[index],
                                                          n1=ns[index],
                                                          n2=ns[index+1],
                                                          S=Q,
                                                          offset=offset,
                                                          signsign=signsign)
                    
                    #Rotate Q2 onto x-axis
                    Q2 = rotation(alpha=-theta, X=Q2)
        
                    theta_tot += theta

                    
                    Q_plot = rs[index]*np.array([np.cos(theta_tot), np.sin(theta_tot)])
                    data[a].append(Q_plot)
                    
                    polar1 = CartToPolar(data[a][-2])
                    polar2 = CartToPolar(data[a][-1])
                    
                    node_iteration = True
                    k=0
                    while node_iteration:
                        k+= 1
                        grid_indices1 = GetIndex(rs = rs,
                                                theta = polar1[0],
                                                r = polar1[1],
                                                n=n_pizza)
                        
                        grid_indices2 = GetIndex(rs = rs,
                                                theta = polar2[0],
                                                r = polar2[1],
                                                n=n_pizza)
                        
       
                        if grid_indices1[0] == grid_indices2[0] or polar1[1] > rs[-1] or k >= intersect_limit :
                            node_iteration = False
                            
                        else:
                       
                            
                            pizza_angle = thetas[grid_indices1[0]]
                            if not abs(pizza_angle) == 0.0:
                                pizza_angle += delta_x/2. * abs(pizza_angle)/pizza_angle*polar1[0]/abs(polar1[0])
                                
                            else:
                                pizza_angle += delta_x/2. *polar1[0]/abs(polar1[0])
                            node = CrossIntersect(Q1=data[a][-1], Q2=data[a][-2], P1=np.zeros([2]),
                                                  P2 = 10*np.array([np.cos(pizza_angle), 
                                                                 np.sin(pizza_angle)]))
            
                            nodes[a].append(node)
                            
                            if not np.isnan(node[0]):
                                polar1 = CartToPolar(node)
                                polar1[0] *= angle_correct
                        
                    #Plot incident point at detector
                    
                Q_add = rotation(alpha=theta, X=Q1)
                
                #Add point at a given distance from sphere in current frame
                Q_add = Q2 + 1.*np.array([np.cos(phi2), np.sin(phi2)])
                
                #Rotate point into original frame
                Q_add = rotation(alpha=theta_tot, X=Q_add)
                data[a].append(Q_add)
   
    
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.set_xlabel('x [a.u.]')
    ax.set_ylabel('y [a.u.]')
    ax.set_xlim(-.6, .6)
    ax.set_ylim(-.6, .6)
    '''
    cmap = plt.get_cmap('jet', 6)
    norm = mpl.colors.Normalize(vmin=.9,vmax=1.3)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ticks=ns, 
                 boundaries=np.linspace(1., 1.2, 4))
    '''
    c = np.linspace(1., 1.2, 7)
    
    colors = [(1, 1, 1), (0, 0, 1)]  # R -> G -> B
    n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins
    cmap_name = 'my_list'
    
    # Create the colormap
    cm = LinearSegmentedColormap.from_list(
        cmap_name, colors, N=7)
    # Fewer bins will result in "coarser" colomap interpolation
    #im = ax.imshow(Z, interpolation='nearest', origin='lower', cmap=cm)
        
    # Make dummie mappable
    dummie_cax = ax.scatter(c, c, c=c, cmap=cm)
    # Clear axis
    #ax.cla()
    cbar = fig.colorbar(dummie_cax, ticks=c, boundaries = np.linspace(1., 1.2, 8))
    cbar.set_label('Refractive Index')
    cbar.ax.get_yaxis().labelpad = 15
    ax.scatter(S[0], S[1], color='k', s=50, marker = 'X')
    ax.plot([-10, 10], [0, 0], color='k', zorder=300)
    ax.plot([0, 0], [-10, 10], color='k', zorder=300)
    
    
    #Plot all rays
    for a in range(len(data)):
        for i in range(len(data[a])-1):
            #ax.scatter(data[a][i][0], data[a][i][1], marker='o', color='k', 
             #          facecolor='None')
            ax.plot([data[a][i][0], data[a][i+1][0]], [data[a][i][1], 
                     data[a][i+1][1]], color = 'y', zorder=200,
                    marker = 'x', markersize=4, markeredgecolor='k')
        
        for i in range(len(nodes[a])):
            ax.scatter(nodes[a][i][0], nodes[a][i][1], color='k', s=5, zorder=300)
    
    for i in range(len(rs)):
        if i == len(rs)-1:
            f = (0., 0.3, 0.3)
        else:
            f = (.9-i/len(rs), .9-i/len(rs), 1.)
   
        circle1 = plt.Circle([0., 0.], rs[-i-1], color=f, fill=f,
                            zorder=0)
        circle2 = plt.Circle([0., 0.], rs[-i-1], color='w', fill=False,
                            zorder=0, linewidth=1)
        ax.add_artist(circle1)
        ax.add_artist(circle2)
    
    
    angle_min = 2*np.pi/2**n_pizza/2
    angle_max = 2*np.pi + angle_min
    angles = np.linspace(angle_min, angle_max, 2**n_pizza+1)
    for i in range(2**n_pizza+1):
        outer = rs[-1]*np.array([np.cos(angles[i]), np.sin(angles[i])])
        ax.plot([0., outer[0]], [0., outer[1]], color='w', linewidth=1)
        
    return data


def a(alpha=None, r=None, R_atmos=None, S=np.array([1., 1.])):
    """
    """
    QR = Qs(alpha=alpha, r=R_atmos, S=S)
    Qr = Qs(alpha=alpha, r=r, S=S)
    
    a =  np.sqrt(sum([(QR[1] - Qr[1])[i]**2 for i in range(len(Qr))]))
    
    b = np.sqrt(sum([(Qr[1] - S)[i]**2 for i in range(len(Qr))]))
   
    if np.isnan(a):
        a = np.sqrt(sum([(QR[1] - Qr[1])[i]**2 for i in range(len(Qr))]))
        b = np.sqrt(sum([(Qr[1] - S)[i]**2 for i in range(len(QR))])) 

    bb = [np.sqrt(sum([(QR[1] - QR[0])[i]**2 for i in range(len(QR))])),
          np.sqrt(sum([(QR[1] - Qr[0])[i]**2 for i in range(len(QR))])),
          np.sqrt(sum([(QR[1] - Qr[1])[i]**2 for i in range(len(QR))]))]
        
    return a, b, bb


def path_length(r=None, theta=None, R_atmos = None, R_planet=None,
                S=np.array([1., 0.])):
    """Compute travelled distance of a ray through the atmosphere at a given
    point defined the polar coordinates r and theta. r and theta can be chosen
    arbitrarily and do not have to lie on a atmospheric layer
    """
    #Compute distance between planet and star
    D = np.sqrt(sum([S[i]**2 for i in range(len(S))]))
    
    alpha_min = np.arcsin(R_planet/D)
    #Compute incident angle of ray passing through point r, theta
    alpha = alpha_incident(r=r, theta=theta, S=S)
    #Compute the path length at this alpha through the outermost atmosphere
    #layer
    t1_atmos, t2_atmos = t12(r=R_atmos, alpha=alpha, S=S)
    Q1_atmos = D_tilde(t=t1_atmos, alpha=alpha, S=S)
    Q2_atmos = D_tilde(t=t2_atmos, alpha=alpha, S=S)

    pl_atmos = np.sqrt(sum([(Q2_atmos[k] - Q1_atmos[k])**2
                        for k in range(len(Q2_atmos))]))    
    #Compute intercection points between incident ray and circle around
    #center of the planet with radius r
    t1, t2 = t12(r=r, alpha=alpha, S=S)
    
    #Compute path length of ray through layer of radius r
    Q1 = D_tilde(t=t1, alpha=alpha, S=S)
    Q2 = D_tilde(t=t2, alpha=alpha, S=S)
    pl = np.sqrt(sum([(Q2[k] - Q1[k])**2
                        for k in range(len(Q1))]))      
    
    #if theta > 3pi/2 always take the closer intersection point
    if theta > np.pi:    
        if 2*np.pi-theta + abs(alpha) >= np.pi/2.:
            ipl = (pl_atmos - pl)/2. + pl   

        else:
            ipl = (pl_atmos - pl)/2.
            
    else:
        if theta + abs(alpha) >= np.pi/2.:
            ipl = (pl_atmos - pl)/2. + pl             

        else:
            ipl = (pl_atmos - pl)/2.
             
    if theta > np.pi/2. and theta < 3.*np.pi/2. and abs(alpha) < abs(alpha_min):
        return 1.
    
    else:
        #print (ipl)
        return ipl


def PlotRays(M=[0.,0.], S=np.array([1., 1.]), rs=[.125, .25, .5], 
             alphas = np.linspace(-np.pi/4, np.pi/4, 5),
             ):
    
    t1s = np.empty([len(alphas), len(rs)])
    t2s = np.empty([len(alphas), len(rs)])
    Q1s = np.empty([len(alphas), len(rs), 2])
    Q2s = np.empty([len(alphas), len(rs), 2])
    rays = np.empty([len(alphas), 2, 2])
    circles = []
    
    R_atmos = min(rs)
    
    #Compute distance between planet and star
    D = np.sqrt(sum([S[i]**2 for i in range(len(S))]))
    
    alpha_max = abs(np.arcsin(R_atmos/D))
    
    for i in range(len(alphas)):
        alpha = alphas[i]
        for j in range(len(rs)):
            r = rs[j]
            t1s[i][j], t2s[i][j] = t12(r=r, alpha=alpha, S=S)
            Q1s[i][j] = D_tilde(t=t1s[i][j], alpha=alpha, S=S)
            Q2s[i][j] = D_tilde(t=t2s[i][j], alpha=alpha, S=S)
        
        Qneg = D_tilde(t=-10, alpha=alpha, S=S)
        Qpos = D_tilde(t=10, alpha=alpha, S=S)
        
        rays[i] = Qneg, Qpos
    
    fig, ax = plt.subplots()
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_aspect('equal')
    ax.set_xlabel('x [a.u.]')
    ax.set_ylabel('y [a.u.]')
    marker_size = 20
    
    #ax.scatter(M[0], M[1], color='k')
    ax.scatter(S[0], S[1], color='y')
    
    for i in range(len(rs)):    
        if i == 0:
            col = (0, 0, .5)
            fill = True
            
        else:
            col ='k'
            fill = False

        circles.append(plt.Circle(M, rs[i], color=col, fill=fill))
        
    for circ in circles:
        ax.add_artist(circ)
    
    for i in range(len(alphas)):
        for j in range(len(rs)):
            print (Q1s[i][j])
            if abs(alphas[i]) < alpha_max*0.:
            
                #ax.scatter(Q1s[i][j][0], Q1s[i][j][1], color=(.75, .75, .75))
                ax.scatter(Q2s[i][j][0], Q2s[i][j][1], color=(.25, .25, .25),
                           s=marker_size)
                
            else:
                ax.scatter(Q1s[i][j][0], Q1s[i][j][1], color=(.75, .75, .75),
                           s=marker_size)
                
                ax.scatter(Q2s[i][j][0], Q2s[i][j][1], color=(.25, .25, .25),
                           s=marker_size)
                pass
            
        ax.plot([rays[i][0][0], rays[i][1][0]], 
                [rays[i][0][1], rays[i][1][1]],
                color='y')
    
    
def PlotMap(M=[0.,0.], S=np.array([1., 0]), R_planet = .25,
         R_atmos=.5):
    
    #Compute distance between planet and star
    D = np.sqrt(sum([(M[i]-S[i])**2 for i in range(len(S))]))
    
    alpha_max = abs(np.arcsin(R_atmos/D))
    
    N = 2
    N_radii = 2
    alphas = np.linspace(-np.pi/4, np.pi/4, N)
    radii = np.linspace(R_planet, R_atmos, N_radii)
    thetas = np.linspace(-np.pi, np.pi, 32)

    t1s = np.empty([len(alphas), len(radii)])
    t2s = np.empty([len(alphas), len(radii)])
    Q1s = np.empty([len(alphas), len(radii), 2])
    Q2s = np.empty([len(alphas), len(radii), 2])
    
    data = np.zeros([len(alphas), len(radii)])
    cdf = np.zeros([len(alphas), len(radii)])
    
    xx, yy = np.meshgrid(alphas, radii)

    #Compute distance between planet and star
    D = np.sqrt(sum([S[i]**2 for i in range(len(S))]))
    
    alpha_max = abs(np.arcsin(R_atmos/D))
    alpha_min = abs(np.arcsin(R_planet/D))
    
    pl = np.zeros([len(alphas), N_radii])
    ipl = np.zeros([len(alphas), len(radii), 2])
    
    for i in range(len(alphas)):
        alpha = alphas[i]
        for j in range(len(radii)):
         
            t1s[i][j], t2s[i][j] = t12(r=radii[j], alpha=alpha, S=S)
            if abs(alpha) > abs(alpha_min)*0:
            
                Q1s[i][j] = D_tilde(t=t1s[i][j], alpha=alpha, S=S)
            
            else:
                Q1s[i][j] = np.array([None, None])

            Q2s[i][j] = D_tilde(t=t2s[i][j], alpha=alpha, S=S)
            pl[i][j] = np.sqrt(sum([(Q2s[i][j][k] - Q1s[i][j][k])**2
                        for k in range(len(Q2s[i][j]))]))
            
    
    #Compute individual path lengths
    for i in range(len(alphas)):
        for j in range(len(radii)):
            for k in range(2):
                ipl[i][j][k] = (pl[i][-1] - pl[i][j])/2. + pl[i][j]*(1.-k)
        
    for i in range(len(alphas)):
        for j in range(len(radii)):
            alpha = alphas[i]
            radius = radii[j]
            
            aa, bb, bbb = a(alpha=alpha, r=radius, R_atmos=R_atmos, S=S)

            cdf[i][j] = aa
            data[i][j] = 1.- aa            
   
    Nx_ticks = 5
    Ny_ticks = 6
    
    x_ticks = np.linspace(0, N, Nx_ticks) 
    y_ticks = np.linspace(0, N, Ny_ticks)
    
    dx = (max(alphas) - min(alphas))/(Nx_ticks-1)
    dy = (max(radii) - min(radii))/(Ny_ticks-1)
    
    x_tick_labels = [(min(alphas)+dx*i)/np.pi for i in range(Nx_ticks)]
    y_tick_labels = [((min(radii)+dy*j)) for j in range(Ny_ticks)]
    
    fig, ax = plt.subplots()
    
    im = ax.imshow(data.T, origin='lower')

    divider=make_axes_locatable(ax)
    cax=divider.append_axes('right', size='5%', pad=.5)
    cbar=fig.colorbar(im, cax=cax)
    cbar.set_label(r'$I(\alpha, r)/I_0$')
    
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    
    ax.set_xticklabels(x_tick_labels)
    ax.set_yticklabels(y_tick_labels)    

    ax.set_xlabel(r'$\alpha \ \rm [rad]$')
    ax.set_ylabel(r'$\rm Radius \ [a.u.]$')
    
    fig, ax2 = plt.subplots()
    ax2.set_xlim(-1, 1)
    ax2.set_ylim(-1, 1)
    
    for i in range(len(alphas)):
        for j in range(len(thetas)):
            if thetas[j]/abs(thetas[j]) == alphas[i]/abs(alphas[i]):
            
                radius = r(theta=thetas[j], alpha=alphas[i], S = S)
                
                x, y = radius*np.array([np.cos(thetas[j]), np.sin(thetas[j])])
            
                #color = (0, data[i][j]/np.nanmax(data), data[i][j]/np.nanmax(data))
                ax2.scatter(x, y, color='k')
                ax2.set_aspect('equal')

    fig, ax3 = plt.subplots()
    ax3.set_xlim(-1, 1)
    ax3.set_ylim(-1, 1)
    ax3.set_aspect('equal')

    ipl_max = np.nanmax(ipl)    
    for i in range(len(alphas)):

        for j in range(len(radii)):
            alpha = alphas[i]
            
            Qs = [Q1s[i][j], Q2s[i][j]]
            for k in range(2):
                if k == 0 and abs(alpha) <= abs(alpha_min):
                    pass
                
                else:
                    color = (0, 1-ipl[i][j][k]/ipl_max, 1-ipl[i][j][k]/ipl_max)
                    ax3.scatter(Qs[k][0], Qs[k][1], color=color)

    for i in range(len(radii)):    
        if i == 0:
            col = (0, 0, .5)
            fill = True
            
        else:
            col ='k'
            fill = False

        circ = plt.Circle(M, radii[i], color=col, fill=fill)
        
        ax3.add_artist(circ)

    #Compute path length for given theta and radius

    radii = np.linspace(R_planet, R_atmos, 4)
    thetas = np.linspace(0, 2*np.pi, 4)
    
    pl = np.zeros([len(radii), len(thetas)])
    x = np.zeros([len(radii), len(thetas)])
    y = np.zeros([len(radii), len(thetas)])

    for i in range(len(radii)):
        for j in range(len(thetas)):
            pl[i][j] = path_length(r=radii[i], theta=thetas[j],
                       R_atmos=R_atmos, R_planet=R_planet, S=S)
    
            x[i][j] = radii[i] * np.cos(thetas[j])
            y[i][j] = radii[i] * np.sin(thetas[j])       
            
    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    for i in range(len(radii)):
        for j in range(len(thetas)):
            col = (0., 1.-pl[i][j]/np.nanmax(pl), 1.-pl[i][j]/np.nanmax(pl))
            
            ax.scatter(x[i][j], y[i][j], color=col)
            
            
    x = np.linspace(-1, 1, 128)
    y = np.linspace(-1, 1, 128)
    
    xx, yy = np.meshgrid(x, y)
    
    data = np.zeros([len(x), len(y)])
    
    for i in range(len(x)):
        for j in range(len(y)):
            radius = np.sqrt(x[i]**2 + y[j]**2)
            theta = np.arccos(x[i]/radius)
            if radius <= R_atmos:
                data[i][j] = np.exp(-path_length(r=radius, 
                                                theta=theta, 
                                                R_atmos=R_atmos,
                                                R_planet=R_planet, S=S))
            
    fig, ax = plt.subplots()
    
    im = ax.imshow(data.T, origin='lower')
    ax.set_aspect('equal')
        
    