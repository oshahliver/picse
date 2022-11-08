#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 12:38:17 2020

@author: oshah
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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


def alpha_incident(r=None, theta=None, S=np.array([1., 0.])):
    """Compute the incident angle alpha as a function of the polar coordinates
    in the planets frame
    """
    sx, sy = S
    D = np.sqrt(sx**2 + sy**2)
    return np.arctan(np.tan(theta)/(D/(r*np.cos(theta))-1.))


def t12(r=None, alpha=None, S=None):
    
    RS = rotation(alpha, S)
    
    SS = np.sqrt(sum(S**2))
    
    RSS = vector_product(X=RS, Y=S)
    
    RSRS = vector_product(X=RS, Y=RS)
    
    #print ('RS:', RS)
    #print ('SS:', SS)
    #print ('RSS:', RSS)
    #print ('RSRS:', RSRS)
    
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


def Qs_refracted(alphas=None, rs=None, ns=None, S=None):
    """Compute for a set of circles around the origin the intersection points 
    with a straight line at angle alpha originating at S while changing the
    the angle of S at each intersection point. Returns an array containing
    the coordinates of the nodes (intersections of ray and circles) 
    for all alphas
    """
    fig, ax = plt.subplots()
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_aspect('equal')
    
    #paths
    data = [[S] for a in range(len(alphas))]
    
    Q1_all = np.zeros([len(alphas), len(rs)+1, 2])
    Q2_all = np.zeros([len(alphas), len(rs)+1, 2])
    
    ax.scatter(S[0], S[1], color='k', s=50, marker = 'X')
    ax.plot([-10, 10], [0, 0], color='grey')
    ax.plot([0, 0], [-10, 10], color='grey')
    
    for i in range(len(rs)):
        circle = plt.Circle([0., 0.], rs[i], color='k', fill=False)
        ax.add_artist(circle)    
    
    for a  in range(len(alphas)):
        
        stop_this = False
        node_count = 0
        alpha = alphas[a]
        print ('\nalpha =', alpha/np.pi)
        #Q1, Q2 = Qs(alpha=alpha, r=rs[-1], S=S)
        #circle = plt.Circle([0., 0.], rs[-1], color='k', fill=False)
        #ax.add_artist(circle)
        #ax.scatter(Q2[0], Q2[1], color=color)
        #data[a][0] = Q1, Q2
        Q2 = S
        Q1 = np.zeros([2])
        phi2 = alpha
        phi1 = alpha
        theta = 0.
        theta_tot = 0.
        
        #Gather front points
        for i in range(len(rs)):
            print ()
            node_count += 1
            iterator = len(rs) - 1 - i
            index = abs(iterator)
            print ('index =', index)
            Q = Q2
    
            #Compute first intersection point with largest circle
            circle = plt.Circle([0., 0.], rs[index], color='k', fill=False)
            ax.add_artist(circle)
            
            #Rotate ray along first rotation angle around origin at Q2 (closer point)
            #in fram of Q2. The direction vector of the rotated ray will be the same
            #in the planet frame but with a shifted origin

            print ('phi1 =', phi1/np.pi)
            print ('phi2 =', phi2/np.pi)                           
            print ('theta before =', theta/np.pi)
            
            #With Q2 as S find now the intersection with the second circle
            Q1, Q2 = Qs(alpha = phi2, r = rs[index], S = Q)
            print ('Q before =', Q1, Q2) 
            Q1_all[a][i] = Q1
            Q2_all[a][i] = Q2
            blabla = True
            
            if np.isnan(Q2[0]) or np.isnan(Q2[1]):
                print ('index =', index)
                blabla = False
                stop_this = True
                stop_i = i-1
                stop_index = index+1
                Q1, Q2 = Qs(alpha = phi2, r = rs[index+1], S = Q)
                Q1_all[a][i] = Q1
                Q2_all[a][i] = Q2
                Q1_plot = rotation(alpha=theta_tot, X=Q1)
                Q2_plot = rotation(alpha=theta_tot, X=Q2)
                ax.scatter(Q1_plot[0], Q1_plot[1], color='b', marker='x')
                ax.scatter(Q2_plot[0], Q2_plot[1], color='k', marker='o')
                data[a].append(Q2_plot)
                data[a].append(Q1_plot)
                print ('Q new =', Q1, Q2) 
                
                #Compute the angular separation between primary and secondary
                #intersection points
                delta_theta = np.arccos(Q1[0]/rs[index+1])*Q1[1]/abs(Q1[1])
                print ('delta theta=', delta_theta/np.pi)
                #The incident angle is positive/negative for negative/positive y
                #coordinates
                sign = abs(phi2)/phi2#phi2/abs(phi2)*(np.pi/2.-abs(theta))/abs(np.pi/2.-abs(theta))
                phi1 = -sign*(np.pi - abs(delta_theta) - abs(phi2))
                
                phi2 = np.arcsin(ns[index+1]/ns[index+2]*np.sin(phi1))
                
                Q1 = rotation(alpha=-delta_theta, X=Q1_all[a][i])
                print ('ns =', ns[index+1], ns[index+2])
                print ('phi1 =', phi1/np.pi)
                print ('phi2 =', phi2/np.pi)
                #Add the net-theta by which it has been rotated with respect to
                #the previous step to the total theta to keep track of all the
                #rotation so far performed
                theta_tot += delta_theta
                #data[a].append(Q2_plot)
                #data[a].append(Q1_plot)
                
                theta = delta_theta
                break
                
            if i == len(rs)-1:
                print ('check1')
                stop_this = True
                stop_i = i
                stop_index = index
            
            theta = np.arccos(Q2[0]/rs[index])*Q2[1]/abs(Q2[1])
            
            #The incident angle is positive/negative for negative/positive y
            #coordinates
            sign = abs(phi2)/phi2#phi2/abs(phi2)*(np.pi/2.-abs(theta))/abs(np.pi/2.-abs(theta))
            phi1 = sign*(abs(theta) + abs(phi2))
            print ('ns=', ns[index+1], ns[index])
            phi2 = np.arcsin(ns[index+1]/ns[index]*np.sin(phi1))
            print ('phi1 =', phi1/np.pi)
            print ('phi2 =', phi2/np.pi)
            print ('theta after =', theta/np.pi)            
            
            #Rotate Q2 onto x-axis
            Q2 = rotation(alpha=-theta, X=Q2)
            print ('Q2 after =', Q2)
            
            #Update theta to track the total rotation that has been performed
            theta_tot += theta
            print ('theta tot=', theta_tot/np.pi)
            
            Q_plot = rs[index]*np.array([np.cos(theta_tot), np.sin(theta_tot)])
            ax.scatter(Q_plot[0], Q_plot[1], color='k')
            
            #Note that Q2_all has not been rotated yet so the rotation angle
            #is theta_tot minus the current theta
            Q_plot = rotation(alpha=theta_tot-theta, X=Q2_all[a][i])
            data[a].append(Q_plot)
            ax.scatter(Q_plot[0], Q_plot[1], color='g', marker='x')
            
            if stop_this and blabla:
                print ('\nbreak after', iterator)
                print ('stop index =', stop_index)
                print ('stop i =', stop_i)
                print ('Q1 before rot =', Q1_all[a][stop_i])
                print ('theta before =', theta/np.pi)  
                node_count += 1
                
                #To plot the secondary intersection point at the right location
                #in the original reference frame, it has to be rotated back
                #by the total theta angle so far minus the current theta
                #because the current theta_tot already took into account the
                #current theta in the previous iteration

                Q1, Q2 = Qs(alpha = phi2, r = rs[index], S = Q2)
                Q1_all[a][i] = Q1
                Q2_all[a][i] = Q2
                Q1_plot = rotation(alpha=theta_tot, X=Q1)
                Q2_plot = rotation(alpha=theta_tot, X=Q2)
                data[a].append(Q2_plot)
                data[a].append(Q1_plot)
                ax.scatter(Q1_plot[0], Q1_plot[1], color='r', marker='x')
                ax.scatter(Q2_plot[0], Q2_plot[1], color='k', marker='o')
                print ('r =', rs[stop_index])
                print ('ns =', ns[stop_index+1], ns[stop_index])
                print ('Q =', Q1, Q2)
                #Compute the angular separation between primary and secondary
                #intersection points
                delta_theta = np.arccos(Q1[0]/rs[index])*Q1[1]/abs(Q1[1])#-theta_old
                
                print ('delta theta =', (delta_theta)/np.pi)
                
                #The incident angle is positive/negative for negative/positive y
                #coordinates
                sign = abs(phi2)/phi2#phi2/abs(phi2)*(np.pi/2.-abs(theta))/abs(np.pi/2.-abs(theta))
                phi1 = -sign*(np.pi - (abs(phi2) + abs(delta_theta)))
                print ('ns=', ns[stop_index+1], ns[stop_index])
                print ('sin(phi1)=', np.sin(phi1))
                phi2 = np.arcsin(ns[stop_index]/ns[stop_index+1]*np.sin(phi1))
                print ('phi1 =', phi1/np.pi)
                print ('phi2 =', phi2/np.pi)
                
                #Rotate Q1 onto x-axis
                #To rotate into the angle of the secondary intersection point,
                #the the angular separation between primary and secondary point
                #must be added to the previous theta. This angle is now the angle
                #by which Q1 must be rotated back in order to bring it onto
                #the x axis -> theta = 0. The rotatio must be backward and thus
                #-(delta_theta + theta)
                Q1 = rotation(alpha=-delta_theta, X=Q1)
                
                #Add the net-theta by which it has been rotated with respect to
                #the previous step to the total theta to keep track of all the
                #rotation so far performed
                theta_tot += delta_theta
                #data[a].append(Q2_plot)
                #data[a].append(Q1_plot)
                
                theta = delta_theta
                print ('theta after =', theta/np.pi)  
                print ('Q1 after rot =', Q1)
                break            
            #Q_plot = rotation(alpha=theta_tot-theta, X=Q1_all[a][i])
            #ax.scatter(Q_plot[0], Q_plot[1], color='y', marker='x')
            
        print ('continuing on the other side...')
        print ('Q1 is:', Q1)
        #Rotate Q1 such that initial theta is zero again
        iterator = stop_index + 1
        Q2 = Q1
        #Gather front points
        for i in range(len(rs) - stop_index - 1):
            node_count += 1
            print ('\ni =', i)
            index_new = iterator + i
            print ('new index =', index_new)
            Q = Q2
            #Rotate ray along first rotation angle around origin at Q2 (closer point)
            #in fram of Q2. The direction vector of the rotated ray will be the same
            #in the planet frame but with a shifted origin
                       
            print ('theta before =', theta/np.pi)
            #With Q2 as S find now the intersection with the second circle
            print ('phi1 =', phi1/np.pi)
            print ('phi2 =', phi2/np.pi)            
            Q1, Q2 = Qs(alpha = phi2, r = rs[index_new], S = Q)
            print ('Q before =', Q)
            
            print ('Q2 =', Q2)
            theta = np.arccos(Q2[0]/rs[index_new])*Q2[1]/abs(Q2[1])
            
            #The incident angle is positive/negative for negative/positive y
            #coordinates
            sign = phi2/abs(phi2)#phi2/abs(phi2)*(np.pi/2.-abs(theta))/abs(np.pi/2.-abs(theta))
            phi1 = sign*(np.pi - (abs(phi2) - abs(theta)))#np.pi/2.-theta#sign*(abs(theta) + abs(phi2))
            print ('ns =', ns[index_new], ns[index_new+1])
            print ('abs(phi2) =', abs(phi2)/np.pi)
            print ('abs(theta) =', abs(theta)/np.pi)

            phi2 = np.arcsin(ns[index_new]/ns[index_new+1]*np.sin(phi1))
            print ('theta after =', theta/np.pi)            
            print ('phi1 =', phi1/np.pi)
            print ('phi2 =', phi2/np.pi)
            
            #Rotate the new Q2 to be aligned with the x-axis for the next iteration
            Q2 = rotation(alpha=-theta, X=Q2)

            theta_tot += theta
            print ('theta tot=', theta_tot/np.pi)
            
            Q_plot = rs[index_new]*np.array([np.cos(theta_tot), np.sin(theta_tot)])
            data[a].append(Q_plot)
            print ('rs[index]=', rs[index_new])
            ax.scatter(Q_plot[0], Q_plot[1], color='k', facecolor='None')
            print ('Q after =', Q2)
            print ('node_count =', node_count)
            #Plot incident point at detector
            if i == len(rs) - stop_index - 2:
                Q_add = rotation(alpha=theta, X=Q1)
                print ('Q_add =', Q_add)
                #Add point at a given distance from sphere in current frame
                Q_add = Q2 + 1.*np.array([np.cos(phi2), np.sin(phi2)])
                print ('Q_add =', Q_add)
                #Rotate point into original frame
                Q_add = rotation(alpha=theta_tot, X=Q_add)
                print ('Q_add =', Q_add)
                ax.scatter(Q_add[0], Q_add[1], color = 'blue')
                data[a].append(Q_add)
                
        for i in range(len(data[a])-1):
            ax.plot([data[a][i][0], data[a][i+1][0]], [data[a][i][1], data[a][i+1][1]],
                    color = 'y')
            
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
        
    