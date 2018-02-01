#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 14:38:59 2018

@author: antariksh
"""
import numpy as np

def axis_main(Ns, d_x, d_y, d_z):
    """
    Define the main axis as a discretized curve.
    
    Args:
        d_z (float): the Z-position of the final point on the axis
        Ns (int) : number of sections on the axis.
        d_x (float) : the x-ordinate of the final point on the axis.
        d_y (float) : the y-cordinate of the final point on the axis.
    
    Returns:
        axis : (Ns+1) floating point numpy array consisting of the vertices
                   of the main axis in the body co-ordinate system.
        ax_length: Length of axis           
                   
    """
    # initialize the size of the axis [xi, yi, zi]
    axis = np.zeros((Ns + 1 , 3), dtype = float)
    # 
    delta_x = np.linspace(0, d_x, num = Ns + 1, endpoint = True)
    #
    delta_y = np.linspace(0, d_y, num = Ns + 1, endpoint = True)
    #
    delta_z = np.linspace(0, d_z, num = Ns + 1, endpoint = True)
    #
    axis[:, 0] = delta_x
    axis[:, 1] = delta_y
    axis[:, 2] = delta_z

    return axis

def circle(Nc, radius = 1, flag = 0):
    """
     Generate a circle in the non-dimensional space. The scaling is the radius.
     
     Args:
         Nc (int) : number of cross-sections
         
     Returns:
         surf_2d : A float numoy array with size Nc x 2 containing the 
                   co-ordinates of the profile in non-dimensional coordinates.
             
    """
    # initialize surface array
    surf_2d = np.zeros((Nc, 2), dtype = float)
    # generate a circular profile in non-dimensional co-ordinate system
    start_theta = 0
    end_theta = 2*np.pi - 0.001*np.pi
    theta_vec= np.linspace(start_theta, end_theta, num = Nc, endpoint = True)

    # obtian the x,y vectors
    x = radius*np.cos(theta_vec)
    y = radius*np.sin(theta_vec)
    
    # surf_2d
    #
    surf_2d[:, 0] = x
    #
    surf_2d[:, 1] = y
    
    # plot check
    if flag:
        from matplotlib import pyplot as plt
        fig = plt.figure('cross_section')
        plt.title('circle')
        plt.plot(x, y, 'xb-')
        plt.xlabel('x')
        plt.ylabel('y')
        
    return surf_2d

def length(axis):
    """
       Returns the length of the axis covered until each of the vertex points
       
       Args:
           axis (float) : main axis vertices in body co-ordinate system
           
       Returns:
           length : Floating point numpy array of size Ns with lengths of axis
                    until the corresponding vertex
    """    
    # 
    length = np.sqrt(np.power(axis[1:,0] - axis[0,0], 2) + 
                     np.power(axis[1:,1] - axis[0,1], 2) +
                     np.power(axis[1:,2] - axis[0,2], 2))
    
    return length