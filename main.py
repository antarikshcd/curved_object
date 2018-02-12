#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build a curved object out of input profile data defined in local co-ordinate
system at locations on a curved main body axis.

Created on Thu Feb  1 14:12:31 2018

@author: antariksh
"""
import numpy as np
from input_surface import axis_main
from input_surface import circle
from input_surface import length
from vector_operations import grad_Fs
from vector_operations import local_cs
from vector_operations import transformation_matrix
from testing import transform_bcs_profile
from testing import normal_vec


# Number of axial sections
Ns = 10
# Number of cross-sectional points
Nc = 100
# initialize the final point of the axis
d_X = 4
d_Y = 5
d_Z = 1

# obtain the axis in blade -coordinate system
axis = axis_main(Ns, d_X, d_Y, d_Z)
# length
ax_length = length(axis)


# get the direction of local z-axes in the body CS at each axis vertex
# gradient dF/ds
Z_bar = grad_Fs(axis)

# obtain the unit vectors of the local blade section cs at the axial vertices
# in the body cooridante system 
X_bar_unit, Y_bar_unit, Z_bar_unit = local_cs(Z_bar)

#
surface = np.zeros((Ns + 1, Nc, 3), dtype = float)

surface_test = np.zeros((Ns + 1, Nc, 3), dtype = float)

Pb = np.zeros((Nc, 3), dtype = float)

norm = np.zeros((Ns+1, 3), dtype = float)

norm_test = np.zeros((Ns+1, 3), dtype = float)
for i in range(0, Ns+1):
    # transformationmatrix
    T = transformation_matrix(X_bar_unit[i, :], Y_bar_unit[i, :],
                              Z_bar_unit[i, :], Nc)

    # construct the profile vector per section
    # obtain the profile in local blade section co-ordinate system 
    profile = circle(Nc, radius = 1, flag = 0) 
    # vector of the profile in format xi,yi,zi where i belongs to [0, Nc)
    p = np.zeros(3*Nc, dtype = float)
    # construct the profile array
    ind_p = np.arange(0, 3*Nc, step = 3, dtype = int)
    #
    p[ind_p] = profile[:, 0]
    p[ind_p + 1] = profile[:, 1]
    
    # multiply and obtain the cross-section in body cs
    P = np.dot(T.toarray(), p)
    
    # vectorially add the location
    
    Pb[:, 0] = P[ind_p] + axis[i, 0]
    Pb[:, 1] = P[ind_p + 1] + axis[i, 1]
    Pb[:, 2] = P[ind_p + 2] + axis[i, 2]
    
    # build the surface
    surface[i , :, 0] = Pb[:, 0]
    surface[i , :, 1] = Pb[:, 1]
    surface[i , :, 2] = Pb[:, 2]
    
    # build the normal vector for each span
    norm[i, :] = -normal_vec(surface[i, :, :])
    
    # test the teansformation by first getting profiles in body cs
    Pb_test = transform_bcs_profile(T, axis[i,:], Nc)
    surface_test[i , :, 0]  = Pb_test[:, 0] 
    surface_test[i , :, 1] = Pb[:, 1]
    surface_test[i , :, 2] = Pb[:, 2]
    
    # build the normal vector for each span
    norm_test[i, :] = -normal_vec(surface_test[i, :, :])
    
plot_flag = 1
if plot_flag:
    # plot each span section and save image
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
       
    fig = plt.figure('surface1')
    ax = fig.add_subplot(111, projection = '3d')
    ax.plot_surface(surface[:,:,0], surface[:,:,1], surface[:,:,2])
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.show()
    
    # figure test
    fig = plt.figure('test_surface')
    ax = fig.add_subplot(111, projection = '3d')
    ax.plot_surface(surface_test[:,:,0], surface_test[:,:,1], surface_test[:,:,2])
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.show()