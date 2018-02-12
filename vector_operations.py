#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 09:56:20 2018

@author: antariksh
"""
import numpy as np

def transformation_matrix(X_bar_unit, Y_bar_unit, Z_bar_unit, Nc):
    """
    Constructs the transformation matrix for every spanwise section to 
    tranform the profile from body co-ordinate to blade co-ordinate system.
    
    Args:
        Nc (int): number of cross-section points
        X_bar_unit (float): Numpy array of size 3, of the X-axis directional 
                            unit vector in body co-ordinate system
        Y_bar_unit (float): Numpy array of size 3, of the Y-axis directional 
                            unit vector in body co-ordinate system
        Z_bar_unit (float): Numpy array of size 3, of the Z-axis directional 
                            unit vector in body co-ordinate system
                            
    Returns:
        T (float): scipy sparse matrix of size (3*Nc X 3*Nc) representing the
                    transformation matrix for the spanwise section
    
    """
    # import scipy sparse matrix
    from scipy.sparse import coo_matrix
    
    # construct the transformation matrix
    data = np.zeros(9*Nc, dtype = float)
    row = np.zeros(9*Nc, dtype = int)
    col = np.zeros(9*Nc, dtype = int)

    #
    ind = np.arange(0, 9*Nc, step = 9, dtype = int)
    # assign the data, row and col 
    row_ind = np.arange(0, 3*Nc, step = 3, dtype = int)
    col_ind = np.arange(0, 3*Nc, step = 3, dtype = int)

    
    # assign the data
    data[ind] = X_bar_unit[0]
    data[ind + 1] = Y_bar_unit[0]
    data[ind + 2] = Z_bar_unit[0]
    #
    data[ind + 3] = X_bar_unit[1]
    data[ind + 4] = Y_bar_unit[1]
    data[ind + 5] = Z_bar_unit[1]
    #
    data[ind + 6] = X_bar_unit[2]
    data[ind + 7] = Y_bar_unit[2]
    data[ind + 8] = Z_bar_unit[2]

    # assign corresponding row
    row[ind] = row_ind 
    row[ind + 1] = row_ind
    row[ind + 2] = row_ind
    #
    row[ind + 3] = row_ind + 1
    row[ind + 4] = row_ind + 1
    row[ind + 5] = row_ind + 1
    #
    row[ind + 6] = row_ind + 2
    row[ind + 7] = row_ind + 2
    row[ind + 8] = row_ind + 2

    # assign corresponding column
    col[ind] = col_ind
    col[ind + 1] = col_ind + 1
    col[ind + 2] = col_ind + 2
    #
    col[ind + 3] = col_ind
    col[ind + 4] = col_ind + 1
    col[ind + 5] = col_ind + 2
    #
    col[ind + 6] = col_ind
    col[ind + 7] = col_ind + 1
    col[ind + 8] = col_ind + 2

    # build the sparse transfromattion matrix
    T = coo_matrix((data,(row, col)), shape= (3*Nc, 3*Nc))

    return T

def local_cs(Z_bar):
    """
    Local co-ordinate systems at axial (spanwise) locations in the body co-
    ordinate system.
    
    Args:
        Z_bar (float): The direction vector representing the local Z-axes at
                       axial locations, tangent to the axis.
    
    Returns:
        X_bar_unit (float): Numpy array of size 3, of the X-axis directional 
                            unit vector in body co-ordinate system
        Y_bar_unit (float): Numpy array of size 3, of the Y-axis directional 
                            unit vector in body co-ordinate system
        Z_bar_unit (float): Numpy array of size 3, of the Z-axis directional 
                            unit vector in body co-ordinate system                   
    """
    
    # Number of sections
    Ns = Z_bar.shape[0]
    
    # Z_bar unit vector
    Z_bar_unit = np.zeros((Ns, 3), dtype = float)
    # X-component
    Z_bar_unit[:, 0]  = np.divide(Z_bar[:, 0], np.linalg.norm(Z_bar, axis = 1)) 
    # Y-component
    Z_bar_unit[:, 1]  = np.divide(Z_bar[:, 1], np.linalg.norm(Z_bar, axis = 1))
    # Z-component
    Z_bar_unit[:, 2]  = np.divide(Z_bar[:, 2], np.linalg.norm(Z_bar, axis = 1))

    # define the Rotation axis unit vector in the body cooridinate system
    # Rotation axis in body system = +ve Y-axis at the root of object
    R_bar_unit = np.zeros((Ns, 3), dtype = float)
    R_bar_unit[: , 1] = 1

    # get the X-axis at the vertices in blade coordinate system
    # X_bar = R_bar_unit (cross) Z_bar_unit
    X_bar = np.cross(R_bar_unit, Z_bar_unit, axisa = 1, axisb = 1)
    # get the unit vector of X_bar
    X_bar_unit = np.zeros((Ns, 3), dtype = float)
    # X-component
    X_bar_unit[:, 0] = np.divide(X_bar[:, 0], np.linalg.norm(X_bar, axis = 1))
    # Y-component
    X_bar_unit[:, 1] = np.divide(X_bar[:, 1], np.linalg.norm(X_bar, axis = 1))
    # Z-component
    X_bar_unit[:, 2] = np.divide(X_bar[:, 2], np.linalg.norm(X_bar, axis = 1))

    # get the Y-axis at the vertices in blade coordiante system
    # Z_bar_unit (cross) X_bar_unit
    Y_bar_unit  = np.cross(Z_bar_unit, X_bar_unit, axisa = 1, axisb = 1)
    
    return X_bar_unit, Y_bar_unit, Z_bar_unit

def grad_Fs(axis):
    """
    Calculates the gradient 
        
    """
    t1 = np.gradient(axis[:,:])[0]
    dFds = np.array([t1[i, :] / np.linalg.norm(t1[i, :]) for i in range(t1.shape[0])])
    
    return dFds