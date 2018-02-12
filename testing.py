#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 13:41:23 2018

@author: antariksh
"""

import numpy as np

def transform_bcs_profile(T, axis, Nc):
    """ Translates the profile to body cs and then transforms it for rotation.
    """
    from input_surface import circle
    profile = circle(Nc, radius = 1, flag = 0) 
    
    Pb_new = np.zeros((Nc, 3), dtype = float)
    ind_p = np.arange(0, 3*Nc, step = 3, dtype = int)
    p_new = np.zeros(3*Nc, dtype = float)
    
    p_new[ind_p] = profile[:, 0] + axis[0]
    p_new[ind_p + 1] = profile[:, 1] + axis[1]
    p_new[ind_p + 2] = axis[2]
    
    P_new = np.dot(T.toarray(), p_new)
    
    Pb_new[:, 0] = P_new[ind_p]
    Pb_new[:, 1] = P_new[ind_p + 1]
    Pb_new[:, 2] = P_new[ind_p + 2]
    
    return Pb_new

def normal_vec(surface):
    """
    Obtain the normal vector of the surface.
    
    Args:
        surface (float): Cross-sectional surface definition of order N x3
        
    Return:
        normal(float) : normal unit vector of shape (1x3)
    """
    # number of points on the surface
    N = surface.shape[0] 
    # take the first and last points on the surface and find the furthest points
    # first point on the lower surface corres. to t=0
    P1 = surface[0, :]
    #last point on the surface corresponding to t = N-1
    Pn = surface[N-1, :]
    
    # find the point whose dif. in distance from P1 and Pn is min
    d1 = np.sqrt(np.power(P1[0] - surface[:, 0], 2) + 
                 np.power(P1[1] - surface[:, 1], 2) + 
                 np.power(P1[2] - surface[:, 2], 2))
    
    dn = np.sqrt(np.power(Pn[0] - surface[:, 0], 2) + 
                 np.power(Pn[1] - surface[:, 1], 2) + 
                 np.power(Pn[2] - surface[:, 2], 2))
    
    ind_min = np.argmin(abs(d1 - dn))
    # the third point on the surface
    Pm = surface[ind_min, :]
    
    # define the vectors
    V1m = Pm - P1 # vector from first point on lower surface to mid point
    Vnm = Pm - Pn # vector from last point on upper surface to mid point
    
    # get the normal
    norm_vec = np.cross(V1m, Vnm)
    norm = norm_vec/np.linalg.norm(norm_vec)
    
    return norm