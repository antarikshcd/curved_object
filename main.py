#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build a curved object out of input profile data defined in local co-ordinate
system at locations on a curved main body axis.

Created on Thu Feb  1 14:12:31 2018

@author: antariksh
"""
import numpy as np
from surface_2d import axis_main
from surface_2d import circle
from surface_2d import length
# Number of axial sections
Ns = 10
# Number of cross-sectional sections
Nc = 100
# initialize the final point of the axis
d_X = 1
d_Y = 1
d_Z = 1
# obtain the axis in blade -coordinate system
axis = axis_main(Ns, d_X, d_Y, d_Z)
# obtain the profile in local blade section co-ordinate system 
profile = circle(Nc, radius = 1, flag = 0) 
# length
ax_length = length(axis)


# get the direction of local z-axes in the body CS at each axis vertex
# gradient dF/ds

