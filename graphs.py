# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 14:41:15 2021

@author: RGHS

Generate graphs for various procedures and find intersection of lines.
"""

import numpy as np

#%% Line intersect
# Find the intersection of two straight lines by specifying two points on each
# Method from Norbu Tsering on Stack Overflow
# https://stackoverflow.com/questions/3252194/numpy-and-line-intersections/42727584#42727584

def get_intersect(a1, a2, b1, b2):
    """ 
    Returns the point of intersection of the lines passing through a2,a1 and b2,b1.
    a1: [x, y] a point on the first line
    a2: [x, y] another point on the first line
    b1: [x, y] a point on the second line
    b2: [x, y] another point on the second line
    """
    s = np.vstack([a1,a2,b1,b2])        # s for stacked
    h = np.hstack((s, np.ones((4, 1)))) # h for homogeneous
    l1 = np.cross(h[0], h[1])           # get first line
    l2 = np.cross(h[2], h[3])           # get second line
    x, y, z = np.cross(l1, l2)          # point of intersection
    if z == 0:                          # lines are parallel
        return (float('inf'), float('inf'))
    
    int_x = x/z
    int_y = y/z
    return int_x, int_y

def lerouxEAD():
    EAD_x = np.arange(0,101,1)
    EAD_y = np.zeros(101)
    
    eq6 = np.where(EAD_x <= 85)
    eq7 = np.where((EAD_x > 85) & (EAD_x < 98))
    eq8 = np.where(EAD_x >= 98)
    
    EAD_y[eq6] = 5e-8*EAD_x[eq6]**4 - 9e-5*EAD_x[eq6]**3 + 1.07e-2*EAD_x[eq6]**2 - 1.7402*EAD_x[eq6] + 179.95
    EAD_y[eq7] = -2.7e-3*(100 - EAD_x[eq7])**4 + 9.53e-2*(100 - EAD_x[eq7])**3 - 1.23*(100 - EAD_x[eq7])**2 + 9.3963*(100 - EAD_x[eq7]) + 5.9
    EAD_y[eq8] = 1.2183*(100 - EAD_x[eq8])**3 - 5.8723*(100 - EAD_x[eq8])**2 + 14.654*(100 - EAD_x[eq8]) + 5
    
    return EAD_x, EAD_y

def cmp(author = 1):
    x_base = np.array([1,2,4,6,8,10,20,40,60,80,100,200,400,600,800,1000])
    if(author == 1): # Van Rijn
        cmp_x = x_base[1:16]
        cmp_y = np.array([0.2, 0.1, 0.08, 0.077, 0.08, 0.1, 0.15, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2])
    elif(author == 2):  # Shields
        cmp_x = x_base
        cmp_y = np.array([0.2, 0.1, 0.06, 0.04, 0.035, 0.03, 0.03, 0.04, 0.045, 0.05, 0.055, 0.055, 0.055, 0.055, 0.055, 0.055])
    elif(author == 3):  # Engelund
        cmp_x = x_base[5:16]
        cmp_y = np.array([0.035, 0.06, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08])
    elif(author == 4):  # Bagnold
        cmp_x = x_base[2:16]
        cmp_y = np.array([0.15, 0.3, 0.4, 0.6, 1, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2])
    
    return cmp_x, cmp_y