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
    EAD_y = np.array([180,177.833,176.3,174.768,173.235,171.703,170.17,
                        168.638,167.105,165.573,164.04,162.508,160.975,159.443,
                        157.972,156.504,155.037,153.569,152.102,150.634,149.167,
                        147.699,146.232,144.764,143.297,141.829,140.362,138.894,
                        137.427,135.959,134.492,133.024,131.557,130.102,128.895,
                        127.688,126.481,125.274,124.067,122.861,121.654,120.447,
                        119.24,118.033,116.826,115.62,114.413,113.206,111.999,
                        110.792,109.585,108.378,107.172,105.965,104.757,103.26,
                        101.763,100.266,98.769,97.272,95.775,94.278,92.781,91.284,
                        89.788,88.291,86.794,85.297,83.8,82.303,80.806,79.303,
                        77.646,75.989,74.333,72.676,71.019,69.362,67.706,66.049,
                        64.392,62.735,61.079,59.422,57.383,55.331,53.279,51.227,
                        49.174,47.122,45.07,43.018,40.966,38.052,35.133,32.214,
                        29.102,23.713,18.324,9.554,0])
    
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