# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 13:58:29 2022

@author: rghs

Collection of small functions applying simple flow dynamics equations. Made for
my convenience only, very little in the way of innovation here.
"""

import numpy as np

def antidune(lenslength):
    g = 9.81
    
    wavelength = np.array([lenslength]*3)/np.array([0.4,0.45,0.5])
    u = np.sqrt((g*wavelength)/(2*np.pi))
    hm = wavelength/(2*np.pi)
    
    return u, hm
    