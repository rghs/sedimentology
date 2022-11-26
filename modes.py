# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 15:50:54 2022

@author: rghs

Small functions (probably massively over-engineered) to calculate modes
of distributions since there aren't really any good ways in base Python
"""

import numpy as np
import scipy.stats as ss
import scipy.optimize as so

#%% Mode calculation methods based on finding peak of a probability density
# function, adapted from here: https://stats.stackexchange.com/questions/19952/computing-the-mode-of-data-sampled-from-a-continuous-distribution

def modePdf(data, bins = 1000):
    '''
    Calculates mode by fitting continuous data to a PDF and finding the peak.
    Accepts an array like as input and returns a float value of mode.
    Also accepts bins as a parameter to determine how many points the
    PDF should be composed of.
    '''
    dist = ss.gaussian_kde(data)
    x = np.linspace(np.min(data), np.max(data), bins)
    y = dist.pdf(x)
    i = np.argmax(y)
    mode = x[i]
    
    return mode

def modeOpt(data, mode = 'min'):
    dist = ss.gaussian_kde(data)
    def objective(x):
        return 1/dist.pdf(x)[0]
    
    bounds = [(np.min(data),np.max(data))]
    
    if mode == 'min':
        solution = so.minimize(objective, [1], bounds = bounds)
    elif mode == 'shgo':
        solution = so.shgo(objective, bounds = bounds, n=100*len(data))
    else:
        raise ValueError('Mode accepts either "min" or "shgo"!')
    
    return solution.x[0]