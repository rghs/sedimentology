# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 13:21:23 2021

@author: RGHS

Functions for calculating paleosinuosity.

Sources:
    - Ghosh, 2000: Estimation of channel sinuosity from paleocurrent data: a
        method using fractal geometry
"""

#%% Import packages

import numpy as np
import pandas as pd
import graphs

#%% Ghosh, 2000

def ghosh(L):
    '''
    Calculates sinuosity after Ghosh (2000).

    Parameters
    ----------
    L : float or array_like
        Consistency ratio of vector mean magnitude, calculated after Curray (1956).

    Returns
    -------
    sinuosity : pandas DataFrame
        Returns calculated sinuosities in data frame with column headings "S",
        "S_max" and "S_abstract".

    '''
    S = np.exp(2.49 - 0.0475 * L + 0.000234 * np.square(L))
    S_max = np.exp(3.68 - 0.0684 * L + 0.00032 * np.square(L))
    S_abstract = np.exp(4.2416 - 0.0801 * L + 0.0004 * np.square(L))
    
    sinuosity = pd.DataFrame({'S': S,
                         'S_max': S_max,
                         'S_abstract': S_abstract})
    
    return sinuosity

#%% Le Roux 1992, 1994

def leroux(L):
    '''
    Calculates sinuosity after Le Roux (1992, 1994).

    Parameters
    ----------
    L : float or array_like
        Consistency ratio of vector mean magnitude, calculated after Curray (1956).

    Returns
    -------
    sinuosity : pandas DataFrame
        Returns calculated sinuosities in data frame with column headings "Phi"
        and "S".

    '''
    eadx,eady = graphs.lerouxEAD()
    
    if(hasattr(L, '__len__') is False):
        L = np.array([L])
    phi = np.zeros(len(L))
    
    for i in range(0,len(L)):
        if(L[i] == 100):
            phi[i] = 0
        elif(L[i] == 0):
            phi[i] = 180
        else:
            a1 = (L[i] , 0)
            a2 = (L[i] , 180)
            b1 = (eadx[np.where(eadx <= L[i])[0][-1]], eady[np.where(eadx <= L[i])[0][-1]])
            b2 = (eadx[np.where(eadx > L[i])[0][0]], eady[np.where(eadx > L[i])[0][0]])
            phi[i] = 2 * (graphs.get_intersect(a1, a2, b1, b2)[1])
    
    r_phi = np.radians(phi)
    S = (r_phi/2)/(np.sin(r_phi/2))
    
    sinuosity = pd.DataFrame({'Phi':phi,
                              'S':S})
    
    return sinuosity