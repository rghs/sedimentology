# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 14:41:15 2021

@author: RGHS

Generate graphs for various procedures and find intersection of lines.
"""

import numpy as np
import matplotlib.pyplot as plt
import warnings

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
    '''
    Generate axes for plotting the Effective Angular Deviation graph from
    Le Roux (1994).

    Returns
    -------
    EAD_x : np.ndarray
        Array of 0:100 for plotting x-axis of graph.
    EAD_y : np.ndarray
        Array of corresponding values of effective angular deviation for plotting
        y-axis.

    '''
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
    '''
    Generates axes of initiation of motion graphs based on van Rijn (1984).

    Parameters
    ----------
    author : int, optional
        Specifies which curve to generate. Options are as follows:
            - 1: Van Rijn, 1984
            - 2: Shields
            - 3: Engelund
            - 4: Bagnold
        The default is 1.

    Returns
    -------
    cmp_x : np.ndarray
        Points on x-axis of selected graph.
    cmp_y : np.ndarray
        Points on y-axis of selected graph.

    '''
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

#%% Stock graph plotting

def plot_leroux():
    '''
    Plots the Le Roux (1994) Effective Angular Deviation graph.

    Returns
    -------
    fig : matplotlib figure
        Plot of EAD vs consistency ratio.

    '''
    fig = plt.figure(figsize = (10,10))
    x,y = lerouxEAD()
    plt.plot(x,y)
    plt.xlim(0,100)
    plt.xticks(np.arange(0,120,20))
    plt.xlabel('Consistency ratio')
    plt.ylim(0, 180)
    plt.yticks(np.arange(0,210,30))
    plt.ylabel('Effective angular deviation / degrees')
    
    return fig

#%% Custom graph plotting

def windrose(x, b_int = 10):
    '''
    Generates windroses of directional (most commonly paleoflow) data.

    Parameters
    ----------
    x : array-like
        List of directional measurements.
    b_int : float, optional
        Interval by which measurements are binned. The default is 10.
        (Note: numbers greater than 10 produce rather ugly graphs)

    Returns
    -------
    fig : matplotlib figure
        Windrose of plotted figure.

    '''
    if(b_int > 10):
        warnings.warn('Bin sizes greater than 10 produce strange looking bars. You have been warned.')
    bins = np.arange(0,360+b_int,b_int)
    binned = np.histogram(x, bins)[0]
    
    xmax = 2 * np.pi
    x_coords = np.linspace(0, xmax, len(binned), endpoint=False)
    width = xmax/len(bins)
    
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(111, polar = True)
    ax.bar(x_coords, binned, width = width, align = 'edge')
    ax.set_theta_direction(-1)
    ax.set_theta_zero_location("N")
    
    return fig