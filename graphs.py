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

def get_y_val(x_series,y_series,x_unknown):
    '''
    Uses get_intersect to get the Y value of a desired x point.
    '''
    x_series = np.asarray(x_series)
    y_series = np.asarray(y_series)
    
    lower = np.max(np.where(x_series <= x_unknown))
    upper = np.max(np.where(x_series >= x_unknown))
    
    a1 = [x_unknown, np.min(y_series)]
    a2 = [x_unknown, np.max(y_series)]
    b1 = [x_series[lower], y_series[lower]]
    b2 = [x_series[upper], y_series[upper]]
    
    _, intersect_y = get_intersect(a1, a2, b1, b2)
    
    return intersect_y

#%% Useful graphs in sedimentology
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

def carling():
    '''
    Generate the dune stability field from Carling, 1999. Returns x and y
    for upper and lower boundaries of field

    Returns
    -------
    lower_x : list
    lower_y : list
    upper_x : list
    upper_y : list
    
    '''
    upper_x = [0.00001151,0.00001243,0.00001362,0.00001545,0.00001803,0.00002074,0.00002403,0.00002689,0.00003051,0.00003461,0.00003955,0.0000455,0.00005236,0.00006152,0.00007229,0.00009241,0.00010339,0.00012064,0.00014077,0.00017132,0.00020274,0.00025198,0.00030453,0.00038927,0.00046064,0.0005413,0.00066808,0.00079058,0.0009225,0.00110707,0.0012213,0.00139543,0.00165128,0.00184742,0.00198166,0.00218615,0.00236152,0.00249783,0.00271719,0.00289427,0.00299758,0.00310458,0.00342494,0.00391324,0.00447115,0.00514456,0.00617385,0.00715371,0.00805975,0.00889143,0.00987799,0.0112074,0.01316963,0.01442705,0.01648391,0.01844181,0.02107105,0.02373974,0.02829029,0.03142927,0.03565912,0.03989458]
    upper_y = [0.66967476,0.70855272,0.75499477,0.82168679,0.90697585,0.98015454,1.04439871,1.10503127,1.16918385,1.22836524,1.29054225,1.36546467,1.42449736,1.4756362,1.51786596,1.59469683,1.60598565,1.64033386,1.64033386,1.64033386,1.61735439,1.57235668,1.51786596,1.42449736,1.34633584,1.28147075,1.17746047,1.08189217,0.98709303,0.86327868,0.79321086,0.71356855,0.59820143,0.53060023,0.48070431,0.43550044,0.39734038,0.35997578,0.32612483,0.29754864,0.27924551,0.26767328,0.26579175,0.26767328,0.26767328,0.26767328,0.26956813,0.27533355,0.27924551,0.28321305,0.2892703,0.29965498,0.31705145,0.32612483,0.33545787,0.33545787,0.34022408,0.34022408,0.33783257,0.33783257,0.34263251,0.34263251]
    
    lower_x = [0.00001167,0.00001333,0.00001492,0.00001692,0.0000192,0.00002103,0.00002472,0.00002804,0.00003094,0.0000361,0.00004096,0.00004582,0.00005163,0.00005695,0.00006599,0.00007593,0.00008615,0.00009438,0.00010708,0.00011896,0.00013497,0.00014995,0.00017374,0.00020132,0.00022841,0.00025198,0.00027411,0.0003024,0.0003336,0.00036036,0.00040317,0.00043857,0.00047709,0.0005413,0.00059298,0.00067278,0.00075799,0.00088448,0.00100351,0.00117921,0.00131927,0.0014251,0.00159437,0.00192683,0.00229618,0.00266061,0.00308288,0.0035472,0.00459838,0.00555724,0.00643924,0.00783674,0.00876756,0.00974038,0.01112906,0.01262684,0.01536723,0.0174354,0.02048804,0.02493453,0.02909533,0.03232363,0.03641749,0.04219739,0.04787644]
    lower_y = [0.11644168,0.11644168,0.11726597,0.11726597,0.11726597,0.11644168,0.11481045,0.11400342,0.11161621,0.11161621,0.109279,0.1077481,0.10623866,0.10401405,0.10255691,0.09900275,0.09692966,0.09357052,0.08969286,0.08537155,0.08183366,0.07844238,0.07361715,0.06812087,0.06438311,0.05999798,0.05630732,0.05173715,0.04654248,0.04461371,0.04216578,0.03985216,0.0387434,0.03793213,0.03713784,0.03687679,0.03687679,0.03687679,0.03687679,0.03687679,0.03713784,0.03687679,0.03661757,0.03687679,0.03713784,0.03740073,0.03793213,0.03820065,0.04041839,0.04306761,0.04524759,0.04959311,0.05284369,0.05591153,0.05999798,0.06348116,0.06908874,0.07361715,0.07899768,0.08658451,0.09032779,0.09357052,0.09624832,0.10040939,0.10183602]
    return lower_x, lower_y, upper_x, upper_y

#%% Stock graph plotting

def plot_leroux(silent = True):
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
    
    if silent is True:
        plt.close()
    
    return fig

#%% Custom graph plotting

def windrose(x, b_int = 10, silent = True):
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
    
    if silent is True:
        plt.close()
    
    return fig