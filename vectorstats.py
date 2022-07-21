# -*- coding: utf-8 -*-
"""
Created on Fri May  7 14:55:16 2021

@author: RGHS

Collection of functions for calculating circular and spherical vector
means. Created for application to palaeoflow analysis, but probably useful
for all kinds of things besides.

Sources:
    - Curray, 1956: The Analysis of Two-Dimensional Orientation Data,
        DOI: 10.1086/626329
    - Allmendinger et al., 2012: Structural geology algorithms: vectors and
        tensors
        DOI: 10.1017/CBO9780511920202
    - Fisher et al., 1987: Statistical Analysis of Spherical Data
    - Jammalamadaka and Sengupta, 2001: Topics in circular statistics
"""

#%% Import packages

import numpy as np
import pandas as pd
import warnings

# Scipy components
from scipy.special import iv
from scipy import integrate

#%% Coordinate conversions - Allmendinger et al., 2012

# Calculate pole to plane measured with dip azimuth/dip
def pole2plane(dip_az, dip):
    '''
    Calculate pole to plane measured with dip azimuth/dip. (Also works
    for converting planes to poles.)

    Parameters
    ----------
    dip_az : float
        Dip azimuth in degrees in range 0 - 360.
    dip : float
        Dip in degrees in range in range 0 - 90.

    Returns
    -------
    pol_az : float
        Trend of calculated pole.
    pol_pln : float
        Plunge of calculated pole.

    '''
    if((hasattr(dip_az, '__len__') is True) and (isinstance(dip_az, np.ndarray) is False)):
        dip_az = np.asarray(dip_az)
        dip = np.asarray(dip)
    elif(hasattr(dip_az, '__len__') is False):
         dip_az = np.array([dip_az])
         dip = np.array([dip])
    
    dip_az[dip_az == 360] = 0

    if(any(dip_az > 360) or any(dip_az < 0)):
        raise Exception('Dip azimuth must fall in range 0 <= x <= 360')
    if(any(dip > 90) or any(dip < 0)):
        raise Exception('Dip must fall in range 0 <= x <= 90')
    
    pol_az = dip_az
    add = np.where(dip_az < 180.0)
    sub = np.where(dip_az >= 180.0)
    
    pol_az[add] += 180
    pol_az[sub] -= 180

    pol_pln = 90 - dip
    
    if(len(dip_az == 1)):
        pol_az = pol_az[0]
        pol_pln = pol_pln[0]
    
    return pol_az, pol_pln

def z2p(x, dp10 = True):
    '''
    Converts circular coordinates in radians to a value between 0 and 2 pi.
    Rounds to 10 d.p. by default to reduce floating point precision errors.
    
    Parameters
    ----------
    x : Float
        Angle in radians.
    dp10 : bool
        Set whether to round to 10 d.p. (True) or to use standard python float
        precision.

    Returns
    -------
    y : Float
        Value between 0 and 2 pi corresponding to provided value.

    '''
    if(hasattr(x, '__len__') is True):
        y = np.array(x)
    else:
        y = np.array([x])
        
    y = y % (2*np.pi)
    
    if(dp10 is True):
        y[np.round(y,10) == np.round(np.pi*2,10)] = 0
        y[np.round(y,10) == 0] = 0
    
    if(len(y) == 1):
        y = y[0]
    
    return y

# Convert line with trend and plunge from spherical to cartesian coordinates
def sph2cart(ln_trend, ln_plunge, degrees = True):
    '''
    Converts spherical coordinates to cartesian coordinates. Accepts inputs
    in degrees or radians.

    Parameters
    ----------
    ln_trend : Float
        Azimuth of line. Plane datasets should be converted to poles first.
    ln_plunge : Float
        Plunge of line. Plane datasets should be converted to poles first.
    degrees : char, optional
        Specifies input as degrees (True) or radians (False). The default is True.

    Returns
    -------
    cd : Float
        Down component of cartesian coordinate.
    ce : Float
        East component of cartesian coordinate.
    cn : Float
        North component of cartesian coordinate.

    '''
    if(type(degrees) is not bool):
        raise Exception('Degrees needs boolean value.')
    
    # Convert degrees to radians
    if(degrees is True):
        ln_trend = np.radians(ln_trend)
        ln_plunge = np.radians(ln_plunge)
    
    cd = np.sin(ln_plunge)
    ce = np.cos(ln_plunge) * np.sin(ln_trend)
    cn = np.cos(ln_plunge) * np.cos(ln_trend)
    
    return cn, ce, cd

# Convert line from cartesian to spherical coordinates
def cart2sph(cn, ce, cd, degrees = True):
    '''
    Converts cartesian coordinates to spherical coordinates. Returns results
    in degrees unless otherwise specified.

    Parameters
    ----------
    cn : Float
        North component of cartesian coordinate.
    ce : Float
        East component of cartesian coordinate.
    cd : Float
        Down component of cartesian coordinate.
    angle : char, optional
        Specifies input as degrees (True) or radians (False). The default is True.

    Returns
    -------
    trend : Float
        Azimuth of line.
    plunge : Float
        Plunge of line.

    '''
    if(type(degrees) is not bool):
        raise Exception('Degrees needs boolean value.')       
    
    # Calculate plunge
    plunge = np.arcsin(cd)
    
    if((hasattr(cn, '__len__') is True) and (isinstance(cn, np.ndarray) is False)):
        cn = np.asarray(cn)
        ce = np.asarray(ce)
        cd = np.asarray(cd)
    elif(hasattr(cn, '__len__') is False):
         cn = np.asarray([cn])
         ce = np.asarray([ce])
         cd = np.asarray([cd])
    
    # Calculate trend
    trend = np.zeros(len(cn))
    trend[np.where((cn == 0.0) & (ce < 0.0))] = 3/2 * np.pi
    trend[np.where((cn == 0.0) & (ce >= 0.0))] = np.pi/2
    
    cn0 = np.where(cn != 0.0)
    cn1 = np.where(cn < 0.0)
    
    trend[cn0] = np.arctan(ce[cn0]/cn[cn0])
    trend[cn1] += np.pi
    trend[cn0] = z2p(trend[cn0])
            
    # Convert to degrees if asked
    if(degrees == True):
        trend = np.degrees(trend)
        plunge = np.degrees(plunge)
    
    if(len(trend) == 1):
        trend = trend[0]
    return trend, plunge

def circ_diff(x, y, deg_in = False, deg_out = False):
    '''
    Calculates distance of point y from point x around a circle, up to 
    a maximum of pi (or 180 degrees). Negative values indicate counter-clockwise,
    positive indicate clockwise.

    Parameters
    ----------
    x : float or array-like
        Azimuth of first point on circle.
    y : float or array-like
        Azimuth of second point on circle.
    deg_in : bool, optional
        Specifies if input is in degrees (True) or radians (False). The default is False.
    deg_out : TYPE, optional
        Specifies if output is in degrees (True) or radians (False). The default is False.

    Returns
    -------
    diff : float or np.ndarray
        Angular distance from point x to point y. Positive indicates clockwise
        distance, negative indicates counter-clockwise distance.
        
    Example
    -------
    >>> circ_diff(10, 30, deg_in = True, deg_out = True)
    20
    >>> circ_diff(10, 270, deg_in = True, deg_out = True)
    -100

    '''
    if hasattr(x, '__len__') is False:
        x = np.array([x])
        y = np.array([y])
    if isinstance(x, np.ndarray) is False:
        x = np.array(x)
        y = np.array(y)
    
    if deg_in:
        x = np.radians(x)
        y = np.radians(y)
    
    diff = (x - y) % (2*np.pi)
    invert = np.where(diff < np.pi)
    loop = np.where(diff >= np.pi)
    diff[invert] *=-1
    diff[loop] = (2*np.pi - diff[loop])
    
    if deg_out:
        diff = np.degrees(diff)
        
    if(len(diff) > 1):
        return diff
    else:
        return diff[0]

#%% Descriptive circular statistics

def circmean(theta, deg_in = False, deg_out = False, verbose = False):
    '''
    Calculates the circular mean after Jammalamadaka and Sengupta, 2001: Topics in circular statistics.
    Takes an array of floats as an input, returning a float value.
    
    Procedure as written in Jammalamadaka and Sengupta, 2001:
    if((s >= 0) and (c > 0)):
        th_bar = np.arctan(s/c)
    elif((c == 0) and (s > 0)):
        th_bar = np.pi/2
    elif(c < 0):
        th_bar = np.arctan(s/c) + np.pi
    elif((s < 0) and (c >= 0)):
        th_bar = np.arctan(s/c) + 2*np.pi
    elif((c == 0) and (s == 0)):
        th_bar = np.nan
        
    Function performs this operation using atan2 to present figure between zero and 2(pi)
    '''
    if(deg_in is True):
        theta = np.radians(theta)
    
    s = np.sum(np.sin(theta))
    c = np.sum(np.cos(theta))
    
    # print(s,c)
    
    if(abs(s) < 1e-15) and (abs(c) < 1e-15):
        th_bar = np.nan
    else:
        th_bar = np.arctan2(s,c)
        
    if(th_bar < 0):
        th_bar = th_bar + 2*np.pi
    
    if(deg_out is True):
        th_bar = np.degrees(th_bar)
    
    if verbose is True:
        return th_bar, s, c
    else:
        return th_bar

def meanveclen(theta, deg_in = False):
    '''
    Calculate R and D_v (mean vector length and circular variance)
    after Jammalamadaka and Sengupta, 2001: Topics in circular statistics.
    Takes an array of floats as an input, returning a float value.
    
    Returns
    -------
    R : float
        Magnitude of mean vector
    D_v : float
        Circular dispersion of dataset
    '''
    if(deg_in is True):
        theta = np.radians(theta)
    
    s = np.sum(np.sin(theta))
    c = np.sum(np.cos(theta))
    
    R = abs(np.sqrt(c**2 + s**2)/len(theta))
    
    D_v = 1-R
    
    return R, D_v

def circsd(R, deg_in = False, deg_out = False):
    '''
    Calculates cicular standard deviation using the mean vector length of
    a dataset. (Designed with unit vectors in mind. Don't really know if it
    makes a difference or not.)

    Parameters
    ----------
    R : float
        Length of the mean vector of the dataset.
    deg_in : bool, optional
        Set whether input is in degrees. The default is False.
    deg_out : bool, optional
        Set whether output is in degrees. The default is False.

    Returns
    -------
    sd : float
        Circular standard deviation of dataset.

    '''
    if(deg_in is True):
        R = np.radians(R)
    
    sd = np.sqrt(np.log(1/R**2))
    
    if(deg_out is True):
        sd = np.degrees(sd)
        
    return sd

def a1inv(x):
    if(hasattr(x, '__len__') is False):
        x = np.array([x])
    indices = np.arange(0, len(x))
    a1 = np.zeros(len(x))
    p1 = np.where((0 <= x) & (x < 0.53))[0]
    p2 = np.where((x >= 0.53) & (x < 0.85))[0]
    p3 = np.where(x == 1)[0]
    p123 = np.append(p1,p2)
    p123 = np.append(p123,p3)
    p4 = np.where(np.isin(indices, p123) == False)[0]
    
    a1[p1] = 2 * x[p1] + x[p1]**3 + (5 * x[p1]**5)/6
    a1[p2] = -0.4 + 1.39 * x[p2] + 0.43/(1-x[p2])
    a1[p3] = np.inf
    a1[p4] = 1/(x[p4]**3 - 4 * x[p4]**2 + 3 * x[p4])
    
    return a1

def est_kappa(theta, deg_in = False, deg_out = False):
    if(deg_in is True):
         theta = np.radians(theta)
    mu = circmean(theta)
    kappa = a1inv(np.mean(np.cos(theta - mu)))[0]
    return kappa

#%% Von Mises distribution generation
# Useful for testing methods relating to circular distributions

def vmdist(k, mu, dist_start = 0, dist_end = np.pi*2, space = int(1e3)):
    '''
    Generates a von Mises (VM) probability distribution function (PDF) and
    cumulative distribution function (CDF).

    Parameters
    ----------
    k : float
        Kappa parameter of VM distribution (concentration factor).
    mu : float
        Mu parameter of VM distribution (centrepoint).
    dist_start : float, optional
        Starting point of distribution. The default is 0.
    dist_end : float, optional
        End point of distribution. The default is 2(pi).
    space : int, optional
        Increments of distribtion. The default is 1000.

    Returns
    -------
    f : np.ndarray, dtype = float64 
        PDF of VM distribution.
    c : np.ndarray, dtype = float64 
        Numerically (trapezoidal) integrated CDF of VM distribution.

    '''
    rng = np.linspace(dist_start, dist_end, space, endpoint = True)
    f = np.exp(k * np.cos(rng - mu))/(2 * np.pi * iv(0, k))
    c = integrate.cumulative_trapezoid(f, rng)
    return f,c

def vmsmp(kappa, mu, numsample, z2p = True, rng_seed = 12345):
    '''
    Generates random samples from a Von Mises distribution using the numpy
    random generator function.

    Parameters
    ----------
    kappa : float
        Spread of the Von Mises distribution.
    mu : float
        Centrepoint of Von Mises distibution. Should fall in range -pi to pi
        or 0 to 2 pi.
    numsample : int
        Number of samples to draw from distribution.
    z2p : bool, optional
        Determines whether range of sampling is from 0 to 2 pi (True) or from
        -pi to pi (False). The default is True.
    rng_seed : int, optional
        Seed for the numpy random generator. The default is 12345.

    Returns
    -------
    vm : np array of floats
        Randomly selected samples from the specified Von Mises distribution.

    '''
    if type(numsample) is not int:
        numsample = int(numsample)
        warnings.warn('Numsample converted to int dtype - make sure its still what you wanted.')
        
    rng = np.random.default_rng(seed = rng_seed)
    if(z2p is True):
        mu -= np.pi
    
    vm = rng.vonmises(mu, kappa, numsample)
    
    if(z2p is True):
        vm += np.pi
        
    return vm


# Method from https://stackoverflow.com/questions/28839246/scipy-gaussian-kde-and-circular-data
def vmkde(data, kappa, n_bins = int(1e3), dist_min = 0):
    '''
    Produces a kernel density estimate for a given kappa of a von Mises
    distribution across a 2 Pi interval.

    Parameters
    ----------
    data : array-like
        Dataset to produce KDE from.
    kappa : float
        Kappa of dataset provided.
    n_bins : int, optional
        Number of bins to calculate KDE across. The default is 1000.
    dist_min : float, optional
        Lower bound of distribution. Common values are 0 and -pi. The default is 0.

    Returns
    -------
    bins : np.ndarray
        Bin boundaries into which data has been split.
    kde : np.ndarray
        KDE results for dataset.

    '''
    if type(n_bins) is not int:
        int(n_bins)
        warnings.warn('n_bins converted to int dtype.')
    dist_max = dist_min + 2*np.pi
    bins = np.linspace(dist_min, dist_max, n_bins)
    
    kde = np.exp(kappa * np.cos(bins[:,None] - data[None,:])).sum(1)/(2 * np.pi * iv(0, kappa))
    kde /= np.trapz(kde, x=bins)
    
    return bins, kde

#%% Curray, 1956 - circular vector mean (with modifications after Sengupta and Rao, 1966)

def currayMean(az, numbins=None, deg_in=False, deg_out=False):
    '''
    Calculates circular vector mean after Curray, 1956 (with modifications
    after Sengupta and Rao, 1966). Superseded by circmean and meanveclen.

    Parameters
    ----------
    az : array-like
        Array containing values of azimuths.
    numbins : integer
        Number of bins to separate values into. Default is None, calculating 
        dataset as ungrouped.
    deg_in : bool, optional
        Specify if input angles are in degrees (True) or radians (False). The default is False.
    deg_out : bool, optional
        Specify if output angles are in degrees (True) or radians (False). The default is False.

    Returns
    -------
    curr : pd.DataFrame
        Dataframe containing calculated properties as follows:
            - n = number of readings
            - W = E-W component of vector mean
            - V = N-S component of vector mean
            - gamma = mean vector
            - r = magnitude of resultant vector
            - L = consistency ratio (i.e., scaled magnitude)

    '''
    if deg_in is True:
        az = np.radians(az)
        
    if numbins is None:
        w_comp = np.sin(az)
        v_comp = np.cos(az)
    else:
        ni, edges = np.histogram(az, numbins, range = (0,2*np.pi))
        binwidth = edges[1] - edges[0]
        midpoints = edges[0:len(edges)-1] + binwidth/2
        
        w_comp = ni * np.sin(midpoints)
        v_comp = ni * np.cos(midpoints)
    
    n = len(az)
    W = np.sum(w_comp)
    V = np.sum(v_comp)
    gamma = np.arctan(W/V)
    
    r = np.sqrt(V**2 + W**2)
    L = r/n * 100
    
    if deg_out is True:
        gamma = np.degrees(gamma)

    curr = pd.DataFrame({'n':n,'W':W,'V':V,'gamma':gamma,'r':r,'L':L}, index=[0])
    
    return curr

#%% Fisher, 1987 - spherical vector mean (code after Allmendinger, 2012)

def fisherMean(trends, plunges, deg_in = True, deg_out = True, lower = True):
    '''
    Calculates spherical vector mean after procedure in Fisher (1987).

    Parameters
    ----------
    trends : array, float
        Array of floats containing input line trends.
    plunges : array, float
        Array of floats containing input line plunges.
    deg_in : bool, optional
        Specifies whether input is in degrees (True) or radians (False).
        The default is True.
    deg_out : bool, optional
        Specifies whether output is in degrees (True) or radians (False).
        The default is True.
    lower : bool, optional
        Specifies whether to convert resultant vector to lower hemisphere.
        The default is True.

    Returns
    -------
    trend_ave : float
        Trend of mean vector.
    plunge_ave : float
        Plunge of mean vector.
    R_ave : float
        Normalised length of resultant vector.
    conc : float
        Fisher (1987) concentration factor. Returns 0 
    d99 : float
        99% uncertainty cone for vector mean.
    d95 : float
        95% uncertainty cone for vector mean.

    '''
    if(hasattr(trends,'__len__') is False):
        trend_ave = trends
        plunge_ave = plunges
        R_ave = np.nan
        conc = np.nan
        d99 = np.nan
        d95 = np.nan
        
        warnings.warn('Only one value provided to fisherMean! Returning input as average.\
                      (Your data probably has some NaNs somewhere now.)')
        return trend_ave, plunge_ave, R_ave, conc, d99, d95
    
    if(isinstance(trends, np.ndarray) is False):
        trends = np.array(trends)
        plunges = np.array(plunges)
    
    if(deg_in is True):
        trends = np.radians(trends)
        plunges = np.radians(plunges)
        
    # Number of lines
    nlines = len(trends)
    
    # Calculate 3 direction cosines
    cn_sum = 0.0
    ce_sum = 0.0
    cd_sum = 0.0
    
    for i in range(0,nlines):
        cn, ce, cd = sph2cart(trends[i], plunges[i], degrees = False)
        cn_sum += cn
        ce_sum += ce
        cd_sum += cd
    
    # R is length of resultant vector, R_ave is R normalised to nlines
    R = np.sqrt(cn_sum**2 + ce_sum**2 + cd_sum**2)
    R_ave = R/nlines
    
    # If R_ave < 0.1, mean vector is insignificant
    if(R_ave < 0.1):
        raise Exception('Mean vector is insignificant (<0.1).')
    else:
        # Divide resultant vector by length to get average unit vector
        cn_sum = cn_sum/R
        ce_sum = ce_sum/R
        cd_sum = cd_sum/R
        # Converts mean vector to lower hemisphere
        if((cd_sum < 0.0) and (lower is True)):
            cn_sum = -cn_sum
            ce_sum = -ce_sum
            cd_sum = -cd_sum
        
        # Convert the mean vector to spherical coordinates
        trend_ave, plunge_ave = cart2sph(cn_sum, ce_sum, cd_sum, degrees = False)
        
        # Calculate Fisher statistics given sufficient lines
        if(R < nlines):
            if(nlines < 16):
                afact = 1.0 - (1.0/nlines)
                conc = (nlines/(nlines - R)) * afact**2
            else:
                conc = (nlines - 1)/(nlines - R)
        else:
            conc = np.nan
                
        if((R_ave >= 0.65) and (R_ave < 1.0)):
            afact = 100.0
            bfact = 1.0/(nlines - 1.0)
        
            # Calculate uncertainty cones
            d99 = np.arccos(1.0 - ((nlines - R)/R) * (afact**bfact - 1.0))
            afact = 20.0
            d95 = np.arccos(1.0 - ((nlines - R)/R) * (afact**bfact - 1.0))
        else:
            d99 = np.nan
            d95 = np.nan
    
    # Convert output to degrees
    if(deg_out is True):
        trend_ave = np.degrees(trend_ave)
        plunge_ave = np.degrees(plunge_ave)
        d99 = np.degrees(d99)
        d95 = np.degrees(d95)
    
    if(np.isnan(conc)):
        warnings.warn('Could not calculate vector concentration. Returning conc as NaN.')
    if(np.isnan(d99)):
        warnings.warn('Could not calculate 99% confidence cone. Returning d99 as NaN.')
    if(np.isnan(conc)):
        warnings.warn('Could not calculate 95% confidence cone. Returning d95 as NaN.')
        
    return trend_ave, plunge_ave, R_ave, conc, d99, d95

#%% Intersection of planes

def intersect(az1, d1, az2, d2, poles = False):
    '''
    Calculates the intersection of planes provided as dip azimuth/dip. Calls
    pole2plane and sph2cart.

    Parameters
    ----------
    az1 : float
        Dip azimuth of plane 1.
    d1 : float
        Dip of plane 1.
    az2 : float
        Dip azimuth of plane 2.
    d2 : float
        Dip of plane 2.
    poles : bool, optional
        Defines whether input data is planes (False) or has already been converted to poles (True). The default is False.

    Returns
    -------
    out_az : float
        Azimuth of line of intersection in spherical coordinates.
    out_plunge : float
        Plunge of plane of intersection in spherical coordinates.

    '''
    if(type(poles) is not bool):
        raise Exception('poles needs boolean value.')
    
    dip1, dip2 = d1, d2
    
    # Convert planes to poles
    if(poles is False):
        az1, d1 = pole2plane(az1, d1)
        az2, d2 = pole2plane(az2, d2)
        
    # Convert lines to cartesian coordinates
    cn1, ce1, cd1 = sph2cart(az1, d1)
    cn2, ce2, cd2 = sph2cart(az2, d2)
    
    # Calculate theta
    theta = np.arccos(ce1*ce2 + cn1*cn2 + cd1*cd2)
    
    if(theta != 0):
        i_ce = (cn1*cd2 - cd1*cn2)/np.sin(theta)
        i_cn = -(ce1*cd2 - cd1*ce2)/np.sin(theta)
        i_cd = (ce1*cn2 - cn1*ce2)/np.sin(theta)
    else:
        i_ce, i_cn, i_cd = 0,0,0
        
    # Convert to lower hemisphere
    if(i_cd < 0):
        i_ce *= -1
        i_cn *= -1
        i_cd *= -1
    
    # Calculate intersection azimuth
    if(dip1 == 90 and dip2 == 90):
        out_az = 0
    elif((i_ce == 0) and (i_cn == 0)):
        out_az = 0
    elif((i_ce < 0) and (i_cn >= 0)):
        out_az = 450 - np.degrees(np.arctan2(i_cn, i_ce))
    else:
        out_az = 90 - np.degrees(np.arctan2(i_cn, i_ce))
    
    # Calculate intersection plunge
    if(dip1 == 90 and dip2 == 90):
        out_plunge = 90
    else:
        out_plunge = 90 - np.degrees(np.arccos(i_cd))
    
    return out_az, out_plunge

#%% Rotate data on arbitrary axis

def rotate(dip_az, dip, rot_ang, ax_az, ax_pln = 0, bed_az_2_ax = True, lower = True):
    '''
    Rotates dip azimuth and dip data about an arbitrary axis. Accepts inputs
    in degrees. Designed with restoring features on inclined bedding to
    horizontal as primary purpose, but can be used for any rotations around a
    single axis. Calls sph2cart and cart2sph.
    
    Parameters
    ----------
    dip_az : float or array like
        Azimuth of feature(s) to be rotated, measured in degrees.
    dip : float or array like
        Dip (or plunge) of feature(s) to be rotated, measured in degrees.
    rot_ang : float
        Angle through which features are to be rotated, measured in degrees.
        For example, if restoring features on folded bedding to horizontal,
        this value should be set to the dip of the bedding.
    ax_az : float
        Azimuth of rotation axis by right-hand-rule. If unfolding bedding,
        entering the value of the bedding's dip azimuth and setting
        bed_az_2_ax to true will automatically convert the dip azimuth to the
        rotation axis azimuth.
    ax_pln : float, optional
        Plunge of rotation axis. The default is 0, the value for restoring
        features on inclined bedding.
    bed_az_2_ax : bool, optional
        Converts dip azimuth of bedding to rotation axis azimuth if set to
        True. The default is True.
    lower : bool, optional
        Converts all results to lower hemisphere measurements if set to True.
        The default is True.

    Returns
    -------
    trend : float or ndarray
        Azimuth of rotated feature(s), measured in degrees.
    plunge : float or ndarray
        Plunge of rotated feature(s), measured in degrees.

    '''
    if((hasattr(dip_az, '__len__') == True) and (isinstance(dip_az, np.ndarray) is False)):
        dip_az = np.asarray(dip_az)
        dip = np.asarray(dip)
        
    # Bed azimuth to axis
    if(bed_az_2_ax is True):
        if((ax_az - 90) > 0):
            ax_az = ax_az - 90
        else:
            ax_az = 360 + (ax_az - 90)
    
    # Calculate data vector in cartesian space
    d_cn, d_ce, d_cd = sph2cart(dip_az, dip)
    
    # Calculate axis vector in cartesian space
    a_cn, a_ce, a_cd = sph2cart(ax_az, ax_pln)
    
    # Calculate dot product
    dot = (d_ce*a_ce + d_cn*a_cn + d_cd*a_cd) * (1 - np.cos(np.radians(rot_ang)))
    
    # Calculate cartesian coordinates of rotated point
    r_ce = np.cos(np.radians(rot_ang))*d_ce + dot*a_ce + np.sin(np.radians(rot_ang))*(a_cn*d_cd - a_cd*d_cn)
    r_cn = np.cos(np.radians(rot_ang))*d_cn + dot*a_cn - np.sin(np.radians(rot_ang))*(a_ce*d_cd - a_cd*d_ce)
    r_cd = np.cos(np.radians(rot_ang))*d_cd + dot*a_cd + np.sin(np.radians(rot_ang))*(a_ce*d_cn - a_cn*d_ce)
    
    # Convert point to lower hemisphere
    if(lower is True):
        if(isinstance(r_cd, np.ndarray) is True):
            convert = np.array((r_cd >= 0), dtype = int)
            convert[convert == 0] = -1
            
            r_ce = r_ce * convert
            r_cn = r_cn * convert
            r_cd = r_cd * convert
        elif(r_cd < 0):
            r_ce *= -1
            r_cn *= -1
            r_cd *= -1
    
    # Convert rotated coordinates to spherical space
    trend, plunge = cart2sph(r_cn, r_ce, r_cd)
    
    return trend, plunge

#%% Geographic midpoint

def geomid(lat,lon):
    '''
    Calculates an approximate midpoint of a set of coordinates in decimal 
    latitude and longitude. Does not account for ellipsoid deviations.

    Parameters
    ----------
    lat : array-like
        List of latitudes.
    lon : array-like
        List of longitudes.

    Returns
    -------
    avlat : float
        Midpoint latitude of the provided coordinates.
    avlong : float
        Midpoint longitude of the provided coordinates.

    '''
    if((hasattr(lat, '__len__') is False) and (hasattr(lat, '__lon__') is False)):
        warnings.warn('Need more than one point to compute average! Returning input lat and long.')
        return lat, lon
    
    if(len(lat) is not len(lon)):
        raise Exception('lat and lon must be of equal length.')
    
    rlat = np.radians(lat)
    rlon = np.radians(lon)
    
    x = np.mean(np.cos(rlat) * np.cos(rlon))
    y = np.mean(np.cos(rlat) * np.sin(rlon))
    z = np.mean(np.sin(rlat))
    
    if((abs(x) < 10**-9) and (abs(y) < 10**-9) and (abs(z) < 10**-9)):
        warnings.warn('Midpoint is centre of the Earth!')
        avlat = np.nan
        avlon = np.nan
    else:
        avlon = np.arctan2(y,x)
        avhyp = np.sqrt(x**2 + y**2)
        avlat = np.arctan2(z, avhyp)
        
        avlon = np.degrees(avlon)
        avlat = np.degrees(avlat)
    
    return avlat, avlon