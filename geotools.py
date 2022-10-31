# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:12:11 2022

@author: RGHS
"""

import numpy as np
import pandas as pd
import warnings

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

#%% Find distance between points

def findDist(point_lat, point_lon, comp_lat, comp_lon, deg=True):
    '''
    Uses the Haversine Equation to give an approximate distance between a single
    point and an arbitrary list of points. Does not account for ellipsoid deviation.

    Parameters
    ----------
    point_lat : float
        Latitude of starting point (in decimal degrees).
    point_lon : float
        Longitude of starting point (in decimal degrees).
    comp_lat : float or array of floats
        Latitude of ending point (or array thereof) (in decimal degrees).
    comp_lon : float or array of floats
        Longitude of ending point (or array thereof) (in decimal degrees).
    deg : bool, optional
        Specifies whether inputted latitude and longitude are in degrees.
        The default is True.

    Returns
    -------
    distance : pandas Series
        Distance in kilometres between provided points.

    '''
    if deg is True:
        point_lat = np.radians(point_lat)
        point_lon = np.radians(point_lon)
        comp_lat = np.radians(comp_lat)
        comp_lon = np.radians(comp_lon)
        
    dlon = comp_lon - point_lon
    dlat = comp_lat - point_lat
    a = np.sin(dlat/2)**2 + np.cos(point_lat) * np.cos(comp_lat) * np.sin(dlon/2)**2
    distance = 2 * np.arcsin(np.sqrt(a)) * 6371
    
    if not isinstance(distance, pd.Series):
        distance = pd.Series(distance)
    
    return distance


#%% Find closest and furthest points

def closest(point_lat, point_lon, comp_lat, comp_lon, deg=True):
    '''
    Calls findDist to find the closest point in a list of points to the starting
    point given. See documentation for findDist for info on inputs.
    
    Returns multiple results if several points are exactly equidistant.
    '''
    dist = findDist(point_lat, point_lon, comp_lat, comp_lon, deg=True)
    idx = dist[dist == np.min(dist)].index
    
    close_lat = comp_lat[idx].reset_index(drop=True)
    close_lon = comp_lon[idx].reset_index(drop=True)
    
    return close_lat, close_lon

def furthest(point_lat, point_lon, comp_lat, comp_lon, deg=True):
    '''
    Calls findDist to find the furthest point in a list of points to the starting
    point given. See documentation for findDist for info on inputs.
    
    Returns multiple results if several points are exactly equidistant.
    '''
    dist = findDist(point_lat, point_lon, comp_lat, comp_lon, deg=True)
    idx = dist[dist == np.max(dist)].index
    
    far_lat = comp_lat[idx].reset_index(drop=True)
    far_lon = comp_lon[idx].reset_index(drop=True)
    
    return far_lat, far_lon