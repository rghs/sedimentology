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
    - Ghosh, 2000: Estimation of channel sinuosity from paleocurrent data: a
        method using fractal geometry
    - Allmendinger et al., 2012: Structural geology algorithms: vectors and
        tensors
"""

#%% Import packages
import numpy as np
import pandas as pd
import math

#%% Curray, 1956 - circular vector mean
### Also includes calculation of sinuosity by relation found in Ghosh, 2000

def currayMean(az, numbins, localities = [], sinuosity = False):
    '''
    Calculates circular vector mean after Curray, 1956.
    Optionally calculates paleosinuosity after Ghosh, 2000

    Parameters
    ----------
    az : float array_like
        Array containing values of azimuths.
    numbins : integer
        Number of bins to separate values into.
    localities : chr array_like, optional
        Array containing values for site measurements were collected.
        Must be of identical length to az.
    sinuosity : bool, optional
        Calculate paleosinuosity after Ghosh, 2000. The default is False.

    Returns
    -------
    results : pandas dataframe
        Statistics calculated during calculation of vector mean +/- sinuosity.

    '''
    az = np.array(az)       ## Convert data to numpy array
    az[az == 360] = 0       ## Convert 360 to 0
    
    # Raise exceptions for bad data
    if any(y > 360 for y in az) or any(y < 0 for y in az):
        raise Exception('Array contains values outside range 0 <= x <= 360.')
    if(type(sinuosity) != bool):
        raise Exception('Argument "sinuosity" must be provided with bool value.')
        
    # Build bins
    bins = np.linspace(0, 360, numbins, endpoint = False)
    bins = np.append(bins, 360)
    
    vec = pd.DataFrame({'bin':bins[0:len(bins)-1],
                        'ni':np.zeros(len(bins)-1),
                        'w_comp':np.zeros(len(bins)-1),
                        'v_comp':np.zeros(len(bins)-1)})
    
    # Caclulate mean on all data or by locality
    if(len(localities) > 0):
        # Check if localities list is the right length
        if(len(az) != len(localities)):
            raise Exception('Localities list must be of same length as azimuth list.')
        
        # Build dataframe for results
        df = pd.DataFrame({'az':az,
                           'locs':localities})
        arealist = localities.unique()
        
        results = pd.DataFrame({'locality':np.zeros(len(arealist)),
                                'n':np.zeros(len(arealist)),
                                'W':np.zeros(len(arealist)),
                                'V':np.zeros(len(arealist)),
                                'gamma':np.zeros(len(arealist)),
                                'R':np.zeros(len(arealist)),
                                'L':np.zeros(len(arealist))})
        
        for j in range(0,len(arealist)):
            trim = df[df.locs == arealist[j]]
            trim = trim.reset_index(drop = True)
            
            for i in range(0, len(bins)-1):
                midpoint = (bins[i] + bins[i+1])/2
                vec.loc[i, 'ni'] = len(trim.az[(trim.az >= bins[i]) & (trim.az < bins[i+1])])
                vec.loc[i, 'w_comp'] = vec.ni[i] * math.sin(math.radians(midpoint))
                vec.loc[i, 'v_comp'] = vec.ni[i] * math.cos(math.radians(midpoint))
            
            results.loc[j, 'locality'] = arealist[j]
            results.loc[j, 'n'] = len(trim.az)
            results.loc[j, 'W'] = sum(vec.w_comp)
            results.loc[j, 'V'] = sum(vec.v_comp)
            if(math.degrees(math.atan(sum(vec.w_comp)/sum(vec.v_comp))) >= 0):
                results.loc[j, 'gamma'] = math.degrees(math.atan(sum(vec.w_comp)/sum(vec.v_comp)))
            else:
                results.loc[j, 'gamma'] = 360 + math.degrees(math.atan(sum(vec.w_comp)/sum(vec.v_comp)))
            results.loc[j, 'R'] = math.sqrt(sum(vec.w_comp) ** 2 + sum(vec.v_comp) ** 2)
            results.loc[j, 'L'] = math.sqrt(sum(vec.w_comp) ** 2 + sum(vec.v_comp) ** 2)/len(trim.az) * 100
                        
    else:    
        for i in range(0, len(bins)-1):
            midpoint = (bins[i] + bins[i+1])/2
            vec.loc[i, 'ni'] = len(az[(az >= bins[i]) & (az < bins[i+1])])
            vec.loc[i, 'w_comp'] = vec.ni[i] * math.sin(math.radians(midpoint))
            vec.loc[i, 'v_comp'] = vec.ni[i] * math.cos(math.radians(midpoint))
            
        results = pd.DataFrame({'n':[0.0],'W':[0.0],'V':[0.0],
                                'gamma':[0.0],'R':[0.0],'L':[0.0]})
        results.loc[0, 'n'] = len(az)
        results.loc[0, 'W'] = sum(vec.w_comp)
        results.loc[0, 'V'] = sum(vec.v_comp)
        if(math.degrees(math.atan(sum(vec.w_comp)/sum(vec.v_comp))) >= 0):
            results.loc[0, 'gamma'] = math.degrees(math.atan(sum(vec.w_comp)/sum(vec.v_comp)))
        else:
            results.loc[0, 'gamma'] = 360 + math.degrees(math.atan(sum(vec.w_comp)/sum(vec.v_comp)))
        results.loc[0, 'R'] = math.sqrt(sum(vec.w_comp) ** 2 + sum(vec.v_comp) ** 2)
        results.loc[0, 'L'] = math.sqrt(sum(vec.w_comp) ** 2 + sum(vec.v_comp) ** 2)/len(az) * 100
    
    if(sinuosity == True):
        sinuosity = []
        for x in results.L:
            def sinu(x): return math.exp(2.49 - 0.0475 * x +
                                         0.000234 * x ** 2)
            sinuosity.append(sinu(x))
        sinuosity = pd.DataFrame(sinuosity)
        sinuosity.rename(columns={0:"sinuosity"}, inplace=True)
        results = results.join(sinuosity)
    
    return results
    
#%% Coordinate conversions - Allmendinger et al., 2012

# Calculate pole to plane measured with dip azimuth/dip
def plane2pole(dip_az, dip):
    '''
    Calculate pole to plane measured with dip azimuth/dip.

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
    # Fix readings of 360 degrees
    if(dip_az == 360.0):
        dip_az = 0.0
        
    # Calculate azimuth of pole to plane
    if(dip_az >= 0.0) and (dip_az < 180.0):
        pol_az = dip_az + 180
    elif(dip_az >= 180.0) and (dip_az < 360.0):
        pol_az = dip_az - 180
    else:
        raise Exception('Dip azimuth must fall in range 0 <= x <= 360')
    
    # Calculate plunge of pole to plane
    if(dip <= 90.0) and (dip >= 0.0):
        pol_pln = 90 - dip
    else:
        raise Exception('Dip must fall in range 0 <= x <= 90')
        
    return pol_az, pol_pln

def z2p(x):
    '''
    Converts circular coordinates in radians to a value between 0 and 2 pi
    provided they are within 2 pi of one of the bounds.
    
    Parameters
    ----------
    x : Float
        Angle in radians.

    Returns
    -------
    y : Float
        Value between 0 and 2 pi corresponding to provided value.

    '''
    y = x
    if(x < 0.0):
        y = x + math.pi * 2
    elif(x >= (math.pi * 2)):
        y = x - math.pi * 2
    
    return y

# Convert line with trend and plunge from spherical to cartesian coordinates
def sph2cart(ln_trend, ln_plunge, degrees = True, plane = False):
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
    if(type(degrees) != bool):
        raise Exception('Degrees needs boolean value.')
    
    # Convert degrees to radians
    if(degrees == True):
        ln_trend = math.radians(ln_trend)
        ln_plunge = math.radians(ln_plunge)
    
    cd = math.sin(ln_plunge)
    ce = math.cos(ln_plunge) * math.sin(ln_trend)
    cn = math.cos(ln_plunge) * math.cos(ln_trend)
    
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
    if(type(degrees) != bool):
        raise Exception('Degrees needs boolean value.')       
    
    # Calculate plunge
    plunge = math.asin(cd)
    
    # Calculate trend
    if(cn == 0.0):
        if(ce < 0.0):
            trend = 3/2 * math.pi
        else:
            trend = math.pi/2
    else:
        trend = math.atan(ce/cn)
        if(cn < 0.0):
            trend += math.pi
        trend = z2p(trend)
        
    # Convert to degrees if asked
    if(degrees == True):
        trend = math.degrees(trend)
        plunge = math.degrees(plunge)
    return trend, plunge

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
        Fisher (1987) concentration factor.
    d99 : float
        99% uncertainty cone for vector mean.
    d95 : float
        95% uncertainty cone for vector mean.

    '''
    trends = np.array(trends)
    plunges = np.array(plunges)
    
    if(deg_in == True):
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
    R = math.sqrt(cn_sum**2 + ce_sum**2 + cd_sum**2)
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
        if((cd_sum < 0.0) and (lower == True)):
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
                
        if((R_ave >= 0.65) and (R_ave < 1.0)):
            afact = 100.0
            bfact = 1.0/(nlines - 1.0)
        
        # Calculate uncertainty cones
        d99 = math.acos(1.0 - ((nlines - R)/R) * (afact**bfact - 1.0))
        afact = 20.0
        d95 = math.acos(1.0 - ((nlines - R)/R) * (afact**bfact - 1.0))
    
    # Convert output to degrees
    if(deg_out == True):
        trend_ave = math.degrees(trend_ave)
        plunge_ave = math.degrees(plunge_ave)
        d99 = math.degrees(d99)
        d95 = math.degrees(d95)
        
    return trend_ave, plunge_ave, R_ave, conc, d99, d95

#%% Intersection of planes

def intersect(az1, d1, az2, d2, poles = False):
    '''
    Calculates the intersection of planes provided as dip azimuth/dip. Calls
    plane2pole and sph2cart.

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
    if(type(poles) != bool):
        raise Exception('poles needs boolean value.')
    
    dip1, dip2 = d1, d2
    
    # Convert planes to poles
    if(poles == False):
        az1, d1 = plane2pole(az1, d1)
        az2, d2 = plane2pole(az2, d2)
        
    # Convert lines to cartesian coordinates
    cn1, ce1, cd1 = sph2cart(az1, d1)
    cn2, ce2, cd2 = sph2cart(az2, d2)
    
    # Calculate theta
    theta = math.acos(ce1*ce2 + cn1*cn2 + cd1*cd2)
    
    if(theta != 0):
        i_ce = (cn1*cd2 - cd1*cn2)/math.sin(theta)
        i_cn = -(ce1*cd2 - cd1*ce2)/math.sin(theta)
        i_cd = (ce1*cn2 - cn1*ce2)/math.sin(theta)
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
        out_az = 450 - math.degrees(math.atan2(i_cn, i_ce))
    else:
        out_az = 90 - math.degrees(math.atan2(i_cn, i_ce))
    
    # Calculate intersection plunge
    if(dip1 == 90 and dip2 == 90):
        out_plunge = 90
    else:
        out_plunge = 90 - math.degrees(math.acos(i_cd))
    
    return out_az, out_plunge