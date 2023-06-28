# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 13:21:23 2021

@author: RGHS

Functions for calculating properties of paleorivers.

Sources:
    - Paleosinuosity:
        - Ghosh, 2000: Estimation of channel sinuosity from paleocurrent data: a
            method using fractal geometry
        - Le Roux, 1992: Determining the Channel Sinuosity of Ancient Fluvial
            Systems from Paleocurrent Data
        - Le Roux, 1994: The Angular Deviation of Paleocurrent Directions as
            Applied to the Calculation of Channel Sinuosities
        - Le Roux, 2001: Estimation of Channel Sinuosity from Paleocurrent Data: 
            A Method Using Fractal Geometry: Discussion
    - Fulcrum approach to source to sink:
        - Holbrook and Wanas, 2014
        - Wright and Parker, 2004
        - Garcia and Parker, 1991
        - Van Rijn, 1984
        - LeClair and Bridge, 2001
        
"""

#%% Import packages

import numpy as np
import pandas as pd
import warnings
import graphs as gr

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
        Returns calculated sinuosities in data frame with column headings "S"
        and "S_max".

    '''
    S = np.exp(2.49 - 0.0475 * L + 0.000234 * np.square(L))
    S_max = np.exp(4.2416 - 0.0801 * L + 0.0004 * np.square(L))
    
    if((hasattr(S, '__len__')) and (len(S) > 1)):
        sinuosity = pd.DataFrame({'S': S,'S_max': S_max})
    else:
        sinuosity = pd.DataFrame({'S': S,'S_max': S_max}, index=[0])
    
    return sinuosity

#%% Le Roux 1992, 1994, 2001

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
    
    if(hasattr(L, '__len__') is False):
        L = np.array([L])
    else:
        L = np.array(L)
    phi = np.zeros(len(L))
    
    # Get indices for each calculation
    eq6 = np.where(L <= 85)
    eq7 = np.where((L > 85) & (L < 98))
    eq8 = np.where(L > 98)
    
    # Do the calculation
    phi[eq6] = 5e-8*L[eq6]**4 - 9e-5*L[eq6]**3 + 1.07e-2*L[eq6]**2 - 1.7402*L[eq6] + 179.95
    phi[eq7] = -2.7e-3*(100 - L[eq7])**4 + 9.53e-2*(100 - L[eq7])**3 - 1.23*(100 - L[eq7])**2 + 9.3963*(100 - L[eq7]) + 5.9
    phi[eq8] = 1.2183*(100 - L[eq8])**3 - 5.8723*(100 - L[eq8])**2 + 14.654*(100 - L[eq8]) + 5
    phi *= 2 # Multiply by two to get operational range
    
    r_phi = np.radians(phi)
    S = (r_phi/2)/(np.sin(r_phi/2))
    
    if((hasattr(S, '__len__')) and (len(S) > 1)):
        sinuosity = pd.DataFrame({'Phi': phi,'P': S})
    else:
        sinuosity = pd.DataFrame({'Phi': phi,'P': S}, index=[0])
    
    return sinuosity

#%% Fulcrum approach calculations, Holbrook and Wanas, 2014

def fulcrum(d16, d50, d84, d90, sm, duration, tbd, b, wc = 'single', depth = 0,
            grainsize_m = False, hm_hbf = 8, ca_method = 'vr',
            crit_mob_param = 1, crit_mob_manual = 0,
            all_parameters=False):
    '''
    UNDER CONSTRUCTION
    
    Performs calculations for source to sink fulcrum method as developed
    by Holbrook and Wanas (2014), using grain size and foreset height.

    Parameters
    ----------
    d16 : float
        Grainsize of the 16th percentile (m or mm).
    d50 : float
        Grainsize of the 50th percentile (m or mm).
    d84 : float
        Grainsize of the 84th percentile (m or mm).
    d90 : float
        Grainsize of the 90th percentile (m or mm).
    sm : float
        Mean foreset height (m).
    duration : float
        Lifetime of the fluvial system (years).
    tbd : float
        Frequency of bankfull discharge events (days).
    b : float
        Multiplier expressing proportion of sediment moved during bankfull events.
    wc : float or string, optional
        Width of the channel. Accepts numerical values for measured channel widths.
        If "single" or "braid" is provided, width will be calculated using H_bf and
        appropriate relationships from Bridge and Mackey (1993) or Leopold and Maddock (1953).
        The default is "single".
    depth : float, optional
        Bankfull depth in m. If left as 0, bankfull depth will be calculated from
        sm.
    grainsize_m : bool, optional
        Designates whether the grainsize is input in m (True) or mm (False). The default is False.
    hm_hbf : float, optional
        The multiplier specifying ratio between Hm and Hbf, constrained by 
        LeClair and Bridge (2001) as falling between 6 and 10. The default is 8.
    ca_method : str, optional
        Specify method used to calculate c_a as follows:
            - 'wp': Wright and Parker, 2004
            - 'gp': Garcia and Parker, 1991
            - 'vr': Van Rijn, 1984
    large_channel : bool, optional
        Specifies whether to calculate ca using the relationship from Wright
        and Parker (2004) for large channels (True) or the standard ca
        calculation (False). The default is False.
    crit_mob_param : int, optional
        Specifies which relationship to use when calculating the critical
        mobility parameter. These are as follows:
            - 1: Van Rijn
            - 2: Shields
            - 3: Engelund
            - 4: Bagnold
        The default is 1. Ignored if crit_mob_manual is specified.
    crit_mob_manual : float, optional
        Manually assigns a value to the critical mobility parameter. If left
        as default (0), critical mobility parameter is calculated by method
        specified by crit_mob_param.

    Returns
    -------
    results : pandas dataframe
        Dataframe containing the following parameters:
            - Bbf (channel width)
            - Hbf (channel depth)
            - S (slope)
            - U-bar (mean flow velocity)
            - m (critical mobility parameter)
            - Qs (Bankfull suspendend discharge)
            - Qtbf (Bankfull bedload discharge)
            - Qmas (Mean annual sediment discharge)
            - Duration (Fluvial system lifespan)
            - Vt (Total sediment discharge across system lifespan)

    '''
    if((hm_hbf < 6) or (hm_hbf > 10)):
        warnings.warn('LeClair and Bridge found Hm/Hbf lies between 6 and 10. Consider choosing a new value.')
    # Generate constants
    tau_bf = 1.86       # Shields number for dimensionless shear stress
    R = 1.65            # Dimensionless submerged density of quartz
    g = 9.81            # Gravitational acceleration (m s^-2)
    tau_c = 0           # Empirical constant
    alpha_EH = 0.05     # Empirical constant
    nt = 2.5            # Empirical constant
    phi_s = 1           # Empirical constant
    
    v = 1.0035e-6       # Water kinematic visc @ 20C, units m^2 s^-1
    s = 2650/1000       # (Density of suspended material)/(Fluid density), units (-)
    c0 = 0.65           # Empirical constant
    k = 0.4             # Von Karman constant
    
    if(grainsize_m is False):
        d16 = d16/1000
        d50 = d50/1000
        d84 = d84/1000
        d90 = d90/1000
    
    # Deal with channel width and depth
    # Assign channel depth
    if(depth == 0):
        beta_1 = sm/1.8
        hm = 5.3 * beta_1 + 0.001 * beta_1**2
        hbf = hm_hbf * hm
    else:
        hbf = depth

    # Assign channel width
    try:
        bbf = wc/1
    except:
        if(wc == 'single'):
            bbf = 8.8 * hbf**1.82
        elif(wc == 'braid'):
            bbf = 42 * hbf**1.11
        else:    
            raise Exception('wc accepts numerical arguments, "single" or "braid" only.')
    
    S = (tau_bf * (R * d50))/hbf
    
    delta = hbf/8
    lam = 7.3 * hbf
    psi = delta / lam
    ks = 3 * d90 + 1.1 * delta * (1 - np.e ** (-25 * psi))
    cf = (8.32 * (hbf/ks)**(1/6))**-2
    
    Q_bf = ((g * hbf * S)/cf * (bbf ** 2 * hbf ** 2)) ** (1/2)
    
    alpha_t = alpha_EH/cf
    Q_tbf = bbf * (R * g * d50) ** (1/2) * d50 * alpha_t * (phi_s * tau_bf - tau_c) ** nt
    
    D_st = d50 * (((s-1)*g)/v**2)**(1/3)
    
    # Construct critical mobility space and select appropriate parameter
    if(crit_mob_manual == 0):
        cmpx, cmpy = gr.cmp(crit_mob_param)    
        try:
            lower = np.max(np.where(cmpx <= D_st))
        except:
            raise Exception(f'D_st too small for chosen critical mobility parameter range. D_st = {D_st}')
        
        try:
            upper = np.max(np.where(cmpx > D_st))
        except:
            raise Exception(f'D_st too small for chosen critical mobility parameter range. D_st = {D_st}')
        a1 = (D_st, 0)
        a2 = (D_st, np.max(cmpy))
        b1 = (cmpx[lower],cmpy[lower])
        b2 = (cmpx[upper],cmpy[upper])
        
        m = gr.get_intersect(a1, a2, b1, b2)[1]

    else:
        m = crit_mob_manual
    
    u_st_cr = (m * ((s-1)*g*d50))**(1/2)
    
    cz = g ** (1/2) * 8.1 * (hbf/ks) ** (1/6)
    Rh = (hbf * bbf)/(2*hbf + bbf)
    u_bar = cz * (Rh * S) ** (1/2)
    u_st_p = (g ** (1/2))/cz * u_bar
    T = (u_st_p ** 2 - u_st_cr ** 2)/u_st_cr ** 2
    
    sig_s = 0.5 * (d84/d50 + d16/d50)
    Ds = (1 + 0.011 * (sig_s - 1) * (T - 25)) * d50
    
    if(Ds < 0.0001):
        w_s = (1/18) * (((s - 1) * g * Ds**2)/v)
    elif(Ds > 0.001):
        w_s = 1.1 * ((s - 1) * g * Ds) ** (1/2)
    else:
        w_s = 10 * v/Ds * ((1 + ((0.01 * (s - 1) * g * Ds ** 3) / (v ** 2)) ** (1/2)) - 1)
        
    u_st = (g * hbf * S) ** (1/2)
    
    a = 0.05 * hbf
    if((ca_method == 'wp') or (ca_method == 'gp')):
        Rep = ((R * g * d50) ** (1/2) * d50) / v
        if(ca_method == 'wp'):
            # Wright and Parker, 2004
            A = 5.7e-7
            Zu = (u_st/w_s) * (Rep**0.6) * (S**0.07)
        else:
            # Garcia and Parker, 1991
            A = 1.3e-7
            Zu = (u_st/w_s) * (Rep**0.6)
        ca = (A * (Zu**5))/(1 + (A/0.3) * (Zu**5))
    elif(ca_method == 'vr'):
        A=np.nan
        Zu=np.nan
        # Van Rijn, 1984
        ca = 0.015 * (d50/a) * T**1.5/D_st**0.3
    else:
        raise Exception('ca_method accepts "wp","gp" or "vr" only.')
    
    beta = 1 + 2 * (w_s/u_st)**2 
    
    phi = 2.5 * (w_s/u_st)**0.8 * (ca/c0)**0.4
    
    Z = w_s/(beta * k * u_st)
    Zp = Z + phi
    
    F = ((a/hbf)**Zp - (a/hbf)**1.2)/((1-(a/hbf))**Zp * (1.2-Zp))
    
    qs = F * u_bar * hbf * ca       # Calculate bankfull suspended specific discharge
    Qs = qs * bbf                   # Calculate bankfull suspended discharge
    
    # Calculate discharge figures
    tbd = tbd * 24 * 60 * 60            # Convert days to seconds
    
    # Qma = 3.102 * Q_bf + 56.581       # Calculate mean annual water discharge (seems wrong, investigate)
    
    Qmas = (Qs + Q_tbf) * tbd * b       # Mean annual sediment discharge (m^3)
    Qmas_km = Qmas * 1e-9               # Mean annual sediment discharge (km^3)
    
    total_km = Qmas_km * duration
    
    # Construct results table
    if all_parameters is False:
        results = pd.DataFrame({'d16':[d16],
                            'd50':[d50],
                            'd84':[d84],
                            'd90':[d90],
                            'mean_foreset_height':[sm],
                            'channel_width':[bbf],
                            'channel_depth':[hbf],
                            'slope':[S],
                            'mean_flow_velocity':[u_bar],
                            'critical_mobility_parameter':[m],
                            'bankfull_discharge':[Q_bf],
                            'bankfull_suspended_discharge':[Qs],
                            'bankfull_bedload_discharge':[Q_tbf],
                            'mean_annual_sediment_discharge':[Qmas],
                            'bankfull_interval':[tbd],
                            'bankfull_multiplier':[b],
                            'duration':[duration],
                            'total_sediment_discharge':[total_km]})
    elif all_parameters is True:
        results = pd.DataFrame({'d16':[d16],
                            'd50':[d50],
                            'd84':[d84],
                            'd90':[d90],
                            'mean_foreset_height':[sm],
                            'channel_width':[bbf],
                            'channel_depth':[hbf],
                            'slope':[S],
                            'mean_flow_velocity':[u_bar],
                            'critical_mobility_parameter':[m],
                            'bankfull_discharge':[Q_bf],
                            'bankfull_suspended_discharge':[Qs],
                            'bankfull_bedload_discharge':[Q_tbf],
                            'mean_annual_sediment_discharge':[Qmas],
                            'bankfull_interval':[tbd],
                            'bankfull_multiplier':[b],
                            'duration':[duration],
                            'total_sediment_discharge':[total_km],
                            'tau_bf':[tau_bf],
                            'R':[R],
                            'g':[g],
                            'tau_c':[tau_c],
                            'alpha_EH':[alpha_EH],
                            'nt':[nt],
                            'phi_s':[phi_s],
                            'v':[v],
                            's':[s],
                            'c0':[c0],
                            'k':[k],
                            'Hbf':[hbf],
                            'Bbf':[bbf],
                            'S':[S],
                            'delta':[delta],
                            'lam':[lam],
                            'psi':[psi],
                            'ks':[ks],
                            'cf':[cf],
                            'alpha_t':[alpha_t],
                            'D_st':[D_st],
                            'm':[m],
                            'u_st_cr':[u_st_cr],
                            'cz':[cz],
                            'Rh':[Rh],
                            'u_bar':[u_bar],
                            'u_st_p':[u_st_p],
                            'T':[T],
                            'sig_s':[sig_s],
                            'Ds':[Ds],
                            'w_s':[w_s],
                            'u_st':[u_st],
                            'a':[a],
                            'Rep':[Rep],
                            'Zu':[Zu],
                            'A':[A],
                            'ca':[ca],
                            'beta':[beta],
                            'phi':[phi],
                            'Z':[Z],
                            'Zp':[Zp],
                            'F':[F]})
    
    return results

def fulcrumEmptyDf(length, all_parameters=False):
    '''
    Produces a blank dataframe for the fulcrum method calculations. Useful for
    creating a placeholder dataframe to put the output of a loop into. See docs 
    for fulcrum for further details.

    Parameters
    ----------
    length : int
        Number of blank rows desired.

    Returns
    -------
    df : pd.DataFrame
        Pandas data frame with zeros for desired number of rows.

    '''
    blank_row = np.zeros(int(length))
    if all_parameters is False:
        df = pd.DataFrame({'d16':blank_row,
                    'd50':blank_row,
                    'd84':blank_row,
                    'd90':blank_row,
                    'mean_foreset_height':blank_row,
                    'channel_width':blank_row,
                    'channel_depth':blank_row,
                    'slope':blank_row,
                    'mean_flow_velocity':blank_row,
                    'critical_mobility_parameter':blank_row,
                    'bankfull_discharge':blank_row,
                    'bankfull_suspended_discharge':blank_row,
                    'bankfull_bedload_discharge':blank_row,
                    'mean_annual_sediment_discharge':blank_row,
                    'bankfull_interval':blank_row,
                    'bankfull_multiplier':blank_row,
                    'duration':blank_row,
                    'total_sediment_discharge':blank_row})
    elif all_parameters is True:
        df = pd.DataFrame({'d16':blank_row,
                            'd50':blank_row,
                            'd84':blank_row,
                            'd90':blank_row,
                            'mean_foreset_height':blank_row,
                            'channel_width':blank_row,
                            'channel_depth':blank_row,
                            'slope':blank_row,
                            'mean_flow_velocity':blank_row,
                            'critical_mobility_parameter':blank_row,
                            'bankfull_discharge':blank_row,
                            'bankfull_suspended_discharge':blank_row,
                            'bankfull_bedload_discharge':blank_row,
                            'mean_annual_sediment_discharge':blank_row,
                            'bankfull_interval':blank_row,
                            'bankfull_multiplier':blank_row,
                            'duration':blank_row,
                            'total_sediment_discharge':blank_row,
                            'tau_bf':blank_row,
                            'R':blank_row,
                            'g':blank_row,
                            'tau_c':blank_row,
                            'alpha_EH':blank_row,
                            'nt':blank_row,
                            'phi_s':blank_row,
                            'v':blank_row,
                            's':blank_row,
                            'c0':blank_row,
                            'k':blank_row,
                            'Hbf':blank_row,
                            'Bbf':blank_row,
                            'S':blank_row,
                            'delta':blank_row,
                            'lam':blank_row,
                            'psi':blank_row,
                            'ks':blank_row,
                            'cf':blank_row,
                            'alpha_t':blank_row,
                            'D_st':blank_row,
                            'm':blank_row,
                            'u_st_cr':blank_row,
                            'cz':blank_row,
                            'Rh':blank_row,
                            'u_bar':blank_row,
                            'u_st_p':blank_row,
                            'T':blank_row,
                            'sig_s':blank_row,
                            'Ds':blank_row,
                            'w_s':blank_row,
                            'u_st':blank_row,
                            'a':blank_row,
                            'Rep':blank_row,
                            'Zu':blank_row,
                            'A':blank_row,
                            'ca':blank_row,
                            'beta':blank_row,
                            'phi':blank_row,
                            'Z':blank_row,
                            'Zp':blank_row,
                            'F':blank_row,
                            })
    
    return df