#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extracting SSWs and analyzing impact on SLP and T2m

Created on Thu Feb 28 09:45:13 2019

@author: semvijverberg
"""
#%%
# =============================================================================
# To get right path connections
# =============================================================================
import inspect, os, sys
curr_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) # script directory
script_dir = os.path.join(curr_dir)
# To link Python modules in current folder to this script
os.chdir(script_dir)
sys.path.append(script_dir)
# =============================================================================
# Load in modules
# =============================================================================
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
# import own python code with functions
import functions

# =============================================================================
# Fill in paths to data
# =============================================================================
path_data = os.path.join(curr_dir, 'data')

filename_T10hpa = 't_10hpa_1979-2017_01okt_31mar_dt-1days_2.5deg_NH.nc'
filename_u10hpa = 'u_10hpa_1979-2017_01okt_31mar_dt-1days_2.5deg_NH.nc'
filename_T2m    = 't2m_1979-2017_01okt_31mar_dt-1days_2.5deg_NH.nc'
filename_SLP    = 'SLP_1979-2017_01okt_31mar_dt-1days_2.5deg_NH.nc'


# =============================================================================
# load data to extract SSW
# =============================================================================
T10hpa = functions.import_dataset(filename_T10hpa, path_data)
u10hpa = functions.import_dataset(filename_u10hpa, path_data)
all_dates = pd.to_datetime(u10hpa.time.values)
#%%
# =============================================================================
# extract the SSW
# =============================================================================

''' Assignment 1:
    Write a function which extracts the dates of the SSW in datetime format.
    See powerpoint slide by Dim for the definition.
    After you accomplished step 1, write it into a function get_u10hpa_reversal()
    and copy it to function.py. '''

# step 1: find where the zonal mean wind is negative (westward, i.e. toward the west)
    
u10hpa_zm, mask_reversals = functions.get_u10hpa_reversal(u10hpa) # output the xarray of 
# zonal mean wind, and a boolean mask of where the reversals occur, True when 
# reversals occur, False everywhere else

mean_gradient, mask_T_grad_rev = functions.get_T10hpa_reversal(T10hpa) # output the xarray of 
# mean zonal mean meridional gradient over time, and the corresponding boolean mask

# step 2: find where the meridional temperature gradient reverses


# When are both of these conditions true?
# mask_SSW = 
dates_SSW = all_dates.where(mask_SSW).dropna() # plug in a numpy mask of 
                # booleans where SSWs are True, and non-SSWs are False
#%%
# =============================================================================
# load data to investigate impact
# =============================================================================
T2m = functions.import_dataset(filename_T2m, path_data)
SLP = functions.import_dataset(filename_SLP, path_data)


''' Assignment 2:
    Start investigating the composites, use the dates_shift_lag function 
    in function.py, use the xarray_plot() function to have a first look '''


#%%
# step 2: investigate and compare the composites at different lags
lags = [-40, -20, 0, 20, 40, 60]

''' Assigignment 3:
    Create a new xarray with the dimensions ['lag','latitude','longitude'].
    Loop over all lags, and store the corresponding composite mean (2D) maps 
    into your new xarray at the right lag. 
    
    The new xarray should look like this:
        xr.DataArray(data=empty_nparray, coords=[lags, lats, lons], 
                          dims=['lag','latitude','longitude'], 
                          name=name,
                          attrs={'units':'units'})
    The Coords contain the actual values, i.e. labels, and dims are simply the 
    dimensions. Also add an appropriate name, you will need that for the plotting
    function.
    
    When you succeed, create a function making_composites() out of you code, such that you can 
    quickly analyze these lags for different variables. 
    
    '''

# I give you the plotting function such that you can already make some nice figures
kwrgs = dict( {'title' : xarray.name, 'clevels' : 'notdefault', 'steps':17,
                'vmin' : -2, 'vmax' : 2, 
               'cmap' : plt.cm.RdBu_r, 'column' : 3, 'subtitles' : None,
               'savefig':True } )

file_name = os.path.join(script_dir,'figures',name+'.pdf')

functions.finalfigure(xarray, file_name, kwrgs)


#%%
# =============================================================================
# Now extract only persistent SSWs and investigate impact on surface temperature
# =============================================================================
''' Assignment 4:
    Now extract SSW events with a certain persistence of choice by using the
    define_events() function and analyze the results. 
    What is your conclusion? '''

    
functions.define_events?



#%%
# =============================================================================
# quantify frequency and trends
# =============================================================================
''' Assignment 5:
    Write a function which return the frequency, e.g. amount of SSWs per year.
    Is there a significant trend in any of the timeseries.
    You can use import statsmodels.api as sm and then 
    OLS = sm.OLS(freq, years).fit() to make the fit.
    
    view all results by typing OLS.summary2()
    
    '''
def get_freq_per_yr(SSW_dates, all_dates):
    '''good luck'''
    
    return freq_per_year

years           = np.unique(all_dates.year)
freq_SSW_normal = get_freq_per_yr(dates_SSW, all_dates)
freq_SSW_7day   = get_freq_per_yr(dates_7day_SSW, all_dates)
freq_SSW_14day  = get_freq_per_yr(dates_14day_SSW, all_dates)
freq_SSW_21day  = get_freq_per_yr(dates_21day_SSW, all_dates)


plt.figure(figsize=(12,8))
plt.plot(years, freq_SSW_normal, label='SSW')
plt.plot(years, freq_SSW_7day, label='7day-SSW')
plt.plot(years, freq_SSW_14day, label='14day-SSW')
plt.plot(years, freq_SSW_21day, label='21day-SSW')
plt.legend()