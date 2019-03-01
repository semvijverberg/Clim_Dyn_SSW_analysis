#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 10:59:34 2019

@author: semvijverberg
"""
import os
# data analysis modules
import xarray as xr
from netCDF4 import num2date
import pandas as pd
import numpy as np
# plotting modules
import matplotlib as plt
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import matplotlib.colors as colors
import matplotlib as mpl


def import_dataset(filename, path_data):

    file_path = os.path.join(path_data, filename) 
    # Load in netcdf lazy, i.e. not storing into working memory yet only 
    # linking to the file
    ds = xr.open_dataset(file_path, decode_cf=True, decode_coords=True, 
                         decode_times=False)
    # Load into working memory, not the smart way
    # marray = np.squeeze(ncdf.to_array(file_path))
    numtime = ds['time']
    dates = num2date(numtime, units=numtime.units, 
                     calendar=numtime.attrs['calendar'])
    dates = pd.to_datetime(dates)
#    print('temporal frequency \'dt\' is: \n{}'.format(dates[1]- dates[0]))
    ds['time'] = dates
    return ds


def get_T10hpa_reversal(T10hpa):
    # select latitudes
    WMO_lats = slice(90,60)
    T10hpa_lats = T10hpa.sel(latitude=WMO_lats)
    # take zonal mean
    T10hpa_zm = T10hpa_lats.mean(dim='longitude')
    # get reversals
    T10hpa_zm.to_array().isel(time=0).plot()
    # calculate gradient with differentiate function of xarray
    mean_gradient = T10hpa_zm.to_array().diff('latitude').mean('latitude') 
    # note, the latitudes are order from 90N towards the equator, thus the 
    # gradient will be positive. To make it more intuitive, i.e. taking the 
    # perspective from the equator, we multiply times - 1
    # .squeeze() is done to remove empty dimensions of 1
    mean_gradient = mean_gradient.squeeze() * -1
    mask_T_grad_rev = mean_gradient.values > 0.
    return mean_gradient, mask_T_grad_rev

def get_u10hpa_reversal(u10hpa):
   
    # select latitudes
    WMO_lats = 60
    u10hpa_lats = u10hpa.sel(latitude=WMO_lats)
    # take zonal mean
    u10hpa_zm = u10hpa_lats.mean(dim='longitude')
    # load in array, .squeeze() is done to remove empty dimensions of 1
    u10hpa_zm = u10hpa_zm.to_array().squeeze()
    # plot one year
    datetime = pd.to_datetime(u10hpa_zm.time.values)
    u10hpa_zm.sel(time=onewinter(datetime)).plot()
    # get reversals
    mask_reversals = u10hpa_zm.values < 0.
    return u10hpa_zm, mask_reversals

def onewinter(datetime, winter=[1979,1980]):
    year1 = datetime.where(datetime.year==winter[0]).dropna()
    year2 = datetime.where(datetime.year==winter[1]).dropna()
    # last months of year1
    winter_p1 = year1.where(np.logical_or(year1.month==11, year1.month==12)).dropna()
    winter_p2 = year2.where(np.logical_or(year2.month==1, year2.month==2,
                                          year2.month==3)).dropna()
    whole_winter = np.concatenate([winter_p1, winter_p2])
    return whole_winter

def dates_shift_lag(all_dates, dates_SSW, lag):
    newdates = dates_SSW + pd.Timedelta('{}D'.format(lag))
    # shift leapdays one day
    leapdays = ((newdates.is_leap_year) & (newdates.month==2) & (newdates.day==29))==True 
    newdates.values[leapdays] = newdates[leapdays] - pd.Timedelta('1D')
    idx_outside_winter_months = [list(newdates).index(d) for d in newdates if d not in all_dates]
    newdates = np.delete(newdates, idx_outside_winter_months)
    
    return newdates


def define_events(all_dates, Ev_dates, min_dur, max_break):
    ''' This function requires:
        all_dates : all dates in your dataset
        Ev_dates  : your event dates
        min_dur   : minimal duration of the event
        max_break : allow for a short break, 0 is no break 
        It returns the central date of an event where min_dur > 2.'''
    event_idx = [list(all_dates).index(E) for E in Ev_dates]
    peak_o_thresh = np.zeros((all_dates.size))
    ev_num = 1
    # group events inter event time less than max_break
    for i in range(Ev_dates.size):
        if i < Ev_dates.size-1:
            curr_ev = event_idx[i]
            next_ev = event_idx[i+1]
        elif i == Ev_dates.size-1:
            curr_ev = event_idx[i]
            next_ev = event_idx[i-1]
                 
        if abs(next_ev - curr_ev) <= max_break:
            peak_o_thresh[curr_ev] = ev_num
        elif abs(next_ev - curr_ev) > max_break:
            peak_o_thresh[curr_ev] = ev_num
            ev_num += 1

    # remove events which are too short
    for i in np.arange(1, max(peak_o_thresh)+1):
        No_ev_ind = np.where(peak_o_thresh==i)[0]
        # if shorter then min_dur, then not counted as event
        if No_ev_ind.size < min_dur:
            peak_o_thresh[No_ev_ind] = 0


    peak_o_thresh[peak_o_thresh == 0] = np.nan
    Ev_labels = xr.DataArray(peak_o_thresh, coords=[all_dates], dims=['time'])
    Ev_labels_dropna = Ev_labels.dropna(how='all', dim='time').values
    new_Ev_dates = Ev_labels.dropna(how='all', dim='time').time.values
    Ev_array = np.array(new_Ev_dates, dtype=str)
    data = np.concatenate([Ev_array[:,None], Ev_labels_dropna[:,None]], axis=1)
    df = pd.DataFrame(data=data, 
                                index=new_Ev_dates,
                                columns=['date', 'event_number'])
    
    df_grouped = df.groupby('event_number')
    mean_idx = []
    for key in df_grouped.indices.keys():
        mean_idx.append( int(df_grouped.indices[key].mean()) )
    mean_idx.sort()
    new_Ev_dates[mean_idx]

    return pd.to_datetime(new_Ev_dates)

# =============================================================================
# Plotting functions
# =============================================================================

def convert_longitude(data):
    import numpy as np
    import xarray as xr
    lon_above = data.longitude[np.where(data.longitude > 180)[0]]
    lon_normal = data.longitude[np.where(data.longitude <= 180)[0]]
    # roll all values to the right for len(lon_above amount of steps)
    data = data.roll(longitude=len(lon_above))
    # adapt longitude values above 180 to negative values
    substract = lambda x, y: (x - y)
    lon_above = xr.apply_ufunc(substract, lon_above, 360)
    if lon_normal[0] == 0.:
        convert_lon = xr.concat([lon_above, lon_normal], dim='longitude')
    else:
        convert_lon = xr.concat([lon_normal, lon_above], dim='longitude')
    data['longitude'] = convert_lon
    return data

def xarray_plot(data, path='default', name = 'default', saving=False):
    # from plotting import save_figure
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import numpy as np
    if type(data) == type(xr.Dataset()):
        data = data.to_array().squeeze()
#    original
    fig = plt.figure( figsize=(10,6) ) 
    if len(data.longitude[np.where(data.longitude > 180)[0]]) != 0:
        if data.longitude.where(data.longitude==0).dropna(dim='longitude', how='all') == 0.:
            print('hoi')   
            data = convert_longitude(data)
    else:
        pass
    if data.ndim != 2:
        print("number of dimension is {}, printing first element of first dimension".format(np.squeeze(data).ndim))
        data = data[0]
    else:
        pass
    if 'mask' in list(data.coords.keys()):
        cen_lon = data.where(data.mask==True, drop=True).longitude.mean()
        data = data.where(data.mask==True, drop=True)
    else:
        cen_lon = data.longitude.mean().values
    proj = ccrs.PlateCarree(central_longitude=cen_lon)
    ax = fig.add_subplot(111, projection=proj)
    ax.coastlines()
    vmin = np.round(float(data.min())-0.01,decimals=2) 
    vmax = np.round(float(data.max())+0.01,decimals=2) 
    vmin = -max(abs(vmin),vmax) ; vmax = max(abs(vmin),vmax)
    # ax.set_global()
    if 'mask' in list(data.coords.keys()):
        data.copy().where(data.mask==True).plot.pcolormesh(ax=ax, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True,
                             vmin=vmin, vmax=vmax,
                             cbar_kwargs={'orientation' : 'horizontal'})
    else:
        data.plot.pcolormesh(ax=ax, cmap=plt.cm.RdBu_r,
                             transform=ccrs.PlateCarree(), add_colorbar=True,
                             vmin=vmin, vmax=vmax, 
                             cbar_kwargs={'orientation' : 'horizontal'})
    if saving == True:
        save_figure(data, path=path)
    plt.show()

def save_figure(data, path):
    import os
#    if 'path' in locals():
#        pass
#    else:
#        path = '/Users/semvijverberg/Downloads'
    if path == 'default':
        path = '/Users/semvijverberg/Downloads'
    else:
        path = path
    import datetime
    today = datetime.datetime.today().strftime("%d-%m-%y_%H'%M")
    if type(data.name) is not type(None):
        name = data.name.replace(' ', '_')
    if 'name' in locals():
        print('input name is: {}'.format(name))
        name = name + '.jpeg'
        pass
    else:
        name = 'fig_' + today + '.jpeg'
    print(('{} to path {}'.format(name, path)))
    plt.savefig(os.path.join(path,name), format='jpeg', dpi=300, bbox_inches='tight')
    

def finalfigure(xrdata, file_name, kwrgs):
    #%%
    map_proj = ccrs.NearsidePerspective(central_longitude=-60, central_latitude=80)  
    lons = xrdata.longitude.values
    lats = xrdata.latitude.values
    strvars = [' {} '.format(var) for var in list(xrdata.dims)]
    var = [var for var in strvars if var not in ' longitude latitude '][0] 
    var = var.replace(' ', '')
    g = xr.plot.FacetGrid(xrdata, col=var, col_wrap=kwrgs['column'], sharex=True,
                      sharey=True, subplot_kws={'projection': map_proj},
                      aspect= 1, size=4)
    figwidth = g.fig.get_figwidth() ; figheight = g.fig.get_figheight()

    lon_tick = xrdata.longitude.values
#    lon_tick[lon_tick > 180] -= 360
    
    longitude_labels = np.linspace(np.min(lon_tick), np.max(lon_tick), 6, dtype=int)
    longitude_labels = np.array(sorted(list(set(np.round(longitude_labels, -1)))))

#    longitude_labels = np.concatenate([ longitude_labels, [longitude_labels[-1]], [180]])
#    longitude_labels = [-150,  -70,    0,   70,  140, 140]
    latitude_labels = np.linspace(xrdata.latitude.min(), xrdata.latitude.max(), 4, dtype=int)
    latitude_labels = sorted(list(set(np.round(latitude_labels, -1))))
    
    g.set_ticks(max_xticks=5, max_yticks=5, fontsize='small')
    g.set_xlabels(label=[str(el) for el in longitude_labels])

    
    if kwrgs['clevels'] == 'default':
        vmin = np.round(float(xrdata.min())-0.01,decimals=2) ; vmax = np.round(float(xrdata.max())+0.01,decimals=2)
        clevels = np.linspace(-max(abs(vmin),vmax),max(abs(vmin),vmax),17) # choose uneven number for # steps
    else:
        vmin=kwrgs['vmin']
        vmax=kwrgs['vmax']
        
        clevels = np.linspace(vmin, vmax,kwrgs['steps'])

    cmap = kwrgs['cmap']
    
    n_plots = xrdata[var].size
    for n_ax in np.arange(0,n_plots):
        ax = g.axes.flatten()[n_ax]
#        print(n_ax)
        plotdata = extend_longitude(xrdata[n_ax]).squeeze().drop('ds')
        im = plotdata.plot.pcolormesh(ax=ax, cmap=cmap,
                               transform=ccrs.PlateCarree(),
                               subplot_kws={'projection': map_proj},
                                levels=clevels, add_colorbar=False)
        ax.coastlines(color='black', alpha=0.5)
        
#        ax.set_extent([lons[0], lons[-1], lats[0], lats[-1]], ccrs.PlateCarree())
        if kwrgs['subtitles'] == None:
            pass
        else:
            fontdict = dict({'fontsize'     : 18,
                             'fontweight'   : 'bold'})
            ax.set_title(kwrgs['subtitles'][n_ax], fontdict=fontdict, loc='center')

        
        if 'ax_text' in kwrgs.keys():
            ax.text(0.0, 1.01, kwrgs['ax_text'][n_ax],
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes,
            color='black', fontsize=15)
            
        if map_proj.proj4_params['proj'] in ['merc', 'eqc']:
#            print(True)
            ax.set_xticks(longitude_labels[:-1], crs=ccrs.PlateCarree())
            ax.set_xticklabels(longitude_labels[:-1], fontsize=12)
            lon_formatter = cticker.LongitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            
            ax.set_yticks(latitude_labels, crs=ccrs.PlateCarree())
            ax.set_yticklabels(latitude_labels, fontsize=12)
            lat_formatter = cticker.LatitudeFormatter()
            ax.yaxis.set_major_formatter(lat_formatter)
            ax.grid(linewidth=1, color='black', alpha=0.3, linestyle='--')
            ax.set_xlabel('')
            ax.set_ylabel('')
            
            
        else:
            pass
        
    g.fig.text(0.5, 0.99, kwrgs['title'], fontsize=20,
               fontweight='heavy', transform=g.fig.transFigure,
               horizontalalignment='center',verticalalignment='top')
    
    if 'adj_fig_h' in kwrgs.keys():
        g.fig.set_figheight(figheight*kwrgs['adj_fig_h'], forward=True)
    if 'adj_fig_w' in kwrgs.keys():
        g.fig.set_figwidth(figwidth*kwrgs['adj_fig_w'], forward=True)

    if 'cbar_vert' in kwrgs.keys():
        cbar_vert = (figheight/40)/(n_plots*2) + kwrgs['cbar_vert']
    else:
        cbar_vert = (figheight/40)/(n_plots*2)
    if 'cbar_hght' in kwrgs.keys():
        cbar_hght = (figheight/40)/(n_plots*2) + kwrgs['cbar_hght']
    else:
        cbar_hght = (figheight/40)/(n_plots*2)
    if 'wspace' in kwrgs.keys():
        g.fig.subplots_adjust(wspace=kwrgs['wspace'])
    if 'hspace' in kwrgs.keys():
        g.fig.subplots_adjust(hspace=kwrgs['hspace'])
    if 'extend' in kwrgs.keys():
        extend = kwrgs['extend'][0]
    else:
        extend = 'neither'

    cbar_ax = g.fig.add_axes([0.25, cbar_vert, 
                                  0.5, cbar_hght], label='xbar')

    if 'clim' in kwrgs.keys(): #adjust the range of colors shown in cbar
        cnorm = np.linspace(kwrgs['clim'][0],kwrgs['clim'][1],11)
        vmin = kwrgs['clim'][0]
    else:
        cnorm = clevels
    norm = colors.BoundaryNorm(boundaries=cnorm, ncolors=256)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, orientation='horizontal', 
                 extend=extend, norm=norm)
    if 'cticks_center' in kwrgs.keys():
        cbar.set_ticks(clevels + 0.5)
        cbar.set_ticklabels(clevels+1, update_ticks=True)
        cbar.update_ticks()
    
    if 'extend' in kwrgs.keys():
        if kwrgs['extend'][0] == 'min':
            cbar.cmap.set_under(cbar.to_rgba(vmin))
    cbar.set_label(xrdata.attrs['units'], fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    #%%
    if kwrgs['savefig'] != False:
        g.fig.savefig(file_name ,dpi=250, frameon=True)
    
    return

def extend_longitude(data):
    import xarray as xr
    import numpy as np
    plottable = xr.concat([data, data.sel(longitude=data.longitude[:1])], dim='longitude').to_dataset(name="ds")
    plottable["longitude"] = np.linspace(0,360, len(plottable.longitude))
    plottable = plottable.to_array(dim='ds')
    return plottable