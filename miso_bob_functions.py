# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 17:34:58 2021

@author: jayes
"""

import numpy as np
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import metpy.calc as mpcalc
import scipy.ndimage as ndimage
from matplotlib.colors import LinearSegmentedColormap
from datetime import datetime

cmap_grey = LinearSegmentedColormap.from_list('cmap',['white','grey'],   N =  64)


def get_hours(t1,t2):
    t1 = datetime(t1[0],t1[1],t1[2],t1[3])
    t2 = datetime(t2[0],t2[1],t2[2],t2[3])
    difference = t2-t1
    return int(difference.total_seconds()/60**2)

def plot_olr_anomaly(lats, lons, data, xlim,ylim,d):
    plt.figure(figsize=(12,3))
                 
    ################### Figure 1 ###############################
    ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.coastlines('10m', linewidth=1.2)
    gl = ax.gridlines(ccrs.PlateCarree(), draw_labels=True, linewidths=0.25, color='white', alpha=0, linestyle='none')
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-180,210,30))
    gl.ylocator = mticker.FixedLocator(np.arange(-90,100,15))
    gl.xlines = False
    gl.ylines = False
    
    levels = np.linspace(-100, 100,10)
    rain_map = plt.contourf(lons, lats, data, levels, cmap = 'bwr', extend = 'both', transform=ccrs.PlateCarree())
    cbar = plt.colorbar(rain_map, ax=ax, extend = 'both')
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    plt.title(str(d) + ' June 2018', fontsize = 16)
    cbar.set_label('OLR anomaly (W/m$^{2}$)',fontsize=14, fontweight = 'bold')
    cbar.ax.tick_params(labelsize = 16, grid_alpha = 0.5)
    
    plt.savefig('D:\\Projects\\MISO-BOB\\Data\\NOAA_OLR\\' + str(d) + 'June2018.png')
    plt.close()
    plt.clf()
    


def get_cape_cin(df):
    """
    Parameters
    ----------
    df : Dataframe of the sounding

    Returns
    -------
    cape : Pint
        CAPE value of the sounding.
    cin  : Pint
        CIN value of the sounding.
    """
    p  = df['Pressure']
    T  = df['Temperature']
    RH = df['RH']
    Td = mpcalc.dewpoint_from_relative_humidity(T, RH/100).to('degC')
    Td = np.round(Td,2)
    cape,cin = mpcalc.most_unstable_cape_cin(p, T, Td)
    
    return cape,cin

def plot_tpw(lats, lons, data, xlim,ylim,d,t):
    plt.figure(figsize=(15,6))
                 
    ################### Figure 1 ###############################
    ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.coastlines('10m', linewidth=1.2)
    gl = ax.gridlines(ccrs.PlateCarree(), draw_labels=True, linewidths=0.25, color='black', alpha=0, linestyle='none')
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-180,210,10))
    gl.ylocator = mticker.FixedLocator(np.arange(-90,100,2))
    gl.xlines = False
    gl.ylines = False
    
    levels = np.arange(30, 71,5)
    tpw_map = plt.pcolormesh(lons, lats, data, vmin = 40, vmax = 70, cmap = 'viridis', alpha = 0.7, transform=ccrs.PlateCarree())
    cbar = plt.colorbar(tpw_map, ax=ax, extend = 'both')
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    plt.title( str(t).zfill(2) + "00 UTC " + str(d) + ' June 2018', fontsize = 16)
    cbar.set_label('TPW (mm)',fontsize=14, fontweight = 'bold')
    cbar.ax.tick_params(labelsize = 16, grid_alpha = 0.5)
    path = "C:\\MISO-BOB\\Data\\TPW\\plots\\" 

    plt.savefig(path + str(d) + str(t) + 'June2018.png', dpi = 300)
    plt.close()
    plt.clf()
    
    
    
def plot_tpw_ir(lats1, lons1, data1,lats2,lons2,data2,u,v,dvu,xlim,ylim, month_str,d,t,locs,df):
    plt.figure(figsize=(15,6))
                 
    ################### Figure 1 ###############################
    ax = plt.subplot(1,1,1, projection=ccrs.PlateCarree())

    ax.coastlines('10m', linewidth=1.2)
    gl = ax.gridlines(ccrs.PlateCarree(), draw_labels=True, linewidths=0.25, color='black', alpha=0, linestyle='none')
    gl.xlabels_top = gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator(np.arange(-180,210,5))
    gl.ylocator = mticker.FixedLocator(np.arange(-90,100,2))
    gl.xlines = False
    gl.ylines = False
    
    tp_levels = np.arange(-10,1,1)
    tpw_map = ax.contourf(lons1, lats1, data1, tp_levels, cmap = 'viridis', extend = 'both',alpha = 0.7, transform=ccrs.PlateCarree())
    cbar = plt.colorbar(tpw_map, ax=ax, orientation='horizontal', shrink=0.35, extend = 'both')
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    plt.title( str(t).zfill(2) + "00 UTC " + str(d) + ' ' + month_str +' 2018', color = 'b', fontsize = 8, loc = 'right')
    plt.title('600 hPa: Dewpoint, Winds', fontsize = 8, loc = 'left')

    cbar.set_label('Dew point  (\N{DEGREE SIGN}C)',fontsize=8)
    cbar.ax.tick_params(labelsize = 8)
    path = "C:\\MISO-BOB\\Data\\TPW\\plots\\"
    
    ir_map = ax.pcolormesh(lons2, lats2, data2, vmin = 200, vmax = 245, cmap = cmap_grey, transform=ccrs.PlateCarree())
    cbar = plt.colorbar(ir_map, ax=ax, orientation='horizontal',shrink=0.35, extend = 'min')
    cbar.set_label('IR temperature (K)',fontsize=8)
    cbar.ax.tick_params(labelsize = 8)
    gl.xlabel_style = {'size': 6, 'color': 'k'}
    gl.ylabel_style = {'size': 6, 'color': 'k'}
    
    #dvs = ndimage.gaussian_filter(dvu, sigma=1.0, order=0)
    #ax.contour(lons1, lats1, -dvs, levels = [-0.3, -0.1], colors ='red', linestyle = 'dashed', linewidths = 0.5, transform=ccrs.PlateCarree())
    
    #dvs = ndimage.gaussian_filter(dvb, sigma=1.0, order=0)
    #ax.contour(lons1, lats1, -dvs, levels = [0.1,0.3], colors ='blue', linestyle = 'dashdotdotted', linewidths = 0.5, transform=ccrs.PlateCarree())

    wf = 7
    Q = ax.quiver(lons1[::wf],lats1[::wf],u[::wf,::wf],v[::wf,::wf],color='black', units='width', scale = 200)
    ax.quiverkey(Q, 0.62, 0.46, 10, r'10 ms$^{-1}$', labelpos='E', coordinates='figure', fontproperties = {'size':10})
    
    ############ Plot flight path ##########
    plt.scatter(df["LON"].values,df["LAT"].values,0.2,'r','.')
    
    ############ Plot dropsonde location ##########
    for loc in locs:
        plt.scatter(loc[0],loc[1],20,'black')
        
    plt.scatter(79.86,6.92,100,'blue','*')
    
    plt.savefig(path + "tpw800-600_" + str(d).zfill(2) + str(t).zfill(2) + month_str + '2018.png',
                dpi = 300, bbox_inches='tight',  pad_inches=0.3)
    plt.close() 
    plt.clf()
    print("Figure saved in " + path)

 

