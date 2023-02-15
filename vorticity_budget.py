# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 15:37:04 2022

@author: jayes
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 19:36:10 2022

@author: jayes
"""

from netCDF4 import Dataset
import miso_bob_functions as mbf
import numpy as np
import time 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import metpy.calc as mpcalc
start_time = time.time()


def crop_data(lats,lons,data, minX,maxX,minY,maxY):
    
    lat_min = np.argmin(np.abs(np.array(lats-minY)))
    lat_max = np.argmin(np.abs(np.array(lats-maxY))) 
    lon_min = np.argmin(np.abs(np.array(lons-minX)))
    lon_max = np.argmin(np.abs(np.array(lons-maxX)))

    lats = lats[lat_min:lat_max]
    lons = lons[lon_min:lon_max]
    data = data[:,lat_min:lat_max,lon_min:lon_max]
    return data,lats,lons


    

minX, maxX, minY, maxY = 70,90,15,8

######################## Load Orography #################################

era_path = 'D:\\Projects\\MISO_BOB\\Data\\'
f_name = "elev.0.25-deg.nc"
ds = Dataset(era_path + f_name)
lats_elv  = ds['lat'][:]
lons_elv  = ds['lon'][:]
elv       = ds['data'][:]


elv_sec,lats_elv_sec,lons_elv_sec  = crop_data(lats_elv, lons_elv, elv , minX , maxX ,minY ,maxY)
elv_sec = elv_sec.squeeze()

elv_sec_m = np.mean(elv_sec, axis = 0)

ph = 1013.25* ((1 - elv_sec_m/44307.7)**5.25)

##################### select date and time ##############################

to = (1900,1,1,0,0,0)

##################### Load ERA-interim datset ##############################

era_path = 'D:\\Projects\\IMONDO\\Data\\'
f_name = "era_interim_112015.nc"
ds = Dataset(era_path + f_name)
lats = ds['latitude'][:]
lons = ds['longitude'][:]
t  = ds['time'][:]
p  = ds['level'][:]
pv = ds['pv'][:]
z = ds['z'][:]
T = ds['t'][:]
w = ds['w'][:]
vo = ds['vo'][:]
dv = ds['d'][:]
u = ds['u'][:]
v = ds['v'][:]

fig, ax = plt.subplots(4,1,figsize=(5,15), sharey = True, sharex = True, dpi = 200)
plt.rcParams['font.size'] = 18

cross_dates = [9,10,11,12]
block_dates = [15,16,17,18]

sel_dates = cross_dates
sel_dates = block_dates

fig = plt.figure(figsize=(20,20))
plt.rcParams['font.size'] = 12

for d,xi in zip(sel_dates, range(4)):
   
    pdesire = 900
    idesire = np.where(p == pdesire)
    idesire = idesire[0][0]
    
    t_sel1 = (2015,11,d,0,0,0)
    t_sel2 = (2015,11,d+1,0,0)

    it1 = np.where(t == mbf.get_hours(to, t_sel1))
    it2 = np.where(t == mbf.get_hours(to, t_sel2))

    it1 = it1[0][0]
    it2 = it2[0][0]
    
    vo_sel0 =  vo[it1:it1+4,idesire,:,:].squeeze()
    dv_sel0 =  dv[it1:it1+4,idesire,:,:].squeeze()
    
    u_sel0 =  u[it1:it1+4,idesire,:,:].squeeze()
    u_sel1 =  u[it2:it2+4,idesire,:,:].squeeze()
    
    v_sel0 =  v[it1:it1+4,idesire,:,:].squeeze()
    v_sel1 =  v[it2:it2+4,idesire,:,:].squeeze()
    
    vo_sel1 =  vo[it2:it2+4,idesire,:,:].squeeze()
    dv_sel1 =  dv[it2:it2+4,idesire,:,:].squeeze()
    
    z_sel0 =  z[it1:it1+4,idesire,:,:].squeeze()
    z_sel1 =  z[it2:it2+4,idesire,:,:].squeeze()


    et_sel0 = np.zeros(vo_sel0.shape)
    et_sel1 = np.zeros(vo_sel1.shape)

    # Calculate planetary vorticity
    for k in range(vo_sel0.shape[0]):
        for i in range(vo_sel0.shape[1]):
            for j in range(vo_sel0.shape[2]):
                et_sel0[k,i,j] =  2*7.2921*10**(-5)*np.sin(lats[i]* np.pi / 180.0) + vo_sel0[k,i,j] 
                et_sel1[k,i,j] =  2*7.2921*10**(-5)*np.sin(lats[i]* np.pi / 180.0) + vo_sel1[k,i,j]
                
    et_sel0 = np.mean(et_sel0, axis = 0)
    et_sel1 = np.mean(et_sel1, axis = 0)

    
    
    ######################### 1. Vorticity tendency ###############################
    ax = plt.subplot(4,4,4*xi+1, projection=ccrs.PlateCarree())
    datacrs = ccrs.PlateCarree()
    
    vo_sel0 = np.mean(vo_sel0, axis = 0)
    vo_sel1 = np.mean(vo_sel1, axis = 0)
    
    # plot vorticity tendnency
    et_levels = [-1,-0.5,-0.25, -0.1, -0.05, 0, 0.05, 0.1, 0.25, 0.5, 1]
    et_map = ax.contourf(lons, lats, (10**4)*(vo_sel1-vo_sel0), levels = et_levels, cmap = "bwr_r", extend = 'both')
    
    # plot geopotential
    z_sel0 = z_sel0/9.81
    z_sel0 = np.mean(z_sel0, axis = 0)
    ax.contour(lons, lats, z_sel0, levels = range(900,1500,10), colors = 'Green')

    # plot orography
    ax.contour(lons_elv, lats_elv, elv.squeeze(), levels = [300,500,750], colors = 'Grey')

    # plot winds
    qf = 1
    u_sel0 = np.mean(u_sel0, axis = 0)
    v_sel0 = np.mean(v_sel0, axis = 0)
    u_sel1 = np.mean(u_sel1, axis = 0)
    v_sel1 = np.mean(v_sel1, axis = 0)
    Q  = ax.quiver(lons[::qf], lats[::qf],u_sel0[::qf, ::qf], v_sel0[::qf,::qf],color='black',units='inches',scale = 100, lw= 2)
    
    ax.set_xticks(np.arange(0,100,5))
    ax.set_yticks(np.arange(0,50,5))
    ax.set_xlim([70,90])         
    ax.set_ylim([0,20])
    ax.coastlines('10m',color='k')
    #ax.set_xlabel('Longitude ($^{o}$E)')
    ax.set_ylabel('Latitude ($^{o}$N)')
    
    
    ######################### 2. Vorticity Stretching ###############################
    ax = plt.subplot(4,4,4*xi+2, projection=ccrs.PlateCarree())
    datacrs = ccrs.PlateCarree()

    # plot vortext stretching term
    vor_stretch = -dv_sel0*et_sel0
    vor_stretch = np.mean(vor_stretch,axis = 0)
    et_levels = [-3,-2,-1,-0.5,-0.25, -0.1, -0.05, 0, 0.05, 0.1, 0.25, 0.5, 1,2,3]
    et_map = ax.contourf(lons, lats, (10**4)*(86400*vor_stretch), levels = et_levels, cmap = "bwr_r", extend = 'both')

    # plot orography
    ax.contour(lons_elv, lats_elv, elv.squeeze(), levels = [300,500,750], colors = 'Grey')
    
    ax.contour(lons, lats, z_sel0, levels = range(900,1500,10), colors = 'Green')

    # plot winds
    qf = 1
    Q  = ax.quiver(lons[::qf], lats[::qf],u_sel0[::qf, ::qf], v_sel0[::qf,::qf],color='black',units='inches',scale = 100)
    
    ax.set_xticks(np.arange(0,100,5))
    ax.set_yticks(np.arange(0,50,5))
    ax.set_xlim([70,90])         
    ax.set_ylim([0,20])
    ax.coastlines('10m',color='k')
    #ax.set_xlabel('Longitude ($^{o}$E)')
    ax.set_ylabel('Latitude ($^{o}$N)')
    
    
    ######################### 3. Vorticity Advection: Zonal ###############################
    # plot vorticity advection
    ax = plt.subplot(4,4,4*xi+3, projection=ccrs.PlateCarree())
    datacrs = ccrs.PlateCarree()
    
    vo_dy, vo_dx = mpcalc.gradient(et_sel0, deltas = (1,1))
    vo_adv =   -( u_sel0*vo_dx)/(110*1000)

    et_map = ax.contourf(lons, lats, (10**4)*(86400*vo_adv), levels = et_levels, cmap = "bwr_r", extend = 'both')

    # plot orography
    ax.contour(lons_elv, lats_elv, elv.squeeze(), levels = [300,500,750], colors = 'Grey')

    ax.contour(lons, lats, z_sel0, levels = range(900,1500,10), colors = 'Green')

    # plot winds
    qf = 1
    Q  = ax.quiver(lons[::qf], lats[::qf],u_sel0[::qf, ::qf], v_sel0[::qf,::qf],color='black',units='inches',scale = 100)
    
    ax.set_xticks(np.arange(0,100,5))
    ax.set_yticks(np.arange(0,50,5))
    ax.set_xlim([70,90])         
    ax.set_ylim([0,20])
    ax.coastlines('10m',color='k')
    #ax.set_xlabel('Longitude ($^{o}$E)')
    ax.set_ylabel('Latitude ($^{o}$N)')
    
    
    ######################### 3. Vorticity Advection: Zonal ###############################

    # plot vorticity advection
    ax = plt.subplot(4,4,4*xi+4, projection=ccrs.PlateCarree())
    datacrs = ccrs.PlateCarree()
    
    wo_sel0 =   w[it1:it1+4,idesire,:,:].squeeze()
    vl2_sel0 =  vo[it1:it1+4,idesire,:,:].squeeze()
    vl1_sel0 =  vo[it1:it1+4,idesire-1,:,:].squeeze()
    vl3_sel0 =  vo[it1:it1+4,idesire+1,:,:].squeeze()
    
    dvodp = (vl3_sel0 - vl1_sel0)/(200*1000) 
    wdvdp = wo_sel0*dvodp
    wdvdp = np.mean(wdvdp, axis = 0)
    
    vo_adv =   -(- v_sel0*vo_dy)/(110*1000)
    
    et_map = ax.contourf(lons, lats, (10**4)*(86400*vo_adv), levels = et_levels, cmap = "bwr_r", extend = 'both')

    # plot orography
    ax.contour(lons_elv, lats_elv, elv.squeeze(), levels = [300,500,750], colors = 'Grey')

    ax.contour(lons, lats, z_sel0, levels = range(900,1500,10), colors = 'Green')

    # plot winds
    qf = 1
    Q  = ax.quiver(lons[::qf], lats[::qf], u_sel0[::qf, ::qf], v_sel0[::qf,::qf], color = 'black',units='inches', scale = 100)
    
    ax.set_xticks(np.arange(0,100,5))
    ax.set_yticks(np.arange(0,50,5))
    ax.set_xlim([70,90])         
    ax.set_ylim([0,20])
    ax.coastlines('10m', color='k')
    #ax.set_xlabel('Longitude ($^{o}$E)')
    ax.set_ylabel('Latitude ($^{o}$N)')

    
    
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.9, 0.25, 0.02, 0.5])
plt.colorbar(et_map, cax = cbar_ax)
cbar_ax.set_ylabel('X 10$^{-4}$ s$^{-1}$', fontsize = 20)

plt.savefig(era_path + "abs_vor_" + str(d) + "_" + str(pdesire) + ".png", dpi = 200, bbox_inches='tight',  pad_inches=0.3)
plt.close() 
plt.clf()
print("Figure saved in " + era_path)   
    












