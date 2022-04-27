flag_run = 1
# ================================================================
# Yu-Chiao @ Fort Lee, NJ Oct 20, 2020
# examination on  surface air temperature in cesm1-cam5 simulations
# ================================================================

# ================================================================
# import functions
# ================================================================
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from math import isnan, radians
#from mpl_toolkits.basemap import Basemap
from IPython import get_ipython
import sys, os
#import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.img_tiles as cimgt
from cartopy.io.img_tiles import StamenTerrain
from scipy import stats
import matplotlib.path as mpath
from matplotlib import rcParams

sys.path.append('/home/yliang/lib/python_functions/data_process/')
import data_process_f

# ================================================================
# define functions 
# ================================================================
def perform_ttest_1d_here(exp1_var,exp2_var,sig_level):
    ttest_map = np.nan
    [xxx, pvalue] = stats.ttest_ind(exp1_var[:],exp2_var[:])
    if pvalue < sig_level:
       ttest_map = 1.
    else:
       ttest_map = 0.3
    pvalue_map = pvalue

    return ttest_map, pvalue_map

def perform_ttest_2d_here(exp1_var,exp2_var,ny,nx,sig_level):
    ttest_map = np.zeros((ny,nx))*np.nan
    pvalue_map = np.zeros((ny,nx))*np.nan
    for JJ in range(ny):
        for II in range(nx):
            [xxx, pvalue] = stats.ttest_ind(exp1_var[:,JJ,II],exp2_var[:,JJ,II])
            if pvalue < sig_level:
               ttest_map[JJ,II] = 1.
            pvalue_map[JJ,II] = pvalue

    return ttest_map, pvalue_map

if flag_run == 1:
# ================================================================
# read files
# ================================================================
   st_index = 0

#   dirname = '/data1/yliang/cesm1_le/data_storage_center/radiative_kernels/'
   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
   filename = 'feedback_parameters_maps_ctr_temp.nc'
#   filename = 'radiation_changes_maps_ctr_temp.nc'
   f = Dataset(dirname + filename, 'r')
   lat = f.variables['lat'][:].data.copy()
   lon = f.variables['lon'][:].data.copy()
   t_feedback_ctr = f.variables['t_feedback_all'][st_index:].data
   t_feedback_map_ctr = f.variables['t_feedback_map_all'][st_index:,:,:].data
   planck_feedback_ctr = f.variables['planck_feedback_all'][st_index:].data
   planck_feedback_map_ctr = f.variables['planck_feedback_map_all'][st_index:,:,:].data
   lapserate_feedback_ctr = f.variables['lapserate_feedback_all'][st_index:].data
   lapserate_feedback_map_ctr = f.variables['lapserate_feedback_map_all'][st_index:,:,:].data
   alb_feedback_ctr = f.variables['alb_feedback_all'][st_index:].data
   alb_feedback_map_ctr = f.variables['alb_feedback_map_all'][st_index:,:,:].data
   q_feedback_ctr = f.variables['q_feedback_all'][st_index:].data
   q_feedback_map_ctr = f.variables['q_feedback_map_all'][st_index:,:,:].data
   lw_cloud_feedback_ctr = f.variables['lw_cloud_feedback_all'][st_index:].data
   lw_cloud_feedback_map_ctr = f.variables['lw_cloud_feedback_map_all'][st_index:,:,:].data
   sw_cloud_feedback_ctr = f.variables['sw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_map_ctr = f.variables['sw_cloud_feedback_map_all'][st_index:,:,:].data
   dts_arctic_ctr = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_ctr = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_ctr = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close() 
   dts_ctr = dts_arctic_ctr.mean()  

   ens_num = len(t_feedback_ctr)

   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
   filename = 'feedback_parameters_maps_xoo_temp.nc'
   f = Dataset(dirname + filename, 'r')
   t_feedback_xco = f.variables['t_feedback_all'][st_index:].data
   t_feedback_map_xco = f.variables['t_feedback_map_all'][st_index:,:,:].data
   planck_feedback_xco = f.variables['planck_feedback_all'][st_index:].data
   planck_feedback_map_xco = f.variables['planck_feedback_map_all'][st_index:,:,:].data
   lapserate_feedback_xco = f.variables['lapserate_feedback_all'][st_index:].data
   lapserate_feedback_map_xco = f.variables['lapserate_feedback_map_all'][st_index:,:,:].data
   alb_feedback_xco = f.variables['alb_feedback_all'][st_index:].data
   alb_feedback_map_xco = f.variables['alb_feedback_map_all'][st_index:,:,:].data
   q_feedback_xco = f.variables['q_feedback_all'][st_index:].data
   q_feedback_map_xco = f.variables['q_feedback_map_all'][st_index:,:,:].data
   lw_cloud_feedback_xco = f.variables['lw_cloud_feedback_all'][st_index:].data
   lw_cloud_feedback_map_xco = f.variables['lw_cloud_feedback_map_all'][st_index:,:,:].data
   sw_cloud_feedback_xco = f.variables['sw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_map_xco = f.variables['sw_cloud_feedback_map_all'][st_index:,:,:].data
   dts_arctic_xco = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_xco = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_xco = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close()
   dts_xco = dts_arctic_xco.mean()

# read xod case
   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel//'
   filename = 'feedback_parameters_maps_xod_temp.nc'
   f = Dataset(dirname + filename, 'r')
   t_feedback_xod = f.variables['t_feedback_all'][st_index:].data
   t_feedback_map_xod = f.variables['t_feedback_map_all'][st_index:,:,:].data
   planck_feedback_xod = f.variables['planck_feedback_all'][st_index:].data
   planck_feedback_map_xod = f.variables['planck_feedback_map_all'][st_index:,:,:].data
   lapserate_feedback_xod = f.variables['lapserate_feedback_all'][st_index:].data
   lapserate_feedback_map_xod = f.variables['lapserate_feedback_map_all'][st_index:,:,:].data
   alb_feedback_xod = f.variables['alb_feedback_all'][st_index:].data
   alb_feedback_map_xod = f.variables['alb_feedback_map_all'][st_index:,:,:].data
   q_feedback_xod = f.variables['q_feedback_all'][st_index:].data
   q_feedback_map_xod = f.variables['q_feedback_map_all'][st_index:,:,:].data
   lw_cloud_feedback_xod = f.variables['lw_cloud_feedback_all'][st_index:].data
   lw_cloud_feedback_map_xod = f.variables['lw_cloud_feedback_map_all'][st_index:,:,:].data
   sw_cloud_feedback_xod = f.variables['sw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_map_xod = f.variables['sw_cloud_feedback_map_all'][st_index:,:,:].data
   dts_arctic_xod = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_xod = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_xod = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close()
   dts_xod = dts_arctic_xod.mean()

# read xc2 case
   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
   filename = 'feedback_parameters_maps_xc2_temp.nc'
   f = Dataset(dirname + filename, 'r')
   t_feedback_xc2 = f.variables['t_feedback_all'][st_index:].data
   t_feedback_map_xc2 = f.variables['t_feedback_map_all'][st_index:,:,:].data
   planck_feedback_xc2 = f.variables['planck_feedback_all'][st_index:].data
   planck_feedback_map_xc2 = f.variables['planck_feedback_map_all'][st_index:,:,:].data
   lapserate_feedback_xc2 = f.variables['lapserate_feedback_all'][st_index:].data
   lapserate_feedback_map_xc2 = f.variables['lapserate_feedback_map_all'][st_index:,:,:].data
   alb_feedback_xc2 = f.variables['alb_feedback_all'][st_index:].data
   alb_feedback_map_xc2 = f.variables['alb_feedback_map_all'][st_index:,:,:].data
   q_feedback_xc2 = f.variables['q_feedback_all'][st_index:].data
   q_feedback_map_xc2 = f.variables['q_feedback_map_all'][st_index:,:,:].data
   lw_cloud_feedback_xc2 = f.variables['lw_cloud_feedback_all'][st_index:].data
   lw_cloud_feedback_map_xc2 = f.variables['lw_cloud_feedback_map_all'][st_index:,:,:].data
   sw_cloud_feedback_xc2 = f.variables['sw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_map_xc2 = f.variables['sw_cloud_feedback_map_all'][st_index:,:,:].data
   dts_arctic_xc2 = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_xc2 = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_xc2 = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close()
   dts_xc2 = dts_arctic_xc2.mean()

# ================================================================
# calculate significance
# ================================================================
   ny = len(lat)
   nx = len(lon)
   sig_level = 0.05
#   [sig_planck_feedback_ctr_xod, pvalue] = perform_ttest_2d_here(planck_feedback_map_ctr,planck_feedback_map_xod,ny,nx,sig_level)
   [sig_lapserate_feedback_ctr_xod, pvalue] = perform_ttest_2d_here(lapserate_feedback_map_ctr,lapserate_feedback_map_xod,ny,nx,sig_level)
   [sig_wv_feedback_ctr_xod, pvalue] = perform_ttest_2d_here(q_feedback_map_ctr,q_feedback_map_xod,ny,nx,sig_level)
   [sig_al_feedback_ctr_xod, pvalue] = perform_ttest_2d_here(alb_feedback_map_ctr,alb_feedback_map_xod,ny,nx,sig_level)
   [sig_cl_feedback_ctr_xod, pvalue] = perform_ttest_2d_here(lw_cloud_feedback_map_ctr+sw_cloud_feedback_map_ctr,lw_cloud_feedback_map_xod+sw_cloud_feedback_map_xod,ny,nx,sig_level)

#   [sig_planck_feedback_ctr_xc2, pvalue] = perform_ttest_2d_here(planck_feedback_map_ctr,planck_feedback_map_xc2,ny,nx,sig_level)
   [sig_lapserate_feedback_ctr_xc2, pvalue] = perform_ttest_2d_here(lapserate_feedback_map_ctr,lapserate_feedback_map_xc2,ny,nx,sig_level)
   [sig_wv_feedback_ctr_xc2, pvalue] = perform_ttest_2d_here(q_feedback_map_ctr,q_feedback_map_xc2,ny,nx,sig_level)
   [sig_al_feedback_ctr_xc2, pvalue] = perform_ttest_2d_here(alb_feedback_map_ctr,alb_feedback_map_xc2,ny,nx,sig_level)
   [sig_cl_feedback_ctr_xc2, pvalue] = perform_ttest_2d_here(lw_cloud_feedback_map_ctr+sw_cloud_feedback_map_ctr,lw_cloud_feedback_map_xc2+sw_cloud_feedback_map_xc2,ny,nx,sig_level)

# ================================================================
# plot figures
# ================================================================
if True:

   plt.close('all')

# plot spatial maps
   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)

   theta = np.linspace(0, 2*np.pi, 100)
   center, radius = [0.5, 0.5], 0.5
   verts = np.vstack([np.sin(theta), np.cos(theta)]).T
   circle = mpath.Path(verts * radius + center)

   lat_sel = 160
#   lat_sel = 0

   rcParams['hatch.linewidth'] = 0.4
# plot sat
   clevel = np.linspace(-10,10,21)

# hist-xc2
   ax1 = fig.add_axes([0.01, 0.66, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
#   map_2d = (np.nanmean(lapserate_feedback_map_ctr[:,lat_sel:,:], axis=0) - np.nanmean(lapserate_feedback_map_xc2[:,lat_sel:,:], axis=0))/(dts_ctr-dts_xc2)
   map_2d_tmp = (lapserate_feedback_map_ctr[:,lat_sel:,:]-lapserate_feedback_map_xc2[:,lat_sel:,:])
   for II in range(10):
       map_2d_tmp[II,:,:] = map_2d_tmp[II,:,:]/(dts_arctic_ctr-dts_arctic_xc2)[II]
   map_2d = np.nanmean(map_2d_tmp, axis=0).copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = sig_lapserate_feedback_ctr_xc2[lat_sel:,:].copy()
#   sig_map[np.isnan(sig_map)] = 1
#   sig_map[sig_lapserate_feedback_ctr_xc2[lat_sel:,:]==1] = np.nan
   ax1.contourf(lon,lat[lat_sel:],sig_map,colors = 'none', hatches=['//'], transform=ccrs.PlateCarree())
#   for II in range(len(lat[lat_sel::3])):
#       ax1.plot(lon[::4],lat[lat_sel+II*3]*sig_map[lat_sel+II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(a) Lapse-rate CO2', fontsize=10)

   ax1 = fig.add_axes([0.26, 0.66, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
#   map_2d = (np.nanmean(alb_feedback_map_ctr[:,lat_sel:,:], axis=0) - np.nanmean(alb_feedback_map_xc2[:,lat_sel:,:], axis=0))/(dts_ctr-dts_xc2)
   map_2d_tmp = (alb_feedback_map_ctr[:,lat_sel:,:]-alb_feedback_map_xc2[:,lat_sel:,:]).copy()
   for II in range(10):
       map_2d_tmp[II,:,:] = map_2d_tmp[II,:,:]/(dts_arctic_ctr-dts_arctic_xc2)[II]
   map_2d = np.nanmean(map_2d_tmp, axis=0).copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = sig_al_feedback_ctr_xc2[lat_sel:,:].copy()
#   sig_map[np.isnan(sig_map)] = 1
#   sig_map[sig_al_feedback_ctr_xc2[lat_sel:,:]==1] = np.nan
   ax1.contourf(lon,lat[lat_sel:],sig_map,colors = 'none', hatches=['//'], transform=ccrs.PlateCarree())
#   for II in range(len(lat[lat_sel::3])):
#       ax1.plot(lon[::4],lat[lat_sel+II*3]*sig_map[lat_sel+II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(b) Albedo CO2', fontsize=10)

   ax1 = fig.add_axes([0.51, 0.66, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
#   map_2d = (np.nanmean(q_feedback_map_ctr[:,lat_sel:,:], axis=0) - np.nanmean(q_feedback_map_xc2[:,lat_sel:,:], axis=0))/(dts_ctr-dts_xc2)
   map_2d_tmp = (q_feedback_map_ctr[:,lat_sel:,:]-q_feedback_map_xc2[:,lat_sel:,:])
   for II in range(10):
       map_2d_tmp[II,:,:] = map_2d_tmp[II,:,:]/(dts_arctic_ctr-dts_arctic_xc2)[II]
   map_2d = np.nanmean(map_2d_tmp, axis=0).copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = sig_wv_feedback_ctr_xc2[lat_sel:,:].copy()
#   sig_map[np.isnan(sig_map)] = 1
#   sig_map[sig_wv_feedback_ctr_xc2[lat_sel:,:]==1] = np.nan
   ax1.contourf(lon,lat[lat_sel:],sig_map,colors = 'none', hatches=['//'], transform=ccrs.PlateCarree())
#   for II in range(len(lat[lat_sel::3])):
#       ax1.plot(lon[::4],lat[lat_sel+II*3]*sig_map[lat_sel+II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(c) Water vapor CO2', fontsize=10)

   ax1 = fig.add_axes([0.76, 0.66, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
#   map_2d = (np.nanmean((lw_cloud_feedback_map_ctr+sw_cloud_feedback_map_ctr)[:,lat_sel:,:], axis=0) - np.nanmean((lw_cloud_feedback_map_xc2+sw_cloud_feedback_map_xc2)[:,lat_sel:,:], axis=0))/(dts_ctr-dts_xc2)
   map_2d_tmp = ((lw_cloud_feedback_map_ctr+sw_cloud_feedback_map_ctr)[:,lat_sel:,:]-(lw_cloud_feedback_map_xc2+sw_cloud_feedback_map_xc2)[:,lat_sel:,:])
   for II in range(10):
       map_2d_tmp[II,:,:] = map_2d_tmp[II,:,:]/(dts_arctic_ctr-dts_arctic_xc2)[II]
   map_2d = np.nanmean(map_2d_tmp, axis=0).copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = sig_cl_feedback_ctr_xc2[lat_sel:,:].copy()
#   sig_map[np.isnan(sig_map)] = 1
#   sig_map[sig_cl_feedback_ctr_xc2[lat_sel:,:]==1] = np.nan
   ax1.contourf(lon,lat[lat_sel:],sig_map,colors = 'none', hatches=['//'], transform=ccrs.PlateCarree())
#   for II in range(len(lat[lat_sel::3])):
#       ax1.plot(lon[::4],lat[lat_sel+II*3]*sig_map[lat_sel+II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(d) Cloud CO2', fontsize=10)

# hist-xod
   ax1 = fig.add_axes([0.01, 0.40, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
#   map_2d = (np.nanmean(lapserate_feedback_map_ctr[:,lat_sel:,:], axis=0) - np.nanmean(lapserate_feedback_map_xod[:,lat_sel:,:], axis=0))/(dts_ctr-dts_xod)
   map_2d_tmp = (lapserate_feedback_map_ctr[:,lat_sel:,:]-lapserate_feedback_map_xod[:,lat_sel:,:])
   for II in range(10):
       map_2d_tmp[II,:,:] = map_2d_tmp[II,:,:]/(dts_arctic_ctr-dts_arctic_xod)[II]
   map_2d = np.nanmean(map_2d_tmp, axis=0).copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = sig_lapserate_feedback_ctr_xod[lat_sel:,:].copy()
#   sig_map[np.isnan(sig_map)] = 1
#   sig_map[sig_lapserate_feedback_ctr_xod[lat_sel:,:]==1] = np.nan
   ax1.contourf(lon,lat[lat_sel:],sig_map,colors = 'none', hatches=['//'], transform=ccrs.PlateCarree())
#   for II in range(len(lat[lat_sel::3])):
#       ax1.plot(lon[::4],lat[lat_sel+II*3]*sig_map[lat_sel+II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(e) Lapse-rate ODS', fontsize=10)

   ax1 = fig.add_axes([0.26, 0.40, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
#   map_2d = (np.nanmean(alb_feedback_map_ctr[:,lat_sel:,:], axis=0) - np.nanmean(alb_feedback_map_xod[:,lat_sel:,:], axis=0))/(dts_ctr-dts_xod)
   map_2d_tmp = (alb_feedback_map_ctr[:,lat_sel:,:]-alb_feedback_map_xod[:,lat_sel:,:])
   for II in range(10):
       map_2d_tmp[II,:,:] = map_2d_tmp[II,:,:]/(dts_arctic_ctr-dts_arctic_xod)[II]
   map_2d = np.nanmean(map_2d_tmp, axis=0).copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = sig_al_feedback_ctr_xod[lat_sel:,:].copy()
#   sig_map[np.isnan(sig_map)] = 1
#   sig_map[sig_al_feedback_ctr_xod[lat_sel:,:]==1] = np.nan
   ax1.contourf(lon,lat[lat_sel:],sig_map,colors = 'none', hatches=['//'], transform=ccrs.PlateCarree())
#   for II in range(len(lat[lat_sel::3])):
#       ax1.plot(lon[::4],lat[lat_sel+II*3]*sig_map[lat_sel+II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(f) Albedo ODS', fontsize=10)

   ax1 = fig.add_axes([0.51, 0.40, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
#   map_2d = (np.nanmean(q_feedback_map_ctr[:,lat_sel:,:], axis=0) - np.nanmean(q_feedback_map_xod[:,lat_sel:,:], axis=0))/(dts_ctr-dts_xod)
   map_2d_tmp = (q_feedback_map_ctr[:,lat_sel:,:]-q_feedback_map_xod[:,lat_sel:,:])
   for II in range(10):
       map_2d_tmp[II,:,:] = map_2d_tmp[II,:,:]/(dts_arctic_ctr-dts_arctic_xod)[II]
   map_2d = np.nanmean(map_2d_tmp, axis=0).copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = sig_wv_feedback_ctr_xod[lat_sel:,:].copy()
#   sig_map[np.isnan(sig_map)] = 1
#   sig_map[sig_wv_feedback_ctr_xod[lat_sel:,:]==1] = np.nan
   ax1.contourf(lon,lat[lat_sel:],sig_map,colors = 'none', hatches=['//'], transform=ccrs.PlateCarree())
#   for II in range(len(lat[lat_sel::3])):
#       ax1.plot(lon[::4],lat[lat_sel+II*3]*sig_map[lat_sel+II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(g) Water vapor ODS', fontsize=10)

   ax1 = fig.add_axes([0.76, 0.4, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
#   map_2d = (np.nanmean((lw_cloud_feedback_map_ctr+sw_cloud_feedback_map_ctr)[:,lat_sel:,:], axis=0) - np.nanmean((lw_cloud_feedback_map_xod+sw_cloud_feedback_map_xod)[:,lat_sel:,:], axis=0))/(dts_ctr-dts_xod)
   map_2d_tmp = ((lw_cloud_feedback_map_ctr+sw_cloud_feedback_map_ctr)[:,lat_sel:,:]-(lw_cloud_feedback_map_xod+sw_cloud_feedback_map_xod)[:,lat_sel:,:])
   for II in range(10):
       map_2d_tmp[II,:,:] = map_2d_tmp[II,:,:]/(dts_arctic_ctr-dts_arctic_xod)[II]
   map_2d = np.nanmean(map_2d_tmp, axis=0).copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = sig_cl_feedback_ctr_xod[lat_sel:,:].copy()
#   sig_map[np.isnan(sig_map)] = 1
#   sig_map[sig_cl_feedback_ctr_xod[lat_sel:,:]==1] = np.nan
   ax1.contourf(lon,lat[lat_sel:],sig_map,colors = 'none', hatches=['//'], transform=ccrs.PlateCarree())
#   for II in range(len(lat[lat_sel::3])):
#       ax1.plot(lon[::4],lat[lat_sel+II*3]*sig_map[lat_sel+II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(h) Cloud ODS', fontsize=10)

   cbaxes = fig.add_axes([0.14, 0.37, 0.7, 0.01])
   cbar = plt.colorbar(im1, cax=cbaxes, orientation='horizontal', ticks=np.linspace(clevel[0],clevel[-1],11))
   cbar.set_label('W/m$^{2}/$K', rotation=0)

   plt.savefig('trend_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()

   sys.exit()



