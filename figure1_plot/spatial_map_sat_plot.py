flag_run = 1
# exec(open('xxx').read())
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

sys.path.append('/home/yliang/lib/python_functions/data_process/')
import data_process_f

# ================================================================
# define functions 
# ================================================================
def perform_ttest_here(exp1_var,exp2_var,ny,nx,sig_level):
    ttest_map = np.zeros((ny,nx))*np.nan
    pvalue_map = np.zeros((ny,nx))*np.nan
    for JJ in range(ny):
        for II in range(nx):
            [xxx, pvalue] = stats.ttest_ind(exp1_var[:,JJ,II],exp2_var[:,JJ,II])
            if pvalue < sig_level:
               ttest_map[JJ,II] = 1.
            pvalue_map[JJ,II] = pvalue

    return ttest_map, pvalue_map

def plot_box(ax,lon1,lon2,lat1,lat2, color_txt):

    ax.plot(np.linspace(lon1,lon1,100), np.linspace(lat1, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=0.6)
    ax.plot(np.linspace(lon2,lon2,100), np.linspace(lat1, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=0.6)
    ax.plot(np.linspace(lon1,lon2,100), np.linspace(lat1, lat1, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=0.6)
    ax.plot(np.linspace(lon1,lon2,100), np.linspace(lat2, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=0.6)

if flag_run == 1:
# ================================================================
# read simulations
# ================================================================
# read sat
   varname = 'TREFHT'
   year1 = 1955
   year2 = 2005
   year_N = int(year2 - year1 + 1)
# read grid basics
   dirname = '/data1/yliang/cesm1_le/data_storage_center/'
   filename = varname + '_annual_mean_temp_output.nc'
   f = Dataset(dirname + filename, 'r')
   lon = f.variables['lon'][:].data
   lat = f.variables['lat'][:].data
   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(0,360,60,89,lon,lat)
   y1 = y1 + 1
   print("region selected lon:" + str(lon[x1]) + '-' + str(lon[x2]) + ' lat:' + str(lat[y1]) + '-' + str(lat[y2]))
   lon_sel = lon.copy()
   lat_sel = lat[y1:].copy()
   ny_sel = len(lat_sel)
   nx_sel = len(lon_sel)
   sat_ctr = f.variables['var_ctr'][:,:,y1:,:].data
   sat_xc2 = f.variables['var_xc2'][:,:,y1:,:].data
   sat_xod = f.variables['var_xod'][:,:,y1:,:].data
   sat_xco = f.variables['var_xco'][:,:,y1:,:].data
   f.close()

   ens_num = sat_ctr.shape[0]
   ny = len(lat)
   nx = len(lon)

# read sic
   varname = 'ICEFRAC'
   year1 = 1955
   year2 = 2005
   year_N = int(year2 - year1 + 1)
# read grid basics
   dirname = '/data1/yliang/cesm1_le/data_storage_center/'
   filename = varname + '_annual_mean_temp_output.nc'
   f = Dataset(dirname + filename, 'r')
   sic_ctr = f.variables['var_ctr'][:,:,y1:,:].data
   sic_xc2 = f.variables['var_xc2'][:,:,y1:,:].data
   sic_xod = f.variables['var_xod'][:,:,y1:,:].data
   sic_xco = f.variables['var_xco'][:,:,y1:,:].data
   f.close()

#   sic_ctr[sic_ctr!=0.] = 1.
#   sic_xc2[sic_xc2!=0.] = 1.
#   sic_xod[sic_xod!=0.] = 1.
#   sic_xco[sic_xco!=0.] = 1.

# ================================================================
# calculate change and trend defined in polvani et al. (2020) ncc
# ================================================================
   sat_ctr_first = np.zeros((ens_num,ny_sel,nx_sel))
   sat_xc2_first = np.zeros((ens_num,ny_sel,nx_sel))
   sat_xod_first = np.zeros((ens_num,ny_sel,nx_sel))
   sat_xco_first = np.zeros((ens_num,ny_sel,nx_sel))
   sat_ctr_last = np.zeros((ens_num,ny_sel,nx_sel))
   sat_xc2_last = np.zeros((ens_num,ny_sel,nx_sel))
   sat_xod_last = np.zeros((ens_num,ny_sel,nx_sel))
   sat_xco_last = np.zeros((ens_num,ny_sel,nx_sel))
   for NENS in range(ens_num):
       print(NENS)
       for JJ in range(ny_sel):
           for II in range(nx_sel):
               sat_ctr_first[NENS,JJ,II] = np.nanmean(sat_ctr[NENS,0:10,JJ,II])
               sat_ctr_last[NENS,JJ,II] = np.nanmean(sat_ctr[NENS,-10:,JJ,II])
               sat_xc2_first[NENS,JJ,II] = np.nanmean(sat_xc2[NENS,0:10,JJ,II])
               sat_xc2_last[NENS,JJ,II] = np.nanmean(sat_xc2[NENS,-10:,JJ,II])
               sat_xod_first[NENS,JJ,II] = np.nanmean(sat_xod[NENS,0:10,JJ,II])
               sat_xod_last[NENS,JJ,II] = np.nanmean(sat_xod[NENS,-10:,JJ,II])
               sat_xco_first[NENS,JJ,II] = np.nanmean(sat_xco[NENS,0:10,JJ,II])
               sat_xco_last[NENS,JJ,II] = np.nanmean(sat_xco[NENS,-10:,JJ,II])


   sic_ctr_first = np.zeros((ens_num,ny_sel,nx_sel))
   sic_xc2_first = np.zeros((ens_num,ny_sel,nx_sel))
   sic_xod_first = np.zeros((ens_num,ny_sel,nx_sel))
   sic_xco_first = np.zeros((ens_num,ny_sel,nx_sel))
   sic_ctr_last = np.zeros((ens_num,ny_sel,nx_sel))
   sic_xc2_last = np.zeros((ens_num,ny_sel,nx_sel))
   sic_xod_last = np.zeros((ens_num,ny_sel,nx_sel))
   sic_xco_last = np.zeros((ens_num,ny_sel,nx_sel))
   for NENS in range(ens_num):
       print(NENS)
       for JJ in range(ny_sel):
           for II in range(nx_sel):
               sic_ctr_first[NENS,JJ,II] = np.nanmean(sic_ctr[NENS,0:10,JJ,II])
               sic_ctr_last[NENS,JJ,II] = np.nanmean(sic_ctr[NENS,-10:,JJ,II])
               sic_xc2_first[NENS,JJ,II] = np.nanmean(sic_xc2[NENS,0:10,JJ,II])
               sic_xc2_last[NENS,JJ,II] = np.nanmean(sic_xc2[NENS,-10:,JJ,II])
               sic_xod_first[NENS,JJ,II] = np.nanmean(sic_xod[NENS,0:10,JJ,II])
               sic_xod_last[NENS,JJ,II] = np.nanmean(sic_xod[NENS,-10:,JJ,II])
               sic_xco_first[NENS,JJ,II] = np.nanmean(sic_xco[NENS,0:10,JJ,II])
               sic_xco_last[NENS,JJ,II] = np.nanmean(sic_xco[NENS,-10:,JJ,II])

# ================================================================
# perform statistical significance
# ================================================================
   sig_level = 0.05

   test_map1 = sat_ctr_first.copy()
   test_map2 = sat_ctr_last.copy()
   [ttest_map_ctr, pvalue_map] = perform_ttest_here(test_map1,test_map2,ny_sel,nx_sel,sig_level)

   test_map1 = sat_xc2_last.copy() - sat_xc2_first.copy()
   test_map2 = sat_ctr_last.copy() - sat_ctr_first.copy()
   [ttest_map_xc2, pvalue_map] = perform_ttest_here(test_map1,test_map2,ny_sel,nx_sel,sig_level)

   test_map1 = sat_xod_last.copy() - sat_xod_first.copy()
   test_map2 = sat_ctr_last.copy() - sat_ctr_first.copy()
   [ttest_map_xod, pvalue_map] = perform_ttest_here(test_map1,test_map2,ny_sel,nx_sel,sig_level)

   test_map1 = sat_xco_last.copy() - sat_xco_first.copy()
   test_map2 = sat_ctr_last.copy() - sat_ctr_first.copy()
   [ttest_map_xco, pvalue_map] = perform_ttest_here(test_map1,test_map2,ny_sel,nx_sel,sig_level)

   test_map1 = sic_ctr_first.copy()
   test_map2 = sic_ctr_last.copy()
   [ttest_sic_ctr, pvalue_map] = perform_ttest_here(test_map1,test_map2,ny_sel,nx_sel,sig_level)

   test_map1 = sic_xc2_last.copy() - sic_xc2_first.copy()
   test_map2 = sic_ctr_last.copy() - sic_ctr_first.copy()
   [ttest_sic_xc2, pvalue_map] = perform_ttest_here(test_map1,test_map2,ny_sel,nx_sel,sig_level)

   test_map1 = sic_xod_last.copy() - sic_xod_first.copy()
   test_map2 = sic_ctr_last.copy() - sic_ctr_first.copy()
   [ttest_sic_xod, pvalue_map] = perform_ttest_here(test_map1,test_map2,ny_sel,nx_sel,sig_level)

   test_map1 = sic_xco_last.copy() - sic_xco_first.copy()
   test_map2 = sic_ctr_last.copy() - sic_ctr_first.copy()
   [ttest_sic_xco, pvalue_map] = perform_ttest_here(test_map1,test_map2,ny_sel,nx_sel,sig_level)

# ================================================================
# plot figures
# ================================================================
if True:
   
   plt.close('all')

   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)

   theta = np.linspace(0, 2*np.pi, 100)
   center, radius = [0.5, 0.5], 0.5
   verts = np.vstack([np.sin(theta), np.cos(theta)]).T
   circle = mpath.Path(verts * radius + center)

# plot sat
   clevel = np.linspace(-0.6,0.6,21)

   ax1 = fig.add_axes([0.01, 0.65, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = (np.nanmean(sat_ctr_last, axis=0)-np.nanmean(sat_ctr_first, axis=0))/51.*10
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat_sel,lon_sel)
   im1 = ax1.contourf(lon_x, lat_sel, x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = ttest_map_ctr.copy()
   for II in range(len(lat_sel[::3])):
       ax1.plot(lon_sel[::4],lat_sel[II*3]*sig_map[II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
#   plot_box(ax1,0,360,41,42, 'w')
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(a) Hist SAT', fontsize=12)  

   ax1 = fig.add_axes([0.26, 0.65, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = (np.nanmean(sat_ctr_last, axis=0)-np.nanmean(sat_ctr_first, axis=0) - (np.nanmean(sat_xc2_last, axis=0)-np.nanmean(sat_xc2_first, axis=0)))/51.*10
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat_sel,lon_sel)
   im1 = ax1.contourf(lon_x, lat_sel, x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = ttest_map_xc2.copy()
   for II in range(len(lat_sel[::3])):
       ax1.plot(lon_sel[::4],lat_sel[II*3]*sig_map[II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(b) Hist-xCO2 SAT', fontsize=12)

   ax1 = fig.add_axes([0.51, 0.65, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = (np.nanmean(sat_ctr_last, axis=0)-np.nanmean(sat_ctr_first, axis=0) - (np.nanmean(sat_xod_last, axis=0)-np.nanmean(sat_xod_first, axis=0)))/51.*10
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat_sel,lon_sel)
   im1 = ax1.contourf(lon_x, lat_sel, x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = ttest_map_xod.copy()
   for II in range(len(lat_sel[::3])):
       ax1.plot(lon_sel[::4],lat_sel[II*3]*sig_map[II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(c) Hist-xODS SAT', fontsize=12)

   ax1 = fig.add_axes([0.76, 0.65, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = (np.nanmean(sat_ctr_last, axis=0)-np.nanmean(sat_ctr_first, axis=0) - (np.nanmean(sat_xco_last, axis=0)-np.nanmean(sat_xco_first, axis=0)))/51.*10
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat_sel,lon_sel)
   im1 = ax1.contourf(lon_x, lat_sel, x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   sig_map = ttest_map_xco.copy()
   for II in range(len(lat_sel[::3])):
       ax1.plot(lon_sel[::4],lat_sel[II*3]*sig_map[II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(d) Hist-xCO SAT', fontsize=12)

   cbaxes = fig.add_axes([0.14, 0.62, 0.7, 0.01])
   cbar = plt.colorbar(im1, cax=cbaxes, orientation='horizontal', ticks=np.linspace(clevel[0],clevel[-1],13))
   cbar.set_label('K per decade', rotation=0)

# plot sic
   clevel = np.linspace(-100,100,21)   

   factor = 1./1000000.
   area = data_process_f.area_calculate_nonuniform(lon_sel,lat_sel)*factor

   ax1 = fig.add_axes([0.01, 0.3, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = (np.nanmean(sic_ctr_last, axis=0)-np.nanmean(sic_ctr_first, axis=0))*area/51.*10
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat_sel,lon_sel)
   im1 = ax1.contourf(lon_x, lat_sel, x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu')
   sig_map = ttest_sic_ctr.copy()
   for II in range(len(lat_sel[::3])):
       ax1.plot(lon_sel[::4],lat_sel[II*3]*sig_map[II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
#   plot_box(ax1,0,360,41,42, 'w')
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(e) Hist SIA', fontsize=12)

   ax1 = fig.add_axes([0.26, 0.3, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = (np.nanmean(sic_ctr_last, axis=0)-np.nanmean(sic_ctr_first, axis=0) - (np.nanmean(sic_xc2_last, axis=0)-np.nanmean(sic_xc2_first, axis=0)))*area/51.*10
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat_sel,lon_sel)
   im1 = ax1.contourf(lon_x, lat_sel, x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu')
   sig_map = ttest_sic_xc2.copy()
   for II in range(len(lat_sel[::3])):
       ax1.plot(lon_sel[::4],lat_sel[II*3]*sig_map[II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(f) Hist-xCO2 SIA', fontsize=12)

   ax1 = fig.add_axes([0.51, 0.3, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = (np.nanmean(sic_ctr_last, axis=0)-np.nanmean(sic_ctr_first, axis=0) - (np.nanmean(sic_xod_last, axis=0)-np.nanmean(sic_xod_first, axis=0)))*area/51.*10
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat_sel,lon_sel)
   im1 = ax1.contourf(lon_x, lat_sel, x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu')
   sig_map = ttest_sic_xod.copy()
   for II in range(len(lat_sel[::3])):
       ax1.plot(lon_sel[::4],lat_sel[II*3]*sig_map[II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(g) Hist-xODS SIA', fontsize=12)

   ax1 = fig.add_axes([0.76, 0.3, 0.22, 0.22], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = (np.nanmean(sic_ctr_last, axis=0)-np.nanmean(sic_ctr_first, axis=0) - (np.nanmean(sic_xco_last, axis=0)-np.nanmean(sic_xco_first, axis=0)))*area/51.*10
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat_sel,lon_sel)
   im1 = ax1.contourf(lon_x, lat_sel, x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu')
   sig_map = ttest_sic_xco.copy()
   for II in range(len(lat_sel[::3])):
       ax1.plot(lon_sel[::4],lat_sel[II*3]*sig_map[II*3,::4],'ko', transform=ccrs.PlateCarree(), markersize=0.3)
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('(h) Hist-xCO SIA', fontsize=12)

   cbaxes = fig.add_axes([0.14, 0.27, 0.7, 0.01])
   cbar = plt.colorbar(im1, cax=cbaxes, orientation='horizontal', ticks=np.linspace(clevel[0],clevel[-1],11))
   cbar.set_label('km$^{2}$ per decade', rotation=0)

   plt.savefig('trend_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()

   sys.exit()






