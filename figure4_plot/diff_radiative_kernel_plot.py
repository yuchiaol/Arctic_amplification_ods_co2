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

if flag_run == 1:
# ================================================================
# read files
# ================================================================
   st_index = 0

   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
   filename = 'feedback_parameters_maps_ctr_temp.nc'
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

   dirname = '/home/yliang/research/co2_ods/data_process/ohc/'
   filename = 'ctr_H_annual_temp_output.nc'
   f = Dataset(dirname + filename, 'r')
   h_ctr = f.variables['H'][:,:].data
   f.close()
   dh_ctr = np.nanmean(h_ctr[:,-10:], axis=1) - np.nanmean(h_ctr[:,0:10], axis=1)

   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
   filename = 'feedback_parameters_maps_xc2_temp.nc'
   f = Dataset(dirname + filename, 'r')
   t_feedback_xoo = f.variables['t_feedback_all'][st_index:].data
   t_feedback_map_xoo = f.variables['t_feedback_map_all'][st_index:,:,:].data
   planck_feedback_xoo = f.variables['planck_feedback_all'][st_index:].data
   planck_feedback_map_xoo = f.variables['planck_feedback_map_all'][st_index:,:,:].data
   lapserate_feedback_xoo = f.variables['lapserate_feedback_all'][st_index:].data
   lapserate_feedback_map_xoo = f.variables['lapserate_feedback_map_all'][st_index:,:,:].data
   alb_feedback_xoo = f.variables['alb_feedback_all'][st_index:].data
   alb_feedback_map_xoo = f.variables['alb_feedback_map_all'][st_index:,:,:].data
   q_feedback_xoo = f.variables['q_feedback_all'][st_index:].data
   q_feedback_map_xoo = f.variables['q_feedback_map_all'][st_index:,:,:].data
   lw_cloud_feedback_xoo = f.variables['lw_cloud_feedback_all'][st_index:].data
   lw_cloud_feedback_map_xoo = f.variables['lw_cloud_feedback_map_all'][st_index:,:,:].data
   sw_cloud_feedback_xoo = f.variables['sw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_map_xoo = f.variables['sw_cloud_feedback_map_all'][st_index:,:,:].data
   dts_arctic_xoo = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_xoo = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_xoo = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close()
   dts_xoo = dts_arctic_xoo.mean()

   dirname = '/home/yliang/research/co2_ods/data_process/ohc/'
   filename = 'xC2_H_annual_temp_output.nc'
   f = Dataset(dirname + filename, 'r')
   h_xoo = f.variables['H'][:,:].data
   f.close()
   dh_xoo = np.nanmean(h_xoo[:,-10:], axis=1) - np.nanmean(h_xoo[:,0:10], axis=1)

# read xod case
   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
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

   dirname = '/home/yliang/research/co2_ods/data_process/ohc/'
   filename = 'xOD_H_annual_temp_output.nc'
   f = Dataset(dirname + filename, 'r')
   h_xod = f.variables['H'][:,:].data
   f.close()
   dh_xod = np.nanmean(h_xod[:,-10:], axis=1) - np.nanmean(h_xod[:,0:10], axis=1)

# ================================================================
# compute feedback parameters
# ================================================================
   p_f_ctr = 1.26
   p_planck_ctr = planck_feedback_ctr.mean()
   p_lr_ctr = lapserate_feedback_ctr.mean()
   p_wv_ctr = q_feedback_ctr.mean()
   p_al_ctr = alb_feedback_ctr.mean()
   p_cllw_ctr = lw_cloud_feedback_ctr.mean()
   p_clsw_ctr = sw_cloud_feedback_ctr.mean()
   p_dh_ctr = dh_ctr.mean()
#   p_sum_ctr = (p_f_ctr+p_planck_ctr+p_lr_ctr+p_wv_ctr+p_al_ctr+p_cllw_ctr+p_clsw_ctr)
   p_sum_ctr = (p_f_ctr+planck_feedback_ctr+lapserate_feedback_ctr+q_feedback_ctr+alb_feedback_ctr+lw_cloud_feedback_ctr+sw_cloud_feedback_ctr).mean()
   p_df_ctr = (dh_ctr - (FLNT_arctic_ctr+FSNT_arctic_ctr)).mean()

# for xoo
   p_f_xoo = 1.01
   p_planck_xoo = planck_feedback_xoo.mean()
   p_lr_xoo = lapserate_feedback_xoo.mean()
   p_wv_xoo = q_feedback_xoo.mean()
   p_al_xoo = alb_feedback_xoo.mean()
   p_cllw_xoo = lw_cloud_feedback_xoo.mean()
   p_clsw_xoo = sw_cloud_feedback_xoo.mean()
   p_dh_xoo = dh_xoo.mean()
#   p_sum_xoo = (p_f_xoo+p_planck_xoo+p_lr_xoo+p_wv_xoo+p_al_xoo+p_cllw_xoo+p_clsw_xoo)
   p_sum_xoo = (p_f_xoo+planck_feedback_xoo+lapserate_feedback_xoo+q_feedback_xoo+alb_feedback_xoo+lw_cloud_feedback_xoo+sw_cloud_feedback_xoo).mean()
   p_df_xoo = (dh_xoo - (FLNT_arctic_xoo+FSNT_arctic_xoo)).mean()

# for xod
   p_f_xod = 1.11
   p_planck_xod = planck_feedback_xod.mean()
   p_lr_xod = lapserate_feedback_xod.mean()
   p_wv_xod = q_feedback_xod.mean()
   p_al_xod = alb_feedback_xod.mean()
   p_cllw_xod = lw_cloud_feedback_xod.mean()
   p_clsw_xod = sw_cloud_feedback_xod.mean()
   p_dh_xod = dh_xod.mean()
#   p_sum_xod = (p_f_xod+p_planck_xod+p_lr_xod+p_wv_xod+p_al_xod+p_cllw_xod+p_clsw_xod)
   p_sum_xod = (p_f_xod+planck_feedback_xod+lapserate_feedback_xod+q_feedback_xod+alb_feedback_xod+lw_cloud_feedback_xod+sw_cloud_feedback_xod).mean()
   p_df_xod = (dh_xod - (FLNT_arctic_xod+FSNT_arctic_xod)).mean()

# ================================================================
# calculate significance
# ================================================================
   sig_level = 0.05
   dts1 = (dts_ctr - dts_xod)
   dts2 = (dts_ctr - dts_xoo)
#   dts3 = (dts_xod - dts_xoo)
   [sig_planck_feedback_ctr_xod, pvalue] = perform_ttest_1d_here(planck_feedback_ctr,planck_feedback_xod,sig_level)
   [sig_planck_feedback_ctr_xoo, pvalue] = perform_ttest_1d_here(planck_feedback_ctr,planck_feedback_xoo,sig_level)
   [sig_planck_feedback_xod_xoo, pvalue] = perform_ttest_1d_here(planck_feedback_xod,planck_feedback_xoo,sig_level)
   [sig_lapserate_feedback_ctr_xod, pvalue] = perform_ttest_1d_here(lapserate_feedback_ctr,lapserate_feedback_xod,sig_level)
   [sig_lapserate_feedback_ctr_xoo, pvalue] = perform_ttest_1d_here(lapserate_feedback_ctr,lapserate_feedback_xoo,sig_level)
   [sig_lapserate_feedback_xod_xoo, pvalue] = perform_ttest_1d_here(lapserate_feedback_xod,lapserate_feedback_xoo,sig_level)
   [sig_wv_feedback_ctr_xod, pvalue] = perform_ttest_1d_here(q_feedback_ctr,q_feedback_xod,sig_level)
   [sig_wv_feedback_ctr_xoo, pvalue] = perform_ttest_1d_here(q_feedback_ctr,q_feedback_xoo,sig_level)
   [sig_wv_feedback_xod_xoo, pvalue] = perform_ttest_1d_here(q_feedback_xod,q_feedback_xoo,sig_level)
   [sig_al_feedback_ctr_xod, pvalue] = perform_ttest_1d_here(alb_feedback_ctr,alb_feedback_xod,sig_level)
   [sig_al_feedback_ctr_xoo, pvalue] = perform_ttest_1d_here(alb_feedback_ctr,alb_feedback_xoo,sig_level)
   [sig_al_feedback_xod_xoo, pvalue] = perform_ttest_1d_here(alb_feedback_xod,alb_feedback_xoo,sig_level)
   [sig_cllw_feedback_ctr_xod, pvalue] = perform_ttest_1d_here(lw_cloud_feedback_ctr,lw_cloud_feedback_xod,sig_level)
   [sig_cllw_feedback_ctr_xoo, pvalue] = perform_ttest_1d_here(lw_cloud_feedback_ctr,lw_cloud_feedback_xoo,sig_level)
   [sig_cllw_feedback_xod_xoo, pvalue] = perform_ttest_1d_here(lw_cloud_feedback_xod,lw_cloud_feedback_xoo,sig_level)
   [sig_clsw_feedback_ctr_xod, pvalue] = perform_ttest_1d_here(sw_cloud_feedback_ctr,sw_cloud_feedback_xod,sig_level)
   [sig_clsw_feedback_ctr_xoo, pvalue] = perform_ttest_1d_here(sw_cloud_feedback_ctr,sw_cloud_feedback_xoo,sig_level)
   [sig_clsw_feedback_xod_xoo, pvalue] = perform_ttest_1d_here(sw_cloud_feedback_xod,sw_cloud_feedback_xoo,sig_level)
   [sig_cl_feedback_ctr_xod, pvalue] = perform_ttest_1d_here(sw_cloud_feedback_ctr+lw_cloud_feedback_ctr,sw_cloud_feedback_xod+lw_cloud_feedback_xod,sig_level)
   [sig_cl_feedback_ctr_xoo, pvalue] = perform_ttest_1d_here(sw_cloud_feedback_ctr+lw_cloud_feedback_ctr,sw_cloud_feedback_xoo+lw_cloud_feedback_xoo,sig_level)
   [sig_cl_feedback_xod_xoo, pvalue] = perform_ttest_1d_here(sw_cloud_feedback_xod+lw_cloud_feedback_xod,sw_cloud_feedback_xoo+lw_cloud_feedback_xoo,sig_level)

   [sig_dh_ctr_xod, pvalue] = perform_ttest_1d_here(dh_ctr,dh_xod,sig_level)
   [sig_dh_ctr_xoo, pvalue] = perform_ttest_1d_here(dh_ctr,dh_xoo,sig_level)
   [sig_dh_xod_xoo, pvalue] = perform_ttest_1d_here(dh_xod,dh_xoo,sig_level)

   ts1 = p_f_ctr+planck_feedback_ctr+lapserate_feedback_ctr+q_feedback_ctr+alb_feedback_ctr+lw_cloud_feedback_ctr+sw_cloud_feedback_ctr
   ts2 = p_f_xod+planck_feedback_xod+lapserate_feedback_xod+q_feedback_xod+alb_feedback_xod+lw_cloud_feedback_xod+sw_cloud_feedback_xod
   [sig_sum_ctr_xod, pvalue] = perform_ttest_1d_here(ts1,ts2,sig_level)
   ts1 = p_f_ctr+planck_feedback_ctr+lapserate_feedback_ctr+q_feedback_ctr+alb_feedback_ctr+lw_cloud_feedback_ctr+sw_cloud_feedback_ctr
   ts2 = p_f_xoo+planck_feedback_xoo+lapserate_feedback_xoo+q_feedback_xoo+alb_feedback_xoo+lw_cloud_feedback_xoo+sw_cloud_feedback_xoo
   [sig_sum_ctr_xoo, pvalue] = perform_ttest_1d_here(ts1,ts2,sig_level)
   ts1 = p_f_xod+planck_feedback_xod+lapserate_feedback_xod+q_feedback_xod+alb_feedback_xod+lw_cloud_feedback_xod+sw_cloud_feedback_xod
   ts2 = p_f_xoo+planck_feedback_xoo+lapserate_feedback_xoo+q_feedback_xoo+alb_feedback_xoo+lw_cloud_feedback_xoo+sw_cloud_feedback_xoo
   [sig_sum_xod_xoo, pvalue] = perform_ttest_1d_here(ts1,ts2,sig_level)

   ts1 = dh_ctr - (FLNT_arctic_ctr+FSNT_arctic_ctr)
   ts2 = dh_xod - (FLNT_arctic_xod+FSNT_arctic_xod)
   [sig_df_ctr_xod, pvalue] = perform_ttest_1d_here(ts1,ts2,sig_level)
   ts1 = dh_ctr - (FLNT_arctic_ctr+FSNT_arctic_ctr)
   ts2 = dh_xoo - (FLNT_arctic_xoo+FSNT_arctic_xoo)
   [sig_df_ctr_xoo, pvalue] = perform_ttest_1d_here(ts1,ts2,sig_level)
   ts1 = dh_xod - (FLNT_arctic_xod+FSNT_arctic_xod)
   ts2 = dh_xoo - (FLNT_arctic_xoo+FSNT_arctic_xoo)
   [sig_df_xod_xoo, pvalue] = perform_ttest_1d_here(ts1,ts2,sig_level)

   ts1 = (FLNT_arctic_ctr+FSNT_arctic_ctr)
   ts2 = (FLNT_arctic_xod+FSNT_arctic_xod)
   [sig_net_ctr_xod, pvalue] = perform_ttest_1d_here(ts1,ts2,sig_level)
   ts1 = (FLNT_arctic_ctr+FSNT_arctic_ctr)
   ts2 = (FLNT_arctic_xoo+FSNT_arctic_xoo)
   [sig_net_ctr_xoo, pvalue] = perform_ttest_1d_here(ts1,ts2,sig_level)
   ts1 = (FLNT_arctic_xod+FSNT_arctic_xod)
   ts2 = (FLNT_arctic_xoo+FSNT_arctic_xoo)
   [sig_net_xod_xoo, pvalue] = perform_ttest_1d_here(ts1,ts2,sig_level)

# ================================================================
# plot figures
# ================================================================
if True:
   
   plt.close('all')
   fig = plt.figure()
   dts1 = (dts_ctr - dts_xod)
   dts2 = (dts_ctr - dts_xoo)
#   fig.set_size_inches(10, 10, forward=True)
#   ax1 = fig.add_axes([0.1, 0.2, 0.75, 0.55])
#   plt.fill_between([0.65, 0.95], 0, (p_f_ctr-p_f_xod)/dts1, color='r',alpha=0.5)
#   plt.fill_between([1.05, 1.35], 0, (p_f_ctr-p_f_xoo)/dts2, color='b',alpha=0.5)
#   plt.fill_between([1.65, 1.95], 0, (p_planck_ctr-p_planck_xod)/dts1, color='r',alpha=sig_planck_feedback_ctr_xod) 
#   plt.fill_between([2.05, 2.35], 0, (p_planck_ctr-p_planck_xoo)/dts2, color='b',alpha=sig_planck_feedback_ctr_xoo)
   plt.fill_between([2.65, 2.95], 0, (p_lr_ctr-p_lr_xod)/dts1, color='b',alpha=sig_lapserate_feedback_ctr_xod, label='Hist-xODS')
   plt.fill_between([3.05, 3.35], 0, (p_lr_ctr-p_lr_xoo)/dts2, color='r',alpha=sig_lapserate_feedback_ctr_xoo, label='Hist-xCO2')
   plt.fill_between([3.65, 3.95], 0, (p_wv_ctr-p_wv_xod)/dts1, color='b',alpha=sig_wv_feedback_ctr_xod)
   plt.fill_between([4.05, 4.35], 0, (p_wv_ctr-p_wv_xoo)/dts2, color='r',alpha=sig_wv_feedback_ctr_xoo)
   plt.fill_between([4.65, 4.95], 0, (p_al_ctr-p_al_xod)/dts1, color='b',alpha=sig_al_feedback_ctr_xod)
   plt.fill_between([5.05, 5.35], 0, (p_al_ctr-p_al_xoo)/dts2, color='r',alpha=sig_al_feedback_ctr_xoo)
   plt.fill_between([5.65, 5.95], 0, (p_cllw_ctr-p_cllw_xod)/dts1, color='b',alpha=sig_cllw_feedback_ctr_xod)
   plt.fill_between([6.05, 6.35], 0, (p_cllw_ctr-p_cllw_xoo)/dts2, color='r',alpha=sig_cllw_feedback_ctr_xoo)
   plt.fill_between([6.65, 6.95], 0, (p_clsw_ctr-p_clsw_xod)/dts1, color='b',alpha=sig_clsw_feedback_ctr_xod)
   plt.fill_between([7.05, 7.35], 0, (p_clsw_ctr-p_clsw_xoo)/dts2, color='r',alpha=sig_clsw_feedback_ctr_xoo)
   plt.fill_between([7.65, 7.95], 0, (p_clsw_ctr+p_cllw_ctr-p_clsw_xod-p_cllw_xod)/dts1, color='b',alpha=sig_cl_feedback_ctr_xod)
   plt.fill_between([8.05, 8.35], 0, (p_clsw_ctr+p_cllw_ctr-p_clsw_xoo-p_cllw_xod)/dts2, color='r',alpha=sig_cl_feedback_ctr_xoo)

#   plt.fill_between([7.65, 7.95], 0, (p_sum_ctr-p_sum_xod)/dts1, color='r',alpha=sig_sum_ctr_xod)
#   plt.fill_between([8.05, 8.35], 0, (p_sum_ctr-p_sum_xoo)/dts2, color='b',alpha=sig_sum_ctr_xoo)
#   plt.fill_between([8.65, 8.95], 0, ((FLNT_arctic_ctr+FSNT_arctic_ctr).mean()-(FLNT_arctic_xod+FSNT_arctic_xod).mean())/dts1, color='r',alpha=sig_net_ctr_xod)
#   plt.fill_between([9.05, 9.35], 0, ((FLNT_arctic_ctr+FSNT_arctic_ctr).mean()-(FLNT_arctic_xoo+FSNT_arctic_xoo).mean())/dts2, color='b',alpha=sig_net_ctr_xod)
#   plt.fill_between([9.65, 9.95], 0, (p_dh_ctr-p_dh_xod)/dts1, color='r',alpha=sig_dh_ctr_xod)
#   plt.fill_between([10.05, 10.35], 0, (p_dh_ctr-p_dh_xoo)/dts2, color='b',alpha=sig_dh_ctr_xoo)
   plt.fill_between([8.65, 8.95], 0, (p_df_ctr-p_df_xod)/dts1, color='b',alpha=sig_df_ctr_xod)
   plt.fill_between([9.05, 9.35], 0, (p_df_ctr-p_df_xoo)/dts2, color='r',alpha=sig_df_ctr_xoo)
   plt.legend(loc='upper left')  

   plt.plot([0,12],[0,0],'k-')
   plt.xticks([3,4,5,6,7,8,9],['lapserate','water\nvapor','albedo','cloud\nlongwave','cloud\nshortwave','cloud net','transport'], rotation=0, fontsize=9)
#   plt.axis([0.5, 3.5, -0.1, 0.8])  
   plt.xlim(2.5,9.5)
   plt.ylabel('W/m$^{2}/K$')

   plt.savefig('trend_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()

   sys.exit()

# plot spatial maps
   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)

   theta = np.linspace(0, 2*np.pi, 100)
   center, radius = [0.5, 0.5], 0.5
   verts = np.vstack([np.sin(theta), np.cos(theta)]).T
   circle = mpath.Path(verts * radius + center)

   lat_sel = 160
#   lat_sel = 0

# plot sat
   clevel = np.linspace(-10,10,21)

   ax1 = fig.add_axes([0.02, 0.65, 0.27, 0.27], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = planck_feedback_map[lat_sel:,:].copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('planck, ' + str(round(planck_feedback, 2)), fontsize=12)

   ax1 = fig.add_axes([0.31, 0.65, 0.27, 0.27], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = lapserate_feedback_map[lat_sel:,:].copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('lapse rate, ' + str(round(lapserate_feedback, 2)), fontsize=12)

   ax1 = fig.add_axes([0.6, 0.65, 0.27, 0.27], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = alb_feedback_map[lat_sel:,:].copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('albedo, ' + str(round(alb_feedback, 2)), fontsize=12)

   cbaxes = fig.add_axes([0.9, 0.65, 0.01, 0.27])
   cbar = plt.colorbar(im1, cax=cbaxes, orientation='vertical', ticks=np.linspace(clevel[0],clevel[-1],11))
   ax1 = fig.add_axes([0.96, 0.33, 0.01, 0.27])
#   ax1.text(0.01,0.05,'Annual SIC (% per decade)',fontsize=9, rotation=270)
   ax1.axis('off')

   clevel = np.linspace(-2,2,10)
   ax1 = fig.add_axes([0.02, 0.3, 0.27, 0.27], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = q_feedback_map[lat_sel:,:].copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('water vapor, ' + str(round(q_feedback, 2)), fontsize=12)

   ax1 = fig.add_axes([0.31, 0.3, 0.27, 0.27], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = sw_cloud_feedback_map[lat_sel:,:].copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('shortwave cloud, ' + str(round(sw_cloud_feedback, 2)), fontsize=12)

   ax1 = fig.add_axes([0.6, 0.3, 0.27, 0.27], projection=ccrs.Orthographic(0, 90))
   ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
   map_2d = lw_cloud_feedback_map[lat_sel:,:].copy()
   [x_out, lon_x] = data_process_f.extend_longitude(map_2d,lat[lat_sel:],lon)
   im1 = ax1.contourf(lon_x, lat[lat_sel:], x_out, levels=clevel, extend='both', transform=ccrs.PlateCarree(), cmap='RdBu_r')
   ax1.coastlines('110m',color='k',linewidth=0.7)
   ax1.set_boundary(circle, transform=ax1.transAxes)
   ax1.set_aspect('auto')
   plt.title('longwave cloud, ' + str(round(lw_cloud_feedback, 2)), fontsize=12)

   cbaxes = fig.add_axes([0.9, 0.3, 0.01, 0.27])
   cbar = plt.colorbar(im1, cax=cbaxes, orientation='vertical', ticks=np.linspace(clevel[0],clevel[-1],11))
#   ax1 = fig.add_axes([0.96, 0.33, 0.01, 0.27])
#   ax1.text(0.01,0.05,'Annual SIC (% per decade)',fontsize=9, rotation=270)
   ax1.axis('off')

   plt.savefig('trend_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()

   sys.exit()






