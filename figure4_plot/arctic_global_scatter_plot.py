flag_run = 0
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
from sklearn.utils import resample
import seaborn as sns

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
   planck_feedback_ctr = f.variables['planck_feedback_all'][st_index:].data
   lapserate_feedback_ctr = f.variables['lapserate_feedback_all'][st_index:].data
   lapserate_feedback_map_ctr = f.variables['lapserate_feedback_map_all'][st_index:,:,:].data
   alb_feedback_ctr = f.variables['alb_feedback_all'][st_index:].data
   q_feedback_ctr = f.variables['q_feedback_all'][st_index:].data
   lw_cloud_feedback_ctr = f.variables['lw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_ctr = f.variables['sw_cloud_feedback_all'][st_index:].data
   dts_arctic_ctr = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_ctr = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_ctr = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close() 

   ens_num = len(t_feedback_ctr)

   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
   filename = 'feedback_parameters_maps_ctr_global_temp.nc'
   f = Dataset(dirname + filename, 'r')
   t_feedback_ctr_g = f.variables['t_feedback_all'][st_index:].data
   planck_feedback_ctr_g = f.variables['planck_feedback_all'][st_index:].data
   lapserate_feedback_ctr_g = f.variables['lapserate_feedback_all'][st_index:].data
   lapserate_feedback_map_ctr_g = f.variables['lapserate_feedback_map_all'][st_index:,:,:].data
   alb_feedback_ctr_g = f.variables['alb_feedback_all'][st_index:].data
   q_feedback_ctr_g = f.variables['q_feedback_all'][st_index:].data
   lw_cloud_feedback_ctr_g = f.variables['lw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_ctr_g = f.variables['sw_cloud_feedback_all'][st_index:].data
   dts_arctic_ctr_g = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_ctr_g = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_ctr_g = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close()

   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
   filename = 'feedback_parameters_maps_xc2_temp.nc'
   f = Dataset(dirname + filename, 'r')
   t_feedback_xoo = f.variables['t_feedback_all'][st_index:].data
   planck_feedback_xoo = f.variables['planck_feedback_all'][st_index:].data
   lapserate_feedback_xoo = f.variables['lapserate_feedback_all'][st_index:].data
   alb_feedback_xoo = f.variables['alb_feedback_all'][st_index:].data
   q_feedback_xoo = f.variables['q_feedback_all'][st_index:].data
   lw_cloud_feedback_xoo = f.variables['lw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_xoo = f.variables['sw_cloud_feedback_all'][st_index:].data
   dts_arctic_xoo = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_xoo = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_xoo = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close()

   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
   filename = 'feedback_parameters_maps_xc2_global_temp.nc'
   f = Dataset(dirname + filename, 'r')
   t_feedback_xoo_g = f.variables['t_feedback_all'][st_index:].data
   planck_feedback_xoo_g = f.variables['planck_feedback_all'][st_index:].data
   lapserate_feedback_xoo_g = f.variables['lapserate_feedback_all'][st_index:].data
   alb_feedback_xoo_g = f.variables['alb_feedback_all'][st_index:].data
   q_feedback_xoo_g = f.variables['q_feedback_all'][st_index:].data
   lw_cloud_feedback_xoo_g = f.variables['lw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_xoo_g = f.variables['sw_cloud_feedback_all'][st_index:].data
   dts_arctic_xoo_g = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_xoo_g = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_xoo_g = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close()


# read xod case
   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
   filename = 'feedback_parameters_maps_xod_temp.nc'
   f = Dataset(dirname + filename, 'r')
   t_feedback_xod = f.variables['t_feedback_all'][st_index:].data
   planck_feedback_xod = f.variables['planck_feedback_all'][st_index:].data
   lapserate_feedback_xod = f.variables['lapserate_feedback_all'][st_index:].data
   alb_feedback_xod = f.variables['alb_feedback_all'][st_index:].data
   q_feedback_xod = f.variables['q_feedback_all'][st_index:].data
   lw_cloud_feedback_xod = f.variables['lw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_xod = f.variables['sw_cloud_feedback_all'][st_index:].data
   dts_arctic_xod = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_xod = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_xod = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close()

   dirname = '/home/yliang/research/co2_ods/change_metric/radiative_kernel/'
   filename = 'feedback_parameters_maps_xod_global_temp.nc'
   f = Dataset(dirname + filename, 'r')
   t_feedback_xod_g = f.variables['t_feedback_all'][st_index:].data
   planck_feedback_xod_g = f.variables['planck_feedback_all'][st_index:].data
   lapserate_feedback_xod_g = f.variables['lapserate_feedback_all'][st_index:].data
   alb_feedback_xod_g = f.variables['alb_feedback_all'][st_index:].data
   q_feedback_xod_g = f.variables['q_feedback_all'][st_index:].data
   lw_cloud_feedback_xod_g = f.variables['lw_cloud_feedback_all'][st_index:].data
   sw_cloud_feedback_xod_g = f.variables['sw_cloud_feedback_all'][st_index:].data
   dts_arctic_xod_g = f.variables['dts_arctic_all'][:].data
   FSNT_arctic_xod_g = f.variables['FSNT_arctic_all'][:].data
   FLNT_arctic_xod_g = f.variables['FLNT_arctic_all'][:].data*-1.
   f.close()

# ================================================================
# compute feedback parameters
# ================================================================
   p_f_ctr = 1.26
   p_f_xoo = 1.01
   p_f_xod = 1.11

   p_pl_ctr_xod = (planck_feedback_ctr-planck_feedback_xod)/(dts_arctic_ctr-dts_arctic_xod)
   p_lr_ctr_xod = (lapserate_feedback_ctr-lapserate_feedback_xod)/(dts_arctic_ctr-dts_arctic_xod)
   p_wv_ctr_xod = (q_feedback_ctr-q_feedback_xod)/(dts_arctic_ctr-dts_arctic_xod)
   p_al_ctr_xod = (alb_feedback_ctr-alb_feedback_xod)/(dts_arctic_ctr-dts_arctic_xod)
   p_cl_ctr_xod = (lw_cloud_feedback_ctr+sw_cloud_feedback_ctr-lw_cloud_feedback_xod-sw_cloud_feedback_xod)/(dts_arctic_ctr-dts_arctic_xod)

   p_pl_ctr_xoo = (planck_feedback_ctr-planck_feedback_xoo)/(dts_arctic_ctr-dts_arctic_xoo)
   p_lr_ctr_xoo = (lapserate_feedback_ctr-lapserate_feedback_xoo)/(dts_arctic_ctr-dts_arctic_xoo)
   p_wv_ctr_xoo = (q_feedback_ctr-q_feedback_xoo)/(dts_arctic_ctr-dts_arctic_xoo)
   p_al_ctr_xoo = (alb_feedback_ctr-alb_feedback_xoo)/(dts_arctic_ctr-dts_arctic_xoo)
   p_cl_ctr_xoo = (lw_cloud_feedback_ctr+sw_cloud_feedback_ctr-lw_cloud_feedback_xoo-sw_cloud_feedback_xoo)/(dts_arctic_ctr-dts_arctic_xoo)

   p_pl_ctr_xod_g = (planck_feedback_ctr_g-planck_feedback_xod_g)/(dts_arctic_ctr_g-dts_arctic_xod_g)
   p_lr_ctr_xod_g = (lapserate_feedback_ctr_g-lapserate_feedback_xod_g)/(dts_arctic_ctr_g-dts_arctic_xod_g)
   p_wv_ctr_xod_g = (q_feedback_ctr_g-q_feedback_xod_g)/(dts_arctic_ctr_g-dts_arctic_xod_g)
   p_al_ctr_xod_g = (alb_feedback_ctr_g-alb_feedback_xod_g)/(dts_arctic_ctr_g-dts_arctic_xod_g)
   p_cl_ctr_xod_g = (lw_cloud_feedback_ctr_g+sw_cloud_feedback_ctr_g-lw_cloud_feedback_xod_g-sw_cloud_feedback_xod_g)/(dts_arctic_ctr_g-dts_arctic_xod_g)

   p_pl_ctr_xoo_g = (planck_feedback_ctr_g-planck_feedback_xoo_g)/(dts_arctic_ctr_g-dts_arctic_xoo_g)
   p_lr_ctr_xoo_g = (lapserate_feedback_ctr_g-lapserate_feedback_xoo_g)/(dts_arctic_ctr_g-dts_arctic_xoo_g)
   p_wv_ctr_xoo_g = (q_feedback_ctr_g-q_feedback_xoo_g)/(dts_arctic_ctr_g-dts_arctic_xoo_g)
   p_al_ctr_xoo_g = (alb_feedback_ctr_g-alb_feedback_xoo_g)/(dts_arctic_ctr_g-dts_arctic_xoo_g)
   p_cl_ctr_xoo_g = (lw_cloud_feedback_ctr_g+sw_cloud_feedback_ctr_g-lw_cloud_feedback_xoo_g-sw_cloud_feedback_xoo_g)/(dts_arctic_ctr_g-dts_arctic_xoo_g)

# ================================================================
# calculate significance
# ================================================================
   sig_level = 0.05
#   dts3 = (dts_xod - dts_xoo)
   [sig_pl_xod_xoo, pvalue_pl_xod_xoo] = perform_ttest_1d_here(p_pl_ctr_xod,p_pl_ctr_xoo,sig_level)
   [sig_lr_xod_xoo, pvalue_lr_xod_xoo] = perform_ttest_1d_here(p_lr_ctr_xod,p_lr_ctr_xoo,sig_level)
   [sig_wv_xod_xoo, pvalue_wv_xod_xoo] = perform_ttest_1d_here(p_wv_ctr_xod,p_wv_ctr_xoo,sig_level)
   [sig_al_xod_xoo, pvalue_al_xod_xoo] = perform_ttest_1d_here(p_al_ctr_xod,p_al_ctr_xoo,sig_level)
   [sig_cl_xod_xoo, pvalue_cl_xod_xoo] = perform_ttest_1d_here(p_cl_ctr_xod,p_cl_ctr_xoo,sig_level)

# ================================================================
# resampling
# ================================================================
   nsample = 10000
   pl_xod_resample = np.zeros((nsample,ens_num))
   lr_xod_resample = np.zeros((nsample,ens_num))
   wv_xod_resample = np.zeros((nsample,ens_num))
   al_xod_resample = np.zeros((nsample,ens_num))
   cl_xod_resample = np.zeros((nsample,ens_num))

   pl_xoo_resample = np.zeros((nsample,ens_num))
   lr_xoo_resample = np.zeros((nsample,ens_num))
   wv_xoo_resample = np.zeros((nsample,ens_num))
   al_xoo_resample = np.zeros((nsample,ens_num))
   cl_xoo_resample = np.zeros((nsample,ens_num))

   pl_xod_resample_g = np.zeros((nsample,ens_num))
   lr_xod_resample_g = np.zeros((nsample,ens_num))
   wv_xod_resample_g = np.zeros((nsample,ens_num))
   al_xod_resample_g = np.zeros((nsample,ens_num))
   cl_xod_resample_g = np.zeros((nsample,ens_num))
   
   pl_xoo_resample_g = np.zeros((nsample,ens_num))
   lr_xoo_resample_g = np.zeros((nsample,ens_num))
   wv_xoo_resample_g = np.zeros((nsample,ens_num))
   al_xoo_resample_g = np.zeros((nsample,ens_num))
   cl_xoo_resample_g = np.zeros((nsample,ens_num))

   for NS in range(nsample):
       print(NS)

       boot = resample(p_pl_ctr_xod, replace=True, n_samples=10, random_state=NS)
       pl_xod_resample[NS,:] = boot.copy()
       boot = resample(p_lr_ctr_xod, replace=True, n_samples=10, random_state=NS)
       lr_xod_resample[NS,:] = boot.copy()
       boot = resample(p_wv_ctr_xod, replace=True, n_samples=10, random_state=NS)
       wv_xod_resample[NS,:] = boot.copy()
       boot = resample(p_al_ctr_xod, replace=True, n_samples=10, random_state=NS)
       al_xod_resample[NS,:] = boot.copy()
       boot = resample(p_cl_ctr_xod, replace=True, n_samples=10, random_state=NS)
       cl_xod_resample[NS,:] = boot.copy()

       boot = resample(p_pl_ctr_xoo, replace=True, n_samples=10, random_state=NS)
       pl_xoo_resample[NS,:] = boot.copy()
       boot = resample(p_lr_ctr_xoo, replace=True, n_samples=10, random_state=NS)
       lr_xoo_resample[NS,:] = boot.copy()
       boot = resample(p_wv_ctr_xoo, replace=True, n_samples=10, random_state=NS)
       wv_xoo_resample[NS,:] = boot.copy()
       boot = resample(p_al_ctr_xoo, replace=True, n_samples=10, random_state=NS)
       al_xoo_resample[NS,:] = boot.copy()
       boot = resample(p_cl_ctr_xoo, replace=True, n_samples=10, random_state=NS)
       cl_xoo_resample[NS,:] = boot.copy()

       boot = resample(p_pl_ctr_xod_g, replace=True, n_samples=10, random_state=NS)
       pl_xod_resample_g[NS,:] = boot.copy()
       boot = resample(p_lr_ctr_xod_g, replace=True, n_samples=10, random_state=NS)
       lr_xod_resample_g[NS,:] = boot.copy()
       boot = resample(p_wv_ctr_xod_g, replace=True, n_samples=10, random_state=NS)
       wv_xod_resample_g[NS,:] = boot.copy()
       boot = resample(p_al_ctr_xod_g, replace=True, n_samples=10, random_state=NS)
       al_xod_resample_g[NS,:] = boot.copy()
       boot = resample(p_cl_ctr_xod_g, replace=True, n_samples=10, random_state=NS)
       cl_xod_resample_g[NS,:] = boot.copy()

       boot = resample(p_pl_ctr_xoo_g, replace=True, n_samples=10, random_state=NS)
       pl_xoo_resample_g[NS,:] = boot.copy()
       boot = resample(p_lr_ctr_xoo_g, replace=True, n_samples=10, random_state=NS)
       lr_xoo_resample_g[NS,:] = boot.copy()
       boot = resample(p_wv_ctr_xoo_g, replace=True, n_samples=10, random_state=NS)
       wv_xoo_resample_g[NS,:] = boot.copy()
       boot = resample(p_al_ctr_xoo_g, replace=True, n_samples=10, random_state=NS)
       al_xoo_resample_g[NS,:] = boot.copy()
       boot = resample(p_cl_ctr_xoo_g, replace=True, n_samples=10, random_state=NS)
       cl_xoo_resample_g[NS,:] = boot.copy()

# ================================================================
# plot figures
# ================================================================
if True:

   factor1 = 3.2

#   To_ctr_xod = ((FSNT_arctic_ctr+FLNT_arctic_ctr).mean()-(FSNT_arctic_xod+FLNT_arctic_xod).mean())/factor1
#   To_ctr_xoo = ((FSNT_arctic_ctr+FLNT_arctic_ctr).mean()-(FSNT_arctic_xoo+FLNT_arctic_xoo).mean())/factor1
#   To_ctr_xod_g = ((FSNT_arctic_ctr_g+FLNT_arctic_ctr_g).mean()-(FSNT_arctic_xod_g+FLNT_arctic_xod_g).mean())/factor1
#   To_ctr_xoo_g = ((FSNT_arctic_ctr_g+FLNT_arctic_ctr_g).mean()-(FSNT_arctic_xoo_g+FLNT_arctic_xoo_g).mean())/factor1

#   dts1 = To_ctr_xoo*(1-(p_lr_ctr_xoo.mean()+p_al_ctr_xoo.mean()+p_cl_ctr_xoo.mean()+p_wv_ctr_xoo.mean())/factor1)
#   dts2 = To_ctr_xod*(1-(p_lr_ctr_xod.mean()+p_al_ctr_xod.mean()+p_cl_ctr_xod.mean()+p_wv_ctr_xod.mean())/factor1)
#   dts1_g = To_ctr_xoo_g*(1-(p_lr_ctr_xoo_g.mean()+p_al_ctr_xoo_g.mean()+p_cl_ctr_xoo_g.mean()+p_wv_ctr_xoo_g.mean())/factor1)
#   dts2_g = To_ctr_xod_g*(1-(p_lr_ctr_xod_g.mean()+p_al_ctr_xod_g.mean()+p_cl_ctr_xod_g.mean()+p_wv_ctr_xod_g.mean())/factor1)

#   dts1 = (dts_arctic_ctr-dts_arctic_xoo).mean()
#   dts2 = (dts_arctic_ctr-dts_arctic_xod).mean()
#   dts1_g = (dts_arctic_ctr_g-dts_arctic_xoo_g).mean()
#   dts2_g = (dts_arctic_ctr_g-dts_arctic_xod_g).mean()

   dts1 = 1.
   dts2 = 1.
   dts1_g = 1.
   dts2_g = 1.

   plt.close('all')
   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)

   ax1 = fig.add_axes([0.09, 0.6, 0.4, 0.35])
   plt.plot(np.nanmean(lr_xoo_resample_g, axis=1)/factor1*dts1_g,np.nanmean(lr_xoo_resample, axis=1)/factor1*dts1, 'm.', alpha=0.05)
   plt.plot(np.nanmean(lr_xod_resample_g, axis=1)/factor1*dts2_g,np.nanmean(lr_xod_resample, axis=1)/factor1*dts2, 'c.', alpha=0.05)
   plt.plot([-2,2],[-2,2],'--', color='darkgray')
   ts_x = np.nanmean(lr_xoo_resample_g, axis=1)/factor1*dts1_g
   ts_y = np.nanmean(lr_xoo_resample, axis=1)/factor1*dts1
   p_xoo_res = stats.theilslopes(ts_y,ts_x)
   p_xoo = np.polyfit(ts_x,ts_y,1)
   ts_x = np.nanmean(lr_xod_resample_g, axis=1)/factor1*dts2_g
   ts_y = np.nanmean(lr_xod_resample, axis=1)/factor1*dts2
   p_xod_res = stats.theilslopes(ts_y,ts_x)
   p_xod = np.polyfit(ts_x,ts_y,1)
   text1 = 'ODS:' + str(round(p_xod_res[0],2)) + '$\pm$' + str(round(abs(p_xod_res[0]-p_xod_res[2]), 2))
   text2 = 'CO2:' + str(round(p_xoo_res[0],2)) + '$\pm$' + str(round(abs(p_xoo_res[0]-p_xoo_res[2]), 2))

   plt.plot(p_lr_ctr_xod_g/factor1*dts2_g,p_lr_ctr_xod/factor1*dts2,'bo',markersize=3)
   plt.plot(p_lr_ctr_xoo_g/factor1*dts1_g,p_lr_ctr_xoo/factor1*dts1,'ro',markersize=3)
   plt.plot(p_lr_ctr_xod_g.mean()/factor1*dts2_g,p_lr_ctr_xod.mean()/factor1*dts2,'bo',markersize=10, label=text1)
   plt.plot(p_lr_ctr_xoo_g.mean()/factor1*dts1_g,p_lr_ctr_xoo.mean()/factor1*dts1,'ro',markersize=10, label=text2)
   plt.plot([-5,5],[p_xoo[0]*-5+p_xoo[1],p_xoo[0]*5+p_xoo[1]],'r-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xoo_res[1] + p_xoo_res[2]*np.array([-5,5]), 'r--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xoo_res[1] + p_xoo_res[3]*np.array([-5,5]), 'r--', linewidth=0.5)
   plt.plot([-5,5],[p_xod[0]*-5+p_xod[1],p_xod[0]*5+p_xod[1]],'b-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xod_res[1] + p_xod_res[2]*np.array([-5,5]), 'b--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xod_res[1] + p_xod_res[3]*np.array([-5,5]), 'b--', linewidth=0.5)
   plt.title('(a) Lapserate')
   plt.axis([-0.15, 0.05, -0.05, 0.45])
   plt.legend(loc='upper left')
#   plt.ylabel('Arctic warming (K)')
#   plt.xticks([-0.12, -0.09, -0.06, -0.03, 0, 0.03])

   ax1 = fig.add_axes([0.56, 0.6, 0.4, 0.35])
   plt.plot(np.nanmean(al_xoo_resample_g, axis=1)/factor1*dts1_g,np.nanmean(al_xoo_resample, axis=1)/factor1*dts1, 'm.', alpha=0.05)
   plt.plot(np.nanmean(al_xod_resample_g, axis=1)/factor1*dts2_g,np.nanmean(al_xod_resample, axis=1)/factor1*dts2, 'c.', alpha=0.05)
   plt.plot([-2,2],[-2,2],'--', color='darkgray')
   ts_x = np.nanmean(al_xoo_resample_g, axis=1)/factor1*dts1_g
   ts_y = np.nanmean(al_xoo_resample, axis=1)/factor1*dts1
   p_xoo_res = stats.theilslopes(ts_y,ts_x)
   p_xoo = np.polyfit(ts_x,ts_y,1)
   ts_x = np.nanmean(al_xod_resample_g, axis=1)/factor1*dts2_g
   ts_y = np.nanmean(al_xod_resample, axis=1)/factor1*dts2
   p_xod_res = stats.theilslopes(ts_y,ts_x)
   p_xod = np.polyfit(ts_x,ts_y,1)
   text1 = 'ODS:' + str(round(p_xod_res[0],2)) + '$\pm$' + str(round(abs(p_xod_res[0]-p_xod_res[2]), 2))
   text2 = 'CO2:' + str(round(p_xoo_res[0],2)) + '$\pm$' + str(round(abs(p_xoo_res[0]-p_xoo_res[2]), 2))

   plt.plot(p_al_ctr_xod_g/factor1*dts2_g,p_al_ctr_xod/factor1*dts2,'bo',markersize=3)
   plt.plot(p_al_ctr_xoo_g/factor1*dts1_g,p_al_ctr_xoo/factor1*dts1,'ro',markersize=3)
   plt.plot(p_al_ctr_xod_g.mean()/factor1*dts2_g,p_al_ctr_xod.mean()/factor1*dts2,'bo',markersize=10, label=text1)
   plt.plot(p_al_ctr_xoo_g.mean()/factor1*dts1_g,p_al_ctr_xoo.mean()/factor1*dts1,'ro',markersize=10, label=text2)
   plt.plot([-5,5],[p_xoo[0]*-5+p_xoo[1],p_xoo[0]*5+p_xoo[1]],'r-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xoo_res[1] + p_xoo_res[2]*np.array([-5,5]), 'r--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xoo_res[1] + p_xoo_res[3]*np.array([-5,5]), 'r--', linewidth=0.5)
   plt.plot([-5,5],[p_xod[0]*-5+p_xod[1],p_xod[0]*5+p_xod[1]],'b-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xod_res[1] + p_xod_res[2]*np.array([-5,5]), 'b--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xod_res[1] + p_xod_res[3]*np.array([-5,5]), 'b--', linewidth=0.5)
   plt.title('(b) Albedo')
#   plt.axis([-0., 0.1, 0.05, 0.7])
   plt.legend(loc='upper left')

   ax1 = fig.add_axes([0.09, 0.17, 0.4, 0.35])
   plt.plot(np.nanmean(wv_xoo_resample_g, axis=1)/factor1*dts1_g,np.nanmean(wv_xoo_resample, axis=1)/factor1*dts1, 'm.', alpha=0.05)
   plt.plot(np.nanmean(wv_xod_resample_g, axis=1)/factor1*dts2_g,np.nanmean(wv_xod_resample, axis=1)/factor1*dts2, 'c.', alpha=0.05)
   plt.plot([-2,2],[-2,2],'--', color='darkgray')
   ts_x = np.nanmean(wv_xoo_resample_g, axis=1)/factor1*dts1_g
   ts_y = np.nanmean(wv_xoo_resample, axis=1)/factor1*dts1
   p_xoo_res = stats.theilslopes(ts_y,ts_x)
   p_xoo = np.polyfit(ts_x,ts_y,1)
   ts_x = np.nanmean(wv_xod_resample_g, axis=1)/factor1*dts2_g
   ts_y = np.nanmean(wv_xod_resample, axis=1)/factor1*dts2
   p_xod_res = stats.theilslopes(ts_y,ts_x)
   p_xod = np.polyfit(ts_x,ts_y,1)
   text1 = 'ODS:' + str(round(p_xod_res[0],2)) + '$\pm$' + str(round(abs(p_xod_res[0]-p_xod_res[2]), 2))
   text2 = 'CO2:' + str(round(p_xoo_res[0],2)) + '$\pm$' + str(round(abs(p_xoo_res[0]-p_xoo_res[2]), 2))

   plt.plot(p_wv_ctr_xod_g/factor1*dts2_g,p_wv_ctr_xod/factor1*dts2,'bo',markersize=3)
   plt.plot(p_wv_ctr_xoo_g/factor1*dts1_g,p_wv_ctr_xoo/factor1*dts1,'ro',markersize=3)
   plt.plot(p_wv_ctr_xod_g.mean()/factor1*dts2_g,p_wv_ctr_xod.mean()/factor1*dts2,'bo',markersize=10, label=text1)
   plt.plot(p_wv_ctr_xoo_g.mean()/factor1*dts1_g,p_wv_ctr_xoo.mean()/factor1*dts1,'ro',markersize=10, label=text2)
   plt.plot([-5,5],[p_xoo[0]*-5+p_xoo[1],p_xoo[0]*5+p_xoo[1]],'r-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xoo_res[1] + p_xoo_res[2]*np.array([-5,5]), 'r--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xoo_res[1] + p_xoo_res[3]*np.array([-5,5]), 'r--', linewidth=0.5)
   plt.plot([-5,5],[p_xod[0]*-5+p_xod[1],p_xod[0]*5+p_xod[1]],'b-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xod_res[1] + p_xod_res[2]*np.array([-5,5]), 'b--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xod_res[1] + p_xod_res[3]*np.array([-5,5]), 'b--', linewidth=0.5)
   plt.title('(c) Water vapor')
#   plt.axis([-0.05, 0.3, -0., 0.25])
#   plt.ylabel('Arctic warming (K)')
#   plt.xlabel('Global warming (K)')
   plt.legend(loc='upper left')

   ax1 = fig.add_axes([0.56, 0.17, 0.4, 0.35])
   plt.plot(np.nanmean(cl_xoo_resample_g, axis=1)/factor1*dts1_g,np.nanmean(cl_xoo_resample, axis=1)/factor1*dts1, 'm.', alpha=0.05)
   plt.plot(np.nanmean(cl_xod_resample_g, axis=1)/factor1*dts2_g,np.nanmean(cl_xod_resample, axis=1)/factor1*dts2, 'c.', alpha=0.05)
   plt.plot([-2,2],[-2,2],'--', color='darkgray')
   ts_x = np.nanmean(cl_xoo_resample_g, axis=1)/factor1*dts1_g
   ts_y = np.nanmean(cl_xoo_resample, axis=1)/factor1*dts1
   p_xoo_res = stats.theilslopes(ts_y,ts_x)
   p_xoo = np.polyfit(ts_x,ts_y,1)
   ts_x = np.nanmean(cl_xod_resample_g, axis=1)/factor1*dts2_g
   ts_y = np.nanmean(cl_xod_resample, axis=1)/factor1*dts2
   p_xod_res = stats.theilslopes(ts_y,ts_x)
   p_xod = np.polyfit(ts_x,ts_y,1)
   text1 = 'ODS:' + str(round(p_xod_res[0],2)) + '$\pm$' + str(round(abs(p_xod_res[0]-p_xod_res[2]), 2))
   text2 = 'CO2:' + str(round(p_xoo_res[0],2)) + '$\pm$' + str(round(abs(p_xoo_res[0]-p_xoo_res[2]), 2))

   plt.plot(p_cl_ctr_xod_g/factor1*dts2_g,p_cl_ctr_xod/factor1*dts2,'bo',markersize=3)
   plt.plot(p_cl_ctr_xoo_g/factor1*dts1_g,p_cl_ctr_xoo/factor1*dts1,'ro',markersize=3)
   plt.plot(p_cl_ctr_xod_g.mean()/factor1*dts2_g,p_cl_ctr_xod.mean()/factor1*dts2,'bo',markersize=10, label=text1)
   plt.plot(p_cl_ctr_xoo_g.mean()/factor1*dts1_g,p_cl_ctr_xoo.mean()/factor1*dts1,'ro',markersize=10, label=text2)
   plt.plot([-5,5],[p_xoo[0]*-5+p_xoo[1],p_xoo[0]*5+p_xoo[1]],'r-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xoo_res[1] + p_xoo_res[2]*np.array([-5,5]), 'r--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xoo_res[1] + p_xoo_res[3]*np.array([-5,5]), 'r--', linewidth=0.5)
   plt.plot([-5,5],[p_xod[0]*-5+p_xod[1],p_xod[0]*5+p_xod[1]],'b-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xod_res[1] + p_xod_res[2]*np.array([-5,5]), 'b--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xod_res[1] + p_xod_res[3]*np.array([-5,5]), 'b--', linewidth=0.5)
   plt.title('(d) Cloud')
#   plt.axis([-0.03, 0.065, -0.2, 0.15])
#   plt.xlabel('Global warming (K)')
   plt.legend(loc='upper left')

   plt.savefig('trend_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()







