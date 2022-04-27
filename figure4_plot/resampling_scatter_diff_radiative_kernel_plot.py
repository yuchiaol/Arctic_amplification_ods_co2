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
   p_f_xoo = 1.01
   p_f_xod = 1.11

# co2
   p_f_ctr_xoo = 0.74
# ods
   p_f_ctr_xod = 0.19

   p_pl_ctr_xod = (planck_feedback_ctr-planck_feedback_xod)/(dts_arctic_ctr-dts_arctic_xod)
   p_lr_ctr_xod = (lapserate_feedback_ctr-lapserate_feedback_xod)/(dts_arctic_ctr-dts_arctic_xod)
   p_wv_ctr_xod = (q_feedback_ctr-q_feedback_xod)/(dts_arctic_ctr-dts_arctic_xod)
   p_al_ctr_xod = (alb_feedback_ctr-alb_feedback_xod)/(dts_arctic_ctr-dts_arctic_xod)
   p_cl_ctr_xod = (lw_cloud_feedback_ctr+sw_cloud_feedback_ctr-lw_cloud_feedback_xod-sw_cloud_feedback_xod)/(dts_arctic_ctr-dts_arctic_xod)
   p_df_ctr_xod = ((dh_ctr-(FLNT_arctic_ctr+FSNT_arctic_ctr))-(dh_xod-(FLNT_arctic_xod+FSNT_arctic_xod)))/(dts_arctic_ctr-dts_arctic_xod)

   p_pl_ctr_xoo = (planck_feedback_ctr-planck_feedback_xoo)/(dts_arctic_ctr-dts_arctic_xoo)
   p_lr_ctr_xoo = (lapserate_feedback_ctr-lapserate_feedback_xoo)/(dts_arctic_ctr-dts_arctic_xoo)
   p_wv_ctr_xoo = (q_feedback_ctr-q_feedback_xoo)/(dts_arctic_ctr-dts_arctic_xoo)
   p_al_ctr_xoo = (alb_feedback_ctr-alb_feedback_xoo)/(dts_arctic_ctr-dts_arctic_xoo)
   p_cl_ctr_xoo = (lw_cloud_feedback_ctr+sw_cloud_feedback_ctr-lw_cloud_feedback_xoo-sw_cloud_feedback_xoo)/(dts_arctic_ctr-dts_arctic_xoo)
   p_df_ctr_xoo = ((dh_ctr-(FLNT_arctic_ctr+FSNT_arctic_ctr))-(dh_xoo-(FLNT_arctic_xoo+FSNT_arctic_xoo)))/(dts_arctic_ctr-dts_arctic_xoo)

   area_ctr_xoo = 1./(dts_arctic_ctr-dts_arctic_xoo) 
   area_ctr_xod = 1./(dts_arctic_ctr-dts_arctic_xod)

# ================================================================
# calculate significance
# ================================================================
   sig_level = 0.05
   dts1 = (dts_ctr - dts_xod)
   dts2 = (dts_ctr - dts_xoo)
#   dts3 = (dts_xod - dts_xoo)
   [sig_pl_xod_xoo, pvalue_pl_xod_xoo] = perform_ttest_1d_here(p_pl_ctr_xod,p_pl_ctr_xoo,sig_level)
   [sig_lr_xod_xoo, pvalue_lr_xod_xoo] = perform_ttest_1d_here(p_lr_ctr_xod,p_lr_ctr_xoo,sig_level)
   [sig_wv_xod_xoo, pvalue_wv_xod_xoo] = perform_ttest_1d_here(p_wv_ctr_xod,p_wv_ctr_xoo,sig_level)
   [sig_al_xod_xoo, pvalue_al_xod_xoo] = perform_ttest_1d_here(p_al_ctr_xod,p_al_ctr_xoo,sig_level)
   [sig_cl_xod_xoo, pvalue_cl_xod_xoo] = perform_ttest_1d_here(p_cl_ctr_xod,p_cl_ctr_xoo,sig_level)
   [sig_df_xod_xoo, pvalue_df_xod_xoo] = perform_ttest_1d_here(p_df_ctr_xod,p_df_ctr_xoo,sig_level)

# ================================================================
# resampling
# ================================================================
   nsample = 10000
   planck_feedback_ctr_resample = np.zeros((nsample,ens_num*2))
   planck_feedback_xod_resample = np.zeros((nsample,ens_num*2))
   planck_feedback_xoo_resample = np.zeros((nsample,ens_num*2))
   lapserate_feedback_ctr_resample = np.zeros((nsample,ens_num*2))
   lapserate_feedback_xod_resample = np.zeros((nsample,ens_num*2))
   lapserate_feedback_xoo_resample = np.zeros((nsample,ens_num*2))
   q_feedback_ctr_resample = np.zeros((nsample,ens_num*2))
   q_feedback_xod_resample = np.zeros((nsample,ens_num*2))
   q_feedback_xoo_resample = np.zeros((nsample,ens_num*2))
   alb_feedback_ctr_resample = np.zeros((nsample,ens_num*2))
   alb_feedback_xod_resample = np.zeros((nsample,ens_num*2))
   alb_feedback_xoo_resample = np.zeros((nsample,ens_num*2))
   cloud_feedback_ctr_resample = np.zeros((nsample,ens_num*2))
   cloud_feedback_xod_resample = np.zeros((nsample,ens_num*2))
   cloud_feedback_xoo_resample = np.zeros((nsample,ens_num*2))
   df_ctr_resample = np.zeros((nsample,ens_num*2))
   df_xod_resample = np.zeros((nsample,ens_num*2))
   df_xoo_resample = np.zeros((nsample,ens_num*2))
   dts_ctr_resample = np.zeros((nsample,ens_num*2))
   dts_xod_resample = np.zeros((nsample,ens_num*2))
   dts_xoo_resample = np.zeros((nsample,ens_num*2))

   for NS in range(nsample):
       print(NS)

       boot = resample(planck_feedback_ctr, replace=True, n_samples=20, random_state=NS)
       planck_feedback_ctr_resample[NS,:] = boot.copy()
       boot = resample(planck_feedback_xod, replace=True, n_samples=20, random_state=NS)
       planck_feedback_xod_resample[NS,:] = boot.copy()
       boot = resample(planck_feedback_xoo, replace=True, n_samples=20, random_state=NS)
       planck_feedback_xoo_resample[NS,:] = boot.copy()

       boot = resample(lapserate_feedback_ctr, replace=True, n_samples=20, random_state=NS)
#       boot = np.random.choice(lapserate_feedback_ctr, size=10)
       lapserate_feedback_ctr_resample[NS,:] = boot.copy()
       boot = resample(lapserate_feedback_xod, replace=True, n_samples=20, random_state=NS)
#       boot = np.random.choice(lapserate_feedback_xod, size=10)
       lapserate_feedback_xod_resample[NS,:] = boot.copy()
       boot = resample(lapserate_feedback_xoo, replace=True, n_samples=20, random_state=NS)
#       boot = np.random.choice(lapserate_feedback_xoo, size=10)
       lapserate_feedback_xoo_resample[NS,:] = boot.copy()

       boot = resample(q_feedback_ctr, replace=True, n_samples=20, random_state=NS)
       q_feedback_ctr_resample[NS,:] = boot.copy()
       boot = resample(q_feedback_xod, replace=True, n_samples=20, random_state=NS)
       q_feedback_xod_resample[NS,:] = boot.copy()
       boot = resample(q_feedback_xoo, replace=True, n_samples=20, random_state=NS)
       q_feedback_xoo_resample[NS,:] = boot.copy()

       boot = resample(alb_feedback_ctr, replace=True, n_samples=20, random_state=NS)
       alb_feedback_ctr_resample[NS,:] = boot.copy()
       boot = resample(alb_feedback_xod, replace=True, n_samples=20, random_state=NS)
       alb_feedback_xod_resample[NS,:] = boot.copy()
       boot = resample(alb_feedback_xoo, replace=True, n_samples=20, random_state=NS)
       alb_feedback_xoo_resample[NS,:] = boot.copy()

       boot = resample((lw_cloud_feedback_ctr+sw_cloud_feedback_ctr), replace=True, n_samples=20, random_state=NS)
       cloud_feedback_ctr_resample[NS,:] = boot.copy()
       boot = resample((lw_cloud_feedback_xod+sw_cloud_feedback_xod), replace=True, n_samples=20, random_state=NS)
       cloud_feedback_xod_resample[NS,:] = boot.copy()
       boot = resample((lw_cloud_feedback_xoo+sw_cloud_feedback_xoo), replace=True, n_samples=20, random_state=NS)
       cloud_feedback_xoo_resample[NS,:] = boot.copy()

       boot = resample((dh_ctr-(FLNT_arctic_ctr+FSNT_arctic_ctr)), replace=True, n_samples=20, random_state=NS)
#       boot = np.random.choice((dh_ctr-(FLNT_arctic_ctr+FSNT_arctic_ctr)), size=10)
       df_ctr_resample[NS,:] = boot.copy()
       boot = resample((dh_xod-(FLNT_arctic_xod+FSNT_arctic_xod)), replace=True, n_samples=20, random_state=NS)
#       boot = np.random.choice((dh_xod-(FLNT_arctic_xod+FSNT_arctic_xod)), size=10)
       df_xod_resample[NS,:] = boot.copy()
       boot = resample((dh_xoo-(FLNT_arctic_xoo+FSNT_arctic_xoo)), replace=True, n_samples=20, random_state=NS)
#       boot = np.random.choice((dh_xoo-(FLNT_arctic_xoo+FSNT_arctic_xoo)), size=10)
       df_xoo_resample[NS,:] = boot.copy()

       boot = resample(dts_arctic_ctr, replace=True, n_samples=20, random_state=NS)
#       boot = np.random.choice(dts_arctic_ctr, size=10)
       dts_ctr_resample[NS,:] = boot.copy()
       boot = resample(dts_arctic_xod, replace=True, n_samples=20, random_state=NS)
#       boot = np.random.choice(dts_arctic_xod, size=10)
       dts_xod_resample[NS,:] = boot.copy()
       boot = resample(dts_arctic_xoo, replace=True, n_samples=20, random_state=NS)
#       boot = np.random.choice(dts_arctic_xoo, size=10)
       dts_xoo_resample[NS,:] = boot.copy()

# ================================================================
# plot figures
# ================================================================
if True:

#   dts1 = (dts_arctic_ctr-dts_arctic_xoo).mean()
#   dts2 = (dts_arctic_ctr-dts_arctic_xod).mean()
   dts1 = 1.
   dts2 = 1.

#   factor1 = 2.64
#   factor1 = 3.2
   factor1 = 1.0

   plt.close('all')
   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)
   ax1 = fig.add_axes([0.25, 0.3, 0.55, 0.45])
   plt.plot([-20,20],[-20,20],'--', color='darkgray')
   plt.plot([-20,20],[0,0],'-',color='darkgray')
   plt.plot([0,0],[-20,20],'-',color='darkgray')

   ts_test1=(np.nanmean((df_ctr_resample-df_xoo_resample)/(dts_ctr_resample-dts_xoo_resample),axis=1))/factor1
   ts_test2=(np.nanmean((df_ctr_resample-df_xod_resample)/(dts_ctr_resample-dts_xod_resample),axis=1))/factor1
   plt.plot(ts_test1,ts_test2,'.', alpha=0.05, color='darkgray')

   ts_test1=(np.nanmean((lapserate_feedback_ctr_resample-lapserate_feedback_xoo_resample)/(dts_ctr_resample-dts_xoo_resample),axis=1))/factor1
   ts_test2=(np.nanmean((lapserate_feedback_ctr_resample-lapserate_feedback_xod_resample)/(dts_ctr_resample-dts_xod_resample),axis=1))/factor1
   plt.plot(ts_test1,ts_test2,'m.', alpha=0.05)

   ts_test1=(np.nanmean((alb_feedback_ctr_resample-alb_feedback_xoo_resample)/(dts_ctr_resample-dts_xoo_resample),axis=1))/factor1
   ts_test2=(np.nanmean((alb_feedback_ctr_resample-alb_feedback_xod_resample)/(dts_ctr_resample-dts_xod_resample),axis=1))/factor1
   plt.plot(ts_test1,ts_test2,'c.', alpha=0.05)

   ts_test1=(np.nanmean((q_feedback_ctr_resample-q_feedback_xoo_resample)/(dts_ctr_resample-dts_xoo_resample),axis=1))/factor1
   ts_test2=(np.nanmean((q_feedback_ctr_resample-q_feedback_xod_resample)/(dts_ctr_resample-dts_xod_resample),axis=1))/factor1
   plt.plot(ts_test1,ts_test2,'.', alpha=0.05, color='gold')

   ts_test1=(np.nanmean((cloud_feedback_ctr_resample-cloud_feedback_xoo_resample)/(dts_ctr_resample-dts_xoo_resample),axis=1))/factor1
   ts_test2=(np.nanmean((cloud_feedback_ctr_resample-cloud_feedback_xod_resample)/(dts_ctr_resample-dts_xod_resample),axis=1))/factor1
   plt.plot(ts_test1,ts_test2,'.', alpha=0.05, color='lime')

   ts_test1 = p_lr_ctr_xoo/factor1*dts1
   ts_test2 = p_lr_ctr_xod/factor1*dts2
   plt.plot(ts_test1,ts_test2,'ro',markersize=4)
#   plt.plot(ts_test1.mean(),ts_test2.mean(),'ro',label='lapserate (p-value=' + str(round(pvalue_lr_xod_xoo,2)) +')', markersize=10)
   plt.plot(ts_test1.mean(),ts_test2.mean(),'ro',label='lapse-rate', markersize=10)

   ts_test1 = p_al_ctr_xoo/factor1*dts1
   ts_test2 = p_al_ctr_xod/factor1*dts2
   plt.plot(ts_test1,ts_test2,'bo',markersize=4)
#   plt.plot(ts_test1.mean(),ts_test2.mean(),'bo',label='albedo (p-value=' + str(round(pvalue_al_xod_xoo,2)) +')', markersize=10)
   plt.plot(ts_test1.mean(),ts_test2.mean(),'bo',label='albedo', markersize=10)

   ts_test1 = p_wv_ctr_xoo/factor1*dts1
   ts_test2 = p_wv_ctr_xod/factor1*dts2
   plt.plot(ts_test1,ts_test2,'o',markersize=4, color='orange')
#   plt.plot(ts_test1.mean(),ts_test2.mean(),'o',label='water vapor (p-avlue=' + str(round(pvalue_wv_xod_xoo,2)) +')', markersize=10, color='orange')
   plt.plot(ts_test1.mean(),ts_test2.mean(),'o',label='water vapor', markersize=10, color='orange')

   ts_test1 = p_cl_ctr_xoo/factor1*dts1
   ts_test2 = p_cl_ctr_xod/factor1*dts2
   plt.plot(ts_test1,ts_test2,'go',markersize=4)
#   plt.plot(ts_test1.mean(),ts_test2.mean(),'go',label='cloud (p-value=' + str(round(pvalue_cl_xod_xoo,2)) +')', markersize=10)
   plt.plot(ts_test1.mean(),ts_test2.mean(),'go',label='cloud', markersize=10)

   ts_test1 = p_df_ctr_xoo/factor1*dts1
   ts_test2 = p_df_ctr_xod/factor1*dts2
   plt.plot(ts_test1,ts_test2,'ko',markersize=4)
#   plt.plot(ts_test1.mean(),ts_test2.mean(),'ko',label='transport (p-value='+ str(round(pvalue_df_xod_xoo,2)) +')', markersize=10)
   plt.plot(ts_test1.mean(),ts_test2.mean(),'ko',label='transport', markersize=10)

   plt.legend(loc='lower right', fontsize='small')
   plt.ylabel('ODS (W/m$^2$/K)')
   plt.xlabel('CO2 (W/m$^2$/K)')
  
#   plt.yticks([-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2])
#   plt.xticks([-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

   plt.axis([-0.3, 2.3, -1.0, 3.3])
#   plt.axis([-0.1, 0.1, -0.1, 0.1])
   plt.title('Feedback parameters')

   plt.savefig('trend_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()






