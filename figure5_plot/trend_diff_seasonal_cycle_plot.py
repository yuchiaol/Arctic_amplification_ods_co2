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
    [xxx, pvalue] = stats.ttest_ind(exp1_var,exp2_var)
    ttest_map = np.nan
    pvalue_map = pvalue
    if pvalue < sig_level:
       ttest_map = 1.
       pvalue_map = pvalue

    return ttest_map, pvalue_map

if flag_run == 1:
# ================================================================
# read simulations
# ================================================================
   year_N = 51
#   dirname = '/data1/yliang/cesm1_le/data_storage_center/'
   dirname = '/home/yliang/research/co2_ods/data_process/'
   filename = 'sat_sie_seasonal_cycle_temp_output.nc'
#   filename = 'seasonal_cycle_temp_output.nc'
   f = Dataset(dirname + filename, 'r')
   global_sat_ctr_change = f.variables['global_sat_ctr_trend'][:,:].data/year_N*10
   global_sat_xc2_change = f.variables['global_sat_xc2_trend'][:,:].data/year_N*10
   global_sat_xod_change = f.variables['global_sat_xod_trend'][:,:].data/year_N*10
   global_sat_xco_change = f.variables['global_sat_xco_trend'][:,:].data/year_N*10
   global_sat_xoo_change = f.variables['global_sat_xoo_trend'][:,:].data/year_N*10

   arctic_sat_ctr_change = f.variables['arctic_sat_ctr_trend'][:,:].data/year_N*10
   arctic_sat_xc2_change = f.variables['arctic_sat_xc2_trend'][:,:].data/year_N*10
   arctic_sat_xod_change = f.variables['arctic_sat_xod_trend'][:,:].data/year_N*10
   arctic_sat_xco_change = f.variables['arctic_sat_xco_trend'][:,:].data/year_N*10
   arctic_sat_xoo_change = f.variables['arctic_sat_xoo_trend'][:,:].data/year_N*10

   arctic_sie_ctr_change = f.variables['arctic_sie_ctr_trend'][:,:].data/year_N*10
   arctic_sie_xc2_change = f.variables['arctic_sie_xc2_trend'][:,:].data/year_N*10
   arctic_sie_xod_change = f.variables['arctic_sie_xod_trend'][:,:].data/year_N*10
   arctic_sie_xco_change = f.variables['arctic_sie_xco_trend'][:,:].data/year_N*10
   arctic_sie_xoo_change = f.variables['arctic_sie_xoo_trend'][:,:].data/year_N*10

   f.close()

   aa_ctr_change = arctic_sat_ctr_change/global_sat_ctr_change
   aa_xc2_change = (arctic_sat_ctr_change-arctic_sat_xc2_change)/(global_sat_ctr_change-global_sat_xc2_change)
   aa_xod_change = (arctic_sat_ctr_change-arctic_sat_xod_change)/(global_sat_ctr_change-global_sat_xod_change)

# ================================================================
# perform statistical test
# ================================================================
   sig_level = 0.05

   sig_arctic_ctr_xc2 = np.zeros((12))
   sig_arctic_ctr_xod = np.zeros((12))
   sig_arctic_ctr_xco = np.zeros((12))
   sig_arctic_xod_xoo = np.zeros((12))
   for NM in range(12):
       [sig_arctic_ctr_xc2[NM], xxx] = perform_ttest_1d_here(arctic_sat_ctr_change[NM,:],arctic_sat_xc2_change[NM,:],sig_level)
       [sig_arctic_ctr_xod[NM], xxx] = perform_ttest_1d_here(arctic_sat_ctr_change[NM,:],arctic_sat_xod_change[NM,:],sig_level)
       [sig_arctic_ctr_xco[NM], xxx] = perform_ttest_1d_here(arctic_sat_ctr_change[NM,:],arctic_sat_xco_change[NM,:],sig_level)
       [sig_arctic_xod_xoo[NM], xxx] = perform_ttest_1d_here(arctic_sat_xod_change[NM,:],arctic_sat_xoo_change[NM,:],sig_level)

   sig_sie_arctic_ctr_xc2 = np.zeros((12))
   sig_sie_arctic_ctr_xod = np.zeros((12))
   sig_sie_arctic_ctr_xco = np.zeros((12))
   sig_sie_arctic_xod_xoo = np.zeros((12))
   for NM in range(12):
       [sig_sie_arctic_ctr_xc2[NM], xxx] = perform_ttest_1d_here(arctic_sie_ctr_change[NM,:],arctic_sie_xc2_change[NM,:],sig_level)
       [sig_sie_arctic_ctr_xod[NM], xxx] = perform_ttest_1d_here(arctic_sie_ctr_change[NM,:],arctic_sie_xod_change[NM,:],sig_level)
       [sig_sie_arctic_ctr_xco[NM], xxx] = perform_ttest_1d_here(arctic_sie_ctr_change[NM,:],arctic_sie_xco_change[NM,:],sig_level)
       [sig_sie_arctic_xod_xoo[NM], xxx] = perform_ttest_1d_here(arctic_sie_xod_change[NM,:],arctic_sie_xoo_change[NM,:],sig_level)

   sig_aa_ctr_xc2 = np.zeros((12))
   sig_aa_ctr_xod = np.zeros((12))
   for NM in range(12):
       [sig_aa_ctr_xc2[NM], xxx] = perform_ttest_1d_here(aa_ctr_change[NM,:],aa_xc2_change[NM,:],sig_level)
       [sig_aa_ctr_xod[NM], xxx] = perform_ttest_1d_here(aa_ctr_change[NM,:],aa_xod_change[NM,:],sig_level)

   sig_arctic_xc2_xod = np.zeros((12))
   for NM in range(12):
       [sig_arctic_xc2_xod[NM], xxx] = perform_ttest_1d_here(arctic_sat_ctr_change[NM,:]-arctic_sat_xc2_change[NM,:],arctic_sat_ctr_change[NM,:]-arctic_sat_xod_change[NM,:],sig_level)

   sig_sie_arctic_xc2_xod = np.zeros((12))
   for NM in range(12):
       [sig_sie_arctic_xc2_xod[NM], xxx] = perform_ttest_1d_here(arctic_sie_ctr_change[NM,:]-arctic_sie_xc2_change[NM,:],arctic_sie_ctr_change[NM,:]-arctic_sie_xod_change[NM,:],sig_level)

   sig_aa_xc2_xod = np.zeros((12))
   for NM in range(12):
       [sig_aa_xc2_xod[NM], xxx] = perform_ttest_1d_here(aa_xc2_change[NM,:],aa_xod_change[NM,:],sig_level)

   sig_level = 0.10

   sig90_arctic_xc2_xod = np.zeros((12))
   for NM in range(12): 
       [sig90_arctic_xc2_xod[NM], xxx] = perform_ttest_1d_here(arctic_sat_ctr_change[NM,:]-arctic_sat_xc2_change[NM,:],arctic_sat_ctr_change[NM,:]-arctic_sat_xod_change[NM,:],sig_level)
   
   sig90_sie_arctic_xc2_xod = np.zeros((12))
   for NM in range(12):
       [sig90_sie_arctic_xc2_xod[NM], xxx] = perform_ttest_1d_here(arctic_sie_ctr_change[NM,:]-arctic_sie_xc2_change[NM,:],arctic_sie_ctr_change[NM,:]-arctic_sie_xod_change[NM,:],sig_level)

   sig90_aa_xc2_xod = np.zeros((12))
   for NM in range(12):
       [sig90_aa_xc2_xod[NM], xxx] = perform_ttest_1d_here(aa_xc2_change[NM,:],aa_xod_change[NM,:],sig_level)

   sig90_arctic_xod_xoo = np.zeros((12))
   for NM in range(12):
       [sig90_arctic_xod_xoo[NM], xxx] = perform_ttest_1d_here(arctic_sat_xoo_change[NM,:],arctic_sat_xod_change[NM,:],sig_level)

   sig90_sie_arctic_xod_xoo = np.zeros((12))
   for NM in range(12):
       [sig90_sie_arctic_xod_xoo[NM], xxx] = perform_ttest_1d_here(arctic_sie_xoo_change[NM,:],arctic_sie_xod_change[NM,:],sig_level)


# ================================================================
# plot figures
# ================================================================
if True:

   tt = np.linspace(1,12,12)   

   plt.close('all')

   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)

# plot dots
   ax1 = fig.add_axes([0.2, 0.6, 0.65, 0.25])
   y_test = np.nanmean(arctic_sat_ctr_change, axis=1) - np.nanmean(arctic_sat_xco_change, axis=1)
   plt.plot(tt,y_test,'g-',label='CO2&ODS')
   y_test = np.nanmean(arctic_sat_ctr_change, axis=1)*2 - np.nanmean(arctic_sat_xc2_change, axis=1) - np.nanmean(arctic_sat_xod_change, axis=1)
   plt.plot(tt,y_test,'-',color='peru',label='CO2+ODS')
   y_test = np.nanmean(arctic_sat_ctr_change, axis=1)
   plt.plot(tt,y_test,'k-',label='ALL')
   y_test = np.nanmean(arctic_sat_ctr_change, axis=1) - np.nanmean(arctic_sat_xc2_change, axis=1)
   plt.plot(tt,y_test,'r-',label='CO2')
   plt.plot(tt,y_test*sig_arctic_xc2_xod,'rx', fillstyle='none', markersize=10)
   plt.plot(tt,y_test*sig90_arctic_xc2_xod,'ro',fillstyle='none', markersize=10)
   y_test = np.nanmean(arctic_sat_ctr_change, axis=1) - np.nanmean(arctic_sat_xod_change, axis=1)
   plt.plot(tt,y_test,'b-',label='ODS')
   plt.plot(tt,y_test*sig_arctic_xc2_xod,'bx', fillstyle='none', markersize=10)
   plt.plot(tt,y_test*sig90_arctic_xc2_xod,'bo',fillstyle='none', markersize=10)
   y_test = np.nanmean(arctic_sat_xod_change, axis=1) - np.nanmean(arctic_sat_xoo_change, axis=1)
   plt.plot(tt,y_test,'y-',label='Strat-O$_3$')
   plt.plot(tt,y_test*sig_arctic_xod_xoo,'yx', fillstyle='none', markersize=10)
   plt.plot(tt,y_test*sig90_arctic_xod_xoo,'yo', fillstyle='none', markersize=10)
   plt.title('(a) Arctic SAT trend', fontsize=11)
   plt.xticks(tt,['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
   plt.ylabel('K per decade', fontsize=11)
#   plt.axis([1, 12, -0., 3.0])
   plt.grid('on')
   plt.legend(ncol=3, loc='upper left', fontsize='small')

   factor = 1./(1.e12)
   ax1 = fig.add_axes([0.2, 0.28, 0.65, 0.25])
   y_test = np.nanmean(arctic_sie_ctr_change, axis=1) - np.nanmean(arctic_sie_xco_change, axis=1)
   plt.plot(tt,y_test*factor,'g-',label='CO2&ODS')
   y_test = np.nanmean(arctic_sie_ctr_change, axis=1)*2 - np.nanmean(arctic_sie_xc2_change, axis=1) - np.nanmean(arctic_sie_xod_change, axis=1)
   plt.plot(tt,y_test*factor,'-',color='peru',label='CO2+ODS')
   y_test = np.nanmean(arctic_sie_ctr_change, axis=1)
   plt.plot(tt,y_test*factor,'k-',label='ALL')
   y_test = np.nanmean(arctic_sie_ctr_change, axis=1) - np.nanmean(arctic_sie_xc2_change, axis=1)
   plt.plot(tt,y_test*factor,'r-',label='CO2')
   plt.plot(tt,y_test*sig_sie_arctic_xc2_xod*factor,'rx', fillstyle='none', markersize=10)
   plt.plot(tt,y_test*sig90_sie_arctic_xc2_xod*factor,'ro',fillstyle='none', markersize=10)
   y_test = np.nanmean(arctic_sie_ctr_change, axis=1) - np.nanmean(arctic_sie_xod_change, axis=1)
   plt.plot(tt,y_test*factor,'b-',label='ODS')
   plt.plot(tt,y_test*sig_sie_arctic_xc2_xod*factor,'bx', fillstyle='none', markersize=10)
   plt.plot(tt,y_test*sig90_sie_arctic_xc2_xod*factor,'bo',fillstyle='none', markersize=10)
   y_test = np.nanmean(arctic_sie_xod_change, axis=1) - np.nanmean(arctic_sie_xoo_change, axis=1)
   plt.plot(tt,y_test*factor,'y-',label='Strat-O$_3$')
   plt.plot(tt,y_test*sig_sie_arctic_xod_xoo*factor,'yx', fillstyle='none', markersize=10)
   plt.plot(tt,y_test*sig90_sie_arctic_xod_xoo*factor,'yo', fillstyle='none', markersize=10)
   plt.title('(b) Arctic SIE trend', fontsize=11)
   plt.xticks(tt,['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
   plt.ylabel('10$^{6}$ km$^{2}$ per decade', fontsize=11)
#   plt.axis([1, 12, -0.5, 0.2])
   plt.grid('on')
   plt.legend(ncol=3, loc='lower left', fontsize='small')

   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)

   ax1 = fig.add_axes([0.1, 0.06, 0.65, 0.25])
   y_test = np.nanmean(aa_ctr_change, axis=1)
   plt.plot(tt,y_test,'k-',label='ALL')
   y_test = np.nanmean(aa_xc2_change, axis=1)
   plt.plot(tt,y_test,'r-',label='CO2')
   plt.plot(tt,y_test*sig90_aa_xc2_xod,'ro', fillstyle='none', markersize=10)
#   plt.plot(tt,y_test*sig99_aa_xc2_xod,'rx', markersize=10)
   y_test = np.nanmean(aa_xod_change, axis=1)
   plt.plot(tt,y_test,'b-',label='ODS')
   plt.plot(tt,y_test*sig90_aa_xc2_xod,'bo', fillstyle='none', markersize=10)
#   plt.plot(tt,y_test*sig99_aa_xc2_xod,'bx', markersize=10)
   plt.title('(c) Arctic amplfication factor', fontsize=11)
   plt.xticks(tt,['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
   plt.ylabel('K per decade', fontsize=11)
#   plt.axis([1, 12, -0., 3.0])
   plt.grid('on')
   plt.legend(ncol=3, loc='upper left')


   plt.savefig('trend_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()

   sys.exit()








