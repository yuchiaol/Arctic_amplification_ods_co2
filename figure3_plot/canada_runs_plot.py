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
    [xxx, pvalue] = stats.ttest_ind(exp1_var,exp2_var)
    ttest_map = np.nan
    pvalue_map = pvalue
    if pvalue < sig_level:
       ttest_map = 1.
       pvalue_map = pvalue

    return ttest_map, pvalue_map

if flag_run == 1:
# ================================================================
# read CanESM5 data
# ================================================================
   ens_num = 20
   arctic_sat_ctr_trend = np.array([2.0266, 0.6527, 1.2231, 1.7442, 2.9844, 2.5921, 0.9912, 1.6062, 2.5874, 2.0722,\
                                    2.9177, 2.3914, 2.0990, 0.9025, 1.5422, 2.2367, 1.9477, 2.2574, 2.3108, 2.8007])/50*10
   global_sat_ctr_trend = np.array([1.1943, 0.9205, 0.9660, 0.9980, 1.1082, 0.9275, 0.9535, 0.8277, 1.1679, 1.0100,\
                                    1.2972, 0.9705, 0.9347, 0.6950, 1.0550, 0.8969, 1.1077, 1.0353, 1.1508, 1.1315])/50*10
   arctic_sat_xod_trend = np.array([1.3665, 0.1396, 0.5783, 1.4017, 1.4807, 0.7484, 1.5574, 1.8920, 1.8295, 1.1142,\
                                    1.5203, 0.2380, 0.4549, 0.9125, 0.5506, 0.8990, 1.0917, 0.9453, 1.2509, 1.3659])/50*10
   global_sat_xod_trend = np.array([0.6533, 0.4439, 0.6154, 0.5674, 0.5611, 0.6489, 0.6785, 0.7275, 0.7749, 0.6316,\
                                    0.7577, 0.3437, 0.3898, 0.5056, 0.6500, 0.4715, 0.8070, 0.5820, 0.7960, 0.6048])/50*10
   arctic_sat_xc2_trend = np.array([-0.7067, 0.1145, -0.5708, -0.9218, 0.5333, 1.4356, 0.5418, 0.6706, 0.4306, -1.1672,\
                                    0.9490, -0.0032, 0.7608, -0.0068, -0.3880, 0.3453, 0.0082, 0.3606, 0.8306, 0.2162])/50*10
   global_sat_xc2_trend = np.array([0.2480, 0.3793, 0.3519, -0.0217, 0.2590, 0.4540, 0.4461, 0.3663, 0.3520, 0.0719,\
                                    0.5379, 0.2803, 0.3530, 0.2097, 0.4420, 0.3486, 0.3789, 0.3278, 0.6039, 0.3240])/50*10

# ================================================================
# resample aa factor
# ================================================================
   aa_sat_ctr_trend = (arctic_sat_ctr_trend)/(global_sat_ctr_trend)
   aa_sat_xod_trend = (arctic_sat_ctr_trend-arctic_sat_xod_trend)/(global_sat_ctr_trend-global_sat_xod_trend)
   aa_sat_xc2_trend = (arctic_sat_ctr_trend-arctic_sat_xc2_trend)/(global_sat_ctr_trend-global_sat_xc2_trend)

   nsample = 10000
   sat_ctr_change_resample_arctic = np.zeros((nsample,ens_num))
   sat_xc2_change_resample_arctic = np.zeros((nsample,ens_num))
   sat_xod_change_resample_arctic = np.zeros((nsample,ens_num))
   sat_ctr_change_resample_global = np.zeros((nsample,ens_num))
   sat_xc2_change_resample_global = np.zeros((nsample,ens_num))
   sat_xod_change_resample_global = np.zeros((nsample,ens_num))

   for NS in range(nsample):
       print(NS)

       boot = resample(arctic_sat_ctr_trend, replace=True, n_samples=20, random_state=NS)
       sat_ctr_change_resample_arctic[NS,:] = boot.copy()
       boot = resample(arctic_sat_xc2_trend, replace=True, n_samples=20, random_state=NS)
       sat_xc2_change_resample_arctic[NS,:] = boot.copy()
       boot = resample(arctic_sat_xod_trend, replace=True, n_samples=20, random_state=NS)
       sat_xod_change_resample_arctic[NS,:] = boot.copy()

       boot = resample(global_sat_ctr_trend, replace=True, n_samples=20, random_state=NS)
       sat_ctr_change_resample_global[NS,:] = boot.copy()
       boot = resample(global_sat_xc2_trend, replace=True, n_samples=20, random_state=NS)
       sat_xc2_change_resample_global[NS,:] = boot.copy()
       boot = resample(global_sat_xod_trend, replace=True, n_samples=20, random_state=NS)
       sat_xod_change_resample_global[NS,:] = boot.copy()

# ================================================================
# plot figures
# ================================================================
if True:
#   aa_sat_ctr_trend = (arctic_sat_ctr_trend)
#   aa_sat_xod_trend = (arctic_sat_ctr_trend-arctic_sat_xod_trend)
#   aa_sat_xc2_trend = (arctic_sat_ctr_trend-arctic_sat_xc2_trend)

   plt.close('all')
   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)

# plot whisker
   ax1 = fig.add_axes([0.58, 0.3, 0.4, 0.37])
   plt.plot([1]*10000,np.nanmean(sat_ctr_change_resample_arctic/sat_ctr_change_resample_global, axis=1),'o',color='darkgray',alpha=0.05)
   plt.plot([2]*10000,np.nanmean((sat_ctr_change_resample_arctic-sat_xc2_change_resample_arctic)/(sat_ctr_change_resample_global-sat_xc2_change_resample_global), axis=1),'o',color='m',alpha=0.05)   
   plt.plot([3]*10000,np.nanmean((sat_ctr_change_resample_arctic-sat_xod_change_resample_arctic)/(sat_ctr_change_resample_global-sat_xod_change_resample_global), axis=1),'o',color='c',alpha=0.05)
   flierprops = dict(marker='o', markerfacecolor='k', markersize=3,linestyle='none', markeredgecolor='k')
   aa_sat_ctr_trend = arctic_sat_ctr_trend/(global_sat_ctr_trend)
   aa_sat_xc2_trend = (arctic_sat_ctr_trend-arctic_sat_xc2_trend)/(global_sat_ctr_trend-global_sat_xc2_trend)
   aa_sat_xod_trend = (arctic_sat_ctr_trend-arctic_sat_xod_trend)/(global_sat_ctr_trend-global_sat_xod_trend)
   data = np.concatenate((aa_sat_ctr_trend,aa_sat_xc2_trend,aa_sat_xod_trend))
   data = data.reshape((int(len(data)/ens_num),ens_num)).T
   plt.boxplot(data, flierprops=flierprops, whis='range')
   for II in range(ens_num):
       plt.plot([1]*2,[aa_sat_ctr_trend[II]]*2,'ko',markersize=3)
       plt.plot([2]*2,[aa_sat_xc2_trend[II]]*2,'ro',markersize=3)
       plt.plot([3]*2,[aa_sat_xod_trend[II]]*2,'bo',markersize=3)
   plt.plot([1]*2,[aa_sat_ctr_trend.mean()]*2,'ko',markersize=10, label='ALL (' + str(round(np.nanmean(arctic_sat_ctr_trend/(global_sat_ctr_trend)), 2))+ ')')
   aa = np.nanmean((arctic_sat_ctr_trend-arctic_sat_xc2_trend)/(global_sat_ctr_trend-global_sat_xc2_trend))
   plt.plot([2]*2,[aa_sat_xc2_trend.mean()]*2,'ro',markersize=10, label='CO2 (' + str(round(aa, 2))  + ')')
   aa = np.nanmean((arctic_sat_ctr_trend-arctic_sat_xod_trend)/(global_sat_ctr_trend-global_sat_xod_trend))
   plt.plot([3]*2,[aa_sat_xod_trend.mean()]*2,'bo',markersize=10, label='ODS (' + str(round(aa, 2))  + ')')
#   plt.plot([1]*2,[aa_sat_obs_trend.mean()]*2,'ro')
   plt.title('(b) Arctic amplification factor', fontsize=12)
   plt.xticks([1,2,3],['ALL','CO2','ODS'], rotation=0)
   plt.ylabel('', fontsize=10)
   plt.grid(axis='y')
   plt.legend()

# plot scatter
   ax1 = fig.add_axes([0.08, 0.3, 0.45, 0.37])
   plt.plot([-5,5],[-5,5],'--', color='darkgray')
   plt.plot(np.nanmean(sat_ctr_change_resample_global, axis=1),np.nanmean(sat_ctr_change_resample_arctic, axis=1), '.', color='darkgray', alpha=0.05)
   plt.plot(np.nanmean(sat_ctr_change_resample_global-sat_xc2_change_resample_global, axis=1),np.nanmean(sat_ctr_change_resample_arctic-sat_xc2_change_resample_arctic, axis=1), 'm.', alpha=0.05)
   plt.plot(np.nanmean(sat_ctr_change_resample_global-sat_xod_change_resample_global, axis=1),np.nanmean(sat_ctr_change_resample_arctic-sat_xod_change_resample_arctic, axis=1), 'c.', alpha=0.05)

   ts_x = np.nanmean(sat_ctr_change_resample_global, axis=1)
   ts_y = np.nanmean(sat_ctr_change_resample_arctic, axis=1)
   p_ctr_res = stats.theilslopes(ts_y,ts_x)
   p_ctr = np.polyfit(ts_x,ts_y,1)
   text1 = 'ALL, slope=' + str(round(p_ctr_res[0],2)) + '$\pm$' + str(round(abs(p_ctr_res[0]-p_ctr_res[2]), 2))
   plt.plot(np.nanmean(global_sat_ctr_trend),np.nanmean(arctic_sat_ctr_trend),'ko',label=text1, markersize=14)
   plt.plot(global_sat_ctr_trend,arctic_sat_ctr_trend,'ko', markersize=6)
   plt.plot([-5,5],[p_ctr[0]*-5+p_ctr[1],p_ctr[0]*5+p_ctr[1]],'k-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_ctr_res[1] + p_ctr_res[2]*np.array([-5,5]), 'k--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_ctr_res[1] + p_ctr_res[3]*np.array([-5,5]), 'k--', linewidth=0.5)

   ts_x = np.nanmean(sat_ctr_change_resample_global-sat_xc2_change_resample_global, axis=1)
   ts_y = np.nanmean(sat_ctr_change_resample_arctic-sat_xc2_change_resample_arctic, axis=1)
   p_xc2_res = stats.theilslopes(ts_y,ts_x)
   p_xc2 = np.polyfit(ts_x,ts_y,1)
   text2 = 'CO2, slope=' + str(round(p_xc2_res[0],2)) + '$\pm$' + str(round(abs(p_xc2_res[0]-p_xc2_res[2]), 2))
   plt.plot(np.nanmean(global_sat_ctr_trend-global_sat_xc2_trend),np.nanmean(arctic_sat_ctr_trend-arctic_sat_xc2_trend),'ro',label=text2, markersize=14)
   plt.plot(global_sat_ctr_trend-global_sat_xc2_trend,arctic_sat_ctr_trend-arctic_sat_xc2_trend,'ro', markersize=6)
#   p_xc2 = np.polyfit(global_sat_ctr_trend-global_sat_xc2_trend,arctic_sat_ctr_trend-arctic_sat_xc2_trend,1)
   plt.plot([-5,5],[p_xc2[0]*-5+p_xc2[1],p_xc2[0]*5+p_xc2[1]],'r-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xc2_res[1] + p_xc2_res[2]*np.array([-5,5]), 'r--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xc2_res[1] + p_xc2_res[3]*np.array([-5,5]), 'r--', linewidth=0.5)


   ts_x = np.nanmean(sat_ctr_change_resample_global-sat_xod_change_resample_global, axis=1)
   ts_y = np.nanmean(sat_ctr_change_resample_arctic-sat_xod_change_resample_arctic, axis=1)
#   ts_x = global_sat_ctr_trend-global_sat_xod_trend
#   ts_y = arctic_sat_ctr_trend-arctic_sat_xod_trend
   p_xod_res = stats.theilslopes(ts_y,ts_x)
   p_xod = np.polyfit(ts_x,ts_y,1)
   text3 = 'ODS, slope=' + str(round(p_xod_res[0],2)) + '$\pm$' + str(round(abs(p_xod_res[0]-p_xod_res[2]), 2))
   plt.plot(np.nanmean(global_sat_ctr_trend-global_sat_xod_trend),np.nanmean(arctic_sat_ctr_trend-arctic_sat_xod_trend),'bo',label=text3, markersize=14)
   plt.plot(global_sat_ctr_trend-global_sat_xod_trend,arctic_sat_ctr_trend-arctic_sat_xod_trend,'bo', markersize=6)
#   p_xod = np.polyfit(global_sat_ctr_trend-global_sat_xod_trend,arctic_sat_ctr_trend-arctic_sat_xod_trend,1)
   plt.plot([-5,5],[p_xod[0]*-5+p_xod[1],p_xod[0]*5+p_xod[1]],'b-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xod_res[1] + p_xod_res[2]*np.array([-5,5]), 'b--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xod_res[1] + p_xod_res[3]*np.array([-5,5]), 'b--', linewidth=0.5)


   plt.axis([0.0,0.30,-0.2,0.80])
   plt.ylabel('Arctic SAT trend (K per decade)', fontsize=12)
   plt.xlabel('Global SAT trend (K per decade)', fontsize=12)
   plt.title('(a) Arctic vs global warming trends', fontsize=12)
   plt.legend(loc='upper left', fontsize='small')

   plt.savefig('trend_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()



