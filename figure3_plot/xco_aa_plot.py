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
# read sat
   varname = 'TREFHT'
   year1 = 1955
   year2 = 2005
   year_N = int(year2 - year1 + 1)
   tt = np.linspace(year1,year2,year_N)

# read grid basics
   dirname = '/data1/yliang/cesm1_le/data_storage_center/'
   filename = varname + '_annual_mean_temp_output.nc'
   f = Dataset(dirname + filename, 'r')
   lon = f.variables['lon'][:].data
   lat = f.variables['lat'][:].data
   area = data_process_f.area_calculate_nonuniform(lon,lat)
   sat_ctr = f.variables['var_ctr'][:,:,:,:].data
   sat_xc2 = f.variables['var_xc2'][:,:,:,:].data
   sat_xod = f.variables['var_xod'][:,:,:,:].data
   sat_xco = f.variables['var_xco'][:,:,:,:].data
   f.close()

   ens_num = sat_ctr.shape[0]
   ny = len(lat)
   nx = len(lon)

# read observation
# sat
   dirname = '/data1/yliang/gistempv4/'
   filename = 'gistemp1200_GHCNv4_ERSSTv5.nc'
   f = Dataset(dirname + filename, 'r')
   lat_obs = f.variables['lat'][:].data
   lon_obs = f.variables['lon'][:].data+180
   area_obs = data_process_f.area_calculate_nonuniform(lon_obs,lat_obs)
   year_ref = 1880
   t0 = int((year1-year_ref)*12)
   t1 = int((year2-year_ref)*12+12)
   sat_tmp = f.variables['tempanomaly'][t0:t1,:,:].data
   f.close()
   sat_tmp[abs(sat_tmp)>30000] = np.nan
   mask_tmp = (sat_tmp/sat_tmp).copy()

   ny_obs = len(lat_obs)
   nx_obs = len(lon_obs)

   sat_obs = np.nanmean(sat_tmp.reshape((year_N, 12, ny_obs, nx_obs)), axis=1)
   mask_tmp = (sat_obs/sat_obs).copy()
   mask_obs = np.nanmean(mask_tmp, axis=0)

# ================================================================
# calculate linear trends
# ================================================================
# simulated sat
   sat_ctr_trend = np.zeros((ens_num,ny,nx))
   sat_xc2_trend = np.zeros((ens_num,ny,nx))
   sat_xod_trend = np.zeros((ens_num,ny,nx))
   sat_xco_trend = np.zeros((ens_num,ny,nx))
   for NENS in range(ens_num):
       print('sat:' + str(NENS))
       for JJ in range(ny):
           for II in range(nx):
               ts_test = (sat_ctr[NENS,:,JJ,II]).copy()
               [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_test)
               sat_ctr_trend[NENS,JJ,II] = slope.copy()*10
               ts_test = (sat_xc2[NENS,:,JJ,II]).copy()
               [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_test)
               sat_xc2_trend[NENS,JJ,II] = slope.copy()*10
               ts_test = (sat_xod[NENS,:,JJ,II]).copy()
               [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_test)
               sat_xod_trend[NENS,JJ,II] = slope.copy()*10
               ts_test = (sat_xco[NENS,:,JJ,II]).copy()
               [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_test)
               sat_xco_trend[NENS,JJ,II] = slope.copy()*10

# observational sat
   sat_obs_trend = np.zeros((ny_obs,nx_obs)) 
   for JJ in range(ny_obs):
       for II in range(nx_obs):
           ts_test = (sat_obs[:,JJ,II]).copy()
           [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_test)
           sat_obs_trend[JJ,II] = slope.copy()*year_N

# ================================================================
# take area-averaged value
# ================================================================
# global domain
   global_sat_ctr_trend = np.zeros((ens_num))
   global_sat_xc2_trend = np.zeros((ens_num))
   global_sat_xod_trend = np.zeros((ens_num))
   global_sat_xco_trend = np.zeros((ens_num))
   for NENS in range(ens_num):
       global_sat_ctr_trend[NENS] = np.nansum(sat_ctr_trend[NENS,:,:]*area)/np.nansum(area)
       global_sat_xc2_trend[NENS] = np.nansum(sat_xc2_trend[NENS,:,:]*area)/np.nansum(area)
       global_sat_xod_trend[NENS] = np.nansum(sat_xod_trend[NENS,:,:]*area)/np.nansum(area)
       global_sat_xco_trend[NENS] = np.nansum(sat_xco_trend[NENS,:,:]*area)/np.nansum(area)

   global_sat_obs_trend = np.nansum(sat_obs_trend*area_obs)/np.nansum(area_obs)

# arctic domain 60N-90N
   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(0,360,60,89,lon,lat)
   print("arctic region selected lon:" + str(lon[x1]) + '-' + str(lon[x2]) + ' lat:' + str(lat[y1]) + '-' + str(lat[y2]))
   arctic_sat_ctr_trend = np.zeros((ens_num))
   arctic_sat_xc2_trend = np.zeros((ens_num))
   arctic_sat_xod_trend = np.zeros((ens_num))
   arctic_sat_xco_trend = np.zeros((ens_num))
   for NENS in range(ens_num):
       arctic_sat_ctr_trend[NENS] = np.nansum(sat_ctr_trend[NENS,y1:,:]*area[y1:,:])/np.nansum(area[y1:,:])
       arctic_sat_xc2_trend[NENS] = np.nansum(sat_xc2_trend[NENS,y1:,:]*area[y1:,:])/np.nansum(area[y1:,:])
       arctic_sat_xod_trend[NENS] = np.nansum(sat_xod_trend[NENS,y1:,:]*area[y1:,:])/np.nansum(area[y1:,:])
       arctic_sat_xco_trend[NENS] = np.nansum(sat_xco_trend[NENS,y1:,:]*area[y1:,:])/np.nansum(area[y1:,:])

   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(0,360,60,88,lon_obs,lat_obs)
   print("arctic region selected lon:" + str(lon_obs[x1]) + '-' + str(lon_obs[x2]) + ' lat:' + str(lat_obs[y1]) + '-' + str(lat_obs[y2]))
   arctic_sat_obs_trend = np.nansum(sat_obs_trend[y1:,:]*area_obs[y1:,:]*mask_obs[y1:,:])/np.nansum(area_obs[y1:,:]*mask_obs[y1:,:])

   aa_sat_obs_trend = arctic_sat_obs_trend/global_sat_obs_trend

# ================================================================
# resample aa factor
# ================================================================
   aa_sat_ctr_trend = (arctic_sat_ctr_trend)/np.nanmean(global_sat_ctr_trend)
   aa_sat_xod_trend = (arctic_sat_ctr_trend-arctic_sat_xod_trend)/np.nanmean(global_sat_ctr_trend-global_sat_xod_trend)
   aa_sat_xc2_trend = (arctic_sat_ctr_trend-arctic_sat_xc2_trend)/np.nanmean(global_sat_ctr_trend-global_sat_xc2_trend)
   aa_sat_xco_trend = (arctic_sat_ctr_trend-arctic_sat_xco_trend)/np.nanmean(global_sat_ctr_trend-global_sat_xco_trend)

   nsample = 10000
   sat_ctr_change_resample_arctic = np.zeros((nsample,ens_num*2))
   sat_xc2_change_resample_arctic = np.zeros((nsample,ens_num*2))
   sat_xod_change_resample_arctic = np.zeros((nsample,ens_num*2))
   sat_xco_change_resample_arctic = np.zeros((nsample,ens_num*2))
   sat_ctr_change_resample_global = np.zeros((nsample,ens_num*2))
   sat_xc2_change_resample_global = np.zeros((nsample,ens_num*2))
   sat_xod_change_resample_global = np.zeros((nsample,ens_num*2))
   sat_xco_change_resample_global = np.zeros((nsample,ens_num*2))

   for NS in range(nsample):
       print(NS)

       boot = resample(arctic_sat_ctr_trend, replace=True, n_samples=20, random_state=NS)
       sat_ctr_change_resample_arctic[NS,:] = boot.copy()
       boot = resample(arctic_sat_xc2_trend, replace=True, n_samples=20, random_state=NS)
       sat_xc2_change_resample_arctic[NS,:] = boot.copy()
       boot = resample(arctic_sat_xod_trend, replace=True, n_samples=20, random_state=NS)
       sat_xod_change_resample_arctic[NS,:] = boot.copy()
       boot = resample(arctic_sat_xco_trend, replace=True, n_samples=20, random_state=NS)
       sat_xco_change_resample_arctic[NS,:] = boot.copy()

       boot = resample(global_sat_ctr_trend, replace=True, n_samples=20, random_state=NS)
       sat_ctr_change_resample_global[NS,:] = boot.copy()
       boot = resample(global_sat_xc2_trend, replace=True, n_samples=20, random_state=NS)
       sat_xc2_change_resample_global[NS,:] = boot.copy()
       boot = resample(global_sat_xod_trend, replace=True, n_samples=20, random_state=NS)
       sat_xod_change_resample_global[NS,:] = boot.copy()
       boot = resample(global_sat_xco_trend, replace=True, n_samples=20, random_state=NS)
       sat_xco_change_resample_global[NS,:] = boot.copy()

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
   plt.plot([4]*10000,np.nanmean((sat_ctr_change_resample_arctic-sat_xco_change_resample_arctic)/(sat_ctr_change_resample_global-sat_xco_change_resample_global), axis=1),'o',color='lawngreen',alpha=0.05)
   flierprops = dict(marker='o', markerfacecolor='k', markersize=3,linestyle='none', markeredgecolor='k')
   aa_sat_ctr_trend = arctic_sat_ctr_trend/(global_sat_ctr_trend)
   aa_sat_xc2_trend = (arctic_sat_ctr_trend-arctic_sat_xc2_trend)/(global_sat_ctr_trend-global_sat_xc2_trend)
   aa_sat_xod_trend = (arctic_sat_ctr_trend-arctic_sat_xod_trend)/(global_sat_ctr_trend-global_sat_xod_trend)
   aa_sat_xco_trend = (arctic_sat_ctr_trend-arctic_sat_xco_trend)/(global_sat_ctr_trend-global_sat_xco_trend)
   data = np.concatenate((aa_sat_ctr_trend,aa_sat_xc2_trend,aa_sat_xod_trend,aa_sat_xco_trend))
   data = data.reshape((int(len(data)/ens_num),ens_num)).T
   plt.boxplot(data, flierprops=flierprops, whis='range')
   for II in range(ens_num):
       plt.plot([1]*2,[aa_sat_ctr_trend[II]]*2,'ko',markersize=3)
       plt.plot([2]*2,[aa_sat_xc2_trend[II]]*2,'ro',markersize=3)
       plt.plot([3]*2,[aa_sat_xod_trend[II]]*2,'bo',markersize=3)
       plt.plot([4]*2,[aa_sat_xco_trend[II]]*2,'go',markersize=3)
   plt.plot([1]*2,[aa_sat_ctr_trend.mean()]*2,'ko',markersize=10, label='ALL (' + str(round(np.nanmean(arctic_sat_ctr_trend/(global_sat_ctr_trend)), 2))+ ')')
   aa = np.nanmean((arctic_sat_ctr_trend-arctic_sat_xc2_trend)/(global_sat_ctr_trend-global_sat_xc2_trend))
   plt.plot([2]*2,[aa_sat_xc2_trend.mean()]*2,'ro',markersize=10, label='CO2 (' + str(round(aa, 2))  + ')')
   aa = np.nanmean((arctic_sat_ctr_trend-arctic_sat_xod_trend)/(global_sat_ctr_trend-global_sat_xod_trend))
   plt.plot([3]*2,[aa_sat_xod_trend.mean()]*2,'bo',markersize=10, label='ODS (' + str(round(aa, 2))  + ')')
   aa = np.nanmean((arctic_sat_ctr_trend-arctic_sat_xco_trend)/(global_sat_ctr_trend-global_sat_xco_trend))
   plt.plot([4]*2,[aa_sat_xco_trend.mean()]*2,'go',markersize=10, label='CO2&ODS (' + str(round(aa, 2))  + ')')
#   plt.plot([1]*2,[aa_sat_obs_trend.mean()]*2,'ro')
   plt.title('(b) Arctic amplification factor', fontsize=12)
   plt.xticks([1,2,3,4],['ALL','CO2','ODS','CO2&ODS'], rotation=0)
   plt.ylabel('', fontsize=10)
   plt.grid(axis='y')
   plt.legend()

# plot scatter
   ax1 = fig.add_axes([0.08, 0.3, 0.45, 0.37])
   plt.plot([-5,5],[-5,5],'--', color='darkgray')
   plt.plot(np.nanmean(sat_ctr_change_resample_global, axis=1),np.nanmean(sat_ctr_change_resample_arctic, axis=1), '.', color='darkgray', alpha=0.05)
   plt.plot(np.nanmean(sat_ctr_change_resample_global-sat_xc2_change_resample_global, axis=1),np.nanmean(sat_ctr_change_resample_arctic-sat_xc2_change_resample_arctic, axis=1), 'm.', alpha=0.05)
   plt.plot(np.nanmean(sat_ctr_change_resample_global-sat_xod_change_resample_global, axis=1),np.nanmean(sat_ctr_change_resample_arctic-sat_xod_change_resample_arctic, axis=1), 'c.', alpha=0.05)
   plt.plot(np.nanmean(sat_ctr_change_resample_global-sat_xco_change_resample_global, axis=1),np.nanmean(sat_ctr_change_resample_arctic-sat_xco_change_resample_arctic, axis=1), '.', color='lawngreen', alpha=0.05)

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

   ts_x = np.nanmean(sat_ctr_change_resample_global-sat_xco_change_resample_global, axis=1)
   ts_y = np.nanmean(sat_ctr_change_resample_arctic-sat_xco_change_resample_arctic, axis=1)
   p_xco_res = stats.theilslopes(ts_y,ts_x)
   p_xco = np.polyfit(ts_x,ts_y,1)
   text3 = 'CO2&ODS, slope=' + str(round(p_xco_res[0],2)) + '$\pm$' + str(round(abs(p_xco_res[0]-p_xco_res[2]), 2))
   plt.plot(np.nanmean(global_sat_ctr_trend-global_sat_xco_trend),np.nanmean(arctic_sat_ctr_trend-arctic_sat_xco_trend),'go',label=text3, markersize=14)
   plt.plot(global_sat_ctr_trend-global_sat_xco_trend,arctic_sat_ctr_trend-arctic_sat_xco_trend,'go', markersize=6)
#   p_xod = np.polyfit(global_sat_ctr_trend-global_sat_xod_trend,arctic_sat_ctr_trend-arctic_sat_xod_trend,1)
   plt.plot([-5,5],[p_xco[0]*-5+p_xco[1],p_xco[0]*5+p_xco[1]],'g-', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xco_res[1] + p_xco_res[2]*np.array([-5,5]), 'g--', linewidth=0.5)
   plt.plot(np.array([-5,5]), p_xco_res[1] + p_xco_res[3]*np.array([-5,5]), 'g--', linewidth=0.5)

   plt.axis([0,0.18,-0.02,0.5])
   plt.ylabel('Arcctic SAT trend (K per decade)', fontsize=12)
   plt.xlabel('Global SAT trend (K per decade)', fontsize=12)
   plt.title('(a) Arctic vs global warming trends', fontsize=12)
   plt.legend(loc='upper left', fontsize='small')

   plt.savefig('trend_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()



