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
   sat_xoo = f.variables['var_xoo'][:,:,:,:].data
   f.close()

   ens_num = sat_ctr.shape[0]
   ny = len(lat)
   nx = len(lon)

# read sic
   varname = 'SIE'
   year1 = 1955
   year2 = 2005
   year_N = int(year2 - year1 + 1)
# read grid basics
   dirname = '/data1/yliang/cesm1_le/data_storage_center/'
#   filename = varname + '_annual_mean_temp_output.nc'
   filename = varname + '_month9_temp_output.nc'
   f = Dataset(dirname + filename, 'r')
   sie_ctr = f.variables['var_ctr'][:,:,:,:].data
   sie_xc2 = f.variables['var_xc2'][:,:,:,:].data
   sie_xod = f.variables['var_xod'][:,:,:,:].data
   sie_xco = f.variables['var_xco'][:,:,:,:].data
   sie_xoo = f.variables['var_xoo'][:,:,:,:].data
   f.close()

# ================================================================
# read observation
# ================================================================
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

#   sat_obs_annual_mean = np.nanmean(sat_tmp.reshape((year_N,12,ny_obs,nx_obs)), axis=1)
   sat_obs_global_mean = np.zeros((year_N*12))
   for NT in range(year_N*12):
       sat_obs_global_mean[NT] = np.nansum(sat_tmp[NT,:,:]*mask_tmp[NT,:,:]*area_obs)/np.nansum(area_obs*mask_tmp[NT,:,:])

   sat_obs_global_annual_mean = np.nanmean(sat_obs_global_mean.reshape((year_N,12)), axis=1)

   ts_test = sat_obs_global_annual_mean[:].copy()
#   global_sat_obs_change = ts_test[-10:].mean() - ts_test[0:10].mean()
   [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_test)
   global_sat_obs_change = slope*year_N

   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(0,360,60,88,lon_obs,lat_obs)
   print("arctic region selected lon:" + str(lon_obs[x1]) + '-' + str(lon_obs[x2]) + ' lat:' + str(lat_obs[y1]) + '-' + str(lat_obs[y2]))
   sat_obs_arctic_mean = np.zeros((year_N*12))
   for NT in range(year_N*12):
       sat_obs_arctic_mean[NT] = np.nansum(sat_tmp[NT,y1:,:]*mask_tmp[NT,y1:,:]*area_obs[y1:,:])/np.nansum(area_obs[y1:,:]*mask_tmp[NT,y1:,:])

   sat_obs_arctic_annual_mean = np.nanmean(sat_obs_arctic_mean.reshape((year_N,12)), axis=1)

   ts_test = sat_obs_arctic_annual_mean[:].copy()
#   arctic_sat_obs_change = ts_test[-10:].mean() - ts_test[0:10].mean()
   [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_test)
   arctic_sat_obs_change = slope*year_N

# SIC
# read nsidc sic
   dirname = '/data1/yliang/nsidc/'
   filename = 'G10010_sibt1850_v2.0.nc'
   f = Dataset(dirname + filename, 'r')
   lat_sic_obs = f.variables['latitude'][::-1].data
   lon_sic_obs = f.variables['longitude'][:].data
   area_sic_obs = data_process_f.area_calculate_nonuniform(lon_sic_obs,lat_sic_obs)
   Gridcell_Area = f.variables['Gridcell_Area'][::-1].data
   year_ref = 1850
   t0 = int((year1-year_ref)*12)
   t1 = int((year2-year_ref)*12+12)
   sic_tmp = f.variables['seaice_conc'][t0:t1,:,:].data*1.
   f.close()
   sic_tmp[sic_tmp==120] = np.nan
   sic_tmp = sic_tmp[:,::-1,:].copy()/100.

# sea-ice extent
   sie_tmp = (sic_tmp/sic_tmp).copy()
   sie_tmp[sic_tmp<0.15] = np.nan

   ny_sic_obs = len(lat_sic_obs)
   nx_sic_obs = len(lon_sic_obs)

   area_default = area_sic_obs.copy()
   for II in range(nx_sic_obs):
       area_default[:,II] = Gridcell_Area[:].copy()*1000000

# sea-ice extent
   sie_tmp = (sic_tmp).copy()
   sie_tmp[sic_tmp<0.15] = 0.
   sie_tmp[sic_tmp>=0.15] = 1.

   ny_sic_obs = len(lat_sic_obs)
   nx_sic_obs = len(lon_sic_obs)

#   sic_obs = sic_tmp.reshape((year_N, 12, ny_sic_obs, nx_sic_obs))[:,8,:,:].copy()
   sie_obs = sie_tmp.reshape((year_N, 12, ny_sic_obs, nx_sic_obs))[:,8,:,:].copy()

   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(0,360,10,88,lon_sic_obs,lat_sic_obs)
   print("arctic region selected lon:" + str(lon_sic_obs[x1]) + '-' + str(lon_sic_obs[x2]) + ' lat:' + str(lat_sic_obs[y1]) + '-' + str(lat_sic_obs[y2]))
   arctic_sie_obs = np.zeros((year_N))
   for NT in range(year_N):
       arctic_sie_obs[NT] = np.nansum(sie_obs[NT,y1:,:]*area_default[y1:,:])

#   arctic_sie_obs_change = np.nanmean(arctic_sie_obs[-10:]) - np.nanmean(arctic_sie_obs[0:10])
   ts_test = arctic_sie_obs.copy()
   [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_test)
   arctic_sie_obs_change = slope*year_N

# read hadi sic
   dirname = '/data1/yliang/hadisst/'
   filename = 'HadISST.2.2.0.0_sea_ice_concentration.nc'
   f = Dataset(dirname + filename, 'r')
   lat_sic_obs = f.variables['latitude'][::-1].data
   lon_sic_obs = f.variables['longitude'][:].data
   area_sic_obs = data_process_f.area_calculate_nonuniform(lon_sic_obs,lat_sic_obs)
   year_ref = 1850
   t0 = int((year1-year_ref)*12)
   t1 = int((year2-year_ref)*12+12)
   sic_tmp = f.variables['sic'][t0:t1,:,:].data
   f.close()
   sic_tmp[abs(sic_tmp)>30000] = np.nan
   sic_tmp = sic_tmp[:,::-1,:].copy()
# sea-ice extent
   sie_tmp = (sic_tmp/sic_tmp).copy()
   sie_tmp[sic_tmp<0.15] = np.nan

   ny_sic_obs = len(lat_sic_obs)
   nx_sic_obs = len(lon_sic_obs)

   sie_obs = sie_tmp.reshape((year_N, 12, ny_sic_obs, nx_sic_obs))[:,8,:,:].copy()

# calculate area average
   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(0,360,30,88,lon_sic_obs,lat_sic_obs)
   print("arctic region selected lon:" + str(lon_sic_obs[x1]) + '-' + str(lon_sic_obs[x2]) + ' lat:' + str(lat_sic_obs[y1]) + '-' + str(lat_sic_obs[y2]))

   ts_sie = np.zeros((year_N))
   for NT in range(year_N):
       ts_sie[NT] = np.nansum(sie_obs[NT,y1:,:]*area_sic_obs[y1:,:])

#   arctic_sie_hadi_change = np.nanmean(ts_sie[-10:]) - np.nanmean(ts_sie[0:10])
   ts_test = ts_sie.copy()
   [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_test)
   arctic_sie_hadi_change = slope*year_N

# ================================================================
# take area-averaged value
# ================================================================
   arctic_sat_ctr = np.zeros((ens_num,year_N))   
   arctic_sat_xc2 = np.zeros((ens_num,year_N))
   arctic_sat_xod = np.zeros((ens_num,year_N))
   arctic_sat_xco = np.zeros((ens_num,year_N))
   arctic_sat_xoo = np.zeros((ens_num,year_N))
   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(0,360,60,89,lon,lat)
   print("arctic region selected lon:" + str(lon[x1]) + '-' + str(lon[x2]) + ' lat:' + str(lat[y1]) + '-' + str(lat[y2]))
   for NENS in range(ens_num):
       print('sat:' + str(NENS))
       for NT in range(year_N):
           arctic_sat_ctr[NENS,NT] = np.nansum(sat_ctr[NENS,NT,y1:,:]*area[y1:,:])/np.nansum(area[y1:,:])
           arctic_sat_xc2[NENS,NT] = np.nansum(sat_xc2[NENS,NT,y1:,:]*area[y1:,:])/np.nansum(area[y1:,:])
           arctic_sat_xod[NENS,NT] = np.nansum(sat_xod[NENS,NT,y1:,:]*area[y1:,:])/np.nansum(area[y1:,:])
#           arctic_sat_xco[NENS,NT] = np.nansum(sat_xco[NENS,NT,y1:,:]*area[y1:,:])/np.nansum(area[y1:,:])
           arctic_sat_xoo[NENS,NT] = np.nansum((sat_xod[NENS,NT,y1:,:]-sat_xoo[NENS,NT,y1:,:])*area[y1:,:])/np.nansum(area[y1:,:])

   arctic_sie_ctr = np.zeros((ens_num,year_N))             
   arctic_sie_xc2 = np.zeros((ens_num,year_N))
   arctic_sie_xod = np.zeros((ens_num,year_N))
   arctic_sie_xco = np.zeros((ens_num,year_N))
   arctic_sie_xoo = np.zeros((ens_num,year_N))
   [x1, x2, y1, y2] = data_process_f.find_lon_lat_index(0,360,10,88,lon,lat)
   print("arctic region selected lon:" + str(lon[x1]) + '-' + str(lon[x2]) + ' lat:' + str(lat[y1]) + '-' + str(lat[y2]))
   for NENS in range(ens_num):
       print('sat:' + str(NENS))
       for NT in range(year_N):
           arctic_sie_ctr[NENS,NT] = np.nansum(sie_ctr[NENS,NT,y1:,:]*area[y1:,:])
           arctic_sie_xc2[NENS,NT] = np.nansum(sie_xc2[NENS,NT,y1:,:]*area[y1:,:])
           arctic_sie_xod[NENS,NT] = np.nansum(sie_xod[NENS,NT,y1:,:]*area[y1:,:])
#           arctic_sie_xco[NENS,NT] = np.nansum(sie_xco[NENS,NT,y1:,:]*area[y1:,:])
           arctic_sie_xoo[NENS,NT] = np.nansum((sie_xod[NENS,NT,y1:,:]-sie_xoo[NENS,NT,y1:,:])*area[y1:,:])

# ================================================================
# calculate change and trend defined in polvani et al. (2020) ncc
# ================================================================
# sat
   arctic_sat_ctr_change = np.zeros((ens_num))
   arctic_sat_xc2_change = np.zeros((ens_num))
   arctic_sat_xod_change = np.zeros((ens_num))
   arctic_sat_xco_change = np.zeros((ens_num))
   arctic_sat_xoo_change = np.zeros((ens_num))
   for NENS in range(ens_num):
       [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,arctic_sat_ctr[NENS,:])
       arctic_sat_ctr_change[NENS] = slope*10.
       [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,arctic_sat_xc2[NENS,:])
       arctic_sat_xc2_change[NENS] = slope*10.
       [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,arctic_sat_xod[NENS,:])
       arctic_sat_xod_change[NENS] = slope*10.
       [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,arctic_sat_xoo[NENS,:])
       arctic_sat_xoo_change[NENS] = slope*10.

# sic
   arctic_sie_ctr_change = np.zeros((ens_num))
   arctic_sie_xc2_change = np.zeros((ens_num))
   arctic_sie_xod_change = np.zeros((ens_num))
   arctic_sie_xco_change = np.zeros((ens_num))
   arctic_sie_xoo_change = np.zeros((ens_num))
   for NENS in range(ens_num):
       [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,arctic_sie_ctr[NENS,:])
       arctic_sie_ctr_change[NENS] = slope*10.
       [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,arctic_sie_xc2[NENS,:])
       arctic_sie_xc2_change[NENS] = slope*10.
       [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,arctic_sie_xod[NENS,:])
       arctic_sie_xod_change[NENS] = slope*10.
       [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,arctic_sie_xoo[NENS,:])
       arctic_sie_xoo_change[NENS] = slope*10.

# ================================================================
# perform resampling 
# ================================================================
   nsample = 10000
   global_sat_ctr_change_resample = np.zeros((nsample,ens_num*2))
   global_sat_xc2_change_resample = np.zeros((nsample,ens_num*2))
   global_sat_xod_change_resample = np.zeros((nsample,ens_num*2))
   global_sat_xco_change_resample = np.zeros((nsample,ens_num*2))
   global_sat_xoo_change_resample = np.zeros((nsample,ens_num*2))
   arctic_sat_ctr_change_resample = np.zeros((nsample,ens_num*2))
   arctic_sat_xc2_change_resample = np.zeros((nsample,ens_num*2))
   arctic_sat_xod_change_resample = np.zeros((nsample,ens_num*2))
   arctic_sat_xco_change_resample = np.zeros((nsample,ens_num*2))
   arctic_sat_xoo_change_resample = np.zeros((nsample,ens_num*2))
   arctic_sie_ctr_change_resample = np.zeros((nsample,ens_num*2))
   arctic_sie_xc2_change_resample = np.zeros((nsample,ens_num*2))
   arctic_sie_xod_change_resample = np.zeros((nsample,ens_num*2))
   arctic_sie_xco_change_resample = np.zeros((nsample,ens_num*2))
   arctic_sie_xoo_change_resample = np.zeros((nsample,ens_num*2))

   for NS in range(nsample):
       print(NS)

       boot = resample(arctic_sat_ctr_change, replace=True, n_samples=20, random_state=NS)
       arctic_sat_ctr_change_resample[NS,:] = boot.copy()
       boot = resample(arctic_sat_xc2_change, replace=True, n_samples=20, random_state=NS)
       arctic_sat_xc2_change_resample[NS,:] = boot.copy()
       boot = resample(arctic_sat_xod_change, replace=True, n_samples=20, random_state=NS)
       arctic_sat_xod_change_resample[NS,:] = boot.copy()
       boot = resample(arctic_sat_xoo_change, replace=True, n_samples=20, random_state=NS)
       arctic_sat_xoo_change_resample[NS,:] = boot.copy()
#       boot = resample(arctic_sat_xco_change, replace=True, n_samples=10, random_state=NS)
#       arctic_sat_xco_change_resample[NS,:] = boot.copy()

       boot = resample(arctic_sie_ctr_change, replace=True, n_samples=20, random_state=NS)
       arctic_sie_ctr_change_resample[NS,:] = boot.copy()
       boot = resample(arctic_sie_xc2_change, replace=True, n_samples=20, random_state=NS)
       arctic_sie_xc2_change_resample[NS,:] = boot.copy()
       boot = resample(arctic_sie_xod_change, replace=True, n_samples=20, random_state=NS)
       arctic_sie_xod_change_resample[NS,:] = boot.copy()
#       boot = resample(arctic_sie_xco_change, replace=True, n_samples=10, random_state=NS)
#       arctic_sie_xco_change_resample[NS,:] = boot.copy()
       boot = resample(arctic_sie_xoo_change, replace=True, n_samples=20, random_state=NS)
       arctic_sie_xoo_change_resample[NS,:] = boot.copy()

# ================================================================
# plot figures
# ================================================================
if True:
 
   flierprops = dict(marker='o', markerfacecolor='k', markersize=3,
                  linestyle='none', markeredgecolor='k')
  
   plt.close('all')

   fig = plt.figure()
   fig.set_size_inches(10, 10, forward=True)

   ax1 = fig.add_axes([0.08, 0.57, 0.3, 0.27])
#   xx = np.nanmean(arctic_sat_ctr_change_resample, axis=1)
#   plt.plot([1]*10000,xx,'.', color='darkgray', alpha=0.05)
#   xx = np.nanmean(arctic_sat_ctr_change_resample, axis=1) - np.nanmean(arctic_sat_xc2_change_resample, axis=1)
#   plt.plot([2]*10000,xx,'.', color='m', alpha=0.05)
#   xx = np.nanmean(arctic_sat_ctr_change_resample, axis=1) - np.nanmean(arctic_sat_xod_change_resample, axis=1)
#   plt.plot([3]*10000,xx,'.', color='c', alpha=0.05)
#   xx = np.nanmean(arctic_sat_xoo_change_resample, axis=1)
#   plt.plot([4]*10000,xx,'.', color='gold', alpha=0.05)
   plt.plot([0.5,4.5],[0,0],'k-')
   data = np.concatenate((arctic_sat_ctr_change,arctic_sat_ctr_change-arctic_sat_xc2_change,arctic_sat_ctr_change-arctic_sat_xod_change,arctic_sat_xoo_change))
   data = data.reshape((int(len(data)/ens_num),ens_num)).T
   plt.boxplot(data, flierprops=flierprops, whis='range')
   for II in range(ens_num):
       plt.plot([1]*2,[arctic_sat_ctr_change[II]]*2,'ko',markersize=3)
       plt.plot([3]*2,[arctic_sat_ctr_change[II]-arctic_sat_xod_change[II]]*2,'bo',markersize=3)
       plt.plot([2]*2,[arctic_sat_ctr_change[II]-arctic_sat_xc2_change[II]]*2,'ro',markersize=3)
       plt.plot([4]*2,[arctic_sat_xoo_change[II]]*2,'yo',markersize=3)
   plt.plot([1]*2,[arctic_sat_ctr_change.mean()]*2,'ko',markersize=10)
   plt.plot([3]*2,[arctic_sat_ctr_change.mean()-arctic_sat_xod_change.mean()]*2,'bo',markersize=10)
   plt.plot([2]*2,[arctic_sat_ctr_change.mean()-arctic_sat_xc2_change.mean()]*2,'ro',markersize=10)
   plt.plot([4]*2,[arctic_sat_xoo_change.mean()]*2,'yo',markersize=10)
   plt.title('(a) Annual Arctic SAT trend', fontsize=12)
   plt.xticks([1,2,3,4],['ALL','CO2','ODS','Strat-O$_3$'],rotation=0)
   plt.ylabel('K per decade', fontsize=10)
   plt.grid(axis='y')
#   plt.axis([0.5, 3.5, -0.6, 1.8])
#   plt.xlim(0.5,5.5)

   ax1 = fig.add_axes([0.445, 0.57, 0.54, 0.27])
   xx = np.nanmean(arctic_sat_ctr_change_resample, axis=1)
   sns.distplot(xx, hist=True, kde=True, color='k', label='ALL') 
   xx = np.nanmean(arctic_sat_ctr_change_resample, axis=1) - np.nanmean(arctic_sat_xod_change_resample, axis=1)
   sns.distplot(xx, hist=True, kde=True, color='c', label='ODS')
   xx = np.nanmean(arctic_sat_ctr_change_resample, axis=1) - np.nanmean(arctic_sat_xc2_change_resample, axis=1)
   sns.distplot(xx, hist=True, kde=True, color='m', label='CO2')
   xx = np.nanmean(arctic_sat_xoo_change_resample, axis=1)
   sns.distplot(xx, hist=True, kde=True, color='y', label='Strat-O$_3$')
   plt.xlabel('K per decade')
   plt.title('(b) Annual Arctic SAT trend distribution', fontsize=12)  
   plt.legend()

   factor = 1./1000000000000
   ax1 = fig.add_axes([0.08, 0.22, 0.3, 0.27])
#   xx = np.nanmean(arctic_sie_ctr_change_resample, axis=1)
#   plt.plot([1]*10000,xx*factor,'.', color='darkgray', alpha=0.05)
#   xx = np.nanmean(arctic_sie_ctr_change_resample, axis=1) - np.nanmean(arctic_sie_xc2_change_resample, axis=1)
#   plt.plot([2]*10000,xx*factor,'.', color='m', alpha=0.05)
#   xx = np.nanmean(arctic_sie_ctr_change_resample, axis=1) - np.nanmean(arctic_sie_xod_change_resample, axis=1)
#   plt.plot([3]*10000,xx*factor,'.', color='c', alpha=0.05)
#   xx = np.nanmean(arctic_sie_xoo_change_resample, axis=1)
#   plt.plot([4]*10000,xx*factor,'.', color='gold', alpha=0.05)
   plt.plot([0.5,4.5],[0,0],'k-')
   data = np.concatenate((arctic_sie_ctr_change*factor,arctic_sie_ctr_change*factor-arctic_sie_xc2_change*factor,arctic_sie_ctr_change*factor-arctic_sie_xod_change*factor,arctic_sie_xoo_change*factor))
   data = data.reshape((int(len(data)/ens_num),ens_num)).T
   plt.boxplot(data, flierprops=flierprops, whis='range')
   for II in range(ens_num):
       plt.plot([1]*2,[arctic_sie_ctr_change[II]*factor]*2,'ko',markersize=3)
       plt.plot([3]*2,[arctic_sie_ctr_change[II]*factor-arctic_sie_xod_change[II]*factor]*2,'bo',markersize=3)
       plt.plot([2]*2,[arctic_sie_ctr_change[II]*factor-arctic_sie_xc2_change[II]*factor]*2,'ro',markersize=3)
       plt.plot([4]*2,[arctic_sie_xoo_change[II]*factor]*2,'yo',markersize=3)
#       plt.plot([4]*2,[arctic_sie_ctr_change[II]*factor-arctic_sie_xco_change[II]*factor]*2,'ko',markersize=3)
   plt.plot([1]*2,[arctic_sie_ctr_change.mean()*factor]*2,'ko',markersize=10)
   plt.plot([3]*2,[arctic_sie_ctr_change.mean()*factor-arctic_sie_xod_change.mean()*factor]*2,'bo',markersize=10)
   plt.plot([2]*2,[arctic_sie_ctr_change.mean()*factor-arctic_sie_xc2_change.mean()*factor]*2,'ro',markersize=10)
   plt.plot([4]*2,[arctic_sie_xoo_change.mean()*factor]*2,'yo',markersize=10)
#   plt.plot([4]*2,[arctic_sie_ctr_change.mean()*factor-arctic_sie_xco_change.mean()*factor]*2,'co')
   plt.title('(c) September Arctic SIE trend', fontsize=12)
   plt.xticks([1,2,3,4],['ALL','CO2','ODS','Strat-O$_3$'],rotation=0)
   plt.ylabel('10$^{6}$ km$^{2}$ per decade', fontsize=10)
   plt.grid(axis='y')
#   plt.axis([0.5, 3.5, -2.5, 0.5])
   plt.gca().invert_yaxis()

   factor = 1./1000000000000
   ax1 = fig.add_axes([0.445, 0.22, 0.54, 0.27])
#   xx = np.nanmean(arctic_sie_ctr_change_resample, axis=1)*factor - np.nanmean(arctic_sie_xco_change_resample, axis=1)*factor
#   sns.distplot(xx, hist=True, kde=True, color='c', label='Hist-xCO')
#   plt.hist(xx, bins=25, color='c',alpha=0.3, rwidth=0.6, label='Hist-xCO')
#   xx = np.nanmean(arctic_sie_ctr_change_resample, axis=1)*factor - np.nanmean(arctic_sie_xoo_change_resample, axis=1)*factor
#   plt.hist(xx, bins=25, color='y',alpha=0.7, rwidth=0.8, label='Hist-xOO')
   xx = np.nanmean(arctic_sie_ctr_change_resample, axis=1)*factor
   sns.distplot(xx, hist=True, kde=True, color='k', label='ALL')
#   plt.hist(xx, bins=25, color='r',alpha=0.3, rwidth=0.6, label='Hist')
   xx = np.nanmean(arctic_sie_ctr_change_resample, axis=1)*factor - np.nanmean(arctic_sie_xod_change_resample, axis=1)*factor
   sns.distplot(xx, hist=True, kde=True, color='c', label='ODS')
#   plt.hist(xx, bins=25, color='b',alpha=0.3, rwidth=0.6, label='ODS')
   xx = np.nanmean(arctic_sie_ctr_change_resample, axis=1)*factor - np.nanmean(arctic_sie_xc2_change_resample, axis=1)*factor
   sns.distplot(xx, hist=True, kde=True, color='m', label='CO2')
   xx = np.nanmean(arctic_sie_xoo_change_resample, axis=1)*factor
   sns.distplot(xx, hist=True, kde=True, color='y', label='Strat-O$_3$')
#   plt.hist(xx, bins=25, color='g',alpha=0.3, rwidth=0.6, label='CO2')
   plt.title('(d) September Arctic SIE trend distribution', fontsize=12)
   plt.xlabel('km$^{2}$ (x10$^{6}$) per decade')
   plt.gca().invert_xaxis()
   plt.legend()

   plt.savefig('trend_tmp_plot.jpg', format='jpeg', dpi=200)

   plt.show()






