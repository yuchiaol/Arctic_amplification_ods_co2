flag_run = 1
# ================================================================
# Yu-Chiao @ NTU
# ERL replies to review comments
# ================================================================
#exec(open('waf_check1.py').read())
#rsync -avzhe ssh —progress

# ================================================================
# Import functions
# ================================================================
import argparse
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from math import isnan, radians
#from mpl_toolkits.basemap import Basemap
from IPython import get_ipython
import sys, os
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from pylab import setp, genfromtxt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata
from datetime import datetime
from sklearn import linear_model
from scipy.io import loadmat
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.img_tiles as cimgt
from cartopy.io.img_tiles import StamenTerrain
import matplotlib.path as mpath
import math

sys.path.append('/work/home/yuchiaol/lib/python_functions/data_process/')
import data_process_f

def plot_box(ax,lon1,lon2,lat1,lat2, color_txt):

    ax.plot(np.linspace(lon1,lon1,100), np.linspace(lat1, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)
    ax.plot(np.linspace(lon2,lon2,100), np.linspace(lat1, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)
    ax.plot(np.linspace(lon1,lon2,100), np.linspace(lat1, lat1, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)
    ax.plot(np.linspace(lon1,lon2,100), np.linspace(lat2, lat2, 100), transform=ccrs.PlateCarree(), color=color_txt, linewidth=1)

def area_calculate_nonuniform(lon,lat):
# https://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2004/msg00023.html

    Re = 6.3781e6
    area = np.zeros((len(lat),len(lon)))*np.nan
    for II in range(len(lon)-1):
       for JJ in range(len(lat)-1):
           area[JJ+1,II+1] = 2*np.pi*Re**2*(abs(np.sin(lat[JJ+1]*np.pi/180.)-np.sin(lat[JJ]*np.pi/180.)))*abs((lon[II+1]-lon[II]))/360.
 
    dlat = abs(lat[1] - lat[0])
    for II in range(len(lon)-1):
        area[0,II+1] = 2*np.pi*Re**2*(abs(np.sin(lat[0]*np.pi/180.)-np.sin((lat[0]-dlat)*np.pi/180.)))*abs((lon[II+1]-lon[II]))/360.

    area[:,0] = area[:,-1].copy()

    return area

def extend_longitude(X,lat,lon):

    X_out = np.zeros((len(lat), len(lon)+1))
    lon_out = np.zeros((len(lon)+1))
    dx = lon[2] - lon[1]

    X_out[:,0:-1] = X[:,:].copy()
    X_out[:,-1] = X[:,-1].copy()
    lon_out[0:-1] = lon[:].copy()
    lon_out[-1] = lon[-1] + dx

    return X_out, lon_out

if flag_run == 1:
# ================================================================
# main starts
# ================================================================
   sig_level = 0.05
   year_ref = 1921 
   year1 = 1955
   year2 = 2005 
   t0 = int((year1-year_ref))
   t1 = int((year2-year_ref)+1)

# read file
   dirname = ''
   filename = 'obs_le_djf_tmp.nc'
   f = Dataset(dirname + filename, 'r')
   ts_global_djf = f.variables['ts_global'][:,t0:t1].data*3.
   ts_arctic_djf = f.variables['ts_arctic'][:,t0:t1].data*3.
   f.close()

   nens = ts_global_djf.shape[0]
   nt = ts_global_djf.shape[1]

   filename = 'obs_le_mam_tmp.nc'
   f = Dataset(dirname + filename, 'r')
   ts_global_mam = f.variables['ts_global'][:,t0:t1].data*3.
   ts_arctic_mam = f.variables['ts_arctic'][:,t0:t1].data*3.
   f.close()

   filename = 'obs_le_jja_tmp.nc'
   f = Dataset(dirname + filename, 'r')
   ts_global_jja = f.variables['ts_global'][:,t0:t1].data*3.
   ts_arctic_jja = f.variables['ts_arctic'][:,t0:t1].data*3.
   f.close()

   filename = 'obs_le_son_tmp.nc'
   f = Dataset(dirname + filename, 'r')
   ts_global_son = f.variables['ts_global'][:,t0:t1].data*3.
   ts_arctic_son = f.variables['ts_arctic'][:,t0:t1].data*3.
   f.close()

   ts_global = (ts_global_djf + ts_global_mam + ts_global_jja + ts_global_son).copy()/12.
   ts_arctic = (ts_arctic_djf + ts_arctic_mam + ts_arctic_jja + ts_arctic_son).copy()/12.

   nens = ts_global.shape[0]
   nt = ts_global.shape[1]

   tt = np.linspace(1,nt,nt)

   ts_global_trend = np.zeros((nens))
   ts_arctic_trend = np.zeros((nens))

   for NENS in range(nens):
       [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_global[NENS,:])
       ts_global_trend[NENS] = slope.copy()*10
       [slope, intercept, r_value, p_value, std_err] = stats.linregress(tt,ts_arctic[NENS,:])
       ts_arctic_trend[NENS] = slope.copy()*10

# ================================================================
# Output variables
# ================================================================
if True:
   filename = 'obs_le_trend_tmp.nc'
   files = os.listdir(os.getcwd())
   for file in files:
       if file == filename:
          print('Delete ' + filename)
          os.system("rm -rf " + file)

   f = Dataset(filename, 'w', format='NETCDF4')
   f.createDimension('samples', nens)

   ts_global_trend_out = f.createVariable('ts_global_trend','f4',('samples'))
   ts_arctic_trend_out = f.createVariable('ts_arctic_trend','f4',('samples'))

   ts_global_trend_out[:] = ts_global_trend[:].copy()
   ts_arctic_trend_out[:] = ts_arctic_trend[:].copy()

   f.close()

   sys.exit()


