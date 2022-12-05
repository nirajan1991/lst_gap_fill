# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 13:55:20 2021
This script is used to fill the gaps in MODIS LST
MODIS LST series consists of all the available data at monthly scale
i.e., MOD11C3 Day, MYD11C3 Day, MOD11C3 Night, MYD11C3 Night in sequence
The sequence in the script is:
1. Decode LST QC flags and separate them into good and bad data based on error
to be tolerated to be good and masking out the bad ones
2. Interpolate ERA5 skin temperature at the pixel centers of MODIS LST;
stack the interpolated data to create a data cube for all the years
3. Create a reference dataset for gap filling using the medain of
difference between good quality data and era-skt
4. Replace bad data in LST by using the referene data
5. Perform smoothing at non-unioform spacing using savitzky-golay filter
6. Test the accuracy of gap filling using rmse and R2

@author: Nirajan
"""

#%%
#from osgeo import gdal
import numpy as np
from netCDF4 import Dataset
from osgeo import osr
import os
import glob
#import scipy # this is needed to detrend timeseries
# import time
# from tqdm import tqdm
from datetime import datetime
#instead of writing all the functions here, I have kept all the functions in a
#separate file so that I can simply import and apply those functions
#%%
from lst_skt_fill_module import interpolate_params, era_modis_difference_median, \
                            apply_trilinear_fill, pad_timeseries, savgol_filter_3d, fill_lst, \
                            timeseries_stack, pixcenters_raster


#%%
proj = osr.SpatialReference() # gdal.GetProjection reads in WKT format so use the same
proj.ImportFromEPSG(4326)


qcbits = np.arange(0,256)
era_file = 'G:\G-Drive\Python_folder\LST_ERA_gapfill\ERA_MODIS_extracted_20210328\ERA_MODIS_time_sample_2003_skt_band1.tif'
modis_file = 'G:\G-Drive\MODIS_LST\MOD11C3.A2003001.006.2015182172300.hdf_Resampled.Clear_sky_days.tif'

pathdir = 'G:\G-Drive\Python_folder\LST_ERA_gapfill\ERA_MODIS_extracted_20210328'
pathdir = pathdir.replace('\\', '/') #use forward slash to avoid the conflict

search_criteria = 'ERA_MODIS_time_sample_*.nc'
search_path = os.path.join(pathdir, search_criteria).replace('\\', '/')
filelist = glob.glob(search_path)
#Convert backwad slash to forward slash
filelist = [i.replace('\\', '/') for i in filelist]

modis_lst_var = 'modis/lst'
era_skt_var = 'era5land/skt'
modis_qc_var = 'modis/qc'
time_var = 'era5land/time'
lst_error_tolerance = 3

window_length = 7
polyorder = 2
#%%
row, col, W1, W2, W3, W4 = interpolate_params(era_file, modis_file)

lst_era_diff_median = era_modis_difference_median(filelist, era_skt_var, modis_lst_var, modis_qc_var, qcbits, lst_error_tolerance, time_var, row, col, W1, W2, W3, W4)
#lons, lats = pixcenters_raster(modis_file)

#apply the function
lst_era_diff_median_filled = apply_trilinear_fill(lst_era_diff_median, 4)
del lst_era_diff_median

timeseries_padded = pad_timeseries(filelist, time_var)
lst_era_diff_median_filled_smoothed = savgol_filter_3d(timeseries_padded, lst_era_diff_median_filled, window_length, polyorder)
del lst_era_diff_median_filled

out = fill_lst(lst_era_diff_median_filled_smoothed, filelist,qcbits, lst_error_tolerance,  row, col, W1, W2, W3, W4, window_length, polyorder)
#%%
yaers = list(range(2003, 2021))
timedata = timeseries_stack(filelist, time_var)
londata, latdata = pixcenters_raster(modis_file)
modis_filled = out*100
modis_filled[np.isnan(out)] = 65535
modis_filled = modis_filled.astype(np.int16)
ds = Dataset(filelist[0])
time_ds = ds[time_var]
timeunit = time_ds.units
outdir = pathdir
#%%
#save as netcdf file
#http://pyhogs.github.io/intro_netcdf4.html
#open a file in write mode
out_filename = 'ERA_MODIS_filled.nc'
outfile = os.path.join(outdir, out_filename)
ds_out = Dataset(outfile,'w', format='NETCDF4') #'w' stands for write

#creat dimensions
ds_out.createDimension('lon', len(londata))
ds_out.createDimension('lat', len(latdata))
ds_out.createDimension('year', len(yaers))
ds_out.createDimension('time', len(timedata))

#create variables
modis_lonvar = ds_out.createVariable('lon', 'f4', 'lon')
modis_latvar = ds_out.createVariable('lat', 'f4', 'lat')
modis_timevar = ds_out.createVariable('time', 'i4', ('time','year'))
modis_year = ds_out.createVariable('year', 'u2', 'year')
modis_tempvar = ds_out.createVariable('lst_filled','u2',('time', 'lat', 'lon', 'year'), fill_value = 65535)

#pass data into variables
modis_lonvar[:] = londata
modis_latvar[:] = latdata
modis_timevar[:] = timedata
modis_year[:] = yaers
modis_tempvar[:] = modis_filled

#adding attributes
modis_lonvar.units = 'degree east'
modis_latvar.units = 'degree north'
modis_timevar.units = timeunit
modis_tempvar.units = 'Kelvin'
modis_tempvar.scale_factor = 0.01

today = datetime.today()
ds_out.history = "Created " + today.strftime("%d/%m/%y")
ds_out.projection = 'EPSG:4326'

ut = ds_out['lst_filled'][:]
#close the file
ds_out.close()
ds_out = None