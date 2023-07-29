# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 17:38:30 2021

@author: Nirajan

This script reads ERA5 skin temperature data and MODIS LST data
Then selects the ERA5 data closest to MODIS data
Builds netcdf file for LST as well as chosen ERA5 data
"""
#%%
import numpy as np
#import netCDF4
from netCDF4 import Dataset
from netCDF4 import date2num #num2date, 
from datetime import datetime, timedelta #timezone, 
from osgeo import gdal
import re
import os
import glob
#from matplotlib import pyplot as plt
import time as time_lib
#%%


#%%
myddir = 'G:\G-Drive\MODIS_LST\myd11C3\Resampled'
moddir = 'G:\G-Drive\MODIS_LST\MOD11C3\Resampled'
era5dir = 'G:\G-Drive\Python_folder\ERA5_data\skin_temp_monthly_mean_by_hour'
outdir = 'G:\G-Drive\Python_folder\LST_ERA_gapfill\ERA_MODIS_extracted_20210328'
yaers = list(range(2003, 2021))

#%%
#define a function to read raster data or import from the list of functions file
from lst_era_fill_functions import pixcenters_raster, read_raster_data, read_raster_nodata

for yr in yaers:
    era5file = 'ERA5_skin_temperature_Nepal_' + str(yr) +'.nc'
    era5file = os.path.join(era5dir, era5file)
    #%%
    #check the time of ERA5 data
    #file_era = 'G:/G-Drive\Python_folder/ERA5_data/skin_temp_monthly_mean_by_hour/ERA5_skin_temperature_Nepal_2003.nc'
    ds = Dataset(era5file)
    time = ds['time']
    sktvar = ds['skt']
    skt_sf = sktvar.scale_factor
    skt_ao = sktvar.add_offset
    skt = ds.variables['skt'][:]
    timedata = time[:]
    londata_era = ds['longitude'][:]
    latdata_era = ds['latitude'][:]
    timeunit = time.units
    fillval = sktvar._FillValue
    ds = None
    #%%
    search_lst = '*C3.A' + str(yr) + '*LST*.tif'
    search_qc = '*C3.A' + str(yr) + '*QC*.tif'
    search_view_time = '*C3.A' + str(yr) + '*view_time.tif'
    
    mod_search_lst_items=os.path.join(moddir,search_lst)
    mod_lst_filelist=glob.glob(mod_search_lst_items)
    
    myd_search_lst_items=os.path.join(myddir,search_lst)
    myd_lst_filelist=glob.glob(myd_search_lst_items)
    
    mod_search_qc_items=os.path.join(moddir,search_qc)
    mod_qc_filelist=glob.glob(mod_search_qc_items)
    
    myd_search_qc_items=os.path.join(myddir,search_qc)
    myd_qc_filelist=glob.glob(myd_search_qc_items)
    
    mod_search_view_time_items=os.path.join(moddir,search_view_time)
    mod_view_time_filelist=glob.glob(mod_search_view_time_items)
    
    myd_search_view_time_items=os.path.join(myddir,search_view_time)
    myd_view_time_filelist=glob.glob(myd_search_view_time_items)
    
    #%%
    
    londata_modis, latdata_modis = pixcenters_raster(mod_lst_filelist[0])
    #%%
    time_list = []
    lst_timeseries = np.zeros((48, len(latdata_modis), len(londata_modis))) #there are four data for each month
    qc_timeseries = np.zeros((48, len(latdata_modis), len(londata_modis))) #there are four data for each month
    
    for fn in range(0,12):
        #It reads day file first
        
        mod_view_time_file = mod_view_time_filelist[2*fn]
        myd_view_time_file = myd_view_time_filelist[2*fn]
        
        pos = re.search('C3.A', mod_view_time_file)
        
        time_nodata = read_raster_nodata(mod_view_time_file)
        time_mod = read_raster_data(mod_view_time_file)
        time_myd = read_raster_data(myd_view_time_file)

        time_mod_mean = time_mod[time_mod != time_nodata].mean()
        time_myd_mean = time_myd[time_myd != time_nodata].mean()
        
        date_modis = mod_view_time_file[pos.span()[1]:pos.span()[1]+7]
        yy = int(date_modis[0:4])
        jd = int(date_modis[4:])
        
        mod_datetime = datetime(yy,1,1,int(round(time_mod_mean*0.2))) + timedelta(days=jd-1)
        mod_datenum = date2num(mod_datetime, timeunit)
        
        myd_datetime = datetime(yy,1,1,int(round(time_myd_mean*0.2))) + timedelta(days=jd-1)
        myd_datenum = date2num(myd_datetime, timeunit)
        
        time_list.append(mod_datenum)
        time_list.append(myd_datenum)
        
        mod_lst_file = mod_lst_filelist[2*fn]
        lst_nodata = read_raster_nodata(mod_lst_file)
        mod_lst_data = read_raster_data(mod_lst_file)

        myd_lst_file = myd_lst_filelist[2*fn]
        myd_lst_data = read_raster_data(myd_lst_file)
        
        mod_qc_file = mod_qc_filelist[2*fn]
        qc_nodata = read_raster_nodata(mod_qc_file)
        mod_qc_data = read_raster_data(mod_qc_file)

        myd_qc_file = myd_qc_filelist[2*fn]
        myd_qc_data = read_raster_data(myd_qc_file)

        lst_timeseries[4*fn, :, :] = mod_lst_data
        lst_timeseries[4*fn+1, :, :] = myd_lst_data
        
        qc_timeseries[4*fn, :, :] = mod_qc_data
        qc_timeseries[4*fn+1, :, :] = myd_qc_data
        
        
        #Then goes to night file
        mod_view_time_file = mod_view_time_filelist[2*fn+1]
        myd_view_time_file = myd_view_time_filelist[2*fn+1]
        
        pos = re.search('C3.A', mod_view_time_file)
        
        time_mod = read_raster_data(mod_view_time_file)
        time_nodata = read_raster_nodata(mod_view_time_file)
        time_myd = read_raster_data(myd_view_time_file)
        
        time_mod_mean = time_mod[time_mod != time_nodata].mean()
        time_myd_mean = time_myd[time_myd != time_nodata].mean()
        
        date_modis = mod_view_time_file[pos.span()[1]:pos.span()[1]+7]
        yy = int(date_modis[0:4])
        jd = int(date_modis[4:])
        
        mod_datetime = datetime(yy,1,1,int(round(time_mod_mean*0.2))) + timedelta(days=jd-1)
        mod_datenum = date2num(mod_datetime, timeunit)
        
        myd_datetime = datetime(yy,1,1,int(round(time_myd_mean*0.2))) + timedelta(days=jd-1)
        myd_datenum = date2num(myd_datetime, timeunit)
        
        time_list.append(mod_datenum)
        time_list.append(myd_datenum)
        
        mod_lst_file = mod_lst_filelist[2*fn+1]
        mod_lst_data = read_raster_data(mod_lst_file)

        myd_lst_file = myd_lst_filelist[2*fn*1]
        myd_lst_data = read_raster_data(myd_lst_file)
        
        mod_qc_file = mod_qc_filelist[2*fn+1]
        mod_qc_data = read_raster_data(mod_qc_file)
        
        myd_qc_file = myd_qc_filelist[2*fn+1]
        myd_qc_data = read_raster_data(myd_qc_file)
        
        lst_timeseries[4*fn+2, :, :] = mod_lst_data
        lst_timeseries[4*fn+3, :, :] = myd_lst_data
        
        qc_timeseries[4*fn+2, :, :] = mod_qc_data
        qc_timeseries[4*fn+3, :, :] = myd_qc_data
        
        fn +=1
        
        #time_lib.sleep(10)
    #%%
    #now test the result of time series
    #lst_series = lst_timeseries[:,0,0]
    #plt.plot(lst_series)
    #plt.plot(time_list)
    
    #%%
    idx = np.isin(timedata,time_list).nonzero() #https://stackoverflow.com/questions/32191029/getting-the-indices-of-several-elements-in-a-numpy-array-at-once
    
    skt_modis = skt[idx[0], :,:].copy()
    modis_datetime = np.ma.core.MaskedArray(time_list)
    skt_modis_compress = (skt_modis - skt_ao)/skt_sf
    #%%
    '''
    https://unidata.github.io/netcdf4-python/
    Notes in datatype of netcdf4 library
    'f4' (32-bit floating point), 
    'f8' (64-bit floating point), 
    'i4' (32-bit signed integer), 
    'i2' (16-bit signed integer), 
    'i8' (64-bit signed integer), 
    'i1' (8-bit signed integer), 
    'u1' (8-bit unsigned integer), 
    'u2' (16-bit unsigned integer), 
    'u4' (32-bit unsigned integer), 
    'u8' (64-bit unsigned integer), or 
    'S1' (single-character string).
    '''
    #http://pyhogs.github.io/intro_netcdf4.html
    #open a file in write mode
    out_filename = 'ERA_MODIS_time_sample_' + str(yr) + '.nc'
    outfile = os.path.join(outdir, out_filename)
    ds_out = Dataset(outfile,'w', format='NETCDF4') #'w' stands for write
    
    #time is common in both so create time variable in root
    
    # create groups to separate era and modis data
    era5_group = ds_out.createGroup('era5land')
    #creat dimensions
    era5_group.createDimension('era_lon', len(londata_era))
    era5_group.createDimension('era_lat', len(latdata_era))
    
    era5_group.createDimension('time', len(modis_datetime))
    #tempgrp.createDimension('z', len(modis_datetime))
    #create variables
    era_lonvar = era5_group.createVariable('lon', 'f4', 'era_lon')
    era_latvar = era5_group.createVariable('lat', 'f4', 'era_lat')
    era_timevar = era5_group.createVariable('time', 'i4', 'time')
    era_tempvar = era5_group.createVariable('skt','i2',('time', 'era_lat', 'era_lon'), fill_value = fillval)
    #pass data into variables
    era_lonvar[:] = londata_era
    era_latvar[:] = latdata_era
    era_timevar[:] = modis_datetime
    era_tempvar[:] = skt_modis_compress
    
    era_lonvar.units = 'degree east'
    era_latvar.units = 'degree north'
    era_timevar.units = timeunit
    era_tempvar.units = 'Kelvin'
    era_tempvar.scale_factor = skt_sf
    era_tempvar.add_offset = skt_ao
    
    #create modis group
    modis_group = ds_out.createGroup('modis')
    #creat dimensions
    modis_group.createDimension('modis_lon', len(londata_modis))
    modis_group.createDimension('modis_lat', len(latdata_modis))
    
    modis_group.createDimension('time', len(modis_datetime))
    #tempgrp.createDimension('z', len(modis_datetime))
    #create variables
    modis_lonvar = modis_group.createVariable('lon', 'f4', 'modis_lon')
    modis_latvar = modis_group.createVariable('lat', 'f4', 'modis_lat')
    modis_timevar = modis_group.createVariable('time', 'i4', 'time')
    modis_qcvar = modis_group.createVariable('qc', 'u1', ('time', 'modis_lat', 'modis_lon'), fill_value = qc_nodata)
    modis_tempvar = modis_group.createVariable('lst','u2',('time', 'modis_lat', 'modis_lon'), fill_value = lst_nodata)
    #pass data into variables
    modis_lonvar[:] = londata_modis
    modis_latvar[:] = latdata_modis
    modis_timevar[:] = modis_datetime
    modis_tempvar[:] = lst_timeseries
    modis_qcvar[:] = qc_timeseries
    #adding attributes
    modis_lonvar.units = 'degree east'
    modis_latvar.units = 'degree north'
    modis_timevar.units = timeunit
    modis_tempvar.units = 'Kelvin'
    modis_tempvar.scale_factor = 0.02
    
    today = datetime.today()
    ds_out.history = "Created " + today.strftime("%d/%m/%y")
    ds_out.projection = 'EPSG:4326'
    
    #close the file
    ds_out.close()
    ds_out = None
    
    time_lib.sleep(30)

