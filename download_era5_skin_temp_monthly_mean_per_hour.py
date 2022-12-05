# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 12:08:48 2021

@author: Nirajan
"""
#%%
# from datetime import datetime, timedelta
import os 
import cdsapi

dirpath = 'G:\G-Drive\Python_folder\ERA5_data\skin_temp_monthly_mean_by_hour' # change folder
yaers=[2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018, 2019, 2020]
#yaers = [2021]
# https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
# this link says it released from 2000 to date
#%%
c = cdsapi.Client()

for yr in yaers:
    out_file="ERA5_skin_temperature_Nepal_" + str(yr) + ".nc" # change filename as per your data
    out_filename=os.path.join(dirpath,out_file).replace("\\","/")
    c.retrieve(
		'reanalysis-era5-land-monthly-means',
		{
			'format':'netcdf',
            'product_type': 'monthly_averaged_reanalysis_by_hour_of_day',
			'variable':'skin_temperature',# change varaible for your data
			'year':[
				str(yr)
			],
            
			'month':[
				
				'01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12', # change the months
			],
            'time':[
                '00:00','01:00','02:00',
                '03:00','04:00','05:00',
                '06:00','07:00','08:00',
                '09:00','10:00','11:00',
                '12:00','13:00','14:00',
                '15:00','16:00','17:00',
                '18:00','19:00','20:00',
                '21:00','22:00','23:00'
    			],
            'area':[
                    31, 79.5, 25.5, 89 # change area
            ],
		},
		out_filename)
    

# TypeError: 'type' object is not subscriptable for date 
##%%
#yr=2003
#out_file="ERA5_skin_temperature_Nepal_" + str(yr) + ".nc"
#outfilename=os.path.join(dirpath,out_file)
#print(outfilename)