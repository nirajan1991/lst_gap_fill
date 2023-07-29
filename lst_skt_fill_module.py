# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 07:44:13 2021
This stores all the functions used in the work
@author: Nirajan
"""
#%%
#import libraries
from osgeo import gdal
import numpy as np
from netCDF4 import Dataset
import time
from tqdm import tqdm

#%%
#First define the functions that are used for gap filling
#Define a functin to convert one sample data from netcdf file to geotiff file
#It is required for interpolation of ERA5 data, I will need a raster file from the era5 data
#So save first layer of the file as raster file

def nc2gtiff(ncfile,ncvar,proj):
    ds = Dataset(ncfile)

    era_lat = ds['era5land/lat'][:]
    era_lon = ds['era5land/lon'][:]
    #print(ds.variables.keys())
    '''
    era_lat = ds['latitude'][:]
    era_lon = ds['longitude'][:]
    '''
    era_data = ds[ncvar][:]
    era_data_layer = era_data[0,:,:]
    pixel_size = 0.1
    xx,yy = era_data_layer.shape
    #gt = [era_lon[0] - (era_lon[1] -era_lon[0]) , era_lon[1] -era_lon[0], 0,\
    #         era_lat[0] - (era_lat[1] -era_lat[0]), 0, era_lat[1] -era_lat[0]] # This gave less satisfactory data
    gt = [era_lon[0] - pixel_size/2 , pixel_size, 0,\
             era_lat[0] + pixel_size/2, 0, -pixel_size]


    #save the layer as tif file
    tif_driver = gdal.GetDriverByName("GTiff")
    tiffile = ncfile[:-3] + '_skt_band1.tif'
    tif_ds = tif_driver.Create(tiffile, yy, xx, 1, gdal.GDT_Float32)
    tif_ds.SetGeoTransform(gt)
    tif_ds.SetProjection(str(proj))
    tif_ds.GetRasterBand(1).SetNoDataValue(-32767)
    tif_ds.GetRasterBand(1).WriteArray(era_data_layer)
    tif_ds.FlushCache()
    ds.close()
    ds = None
    tif_ds = None
    return tiffile
#%%
# define a function to read lat and lon values as center of pixels
def pixcenters_raster(rasterfile):
    ds = gdal.Open(rasterfile)
    width = ds.RasterXSize
    height = ds.RasterYSize
    gt = ds.GetGeoTransform()
    gt = np.asarray(gt)

    lat_rows = np.arange(1,height+1).reshape(-1,1)
    lon_cols = np.arange(1,width+1).reshape(-1,1)

    cols_mat = np.concatenate((np.tile(1, (width, 1)), lon_cols, np.tile(1, (width,1))),axis=1)
    lons = cols_mat @ gt[:3].reshape(-1,1)
    lons_adj = lons - gt[1]/2

    rows_mat = np.concatenate((np.tile(1, (height, 1)), np.tile(1, (height,1)), lat_rows),axis=1)
    lats = rows_mat @ gt[3:].reshape(-1,1)
    lats_adj = lats - gt[5]/2

    ds = None

    return lons_adj, lats_adj

#%%
# define a function to read only data from raster
def read_raster_data(file):
    ds = gdal.Open(file)
    data = ds.ReadAsArray()
    ds = None
    return data

def read_raster_nodata(file):
    ds = gdal.Open(file)
    band = ds.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    ss = None
    return nodata

#%%
#define distance function
#The function was obtained from Martin Thomas' answer in the link and replaced math by numpy
#https://stackoverflow.com/questions/19412462/getting-distance-between-two-points-based-on-latitude-longitude

def distance_np(lat1,lon1,lat2,lon2):
    """
    Calculate the Haversine distance.
    The function was actually in math but I changed math to numpy
    It used coordinates of two points which I chaged to lats and lons of two points
    This is because my lats and lons are in 2D
    If places of lat and lon for both points are changed than it will be same
    """
    radius = 6371  # km

    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = (np.sin(dlat / 2) * np.sin(dlat / 2) +
         np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) *
         np.sin(dlon / 2) * np.sin(dlon / 2))
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    d = radius * c

    return d
#%%
#Now define a function to calculate the lists of rows and columns of the indices for given lat and lon
def interpolate_params(era_file, modis_file):
    #read era raster file
    era_ds = gdal.Open(era_file)
    era_gt = era_ds.GetGeoTransform()
    #era_gt = np.asarray(era_gt)

    xOrigin = era_gt[0]
    yOrigin = era_gt[3]
    pixelWidth = era_gt[1]
    pixelHeight = -era_gt[5]

    #get the lat lon of modis data file
    lon_modis, lat_modis = pixcenters_raster(modis_file)
    lon_era, lat_era = pixcenters_raster(era_file)
    #%%
    
    #Vectorizing the process created hussle as the original pixel retrival was nontrivial
    #So create function to extract the row col value of each pixel at the pixel (2D)

    lon_modis_mesh, lat_modis_mesh = np.meshgrid(lon_modis,lat_modis)

    col = (lon_modis_mesh - xOrigin)/pixelWidth
    row = (yOrigin - lat_modis_mesh)/pixelHeight

    # create coordinate pairs for each corners era pixel for each modis pixel center
    c1_lat, c1_lon = lat_era[np.floor(row).astype(int)],lon_era[np.floor(col).astype(int)]
    c2_lat, c2_lon = lat_era[np.floor(row).astype(int)],lon_era[np.ceil(col).astype(int)]
    c3_lat, c3_lon = lat_era[np.ceil(row).astype(int)],lon_era[np.floor(col).astype(int)]
    c4_lat, c4_lon = lat_era[np.ceil(row).astype(int)],lon_era[np.ceil(col).astype(int)]

    #Since these values are in 3D need to squeeze them
    c1_lat,c1_lon = np.squeeze(c1_lat,axis = 2),np.squeeze(c1_lon,axis = 2)
    c2_lat,c2_lon = np.squeeze(c2_lat,axis = 2),np.squeeze(c2_lon,axis = 2)
    c3_lat,c3_lon = np.squeeze(c3_lat,axis = 2),np.squeeze(c3_lon,axis = 2)
    c4_lat,c4_lon = np.squeeze(c4_lat,axis = 2),np.squeeze(c4_lon,axis = 2)

    #Claculate distance between points
    d1 = distance_np(c1_lat, c1_lon, lat_modis_mesh, lon_modis_mesh)
    d2 = distance_np(c2_lat, c2_lon, lat_modis_mesh, lon_modis_mesh)
    d3 = distance_np(c3_lat, c3_lon, lat_modis_mesh, lon_modis_mesh)
    d4 = distance_np(c4_lat, c4_lon, lat_modis_mesh, lon_modis_mesh)
    dmax = distance_np(c1_lat, c1_lon,c4_lat, c4_lon)


    #weighted distance
    D1 = np.power(np.cos((np.pi/2) * (d1/dmax)),4)
    D2 = np.power(np.cos((np.pi/2) * (d2/dmax)),4)
    D3 = np.power(np.cos((np.pi/2) * (d3/dmax)),4)
    D4 = np.power(np.cos((np.pi/2) * (d4/dmax)),4)

    Di_sum = D1 + D2 + D3 + D4

    #calculate weights
    W1 = D1/Di_sum
    W2 = D2/Di_sum
    W3 = D3/Di_sum
    W4 = D4/Di_sum

    era_ds = None

    #need to return indices_np_2d as well because I need to interpolate to the points as reshaping does not do the job
    return row, col, W1, W2, W3, W4 #indices_np_2d,lon_modis, lat_modis#, c1,lon_era,lat_era

#%%

#define function to interpolate for a single point
import numpy as np
def interpolate_data(era_data, row, col, W1, W2, W3, W4):
    V1=era_data[np.floor(row).astype(int),np.floor(col).astype(int)]#.reshape(-1,1)
    V2=era_data[np.floor(row).astype(int),np.ceil(col).astype(int)]#.reshape(-1,1)
    V3=era_data[np.ceil(row).astype(int),np.floor(col).astype(int)]#.reshape(-1,1)
    V4=era_data[np.ceil(row).astype(int),np.ceil(col).astype(int)]#.reshape(-1,1)
    V=W1*V1+W2*V2+W3*V3+W4*V4

    return V #interpolated_data

#define a function to get interpolated data for all the time period in the file
def interpolate_timeseries(era_timeseries,era_skt_var, row, col, W1, W2, W3, W4):
    timelength, _, _ = era_timeseries.shape
    xsize, ysize = row.shape
    interpolated_timeseries = np.zeros((timelength, xsize, ysize))
    for ii in range(timelength):
        era_data = era_timeseries[ii,:,:]
        interpolated_data = interpolate_data(era_data, row, col, W1, W2, W3, W4)
        interpolated_timeseries[ii,:,:] = interpolated_data
    return interpolated_timeseries

#%%

#First write the function to decode the quality flag and separate them into lists different categories

#Make a vectorized function for converting decimel to binary
de2bi = np.vectorize(np.binary_repr)

#Then create function to slice an array of strings
def slicer_vectorized(a,start,end):
    b = a.view((str,1)).reshape(len(a),-1)[:,start:end]
    return np.frombuffer(b.tobytes(),dtype=(str,end-start))

#%%
def decode_qc(qcbits):
    qc_big_endian = de2bi(qcbits,8) # 8 for each of length 8 by padding zeros
    #get modland qc and lst error flags
    modland = slicer_vectorized(qc_big_endian, 6, 8)
    lsterr = slicer_vectorized(qc_big_endian, 0, 2)
    # determine qc flags according to modland qc bits
    mod1 = modland == '00'
    mod2 = modland == '01'
    mod3 = modland == '10'
    mod4 = modland == '11'
    # determine qc flags according to lst error qc bits
    lsterr1 = lsterr == '00'
    lsterr2 = lsterr == '01'
    lsterr3 = lsterr == '10'
    lsterr3plus = lsterr == '11'
    # combine modland qc bits and lst error tolerance
    lst_error_lessthan1 = qcbits[np.logical_and(np.logical_or(mod1,mod2), lsterr1)]
    lst_error_lessthan2 = qcbits[np.logical_and(np.logical_or(mod1,mod2), lsterr2)]
    lst_error_lessthan3 = qcbits[np.logical_and(np.logical_or(mod1,mod2), lsterr3)]
    lst_error_morethan3 = qcbits[np.logical_and(np.logical_or(mod1,mod2), lsterr3plus)]
    lst_cloud = qcbits[mod3];
    lst_other = qcbits[mod4];
    # keep every thing in a single dictinary
    lst_error_qcbits = {'err1': list(lst_error_lessthan1),
                     'err2': list(lst_error_lessthan2),
                     'err3': list(lst_error_lessthan3),
                     'err3plus': list(lst_error_morethan3),
                     'cloud': list(lst_cloud),
                     'other': list(lst_other)}
    #decoding qcbits can be done more efficiently with leftshift which I was not aware of before writing the function
    return lst_error_qcbits


#%%
#Define a function to convert qc_data into qc_rank
def qc_data2rank(qc_data, qcbits):
    #first separate the qcbits according to lst error
    lst_err_qcbits = decode_qc(qcbits)
    #get the masks for each error level
    qc_1 = np.isin(qc_data,lst_err_qcbits['err1']);
    qc_2 = np.isin(qc_data,lst_err_qcbits['err2']);
    qc_3 = np.isin(qc_data,lst_err_qcbits['err3']);
    qc_3plus = np.isin(qc_data,lst_err_qcbits['err3plus']);
    qc_cloud = np.isin(qc_data,lst_err_qcbits['cloud']);
    qc_other = np.isin(qc_data,lst_err_qcbits['other']);

    #convert them into numbers
    qc_rank_1 = np.multiply(qc_1, 1)
    qc_rank_2 = np.multiply(qc_2, 2)
    qc_rank_3 = np.multiply(qc_3, 3)
    qc_rank_3plus = np.multiply(qc_3plus, 4)
    qc_rank_cloud = np.multiply(qc_cloud, 5)
    qc_rank_other = np.multiply(qc_other, 6)

    # Adding all qc ranks means keeping one value in one pixel as other pixel vlaues for each ranke are zeros
    qc_rank_all = qc_rank_1 + qc_rank_2 + qc_rank_3 + qc_rank_3plus + qc_rank_cloud + qc_rank_other

    return qc_rank_all

#%%
#Define a function to stack data for all the years for that variable
def timeseries_stack(filelist, input_var):
    ncfilename = filelist[0]
    nc_dataset = Dataset(ncfilename)
    data_timeseries = nc_dataset[input_var][:] #there are no variables in the root
    nc_dataset = None
    #check the number of dimensions to separate the type of output required
    ndim = data_timeseries.ndim

    if ndim == 3:

        timesize, xsize, ysize = data_timeseries.shape
        data_stacked = np.zeros((timesize, xsize, ysize, len(filelist)))


        for i,ncfilename in enumerate(filelist):
            nc_dataset = Dataset(ncfilename)
            data_timeseries = nc_dataset[input_var][:]
            data_stacked[:,:,:,i] = data_timeseries
            nc_dataset = None

    if ndim == 1:
        timesize = len(data_timeseries)
        data_stacked = np.zeros((timesize, len(filelist)))

        for i,ncfilename in enumerate(filelist):
            nc_dataset = Dataset(ncfilename)
            data_timeseries = nc_dataset[input_var][:]
            data_stacked[:,i] = data_timeseries
            nc_dataset = None
    return data_stacked


#%%

# Define a function to create a median reference
def era_modis_difference_median(filelist, era_skt_var, modis_lst_var, modis_qc_var, qcbits, lst_error_tolerance, time_var, row, col, W1, W2, W3, W4):
    #Stack data for all the variables in the file
    era_skt_stacked = timeseries_stack(filelist, era_skt_var)
    modis_lst_stacked = timeseries_stack(filelist, modis_lst_var)
    modis_qc_stacked = timeseries_stack(filelist, modis_qc_var)
    #time_stacked = timeseries_stack(filelist, time_var)
    #Then interpolate era5 data to the size of modis data
    era_interpolated_stacked = np.zeros(modis_lst_stacked.shape)
    for i in range(len(filelist)):
        era_skt_timeseries = era_skt_stacked[:,:,:,i]
        era_interpolated_timeseries = interpolate_timeseries(era_skt_timeseries,era_skt_var, row, col, W1, W2, W3, W4)
        era_interpolated_stacked[:,:,:,i] = era_interpolated_timeseries

    #Get qc ranks from qc data for applying fiter
    modis_qc_rank = qc_data2rank(modis_qc_stacked, qcbits)
    modis_qc_bad = np.logical_and(modis_qc_rank > lst_error_tolerance,\
                                  modis_lst_stacked == 0) #Should be 1K, 2K or 3K because anything morethan 3K error is ambigious
    modis_lst_stacked[modis_qc_bad] = np.nan
    #reading lst data with netCDF4.Dataset automatically converts fill value to nan
    #so need to use lst fill value to create void for lst missing data

    #calculate differnce between modis_lst_stacked with voids and era_skt_data
    lst_era_diff = modis_lst_stacked - era_interpolated_stacked
    #This creates voids in the difference matrix corresponding to missing values in lst

    #then caculate the median of difference
    lst_era_diff_median_3d = np.median(lst_era_diff, axis = 3, keepdims = False) #If keepdims is set to True it will become 4D

    #Also calculate median of time_stacked
    #time_stacked_median = np.median(time_stacked, axis =1)
    return lst_era_diff_median_3d#, time_stacked

#%%
#Write a function to fill the missing values in loop
def trilinear_interpolate_fill(lst_era_diff_median, min_neighbor):
    #find the places where nan is present
    nan_pos = np.where(np.isnan(lst_era_diff_median))
    #read shape to get the last item in each dimension
    xx,yy,zz = lst_era_diff_median.shape
    #get a filter that selects nan position only for non edge positions
    nan_pos_reduced_idx = np.where(np.logical_and.reduce(\
                                  (nan_pos[0]!=0,nan_pos[1]!=0,nan_pos[2]!=0,\
                                   nan_pos[0]!=xx-1,nan_pos[1]!=yy-1,nan_pos[2]!=zz-1)))

    #apply the filter to remove edge positions
    nan_pos_reduced = (nan_pos[0][nan_pos_reduced_idx],nan_pos[1][nan_pos_reduced_idx],nan_pos[2][nan_pos_reduced_idx])
    #separate them into x,y,z cordinate points
    xnan, ynan, znan = nan_pos_reduced
    #I could have used nan_pos_reduced but it made me easier to understand

    #initialize a vector to collect all the fill values
    fill_values = np.zeros(xnan.shape)
    valid_buffer_data = np.zeros(xnan.shape)
    for ii in range(len(xnan)):
        #get the iith nan pixel's neighbors
        nan_buffer_x = lst_era_diff_median[[xnan[ii]-1,xnan[ii]+1], ynan[ii], znan[ii]]
        nan_buffer_y = lst_era_diff_median[xnan[ii], [ynan[ii]-1,ynan[ii]+1], znan[ii]]
        nan_buffer_z = lst_era_diff_median[xnan[ii], ynan[ii], [znan[ii]-1,znan[ii]+1]]
        #create an array of all the neighbors
        nan_buffer = np.array([nan_buffer_x,nan_buffer_y,nan_buffer_z])
        #to avoid over reliance limited data check how may neighbor values are valid
        no_valid = len(np.where(~np.isnan(nan_buffer.ravel()))[0])
        #find the mean value of surrounding pixels
        nan_buffer_mean = np.nanmean(nan_buffer) #although it creates 2 d array it calculates mean using all the values

        # #put the value in the position so that next item will make the benefit of this filled value
        # #I could have collected only the fill values and replaced outside loop
        # #but continuously filling the data would increase the
        # lst_era_diff_median[xnan[ii], ynan[jj], znan[kk]] = nan_buffer_mean

        #I abandoned that idea of filling continuously using the same filling values because it would create serious bias
        #especially around that area where one value will be continuously used to fill all the nans in extended voids
        #Although the function needs to be iterated many times in the second case

        fill_values[ii] = nan_buffer_mean
        valid_buffer_data[ii] = no_valid

    #remove those filled values that were generated with less than half of neighbors
    invalid_fill = valid_buffer_data<min_neighbor #could not fill the voids if at least 3 neighbors were required, so settled down for 2 values
    fill_values[invalid_fill] = np.nan
    #now use these mean values to fill the gaps in their respective positions
    lst_era_diff_median[xnan, ynan, znan] = fill_values
    #Viola this works superbly, from 46k to 39k

    return lst_era_diff_median

#define another function to calculate the number of missing value and
#apply the trilinear interpolation iteratively till all the values are not filled
def apply_trilinear_fill(lst_era_diff_median, min_neighbor):
    #pad a days data before
    lst_era_diff_median_fill = np.concatenate((lst_era_diff_median[-4:,:,:].copy(), lst_era_diff_median.copy(),lst_era_diff_median[:4,:,:].copy()),axis = 0)
    #count nans before applying the function
    diff_nans, _, _ = np.where(np.isnan(lst_era_diff_median_fill))
    no_of_nans = len(diff_nans)
    #get the number of pixels in border of each dimension to set as maximum number of voids that cannot be filled
    xx,yy,zz = lst_era_diff_median_fill.shape
    nan_pos = np.where(np.isnan(lst_era_diff_median_fill))
    nan_pos_border_idx = np.where(~np.logical_and.reduce(\
                                  (nan_pos[0]!=0,nan_pos[1]!=0,nan_pos[2]!=0,\
                                    nan_pos[0]!=xx-1,nan_pos[1]!=yy-1,nan_pos[2]!=zz-1)))
    #I ignored the borders in time dimension as I padded one months data in start and end

    no_border_void = len(nan_pos_border_idx[0])
    print('Number of voids in border: ', no_border_void)
    while no_of_nans > no_border_void:
        time.sleep(5)
        lst_era_diff_median_fill = trilinear_interpolate_fill(lst_era_diff_median_fill, min_neighbor)
        #count nans before applying the function
        diff_nans, _, _ = np.where(np.isnan(lst_era_diff_median_fill))
        no_of_nans_after = len(diff_nans)
        if no_of_nans_after == no_of_nans:
            min_neighbor = min_neighbor - 1
        no_of_nans = no_of_nans_after
        print(no_of_nans_after)

    return lst_era_diff_median_fill #do not restore the original size because I need to do smoothening


#%%
#define a function to pad the timeseries to use as input to the smoothening function
def pad_timeseries(filelist, time_var):

    time_stacked = timeseries_stack(filelist, time_var);
    time_stacked_median = np.median(time_stacked, axis =1)
    time_stacked_median_monthly = time_stacked_median.reshape(12,4)
    time_difference_monthly = np.diff(time_stacked_median_monthly, axis=0)
    time_stacked_median_monthly_modified = time_stacked_median_monthly.copy()
    for i in range(time_stacked_median_monthly.shape[0]-1):
        j = time_stacked_median_monthly.shape[0]-1-i
        time_stacked_median_monthly_modified[j:,:] = time_stacked_median_monthly_modified[j:,:] - time_difference_monthly[j-1,:] + 24

    timeseries = time_stacked_median_monthly_modified.flatten()
    timeseries_padded = np.concatenate((timeseries[0:4]-24, timeseries, timeseries[-4:]+24))

    return timeseries_padded

#%%
#apply smoothing in timeseries
import numpy as np

def non_uniform_savgol(x, y, window, polynom):
    """
    Applies a Savitzky-Golay filter to y with non-uniform spacing
    as defined in x

    This is based on https://dsp.stackexchange.com/questions/1676/savitzky-golay-smoothing-filter-for-not-equally-spaced-data
    The borders are interpolated like scipy.signal.savgol_filter would do

    Parameters
    ----------
    x : array_like
        List of floats representing the x values of the data
    y : array_like
        List of floats representing the y values. Must have same length
        as x
    window : int (odd)
        Window length of datapoints. Must be odd and smaller than x
    polynom : int
        The order of polynom used. Must be smaller than the window size

    Returns
    -------
    np.array of float
        The smoothed y values
    """
    if len(x) != len(y):
        raise ValueError('"x" and "y" must be of the same size')

    if len(x) < window:
        raise ValueError('The data size must be larger than the window size')

    if type(window) is not int:
        raise TypeError('"window" must be an integer')

    if window % 2 == 0:
        raise ValueError('The "window" must be an odd integer')

    if type(polynom) is not int:
        raise TypeError('"polynom" must be an integer')

    if polynom >= window:
        raise ValueError('"polynom" must be less than "window"')

    half_window = window // 2
    polynom += 1

    # Initialize variables
    A = np.empty((window, polynom))     # Matrix
    tA = np.empty((polynom, window))    # Transposed matrix
    t = np.empty(window)                # Local x variables
    y_smoothed = np.full(len(y), np.nan)

    # Start smoothing
    for i in range(half_window, len(x) - half_window, 1):
        # Center a window of x values on x[i]
        for j in range(0, window, 1):
            t[j] = x[i + j - half_window] - x[i]

        # Create the initial matrix A and its transposed form tA
        for j in range(0, window, 1):
            r = 1.0
            for k in range(0, polynom, 1):
                A[j, k] = r
                tA[k, j] = r
                r *= t[j]

        # Multiply the two matrices
        tAA = np.matmul(tA, A)

        # Invert the product of the matrices
        tAA = np.linalg.inv(tAA)

        # Calculate the pseudoinverse of the design matrix
        coeffs = np.matmul(tAA, tA)

        # Calculate c0 which is also the y value for y[i]
        y_smoothed[i] = 0
        for j in range(0, window, 1):
            y_smoothed[i] += coeffs[0, j] * y[i + j - half_window]

        # If at the end or beginning, store all coefficients for the polynom
        if i == half_window:
            first_coeffs = np.zeros(polynom)
            for j in range(0, window, 1):
                for k in range(polynom):
                    first_coeffs[k] += coeffs[k, j] * y[j]
        elif i == len(x) - half_window - 1:
            last_coeffs = np.zeros(polynom)
            for j in range(0, window, 1):
                for k in range(polynom):
                    last_coeffs[k] += coeffs[k, j] * y[len(y) - window + j]

    # Interpolate the result at the left border
    for i in range(0, half_window, 1):
        y_smoothed[i] = 0
        x_i = 1
        for j in range(0, polynom, 1):
            y_smoothed[i] += first_coeffs[j] * x_i
            x_i *= x[i] - x[half_window]

    # Interpolate the result at the right border
    for i in range(len(x) - half_window, len(x), 1):
        y_smoothed[i] = 0
        x_i = 1
        for j in range(0, polynom, 1):
            y_smoothed[i] += last_coeffs[j] * x_i
            x_i *= x[i] - x[-half_window - 1]

    return y_smoothed
#%%
#Also apply smoothening in the temporal domain in filled median array
#from scipy.signal import savgol_filter #it uses only in uniform spaced
def savgol_filter_3d(timeseries_padded, lst_era_diff_median_filled, window_length, polyorder):
    lst_era_diff_median_filled_smoothed = lst_era_diff_median_filled.copy()
    xx, yy, zz = lst_era_diff_median_filled.shape
    for ii in tqdm(range(yy)):
        for jj in range(zz):
            lst_median_timeseries = lst_era_diff_median_filled[:,ii,jj]
            lst_median_sg = non_uniform_savgol(timeseries_padded, lst_median_timeseries, window_length, polyorder)
            lst_era_diff_median_filled_smoothed[:,ii,jj] = lst_median_sg

        #pause to prevent over heating
        time.sleep(5)
    #after savitzky golay filter has been applied clip out the padded data
    lst_era_diff_median_filled_smoothed_clipped = lst_era_diff_median_filled_smoothed[4:-4,:,:].copy()
    return lst_era_diff_median_filled_smoothed_clipped

#%%

def fill_lst(lst_era_diff_median_filled_smoothed, filelist, qcbits, lst_error_tolerance, row, col, W1, W2, W3, W4, window_length, polyorder):
    modis_lst_var = 'modis/lst'
    era_skt_var = 'era5land/skt'
    modis_qc_var = 'modis/qc'
    time_var = 'era5land/time'

    #Stack data for all the variables in the file
    era_skt_stacked = timeseries_stack(filelist, era_skt_var)
    modis_lst_stacked = timeseries_stack(filelist, modis_lst_var)
    modis_qc_stacked = timeseries_stack(filelist, modis_qc_var)
    time_stacked = timeseries_stack(filelist, time_var)
    data_shape = modis_lst_stacked.shape

    #Get qc ranks from qc data for applying fiter
    modis_qc_rank = qc_data2rank(modis_qc_stacked, qcbits)
    del modis_qc_stacked
    modis_qc_bad = np.logical_and(modis_qc_rank > lst_error_tolerance,\
                                  modis_lst_stacked == 1)#Should be 1K, 2K or 3K because anything morethan 3K error is ambigious
    del modis_qc_rank
    modis_lst_stacked[modis_qc_bad] = np.nan

    #Then interpolate era5 data to the size of modis data
    era_interpolated_stacked = np.zeros(data_shape, dtype=np.float32)
    for i in range(len(filelist)):
        era_skt_timeseries = era_skt_stacked[:,:,:,i]
        era_interpolated_timeseries = interpolate_timeseries(era_skt_timeseries,era_skt_var, row, col, W1, W2, W3, W4)
        era_interpolated_stacked[:,:,:,i] = era_interpolated_timeseries

    del era_skt_stacked
    #era_interpolated_stacked_3d = np.concatenate(era_interpolated_stacked, axis = 2) #This concatenates in 3rd dimension with 3rd dimension as time

    #create a timeseries at daily basis
    #flatten in columns major to get continuous time series
    time_stacked_flat = time_stacked.flatten('F')
    time_stacked_monthly = time_stacked_flat.reshape(-1, 4)
    time_difference_monthly = np.diff(time_stacked_monthly, axis=0)
    time_stacked_monthly_modified = time_stacked_monthly.copy()
    for i in range(time_stacked_monthly.shape[0]-1):
        j = time_stacked_monthly.shape[0]-1-i
        time_stacked_monthly_modified[j:,:] = time_stacked_monthly_modified[j:,:] - time_difference_monthly[j-1,:] + 24
    timeseries = time_stacked_monthly_modified.flatten('F')

    #expand the dimension of reference difference to match the dimension of era data
    diff_ref_newaxis = np.expand_dims(lst_era_diff_median_filled_smoothed,axis = 3)
    del lst_era_diff_median_filled_smoothed
    diff_ref = np.concatenate([diff_ref_newaxis]*18, axis = 3)
    del diff_ref_newaxis
    modis_lst_reference  = era_interpolated_stacked + diff_ref
    del diff_ref
    modis_lst_filled = modis_lst_stacked.copy()
    del modis_lst_stacked
    modis_lst_filled [modis_qc_bad] = modis_lst_reference[modis_qc_bad]
    del modis_lst_reference

    #smoothout the filled data to eliminate the sudden change in lsts
    #for that is is better to change the data to 3d
    modis_lst_filled_3d = np.concatenate(modis_lst_filled, axis = 2)

    modis_lst_filled_smoothed = np.zeros(modis_lst_filled_3d.shape)
    xx, yy, zz = modis_lst_filled_3d.shape
    for ii in tqdm(range(xx)):
        for jj in range(yy):
            lst_timeseries = modis_lst_filled_3d[ii,jj,:]
            lst_sg = non_uniform_savgol(timeseries, lst_timeseries, window_length, polyorder)
            modis_lst_filled_smoothed[ii,jj,:] = lst_sg

    #checking matlab script i realized tht ther may be some anamolous values undetected by quality flag
    modis_lst_filled_std = np.nanstd(modis_lst_filled_smoothed, axis =2, keepdims = True)
    modis_lst_filled_mean = np.nanmean(modis_lst_filled_smoothed, axis =2, keepdims = True)

    #define the criteria to treat the values as anamolous based on confidence interval
    #95% = 1.96, 98% = 2.326, 99% = 2.576
    cri_low = modis_lst_filled_mean-1.96*modis_lst_filled_std
    cri_high = modis_lst_filled_mean+1.96*modis_lst_filled_std

    #get the mask of positions which are invalid statistically
    invalid_stat = np.logical_or(modis_lst_filled_smoothed<cri_low,\
                                 modis_lst_filled_smoothed>cri_high)
    #convert the filled smoothed data to the same 4 d array
    #remember that the shape position changes while concatenating to make the time series
    modis_lst_filled_smoothed = modis_lst_filled_smoothed.reshape(\
                                data_shape[1], data_shape[2], data_shape[0], data_shape[3])
    #transpose it to get back to origina shape of input data
    modis_lst_filled_smoothed = np.transpose(modis_lst_filled_smoothed, (2,0,1,3))
    #do that for invalid dtat too
    invalid_stat = invalid_stat.reshape(\
                        data_shape[1], data_shape[2], data_shape[0], data_shape[3])
    invalid_stat = np.transpose(invalid_stat,(2,0,1,3))

    #update the modis bad msk
    modis_bad = np.logical_or(modis_qc_bad, invalid_stat)
    #now replace the values void with the smoothed data keeping others same
    modis_lst_filled [modis_bad] = modis_lst_filled_smoothed[modis_bad]
    return modis_lst_filled

#%%
#Now write the functions to be used for validation of the result using artificial gaps

def modis_bad_filter(filelist, qcbits, modis_lst_var, modis_qc_var, lst_error_tolerance):
    modis_qc_stacked = timeseries_stack(filelist, modis_qc_var)
    modis_lst_stacked = timeseries_stack(filelist, modis_lst_var)
    modis_qc_rank = qc_data2rank(modis_qc_stacked, qcbits)
    modis_qc_bad = np.logical_and(modis_qc_rank > lst_error_tolerance,\
                                  modis_lst_stacked == 0)
    return modis_qc_bad

def read_lst_filled(lst_filled_file):
    ds = Dataset(lst_filled_file)
    lst_data = ds['lst_filled'][:]
    ds = None
    return lst_data

#define a function to generate artificial gaps in the fille lst
#make sure that the voids are not in the border in each dimension or overlap real gaps
def generate_artificial_gap(lst_data, modis_qc_bad, frac):
    ww,xx,yy,zz = lst_data.shape
    size = ww*xx*yy*zz
    widx = np.random.randint(1, ww -1, round(frac*size)) #get the idxs along 0 axis
    xidx = np.random.randint(1, xx -1, round(frac*size))
    yidx = np.random.randint(1, yy -1, round(frac*size))
    zidx = np.random.randint(1, zz -1, round(frac*size))
    lst_data_void = lst_data.copy() #copy to make sure that the changes are not trasnferred to original data
    lst_data_void[widx,xidx, yidx, zidx] = np.nan
    lst_data_void[modis_qc_bad] = lst_data[modis_qc_bad] #restore the values in real voids
    #return widx,xidx, yidx, zidx
    return lst_data_void, widx, xidx, yidx, zidx

def era_modisvoid_difference_median(filelist, lst_data_void, era_file, modis_file):
    #skt_shape = era_skt_stacked.shape
    era_skt_var = 'era5land/skt'
    era_skt_stacked = timeseries_stack(filelist, era_skt_var)
    lst_shape = lst_data_void.shape
    row,col,W1,W2,W3,W4 = interpolate_params(era_file, modis_file)

    era_interpolated_stacked = np.zeros(lst_data_void.shape)
    for i in range(lst_shape[3]):
        era_skt_timeseries = era_skt_stacked[:,:,:,i]
        era_interpolated_timeseries = interpolate_timeseries(era_skt_timeseries,era_skt_var, row, col, W1, W2, W3, W4)
        era_interpolated_stacked[:,:,:,i] = era_interpolated_timeseries

    #calculate differnce between modis_lst_stacked with voids and era_skt_data
    lst_era_diff = lst_data_void - era_interpolated_stacked
    #then caculate the median of difference
    lst_era_diff_median_3d = np.median(lst_era_diff, axis = 3, keepdims = False)
    return lst_era_diff_median_3d
#%%
def fill_lst_void(lst_era_diff_median_filled_smoothed, lst_data_void, filelist, row, col, W1, W2, W3, W4, window_length, polyorder):
    era_skt_var = 'era5land/skt'
    time_var = 'era5land/time'

    era_skt_stacked = timeseries_stack(filelist, era_skt_var)
    time_stacked = timeseries_stack(filelist, time_var)
    data_shape = lst_data_void.shape

    #Then interpolate era5 data to the size of modis data
    era_interpolated_stacked = np.zeros(data_shape, dtype=np.float32)
    for i in range(len(filelist)):
        era_skt_timeseries = era_skt_stacked[:,:,:,i]
        era_interpolated_timeseries = interpolate_timeseries(era_skt_timeseries,era_skt_var, row, col, W1, W2, W3, W4)
        era_interpolated_stacked[:,:,:,i] = era_interpolated_timeseries

    del era_skt_stacked
    #era_interpolated_stacked_3d = np.concatenate(era_interpolated_stacked, axis = 2) #This concatenates in 3rd dimension with 3rd dimension as time

    #create a timeseries at daily basis
    #flatten in columns major to get continuous time series
    time_stacked_flat = time_stacked.flatten('F')
    time_stacked_monthly = time_stacked_flat.reshape(-1, 4)
    time_difference_monthly = np.diff(time_stacked_monthly, axis=0)
    time_stacked_monthly_modified = time_stacked_monthly.copy()
    for i in range(time_stacked_monthly.shape[0]-1):
        j = time_stacked_monthly.shape[0]-1-i
        time_stacked_monthly_modified[j:,:] = time_stacked_monthly_modified[j:,:] - time_difference_monthly[j-1,:] + 24
    timeseries = time_stacked_monthly_modified.flatten('F')

    #expand the dimension of reference difference to match the dimension of era data
    diff_ref_newaxis = np.expand_dims(lst_era_diff_median_filled_smoothed,axis = 3)
    del lst_era_diff_median_filled_smoothed
    diff_ref = np.concatenate([diff_ref_newaxis]*18, axis = 3)
    del diff_ref_newaxis
    modis_lst_reference  = era_interpolated_stacked + diff_ref
    del diff_ref
    modis_lst_filled = lst_data_void.copy()

    modis_lst_filled [np.isnan(lst_data_void)] = modis_lst_reference[np.isnan(lst_data_void)]
    del modis_lst_reference

    #smoothout the filled data to eliminate the sudden change in lsts
    #for that is is better to change the data to 3d
    modis_lst_filled_3d = np.concatenate(modis_lst_filled, axis = 2)

    modis_lst_filled_smoothed = np.zeros(modis_lst_filled_3d.shape)
    xx, yy, zz = modis_lst_filled_3d.shape
    for ii in tqdm(range(xx)):
        for jj in range(yy):
            lst_timeseries = modis_lst_filled_3d[ii,jj,:]
            lst_sg = non_uniform_savgol(timeseries, lst_timeseries, window_length, polyorder)
            modis_lst_filled_smoothed[ii,jj,:] = lst_sg

    #checking matlab script i realized tht ther may be some anamolous values undetected by quality flag
    modis_lst_filled_std = np.nanstd(modis_lst_filled_smoothed, axis =2, keepdims = True)
    modis_lst_filled_mean = np.nanmean(modis_lst_filled_smoothed, axis =2, keepdims = True)

    #define the criteria to treat the values as anamolous based on confidence interval
    #95% = 1.96, 98% = 2.326, 99% = 2.576
    cri_low = modis_lst_filled_mean-1.96*modis_lst_filled_std
    cri_high = modis_lst_filled_mean+1.96*modis_lst_filled_std

    #get the mask of positions which are invalid statistically
    invalid_stat = np.logical_or(modis_lst_filled_smoothed<cri_low,\
                                 modis_lst_filled_smoothed>cri_high)
    #convert the filled smoothed data to the same 4 d array
    #remember that the shape position changes while concatenating to make the time series
    modis_lst_filled_smoothed = modis_lst_filled_smoothed.reshape(\
                                data_shape[1], data_shape[2], data_shape[0], data_shape[3])
    #transpose it to get back to origina shape of input data
    modis_lst_filled_smoothed = np.transpose(modis_lst_filled_smoothed, (2,0,1,3))
    #do that for invalid dtat too
    invalid_stat = invalid_stat.reshape(\
                        data_shape[1], data_shape[2], data_shape[0], data_shape[3])
    invalid_stat = np.transpose(invalid_stat,(2,0,1,3))

    #update the modis bad msk
    modis_bad = np.logical_or(np.isnan(lst_data_void), invalid_stat)
    #now replace the values void with the smoothed data keeping others same
    modis_lst_filled [modis_bad] = modis_lst_filled_smoothed[modis_bad]
    return modis_lst_filled
#%%

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def me(predictions, targets):
    return (predictions - targets).mean()

def mae (predictions, targets):
    return np.mean(np.abs(predictions - targets))

def validation(lst_data, lst_data_filled, widx, xidx, yidx, zidx):
    lst_data_org = lst_data[widx, xidx, yidx, zidx]
    lst_data_regenerate = lst_data_filled[widx, xidx, yidx, zidx]

    value_rmse = rmse(lst_data_regenerate, lst_data_org)
    value_me = me(lst_data_regenerate, lst_data_org)
    value_mae = mae(lst_data_regenerate, lst_data_org)

    return value_rmse, value_me, value_mae
