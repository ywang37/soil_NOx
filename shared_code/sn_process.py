"""
Created on May 31, 2020

@author: Yi Wang
"""

from calendar import monthrange
from copy import deepcopy
import datetime
from netCDF4 import Dataset
import numpy as np
import os

from mylib.grid_utility import generate_grid_gc_2x25, get_2D_center_index
from mylib.io import read_nc, write_nc

#
#------------------------------------------------------------------------------
#
def daily_to_monthly_emissions(inDir, outDir, year, month, res='2x2.5'):
    """
    """

    if res == '2x2.5':
        lat_e, lon_e, lat_c, lon_c = generate_grid_gc_2x25()
    lon_e, lat_e = np.meshgrid(lon_e, lat_e)
    lon_c, lat_c = np.meshgrid(lon_c, lat_c)

    emi_3D_varnames = ['EmisNO2_Ship', 'EmisNO_BioBurn', 'EmisNO_Ship', \
            'EmisNO_Soil']

    emi_4D_varnames = ['EmisNO2_Anthro', 'EmisNO_Aircraft', \
            'EmisNO_Anthro', 'EmisNO_Lightning', 'EmisNO_Total']

    if inDir[-1] != '/':
        inDir = inDir + '/'

    if outDir[-1] != '/':
        outDir = outDir + '/'

    if not os.path.isdir(outDir):
        os.system('mkdir -p ' + outDir)

    yyyy   = str(year)
    mm     = str(month).zfill(2)
    yyyymm = yyyy + mm
    currDate_D = datetime.datetime.strptime(yyyymm + '01', '%Y%m%d')

    data_3D = {}
    for varn in emi_3D_varnames:
        data_3D[varn] = {}
        data_3D[varn]['data'] = []
        data_3D[varn]['units'] = ''

    data_4D = {}
    for varn in emi_4D_varnames:
        data_4D[varn] = {}
        data_4D[varn]['data'] = []
        data_4D[varn]['units'] = ''

    ##################
    # read daily data
    ##################
    first = True
    while True:

        # Date
        currDate = str(currDate_D)
        yyyy = currDate[0:4]
        mm   = currDate[5:7]
        dd   = currDate[8:10]
        yyyymmdd = yyyy + mm + dd

        infile = inDir + 'merra2_2x25_tropchem_' + yyyymm + '/OutputDir/' + \
                'HEMCO_diagnostics.' + yyyymmdd + '0000.nc'

        print(' - daily_to_monthly_emissions: reading ' + infile)

        # open infile
        nc_in = Dataset(infile, 'r')

        # 3D variables
        for varn in emi_3D_varnames:

            var_3D = nc_in[varn]

            data_3D[varn]['data'].append(var_3D[:])

            if first:
                data_3D[varn]['units'] = getattr(var_3D, 'units')

        # 4D variables
        for varn in emi_4D_varnames:

            var_4D = nc_in[varn]

            data_4D[varn]['data'].append(var_4D[:])

            if first:
                data_4D[varn]['units'] = getattr(var_4D, 'units')
        
        # close infile
        nc_in.close()

        first = False

        # Next day
        nextDate_D = currDate_D + datetime.timedelta(days=1)
        if (str(nextDate_D)[5:7] != mm):
            break
        currDate_D = nextDate_D

    #########################
    # calculate monthly mean
    #########################

    # 3D variables
    for varn in emi_3D_varnames:

        data_3D[varn]['data'] = np.array(data_3D[varn]['data'])
        data_3D[varn]['data'] = np.mean(data_3D[varn]['data'], axis=0)

    # 4D variables
    for varn in emi_4D_varnames:

        data_4D[varn]['data'] = np.array(data_4D[varn]['data'])
        data_4D[varn]['data'] = np.mean(data_4D[varn]['data'], axis=0)

    ##############
    # output data
    ##############

    # output file
    outfile = outDir + 'HEMCO_diagnostics.' + yyyymm + '.nc'

    print(' - daily_to_monthly_emissions: writing ' + outfile)

    # open infile
    nc_out = Dataset(outfile, 'w')

    # Dimensions of a netCDF file
    dim_lat = nc_out.createDimension('Latitude',  lat_c.shape[0])
    dim_lon = nc_out.createDimension('Longitude', lat_c.shape[1])
    dim_lat_e = nc_out.createDimension('Latitude_e',  lat_e.shape[0])
    dim_lon_e = nc_out.createDimension('Longitude_e', lat_e.shape[1])
    dim_time = nc_out.createDimension('time', 1)
    nlev = data_4D[emi_4D_varnames[0]]['data'].shape[1]
    dim_lev = nc_out.createDimension('lev', nlev)

    # create variables in a netCDF file

    # lat and lon
    Latitude_v = nc_out.createVariable('Latitude', 'f4',
            ('Latitude', 'Longitude'))
    Longitude_v = nc_out.createVariable('Longitude', 'f4',
            ('Latitude', 'Longitude'))
    Latitude_e_v = nc_out.createVariable('Latitude_e', 'f4',
            ('Latitude_e', 'Longitude_e'))
    Longitude_e_v = nc_out.createVariable('Longitude_e', 'f4',
            ('Latitude_e', 'Longitude_e'))

    # 3D variables
    out_var_3D_dict = {}
    for varn in emi_3D_varnames:
        out_var_3D = nc_out.createVariable(varn, 'f4', 
                ('time', 'Latitude', 'Longitude'))
        out_var_3D_dict[varn] = out_var_3D

    # 4D variables
    out_var_4D_dict = {}
    for varn in emi_4D_varnames:
        out_var_4D = nc_out.createVariable(varn, 'f4', 
                ('time', 'lev', 'Latitude', 'Longitude'))
        out_var_4D_dict[varn] = out_var_4D

    # write variables

    # lat and lon
    Latitude_v[:]    = lat_c
    Longitude_v[:]   = lon_c
    Latitude_e_v[:]  = lat_e
    Longitude_e_v[:] = lon_e

    # 3D variables
    for varn in emi_3D_varnames:
        out_var_3D_dict[varn][:] = data_3D[varn]['data']
        out_var_3D_dict[varn].units = data_3D[varn]['units']

    # 4D variables
    for varn in emi_4D_varnames:
        out_var_4D_dict[varn][:] = data_4D[varn]['data']
        out_var_4D_dict[varn].units = data_4D[varn]['units']

    nc_out.source = inDir

    # close file
    nc_out.close()
#
#------------------------------------------------------------------------------
#
def average_of_multi_month_emissions(inDir, outDir, 
        start_year, end_year,
        start_month, end_month, res='2x2.5',
        month_day_flag=True):
    """
    if month_day_flag is True, month_day_weight is considered.
    """

    if res == '2x2.5':
        lat_e, lon_e, lat_c, lon_c = generate_grid_gc_2x25()
    lon_e, lat_e = np.meshgrid(lon_e, lat_e)
    lon_c, lat_c = np.meshgrid(lon_c, lat_c)

    emi_3D_varnames = ['EmisNO2_Ship', 'EmisNO_BioBurn', 'EmisNO_Ship', \
            'EmisNO_Soil']

    emi_4D_varnames = ['EmisNO2_Anthro', 'EmisNO_Aircraft', \
            'EmisNO_Anthro', 'EmisNO_Lightning', 'EmisNO_Total']

    if inDir[-1] != '/':
        inDir = inDir + '/'

    if outDir[-1] != '/':
        outDir = outDir + '/'

    if not os.path.isdir(outDir):
        os.system('mkdir -p ' + outDir)

#    yyyy   = str(year)
#    mm     = str(month).zfill(2)
#    yyyymm = yyyy + mm
#    currDate_D = datetime.datetime.strptime(yyyymm + '01', '%Y%m%d')

    data_3D = {}
    for varn in emi_3D_varnames:
        data_3D[varn] = {}
        data_3D[varn]['data'] = []
        data_3D[varn]['data_yearly'] = []
        data_3D[varn]['units'] = ''

    data_4D = {}
    for varn in emi_4D_varnames:
        data_4D[varn] = {}
        data_4D[varn]['data'] = []
        data_4D[varn]['data_yearly'] = []
        data_4D[varn]['units'] = ''

    ##################
    # read daily data
    ##################
    first = True
    month_day_weight = [] # weight of every month every year
    month_day_weight_one_year = [] # weight of every month in one year
    first_year = True
    for year in range(start_year, end_year+1):

        # save monthly data in a year
        data_tmp = {}
        for varn in emi_3D_varnames:
            data_tmp[varn] = []
        for varn in emi_4D_varnames:
            data_tmp[varn] = []

        for month in range(start_month, end_month+1):

            # yyyymm
            yyyy = str(year)
            mm   = str(month).zfill(2)
            yyyymm = yyyy + mm

            # month_day_weight
            if month_day_flag:
                days = monthrange(year, month)[1]
                weight = float(days)
            else:
                weight = 1.0
            month_day_weight.append(weight)
            if first_year:
                month_day_weight_one_year.append(weight)

            infile = inDir + 'HEMCO_diagnostics.' + yyyymm + '.nc'

            print(' - average_of_multi_month_emissions: reading ' + infile)

            # open infile
            nc_in = Dataset(infile, 'r')

            # 3D variables
            for varn in emi_3D_varnames:

                var_3D = nc_in[varn]

                data_3D[varn]['data'].append(var_3D[:] * weight)

                data_tmp[varn].append(var_3D[:])

                if first:
                    data_3D[varn]['units'] = getattr(var_3D, 'units')

            # 4D variables
            for varn in emi_4D_varnames:

                var_4D = nc_in[varn]

                data_4D[varn]['data'].append(var_4D[:] * weight)

                data_tmp[varn].append(var_4D[:])

                if first:
                    data_4D[varn]['units'] = getattr(var_4D, 'units')
        
            # close infile
            nc_in.close()

            first = False

        first_year = False
        month_day_weight_one_year = np.array(month_day_weight_one_year)

        for varn in emi_3D_varnames + emi_4D_varnames:
            data_tmp[varn] = np.array(data_tmp[varn])
            for i in range(len(month_day_weight_one_year)):
                data_tmp[varn][i,:,:,:] *= month_day_weight_one_year[i]
            data_tmp[varn] = np.nansum(data_tmp[varn], axis=0) \
                    / np.sum(month_day_weight_one_year)

        # 3D variables
        for varn in emi_3D_varnames:
            data_3D[varn]['data_yearly'].append(data_tmp[varn])

        # 4D variables
        for varn in emi_4D_varnames:
            data_4D[varn]['data_yearly'].append(data_tmp[varn])


    ####################
    # yearly data
    ####################
    for varn in emi_3D_varnames:
        data_3D[varn]['data_yearly'] = np.array(data_3D[varn]['data_yearly'])
    for varn in emi_4D_varnames:
        data_4D[varn]['data_yearly'] = np.array(data_4D[varn]['data_yearly'])
    

    #########################
    # calculate monthly mean
    #########################

    # weight
    month_day_weight = np.array(month_day_weight)
    month_day_weight_mean = np.mean(month_day_weight)

    # 3D variables
    for varn in emi_3D_varnames:

        data_3D[varn]['data'] = np.array(data_3D[varn]['data'])
        data_3D[varn]['data'] = \
                np.mean(data_3D[varn]['data'], axis=0) \
                / month_day_weight_mean

    # 4D variables
    for varn in emi_4D_varnames:

        data_4D[varn]['data'] = np.array(data_4D[varn]['data'])
        data_4D[varn]['data'] = \
                np.mean(data_4D[varn]['data'], axis=0) \
                / month_day_weight_mean

    ##############
    # output data
    ##############

    # output file
    year_range = str(start_year) + '-' + str(end_year)
    month_range = str(start_month).zfill(2) + '-' \
            + str(end_month).zfill(2)
    time_range = year_range + '_' + month_range
    outfile = outDir + 'HEMCO_diagnostics.' + time_range + '.nc'

    print(' - average_of_multi_month_emissions: writing ' + outfile)

    # open infile
    nc_out = Dataset(outfile, 'w')

    # Dimensions of a netCDF file
    dim_lat = nc_out.createDimension('Latitude',  lat_c.shape[0])
    dim_lon = nc_out.createDimension('Longitude', lat_c.shape[1])
    dim_lat_e = nc_out.createDimension('Latitude_e',  lat_e.shape[0])
    dim_lon_e = nc_out.createDimension('Longitude_e', lat_e.shape[1])
    dim_time = nc_out.createDimension('time', 1)
    dim_time_year = nc_out.createDimension('time_year', 
            end_year - start_year + 1)
    nlev = data_4D[emi_4D_varnames[0]]['data'].shape[1]
    dim_lev = nc_out.createDimension('lev', nlev)

    # create variables in a netCDF file

    # lat and lon
    Latitude_v = nc_out.createVariable('Latitude', 'f4',
            ('Latitude', 'Longitude'))
    Longitude_v = nc_out.createVariable('Longitude', 'f4',
            ('Latitude', 'Longitude'))
    Latitude_e_v = nc_out.createVariable('Latitude_e', 'f4',
            ('Latitude_e', 'Longitude_e'))
    Longitude_e_v = nc_out.createVariable('Longitude_e', 'f4',
            ('Latitude_e', 'Longitude_e'))

    # 3D variables
    out_var_3D_dict = {}
    for varn in emi_3D_varnames:
        out_var_3D = nc_out.createVariable(varn, 'f4', 
                ('time', 'Latitude', 'Longitude'))
        out_var_3D_dict[varn] = out_var_3D
        out_var_3D_yearly = nc_out.createVariable(varn+'_yearly', 'f4',
                ('time_year', 'Latitude', 'Longitude'))
        out_var_3D_dict[varn+'_yearly'] = out_var_3D_yearly

    # 4D variables
    out_var_4D_dict = {}
    for varn in emi_4D_varnames:
        out_var_4D = nc_out.createVariable(varn, 'f4', 
                ('time', 'lev', 'Latitude', 'Longitude'))
        out_var_4D_dict[varn] = out_var_4D
        out_var_4D_yearly = nc_out.createVariable(varn+'_yearly', 'f4',
                ('time_year', 'lev', 'Latitude', 'Longitude'))
        out_var_4D_dict[varn+'_yearly'] = out_var_4D_yearly

    # write variables

    # lat and lon
    Latitude_v[:]    = lat_c
    Longitude_v[:]   = lon_c
    Latitude_e_v[:]  = lat_e
    Longitude_e_v[:] = lon_e

    # 3D variables
    for varn in emi_3D_varnames:
        out_var_3D_dict[varn][:] = data_3D[varn]['data']
        out_var_3D_dict[varn].units = data_3D[varn]['units']
        out_var_3D_dict[varn+'_yearly'][:] = data_3D[varn]['data_yearly']
        out_var_3D_dict[varn+'_yearly'].units = data_3D[varn]['units']

    # 4D variables
    for varn in emi_4D_varnames:
        out_var_4D_dict[varn][:] = data_4D[varn]['data']
        out_var_4D_dict[varn].units = data_4D[varn]['units']
        out_var_4D_dict[varn+'_yearly'][:] = data_4D[varn]['data_yearly']
        out_var_4D_dict[varn+'_yearly'].units = data_4D[varn]['units']

    nc_out.source = inDir

    # close file
    nc_out.close()
#
#------------------------------------------------------------------------------
#
def combine_NH_SH(north_file, south_file, combine_file):
    """ Put data of Northern Hemisphere from *north_file* and 
    data of Southern Hemisphere from *south_file* to 
    *combine_file*.
    (ywang, 06/03/2020)
    """

    month_list = ['06', '07', '08']

    month_NH_SH_dict = {}
    month_NH_SH_dict['06'] = '12'
    month_NH_SH_dict['07'] = '01'
    month_NH_SH_dict['08'] = '02'


    data_dict = {}
    units_dict = {}

    coord_varnames = ['Latitude', 'Latitude_e', 'Longitude', 'Longitude_e']

    #
    print(' - combine_NH_SH: read ' + north_file)
    nc_north = Dataset(north_file, 'r') 

    # get all variable name in *north_file*
    all_varnames = list(nc_north.variables.keys())
    varnames = deepcopy(all_varnames)
    for coord_var in coord_varnames:
        varnames.pop( varnames.index(coord_var) )

    # get units
    for varn in varnames:
        units_dict[varn] = getattr(nc_north.variables[varn], 'units')

    # get coordinates
    for coord_var in coord_varnames:
        data_dict[coord_var] = nc_north.variables[coord_var][:]

    # get *north_file* data
    north_data_dict = {}
    for varn in varnames:
        north_data_dict[varn] = nc_north.variables[varn][:]

    #
    nc_north.close()

    # get *south_file* data
    print(' - combine_NH_SH: read ' + south_file)
    nc_south = Dataset(south_file, 'r')
    south_data_dict = {}
    for varn in varnames:
        if varn[-2:] in month_list:
            south_varn = varn[:-2] + month_NH_SH_dict[varn[-2:]]
            south_data_dict[varn] = nc_south.variables[south_varn][:]
        else:
            south_data_dict[varn] = nc_south.variables[varn][:]
    nc_south.close()

    # index that seperate NH and SH
    ind = data_dict['Latitude'].shape[0] // 2

    # extract and combine data
    for varn in varnames:
        data_dict[varn] = np.full_like(data_dict['Latitude'], np.nan)
        # SH
        data_dict[varn][:ind,:] = south_data_dict[varn][:ind,:]
        # NH
        data_dict[varn][ind:,:] = north_data_dict[varn][ind:,:]

    # output data
    write_nc(combine_file, data_dict, units_dict)
#
#------------------------------------------------------------------------------
#
def read_file_list(filename_list, varname_list):
    """
    """

    data_dict = {}
    for varname in varname_list:
        data_dict[varname] = []

    for filename in filename_list:

        indata = read_nc(filename, varname_list, verbose=True, squeeze=True)

        for varname in varname_list:
            data_dict[varname].append(indata[varname])

    for varname in varname_list:
        data_dict[varname] = np.array()
#
#------------------------------------------------------------------------------
#
def get_grid_data(indata_dict, varname_list, grid_list, 
        lon_edge, lat_edge, verbose=True):
    """ get data in grids

    indata_dict:
        Elements are 3-D array (time, lat, lon)
    varname_list :
        Elements are keys for indata_dict
    grid_list :
        Elements in the grid_list is tuple. The tuple is 
        like (latitude, longitude)
    """

    # find all indexes
    ind_list = []
    for grid in grid_list:
        
        lat_value = grid[0]
        lon_value = grid[1]

        i, j = get_2D_center_index(lat_value, lat_edge, lon_value, lon_edge)
        ind_list.append( (i, j) )

        if verbose:

            print('--------------------------')
            print('lat={:}, lon={:}'.format(lat_value, lon_value))
            print('i={:}, j={:}'.format(i, j))
            print('lat_e1={:}, lat_e2={:}'.format(lat_edge[i], lat_edge[i+1]))
            print('lon_e1={:}, lon_e2={:}'.format(lon_edge[i], lon_edge[i+1]))

    # get data
    out_dict = {}
    for varn in varname_list:

        out_dict[varn] = []

        for ind in ind_list:

            i = ind[0]
            j = ind[1]
            out_dict[varn].append( indata_dict[varn][:,i,j] )

    return out_dict
#
#------------------------------------------------------------------------------
#
