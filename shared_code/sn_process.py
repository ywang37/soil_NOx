"""
Created on May 31, 2020

@author: Yi Wang
"""

from copy import deepcopy
import datetime
from netCDF4 import Dataset
import numpy as np
import os

from mylib.grid_utility import generate_grid_gc_2x25
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
