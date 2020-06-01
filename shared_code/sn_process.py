"""
Created on May 31, 2020

@author: Yi Wang
"""

import datetime
from netCDF4 import Dataset
import numpy as np
import os

from mylib.grid_utility import generate_grid_gc_2x25

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
