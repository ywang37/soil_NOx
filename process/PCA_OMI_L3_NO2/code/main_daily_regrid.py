"""
Created on April 11, 2020

@author: Yi Wang
"""

import datetime
import glob
import numpy as np
import os
import sys

from mylib.grid_utility import generate_grid
from mylib.io import write_nc
from mylib.pro_omi_no2_l3.io_omi_no2_l3 import read_OMI_NO2_L3
from mylib.pro_satellite.regrid import drop_in_the_box

#######################
# Start user parameters
#

# usage:
# python main_daily_regrid.py startDate endDate
# startDate and endDate formats are like YYYY-MM-DD

startDate = sys.argv[1]
endDate   = sys.argv[2]

inRootDir = '/Dedicated/jwang-data/ywang/OMI_NO2/level3/'
outRootDir = '/Dedicated/jwang-data/ywang/soil_NOx/process/PCA/data/daily/'

min_val = -1e-30

res = '2x25'

# max longitude edge for GC 2x2.5 grid
max_lat = 178.75

# half grid box for GC 2x2.5 grid at polar
half_grid = True

#
# End user parameters
#####################

grid_dict = {}
grid_dict['proj']  = 'PlateCarree'
if res == '2x25':
    grid_dict['lat_start'] =  -91.0
    grid_dict['lat_end']   =   91.0
    grid_dict['lat_step']  =    2.0
    grid_dict['lon_start'] = -181.25
    grid_dict['lon_end']   =  178.75
    grid_dict['lon_step']  =    2.5

currDate  = startDate
currDate_D = datetime.datetime.strptime(currDate, '%Y-%m-%d')
endDate_D  = datetime.datetime.strptime(endDate,  '%Y-%m-%d')


nlat_c, nlon_c = 720, 1440
lat_e, lon_e, lat_c, lon_c = generate_grid(nlat_c, nlon_c)
lon, lat = np.meshgrid(lon_c, lat_c)
lon_flag = np.logical_and(lon>max_lat, lon<=180.0)
lon[lon_flag] = lon[lon_flag] - 360.0

while currDate_D <= endDate_D:

    currDate = str(currDate_D)[0:10]
    print('----------------------------------------------------')
    print('processing ' + currDate)

    currDir = inRootDir + '/' + currDate[0:4] + '_ori'

    curr_files = glob.glob(currDir + '/OMI-Aura_L3-OMNO2d_' \
            + currDate[0:4] + 'm' + currDate[5:7] + currDate[8:10] \
            + '_v003*.he5')

    curr_files.sort()

    all_files = list(curr_files)
    if len(all_files) > 1:
        print('ERROR')
        exit()

    if (len(all_files) == 1):
        filename = all_files[0]

        # read data
        OMI_data = read_OMI_NO2_L3(filename)
        NO2_Trop_CS = OMI_data['ColumnAmountNO2TropCloudScreened']

        # regrid data
        sat_flag = NO2_Trop_CS > min_val
        regrid_data = drop_in_the_box(grid_dict, lat, lon, NO2_Trop_CS,
                sat_flag=sat_flag)

        # change variable names
        regrid_data['NO2_Trop_CS'] = regrid_data.pop('sat_grid')

        # units
        units_dict = {}
        units_dict['NO2_Trop_CS'] = 'molec/cm2'

        # half grid
        if half_grid:
            regrid_data['Latitude'][0,:]    = -89.5
            regrid_data['Latitude'][-1,:]   = 89.5
            regrid_data['Latitude_e'][0,:]  = -90.0
            regrid_data['Latitude_e'][-1,:] = 90.0

        # output data
        out_dir = outRootDir + currDate[0:4] + '/'
        if not os.path.isdir(out_dir):
            os.system('mkdir -p ' + out_dir)
        out_file = out_dir + 'OMI_NO2_' + str(currDate)[0:10] \
                + '_' + res + '.nc'
        write_nc(out_file, regrid_data, units_dict=units_dict)
        

    else:
        print('No files at current day')


    # go to next day
    currDate_D = currDate_D + datetime.timedelta(days=1)









