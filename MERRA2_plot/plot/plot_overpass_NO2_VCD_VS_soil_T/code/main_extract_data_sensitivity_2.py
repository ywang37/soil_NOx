"""
Created on January 20, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import sys, getopt

from mylib.grid_utility import get_center_index_latlon
from mylib.io import read_nc, write_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
import sn_p

#######################
# Start user parameters
#

in_root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/\
MERRA2_plot/process/resample/data/sensitivity_2/daily/'

out_dir = '../data/sensitivity_2/'


verbose = True


# NO2 variable name
varn_list = [ \
        'mod_NO2Trop_AK_tp_sat_ori',
        'mod_NO2Trop_AK_tp_sat_ori_all_scale_0.50',
        'mod_NO2Trop_AK_tp_sat_soil_T_obs',
        'mod_NO2Trop_AK_tp_sat_soil_T_obs_all_scale_0.30',
        'mod_NO2Trop_AK_tp_sat_soil_T_obs_all_scale_0.40',
        'mod_NO2Trop_AK_tp_sat_soil_T_obs_all_scale_0.50',
        'sat_ColumnAmountNO2Trop',
        'mod_Met_TSOIL1',
        ]


#
# End user parameters
#####################


# agruements
start_year = 2005
end_year = 2019

start_month = 6
end_month = 8

region_name = 'Central_US'

argv = sys.argv[1:]

opts = 'sy:ey:sm:em:stt:ts:rn'
long_opts = ['start_year=', 'end_year=', 'start_month=', 'end_month=',
        'region_name=']

try:
    opts, args = getopt.getopt(argv, opts, long_opts)
except getopt.GetoptError:
    print('python main_calc_ave_emi.py -sy <start_year> ' +
            '-ey <end_year> -sm <start_month> -em <end_month> ')
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-sy', '--start_year'):
        start_year = int(arg)
    if opt in ('-ey', '--end_year'):
        end_year = int(arg)
    if opt in ('-sm', '--start_month'):
        start_month = int(arg)
    if opt in ('-em', '--end_month'):
        end_month = int(arg)
    if opt in ('-rn', '--region_name'):
        region_name = arg

# year and month
sy_c = str(start_year)
ey_c = str(end_year)
sm_c = str(start_month).zfill(2)
em_c = str(end_month).zfill(2)

month_days_dict = {}
month_days_dict['06'] = 30
month_days_dict['07'] = 31
month_days_dict['08'] = 31


# read data
file_count = 0
file_no_count = 0
file_no = []
coor_flag = True
coor_var = ['Latitude', 'Longitude', 'Latitude_e', 'Longitude_e']
data_dict = {}
for varn in varn_list:
    data_dict[varn] = []
for iyr in range(start_year, end_year+1):
    for imo in range(start_month, end_month+1):

        yyyy   = str(iyr)
        mm     = str(imo).zfill(2)
        yyyymm = yyyy + mm

        print('---- process ' + yyyymm  + ' ----')

        ndays = month_days_dict[mm]

        for iday in range(1, ndays+1):

            dd = str(iday).zfill(2)

            infile = in_root_dir + yyyy + '/' + \
                    'model_satellite_' + yyyy + '-' + mm + '-' + dd + '.nc'

            if os.path.exists(infile):

                file_count += 1

                if coor_flag:
                    region_limit = sn_p.region_dict[region_name]
                    coor_flag = False
                    indata = read_nc(infile, varnames=varn_list+coor_var,
                            verbose=True)
                    lat   = indata['Latitude']
                    lon   = indata['Longitude']
                    lat_e = indata['Latitude_e']
                    lon_e = indata['Longitude_e']
                    i1, i2, j1, j2 = \
                            get_center_index_latlon(lat_e[:,0], lon_e[0,:],
                                    region_limit)
                    i1 = i1 + 1
                    j1 = j1 + 1
                    lat = lat[i1:i2+1,j1:j2+1]
                    lon = lon[i1:i2+1,j1:j2+1]
                    lat_e = lat_e[i1:i2+2,j1:j2+2]
                    lon_e = lon_e[i1:i2+2,j1:j2+2]
                else:
                    indata = read_nc(infile, varnames=varn_list,
                            verbose=True)

                # save data
                for varn in varn_list:
                    data_dict[varn].append(indata[varn][i1:i2+1,j1:j2+1])

            else:

                file_no_count += 1
                file_no.append(infile)

for varn in varn_list:
    data_dict[varn] = np.array(data_dict[varn])


# prepare to output data
data_2D_dict = {}
data_2D_dict['Latitude']    = lat
data_2D_dict['Longitude']   = lon
data_2D_dict['Latitude_e']  = lat_e
data_2D_dict['Longitude_e'] = lon_e
outfile = out_dir + region_name + '_model_satellite_' + \
        str(start_year) + '-' + str(end_year) + '_' + \
        str(start_month).zfill(2) + '-' + str(end_month).zfill(2) +'.nc'
units_dict = {}
units_dict['mod_NO2Trop_AK_tp_sat_ori'] = 'molec/cm2'
#units_dict['mod_NO2Trop_AK_tp_sat_soil_T_ori'] = 'molec/cm2'
#units_dict['mod_NO2Trop_AK_tp_sat_surf_T_obs'] = 'molec/cm2'
units_dict['mod_NO2Trop_AK_tp_sat_soil_T_obs'] = 'molec/cm2'
units_dict['sat_ColumnAmountNO2Trop'] = 'molec/cm2'
units_dict['mod_Met_TSOIL1'] = 'K'
data_3D_time_dict = data_dict
write_nc(outfile, data_2D_dict, units_dict=units_dict,
        data_3D_time_dict=data_3D_time_dict)


print('--------- Done ----------')
print('start_year: {}'.format(start_year))
print('end_year: {}'.format(end_year))
print('start_month: {}'.format(start_month))
print('end_month: {}'.format(end_month))
print('{} files were read.'.format(file_count))
print('{} files were not found. They were'.format(file_no_count))
for i in range(len(file_no)):
    print(file_no[i] + '\n')
