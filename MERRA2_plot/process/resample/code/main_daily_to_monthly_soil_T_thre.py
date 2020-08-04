"""
Created on January 13, 2020

@author: Yi Wang
"""

import datetime
from copy import deepcopy
import glob
import numpy as np
import os 
import sys

from mylib.io import read_nc

sys.path.append("/Dedicated/jwang-data/ywang/soil_NOx/shared_code/")
from sn_io import save_ave

#######################
# Start user parameters
#

# usage:
# python main_daily_to_monthly_soil_T.py startDate endDate
# startDate and endDate formats are like YYYY-MM-DD
# startDate and endDate must be in the same month

startDate = sys.argv[1]
endDate   = sys.argv[2]

# celsius degree (int)
soil_T_thre = int(sys.argv[3])

#startDate = '2018-06-01'
#endDate   = '2018-06-30'

scene_list = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

tp_list = ['tp_sat']

root_daily_dir = '../data/daily/'

above_monthly_dir = '../data/soil_T_larger_{:}/monthly/'.format(soil_T_thre)

below_monthly_dir = '../data/soil_T_smaller_{:}/monthly/'.format(soil_T_thre)

# model species name
mod_spe_name = 'mod_NO2Trop'

# model AMF name
mod_amf_name = 'mod_AmfTrop'

# number var name
num_name = 'count'

# some variables
varname_list = [ \
        num_name, \
        'sat_AmfTrop', \
        'sat_ColumnAmountNO2Trop', \
        'sat_TropopausePressure', \
        'mod_Met_GWETTOP', \
        'mod_Met_TS', \
        'mod_Met_TSOIL1', \
        ]

coord_name_list = [ 
        'Latitude',
        'Latitude_e',
        'Longitude',
        'Longitude_e'
        ]


count_thre = 1

#
# End user parameters
#####################

# add model variable names of difference scenarios.
for scene in scene_list:
    for tp in tp_list:
        
        tmp_name = '_' + tp + '_' + scene
        # species
        varname_list.append( mod_spe_name + tmp_name )
        varname_list.append( mod_spe_name + '_AK' + tmp_name )
        # amf
        varname_list.append( mod_amf_name + tmp_name )

# Date
currDate   = startDate
currDate_D = datetime.datetime.strptime(currDate, '%Y-%m-%d')
endDate_D  = datetime.datetime.strptime(endDate,  '%Y-%m-%d')

above_out_dict = {}
below_out_dict = {}
latlon_flag = True
above_in_data_list = []
below_in_data_list = []
above_in_count_list = []
below_in_count_list = []
while currDate_D <= endDate_D:

    # current date
    currDate = str(currDate_D)[0:10]
    print(''.join(np.full((79,), '-')))
    print('processing ' + currDate)

    # read data
    daily_dir = root_daily_dir + currDate[0:4] + '/'
    in_filename = daily_dir + 'model_satellite_' + \
            currDate + '.nc'
    if not os.path.exists(in_filename):
        # go to next day
        print('  File does not exist.')
        currDate_D = currDate_D + datetime.timedelta(days=1)
        continue
    if latlon_flag:
        latlon_flag = False
        in_data = read_nc(in_filename, varname_list + coord_name_list, 
                verbose=True)
        for varname in coord_name_list:
            above_out_dict[varname] = in_data.pop(varname)
            below_out_dict[varname] = deepcopy(above_out_dict[varname])
    else:
        in_data = read_nc(in_filename, varname_list, verbose=True)

    # copy data
    above_in_data = deepcopy(in_data)
    below_in_data = deepcopy(in_data)

    # copy mod_Met_TSOIL1 data
    mod_Met_TSOIL1 = in_data['mod_Met_TSOIL1'] - 273.15

    # flag for mod_Met_TSOIL1 less than soil_T_thre
    # if mod_Met_TSOIL1 is NAN, below_flag is set True
    below_flag = (mod_Met_TSOIL1 < soil_T_thre)
    below_flag = np.logical_or(below_flag, np.isnan(mod_Met_TSOIL1))
    below_flag = np.logical_or(below_flag, in_data['mod_Met_TSOIL1']<1e-6)

    # flag for mod_Met_TSOIL1 larger than or equal soil_T_thre
    # if mod_Met_TSOIL1 is NAN, above_flag is set True
    above_flag = (mod_Met_TSOIL1 >= soil_T_thre)
    above_flag = np.logical_or(above_flag, np.isnan(mod_Met_TSOIL1))
    above_flag = np.logical_or(above_flag, in_data['mod_Met_TSOIL1']<1e-6)

    for varname in varname_list:
        
        # For above data, NAN or zero if below_flag is True
        # For below data, NAN or zero if above_flag is True
        if (varname != num_name):

            above_in_data[varname][below_flag] = np.nan

            below_in_data[varname][above_flag] = np.nan

        else:

            above_in_data[varname][below_flag] = 0

            below_in_data[varname][above_flag] = 0

    # for above data
    above_in_data_list.append(above_in_data)
    above_in_count_list.append(above_in_data[num_name])

    # for below data
    below_in_data_list.append(below_in_data)
    below_in_count_list.append(below_in_data[num_name])


    # go to next day
    currDate_D = currDate_D + datetime.timedelta(days=1)

# for above data
above_in_count_sum = np.array(above_in_count_list, dtype=int)
above_in_count_sum = (above_in_count_sum > 0)
above_in_count_sum = np.sum(above_in_count_sum, axis=0)
above_out_dict['count'] = above_in_count_sum

# for below data
below_in_count_sum = np.array(below_in_count_list, dtype=int)
below_in_count_sum = (below_in_count_sum > 0)
below_in_count_sum = np.sum(below_in_count_sum, axis=0)
below_out_dict['count'] = below_in_count_sum

# process daily data to month data
for varname in varname_list:

    print('  processing ' + varname)
    above_in_var_list = []
    below_in_var_list = []

    if (varname != num_name):

        for i in range(len(above_in_data_list)):

            above_in_var = above_in_data_list[i][varname]
            above_in_var[above_in_count_sum<count_thre] = np.nan

            below_in_var = below_in_data_list[i][varname]
            below_in_var[below_in_count_sum<count_thre] = np.nan

            # only use data when there are sufficient
            # pixels in grid

            above_in_var_list.append(above_in_var)

            below_in_var_list.append(below_in_var)

        above_in_var_all = np.array(above_in_var_list)

        below_in_var_all = np.array(below_in_var_list)

        # calculate monthly

        above_in_var_ave = np.nanmean(above_in_var_all, axis=0)
        above_out_dict[varname] = above_in_var_ave 

        below_in_var_ave = np.nanmean(below_in_var_all, axis=0)
        below_out_dict[varname] = below_in_var_ave

# output data
if not os.path.exists(above_monthly_dir):
    os.system('mkdir -p ' + above_monthly_dir)
above_out_filename = above_monthly_dir + \
        'above_{:}_model_satellite_'.format(soil_T_thre) + \
        startDate + '_' + endDate + '.nc'
save_ave(above_out_filename, above_out_dict)
if not os.path.exists(below_monthly_dir):
    os.system('mkdir -p ' + below_monthly_dir)
below_out_filename = below_monthly_dir + \
        'below_{:}_model_satellite_'.format(soil_T_thre) + \
        startDate + '_' + endDate + '.nc'
save_ave(below_out_filename, below_out_dict)



