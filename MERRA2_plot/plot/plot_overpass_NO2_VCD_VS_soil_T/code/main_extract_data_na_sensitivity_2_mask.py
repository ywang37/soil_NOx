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
#import sn_p

#######################
# Start user parameters
#

in_root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/\
MERRA2_plot/process/resample/data/na_sensitivity_2/daily/'

out_dir = '../data/na_sensitivity_2/'


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

region_mask_file_dict = {}
region_mask_file_dict['CA'] = '/Dedicated/jwang-data/ywang/soil_NOx/\
MERRA2_plot/process/select_grids/data/CA_flag_05x0625.nc'

soil_emi_ratio_file = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
plot/plot_soil_NOx_ratio/data/soil_ratio_2009_and_2012.nc'

#
# End user parameters
#####################


# agruements
start_year = 2005
end_year = 2019

start_month = 6
end_month = 8

region_name = 'CA'

soil_emi_ratio_thre = '0.0'
soil_emi_abs_thre = '0' # ng N m-2 s-1

argv = sys.argv[1:]


opts = 'sy:ey:sm:em:stt:ts:rn:ro:ab'
long_opts = ['start_year=', 'end_year=', 'start_month=', 'end_month=',
        'region_name=', 'ratio=', 'abs=']

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
    if opt in ('-ro', '--ratio'):
        soil_emi_ratio_thre = arg
    if opt in ('-ro', '--abs'):
        soil_emi_abs_thre = arg

# year and month
sy_c = str(start_year)
ey_c = str(end_year)
sm_c = str(start_month).zfill(2)
em_c = str(end_month).zfill(2)

month_days_dict = {}
month_days_dict['06'] = 30
month_days_dict['07'] = 31
month_days_dict['08'] = 31

# flag mask
region_mask_file = region_mask_file_dict[region_name]
flag_mask = read_nc(region_mask_file, ['flag'], verbose=verbose)
flag_mask = flag_mask['flag']
flag_mask = flag_mask < 0.5

# soil emission ratio
soil_emi = read_nc(soil_emi_ratio_file, ['EmisNO_Soil_ratio',
    'EmisNO_Soil'],
        verbose=verbose)
soil_emi_ratio = soil_emi['EmisNO_Soil_ratio']
soil_emi_abs = soil_emi['EmisNO_Soil']
flag_mask = np.logical_or(flag_mask, 
        soil_emi_ratio < float(soil_emi_ratio_thre))
flag_mask = np.logical_or(flag_mask, 
        soil_emi_abs < float(soil_emi_abs_thre))

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
                    coor_flag = False
                    indata = read_nc(infile, varnames=varn_list+coor_var,
                            verbose=True)
                    lat   = indata['Latitude']
                    lon   = indata['Longitude']
                    lat_e = indata['Latitude_e']
                    lon_e = indata['Longitude_e']
                else:
                    indata = read_nc(infile, varnames=varn_list,
                            verbose=True)

                # save data
                for varn in varn_list:
                    tmp = indata[varn]
                    tmp[flag_mask] = np.nan
                    data_dict[varn].append(tmp)

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
        str(start_month).zfill(2) + '-' + str(end_month).zfill(2) + \
        '_ratio_' + soil_emi_ratio_thre + '_abs_' + soil_emi_abs_thre + '.nc'
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
