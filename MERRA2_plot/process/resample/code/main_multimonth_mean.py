"""
Created on June 28, 2020

@author: Yi Wang
"""

from calendar import monthrange
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import sys, getopt
from tqdm import tqdm

from mylib.grid_utility import generate_grid_gc_2x25, get_center_index_latlon
from mylib.io import read_nc, write_nc
from mylib.pca.plot import plot_2D_components_and_coeffs
from mylib.trend_analysis import trend_analysis

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_io import read_NO_emissions
import sn_p
from sn_plot import plot_ave_series

#######################
# Start user parameters
#

in_root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/\
MERRA2_plot/process/resample/data/'

out_root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/\
MERRA2_plot/process/resample/data/'

verbose = True


# NO2 variable name
varn_list = [ 
             'mod_NO2Trop_AK_tp_sat_ori',
             'mod_NO2Trop_AK_tp_sat_soil_T_ori',
             'mod_NO2Trop_AK_tp_sat_surf_T_obs',
             'mod_NO2Trop_AK_tp_sat_soil_T_obs',
             'mod_NO2Trop_tp_sat_ori',
             'mod_NO2Trop_tp_sat_soil_T_ori',
             'mod_NO2Trop_tp_sat_surf_T_obs',
             'mod_NO2Trop_tp_sat_soil_T_obs',
             'sat_ColumnAmountNO2Trop',
             'mod_Met_GWETTOP', 
             'mod_Met_TS', 
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



argv = sys.argv[1:]

opts = 'sy:ey:sm:em:stt,ts'
long_opts = ['start_year=', 'end_year=', 'start_month=', 'end_month=',
        'soil_T_thre', 'thre_scene']

try:
    opts, args = getopt.getopt(argv, opts, long_opts)
except getopt.GetoptError:
    print('python main_calc_ave_emi.py -sy <start_year> ' + 
            '-ey <end_year> -sm <start_month> -em <end_month> ' +
            '-stt <soil_T_thre> -ts <thre_scene>')
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
    if opt in ('-stt', '--soil_T_thre'):
        soil_T_thre = arg
    if opt in ('-ts', '--thre_scene'):
        thre_scene = arg

# year and month
sy_c = str(start_year)
ey_c = str(end_year)
sm_c = str(start_month).zfill(2)
em_c = str(end_month).zfill(2)



month_days_dict = {}
month_days_dict['06'] = '30'
month_days_dict['07'] = '31'
month_days_dict['08'] = '31'


all_var_data = {}
for varn in varn_list:
    all_var_data[varn] = {}
    all_var_data[varn]['multimonth_mean'] = []
    for imo in range(start_month, end_month+1):
        mm = str(imo).zfill(2)
        all_var_data[varn]['month_'+mm] = []

coor_flag = True
coor_var = ['Latitude', 'Longitude', 'Latitude_e', 'Longitude_e']
for iyr in range(start_year, end_year+1):

    tmp_data = {}
    for varn in varn_list:
        tmp_data[varn] = []
    for imo in range(start_month, end_month+1):

        yyyy   = str(iyr)
        mm     = str(imo).zfill(2)
        yyyymm = yyyy + mm

        print('---- process ' + yyyymm  + ' ----')

        infile = in_dir + 'model_satellite_' + \
                yyyy + '-' + mm + '-01_' + \
                yyyy + '-' + mm + '-' + month_days_dict[mm] + '.nc'
        if coor_flag:
            coor_flag = False
            indata = read_nc(infile, varnames=varn_list+coor_var)
            lat   = indata['Latitude']
            lon   = indata['Longitude']
            lat_e = indata['Latitude_e']
            lon_e = indata['Longitude_e']
        else:
            indata = read_nc(infile, varnames=varn_list)

        for varn in varn_list:

            # varible value
            value = indata[varn]

            # variable dict
            var_data = all_var_data[varn]

            # save one-month data of every year
            var_data['month_'+mm].append(value)

            # save all-month data in a year
            tmp_data[varn].append(value)

    # multimonth mean in a year
    for varn in varn_list:
        multimonth_mean = np.mean(np.array(tmp_data[varn]), axis=0)
        all_var_data[varn]['multimonth_mean'].append(multimonth_mean) 


# multiyear mean
for varn in varn_list:

    # for every month
    for imo in range(start_month, end_month+1):

        mm     = str(imo).zfill(2)

        # Every month data
        all_var_data[varn]['month_'+mm] = \
                np.array(all_var_data[varn]['month_'+mm])

        # mean
        all_var_data[varn]['month_'+mm+'_yearly_mean'] = \
                np.mean(all_var_data[varn]['month_'+mm], axis=0)

    # for multimonth
    all_var_data[varn]['multimonth_mean'] = \
            np.array(all_var_data[varn]['multimonth_mean'])
    all_var_data[varn]['multimonth_mean_yearly_mean'] = \
            np.mean(all_var_data[varn]['multimonth_mean'], axis=0)



####################
# Start saving data
####################

# prepare coordinate data
coor_dict = {}
coor_dict['Latitude']    = lat
coor_dict['Longitude']   = lon
coor_dict['Latitude_e']  = lat_e
coor_dict['Longitude_e'] = lon_e

# prepare output data
out_dict = {}
data_3D_time_dict = {}
for varn in varn_list:

    out_dict[varn+'_multimonth_mean_yearly_mean'] = \
            all_var_data[varn]['multimonth_mean_yearly_mean']

    data_3D_time_dict[varn+'_multimonth_mean'] = \
            all_var_data[varn]['multimonth_mean']

    for imo in range(start_month, end_month+1):

        mm = str(imo).zfill(2)

        out_dict[varn+'_month_'+mm+'_yearly_mean'] = \
                all_var_data[varn]['month_'+mm+'_yearly_mean']

        data_3D_time_dict[varn+'_month_'+mm] = \
                all_var_data[varn]['month_'+mm]
data_3D_time_dict['time'] = np.array(range(start_year, end_year+1), 
        dtype='int')

# add coordinate
out_dict.update(coor_dict)

units_dict = {}
all_nc_varn = list(out_dict.keys()) + list(data_3D_time_dict.keys())
for varn in all_nc_varn:
    if ('mod_NO2Trop' in varn) or ('sat_ColumnAmountNO2Trop' in varn):
        units_dict[varn] = 'molec/cm2'
    if 'mod_Met_GWETTOP' in varn:
        units_dict[varn] = 'unitless'
    if ('mod_Met_TS' in varn) or ('mod_Met_TSOIL1' in varn):
        units_dict[varn] = 'K'

# output
outfile = out_dir + 'model_satellite_' + \
        sy_c + '-' + ey_c + '_' + sm_c + '-' + em_c + '.nc'
write_nc(outfile, out_dict, units_dict=units_dict, 
        data_3D_time_dict=data_3D_time_dict, verbose=True)
