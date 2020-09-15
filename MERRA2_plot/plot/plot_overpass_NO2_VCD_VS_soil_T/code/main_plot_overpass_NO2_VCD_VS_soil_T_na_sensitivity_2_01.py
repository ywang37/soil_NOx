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
from sn_plot import plot_NO2_VS_T

#######################
# Start user parameters
#


data_dir = '../data/na_sensitivity_2/'

fig_dir = '../figure/na_sensitivity_2/'

sensitivity = 'sensitivity_01'

verbose = True

old_labels = [ \
        'OMI', 
        'Control',
        r'T$_{air}$_S$_{old}$_s0.5',
        r'T$_{soil}$_S$_{new}$',
        r'T$_{soil}$_S$_{new}$_s0.50',
        r'T$_{soil}$_S$_{new}$_s0.40',
        r'T$_{soil}$_S$_{new}$_s0.30',
        ]

# NO2 variable name
varn_list = [ \
        'mod_NO2Trop_AK_tp_sat_ori',
        'mod_NO2Trop_AK_tp_sat_ori_all_scale_0.50',
        'mod_NO2Trop_AK_tp_sat_soil_T_obs',
        'mod_NO2Trop_AK_tp_sat_soil_T_obs_all_scale_0.50',
        'mod_NO2Trop_AK_tp_sat_soil_T_obs_all_scale_0.40',
        'mod_NO2Trop_AK_tp_sat_soil_T_obs_all_scale_0.30',
        'sat_ColumnAmountNO2Trop',
        'mod_Met_TSOIL1',
        ]

sat_var_list = ['sat_ColumnAmountNO2Trop']

mod_varname_str = 'mod_NO2Trop_AK_tp_sat_'

scene_tup = [ \
        'ori',
        'ori_all_scale_0.50',
        'soil_T_obs',
        'soil_T_obs_all_scale_0.50',
        'soil_T_obs_all_scale_0.40',
        'soil_T_obs_all_scale_0.30',
        ]

T_name = 'mod_Met_TSOIL1'

xticks = [-120.0, -100.0, -80.0] 
yticks = [30.0, 40.0, 50.0]

#
# End user parameters
#####################


# agruements
start_year = 2005
end_year = 2019

start_month = 6
end_month = 8

vmax = 3e15

region_name = 'Central_US'

soil_emi_ratio_thre = None
soil_emi_abs_thre = None

argv = sys.argv[1:]

opts = 'sy:ey:sm:em:vm:rn:ro:ab'
long_opts = ['start_year=', 'end_year=', 'start_month=', 'end_month=',
        'vmax=', 'region_name=', 'ratio=', 'abs=']

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
    if opt in ('-vm', '--vmax'):
        vmax = float(arg)
    if opt in ('-rn', '--region_name'):
        region_name = arg
    if opt in ('-ro', '--ratio'):
        soil_emi_ratio_thre = arg
    if opt in ('-ab', '--abs'):
        soil_emi_abs_thre = arg

# year and month
sy_c = str(start_year)
ey_c = str(end_year)
sm_c = str(start_month).zfill(2)
em_c = str(end_month).zfill(2)


title =  region_name + '_model_satellite_' + \
        str(start_year) + '-' + str(end_year) + '_' + \
        str(start_month).zfill(2) + '-' + str(end_month).zfill(2)
if soil_emi_ratio_thre is not None:
    title = title + '_ratio_' + soil_emi_ratio_thre
if soil_emi_abs_thre is not None:
    title = title + '_abs_' + soil_emi_abs_thre

# read data
infile = data_dir + title + '.nc'
coor_var = ['Latitude', 'Longitude', 'Latitude_e', 'Longitude_e']
data_dict = read_nc(infile, varnames=varn_list+coor_var, verbose=True)


# plot
map_region_limit = sn_p.region_dict['US']
if 'CA' in region_name:
    map_region_limit = [31.0, -126, 43.0, -113]
    xticks = np.arange(map_region_limit[1], map_region_limit[3]+0.1, 4)
    yticks = np.arange(map_region_limit[0], map_region_limit[2]+0.1, 2)
out_dict = plot_NO2_VS_T(data_dict, sat_var_list[0],
        mod_varname_str, scene_tup, T_name=T_name,
        xticks=xticks, yticks=yticks, 
        ax1_ylim=[0, vmax],
        old_labels=old_labels,
        map_region_limit=map_region_limit)
figname =  fig_dir +  title + '_NO2_VS_soil_T_' + sensitivity + '.png'
plt.savefig(figname, format='png', dpi=300)

# output data
df = out_dict['df']
df.to_csv(data_dir + 'sat_model_NO2_' + title + \
        '_NO2_VS_soil_T_' + sensitivity  + '.csv')
