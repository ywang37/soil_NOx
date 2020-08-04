"""
Created on May 21, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import sys, getopt

from mylib.colormap.colormap_utility import truncate_colormap
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.io import read_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
import sn_p
from sn_plot import plot_panel_variables, layout_2

#######################
# Start user parameters
#

data_dir = '../data/'

fig_dir = '../figure/'

p_thre = 0.05


r_signi_vmin = 0.5
r_signi_vmax = 1.0
r_signi_min_valid = r_signi_vmin

vmax_diff = 0.3
vmin_diff = -vmax_diff

region_name = 'US'
region_limit = sn_p.region_dict[region_name]
xtick = (-120, -100, -80)
ytick = (30, 40, 50)

ocean_color = 'darkgrey'

varname_list = ['multimonth_mean', 'monthly',\
        'month_06', 'month_07', 'month_08']

#r_signi_cmap = deepcopy(gbcwpry_map)
#r_signi_cmap = plt.get_cmap('rainbow')
r_signi_cmap = truncate_colormap(WhGrYlRd_map,
        minval=0.05, maxval=1.0, n=200)

varn_list = ['mod_NO2Trop_AK_tp_sat_ori',
             'mod_NO2Trop_AK_tp_sat_soil_T_ori',
             'mod_NO2Trop_AK_tp_sat_surf_T_obs',
             'mod_NO2Trop_AK_tp_sat_soil_T_obs',
             'mod_NO2Trop_tp_sat_ori',
             'mod_NO2Trop_tp_sat_soil_T_ori',
             'mod_NO2Trop_tp_sat_surf_T_obs',
             'mod_NO2Trop_tp_sat_soil_T_obs']

ak_list = ['_AK', '']

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

title_list = \
        [
        'Control', \
        r'T$_{soil}$_S$_{old}$', \
        r'T$_{air}$_S$_{new}$', \
        r'T$_{soil}$_S$_{new}$', \
        ]

diff_title_list = \
        [
        'Control', \
        r'T$_{soil}$_S$_{old}$ - Control', \
        r'T$_{air}$_S$_{new}$ - Control', \
        r'T$_{soil}$_S$_{new}$ - Control', \
        ]

lw = 0.5

#
# End user parameters
#####################

time_range = '2005-2019_06-08'

#
argv = sys.argv[1:]
opts = 'tr'
long_opts = ['time_range=']
try:
    opts, args = getopt.getopt(argv, opts, long_opts)
except getopt.GetoptError:
    print('python main_plot_soil_T_OMI_L3_NO2_corr.py -tr <time_range> ')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-tr', '--time_range'):
        time_range = arg
#
for varname in varname_list:

    for ak in ak_list:

        first = True
        data_dict = {}
        data_dict[varname] = {}
        data_dict[varname+'_nan'] = {}
        data_dict[varname+'_diff_nan'] = {}
        data_dict['p'] = {}
        for scene in scene_tup:

            name = 'corr_overpass_sat_ColumnAmountNO2Trop_' + \
                    'mod_NO2Trop' + ak + '_tp_sat_' + scene  \
                    + '_' + time_range

            corr_file = data_dir + name + '.nc'

            if first:
 
                tmp_data = read_nc(corr_file, ['r_' + varname, 
                    'p_' + varname, 'Latitude_e', 'Longitude_e'])

                data_dict['lat_e'] = tmp_data['Latitude_e']
                data_dict['lon_e'] = tmp_data['Longitude_e']

            else:

               tmp_data = read_nc(corr_file, ['r_' + varname, 
                    'p_' + varname])

            first = False

            data_dict['p'][scene] = tmp_data['p_' + varname]
            data_dict[varname][scene] = tmp_data['r_' + varname]
            data_dict[varname+'_nan'][scene] = \
                    tmp_data['p_' + varname] >= p_thre

        # correlation
        plot_panel_variables(data_dict, scene_tup, '',
                varname=varname,
                read_func_varname=varname,
                vmin=r_signi_vmin,
                vmax=r_signi_vmax,
                valid_min=r_signi_min_valid,
                flag_nan=True,
                hspace=-0.3,
                ocean_color=ocean_color,
                xtick=xtick, ytick=ytick,
                region_limit=region_limit,
                title_list=title_list,
                cb1_extend='neither',
                lw=lw)
        figname = fig_dir + region_name + '_' + 'panel_' + varname  + '_'\
                'corr_overpass_sat_ColumnAmountNO2Trop_' + \
                'mod_NO2Trop' + ak + '_' + time_range[0:9] + '.png'
        plt.savefig(figname, format='png', dpi=300)

        # correlation difference
        diff_avail_flag = np.full_like(data_dict['p'][scene_tup[0]], False)
        for scene in scene_tup:
            diff_avail_flag = np.logical_or(diff_avail_flag,
                    data_dict['p'][scene] < p_thre)
        diff_nan_flag = np.logical_not(diff_avail_flag)
        for scene in scene_tup:
            data_dict[varname+'_diff_nan'][scene] = diff_nan_flag
        plot_panel_variables(data_dict, scene_tup, '',
                varname=varname,
                read_func_varname=varname,
                vmin_diff=vmin_diff, vmax_diff=vmax_diff,
                vmin=r_signi_vmin,
                vmax=r_signi_vmax,
                flag_nan=True,
                flag_diff_nan=True,
                hspace=-0,
                ocean_color=ocean_color,
                xtick=xtick, ytick=ytick,
                region_limit=region_limit,
                title_list=diff_title_list,
                cb2_extend='neither',
                diff=True,
                lw=lw)
        figname = fig_dir + region_name + '_' + 'panel_diff_' + varname  + '_'\
                'corr_overpass_sat_ColumnAmountNO2Trop_' + \
                'mod_NO2Trop' + ak + '_' + time_range[0:9] + '.png'
        plt.savefig(figname, format='png', dpi=300)




