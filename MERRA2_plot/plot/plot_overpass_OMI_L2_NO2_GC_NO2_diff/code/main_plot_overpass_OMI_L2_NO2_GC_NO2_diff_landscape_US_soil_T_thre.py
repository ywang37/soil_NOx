"""
Created on May 21, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import os
import sys, getopt

from mylib.colormap.colormap_utility import truncate_colormap
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.colormap.emphasize_small_map import emphasize_small_map
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.io import read_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
import sn_p
from sn_plot import plot_compare_4_to_1_new, layout_2

#######################
# Start user parameters
#

root_data_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/resample/data/'

root_fig_dir = '../figure/'

portrait = False

ocean_color = 'darkgrey'

month_list = ['multimonth_mean', 'month_06', 'month_07', 'month_08']

cb_ratio = 0.10
cb1_y_off = -0.08
cb_diff_y_off = -0.08

hspace = -0.3

cmap1 = deepcopy(WhGrYlRd_map)
vmin1 = 0.0
vmax1 = 3
units1 = r'[10$^{15}$ molec cm$^{-2}$]'

cmap_diff = deepcopy(gbcwpry_map)
vmin_diff = -0.8
vmax_diff = -vmin_diff
units_diff = r'[10$^{15}$ molec cm$^{-2}$]'
cb_ticks_diff = [-0.8, -0.4, 0.0, 0.4, 0.8]

verbose = True

ak_list = ['_AK', '']

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

title_dict = {}
title_dict['OMI'] = 'OMI'
title_dict['ori'] = 'Control'
title_dict['diff_ori'] = 'Control - OMI'
title_dict['diff_soil_T_ori'] = r'T$_{soil}$_S$_{old}$ - OMI'
title_dict['diff_surf_T_obs'] = r'T$_{air}$_S$_{new}$ - OMI'
title_dict['diff_soil_T_obs'] = r'T$_{soil}$_S$_{new}$ - OMI'


region_name = 'US'
region_limit = sn_p.region_dict[region_name]
xticks = (-120, -100, -80)
yticks = (30, 40, 50)

scale = 1e-15

lw = 0.5

#
# End user parameters
#####################

secne_dict = {}
secne_dict['above'] = 'larger'
secne_dict['below'] = 'smaller'

time_range = '2005-2019_06-08'

soil_T_thre = '27'

thre_scene = 'above'

#
argv = sys.argv[1:]
opts = 'tr'
long_opts = ['time_range=', 'soil_T_thre=', 'thre_scene=']
try:
    opts, args = getopt.getopt(argv, opts, long_opts)
except getopt.GetoptError:
    print('python main_plot_soil_T_OMI_L3_NO2_corr.py -tr <time_range> ' + 
            '-stt <soil_T_thre> -ts <thre_scene>')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-tr', '--time_range'):
        time_range = arg
    if opt in ('-stt', '--soil_T_thre'):
        soil_T_thre = arg
    if opt in ('-ts', '--thre_scene'):
        thre_scene = arg
#
for month in month_list:

    for ak in ak_list:

        sat_varname = 'sat_ColumnAmountNO2Trop_'+month+'_yearly_mean'
        varname = ['Latitude_e', 'Longitude_e', sat_varname]
        mod_varname_list = []
        for scene in scene_tup:

            mod_varname = 'mod_NO2Trop' + ak + '_tp_sat_' + scene \
                    + '_' + month + '_yearly_mean'
            varname.append(mod_varname)
            mod_varname_list.append(mod_varname)

        infile = root_data_dir + 'soil_T_' + secne_dict[thre_scene] + \
                '_' + soil_T_thre + '/multimonth_mean/' + \
                thre_scene + '_' + soil_T_thre + \
                '_model_satellite_' + time_range + '.nc'
 
        tmp_data = read_nc(infile, varname, verbose=verbose)

        data_dict = {}
        data_dict['Latitude_e'] = tmp_data['Latitude_e']
        data_dict['Longitude_e'] = tmp_data['Longitude_e']
        data_dict['OMI'] = tmp_data[sat_varname] * scale
        data_dict['ori'] = tmp_data[mod_varname_list[0]] * scale
        for i in range(len(scene_tup)):
            scene = scene_tup[i]
            mod_varname = mod_varname_list[i]
            data_dict['diff_'+scene] = (tmp_data[mod_varname] \
                    - tmp_data[sat_varname]) * scale

        plot_compare_4_to_1_new(data_dict,
                cmap1=cmap1, cmap_diff=cmap_diff,
                vmin1=vmin1, vmax1=vmax1, units1=units1,
                vmin_diff=vmin_diff, vmax_diff=vmax_diff, 
                units_diff=units_diff,
                title_dict=title_dict,
                portrait=portrait,
                cb_ratio=cb_ratio,
                cb1_y_off=cb1_y_off,
                cb_diff_y_off=cb_diff_y_off,
                hspace=hspace,
                xticks=xticks, yticks=yticks,
                region_limit=region_limit,
                )

        fig_dir = root_fig_dir + 'soil_T_' + secne_dict[thre_scene] + \
                '_' + soil_T_thre + '/'
        if not os.path.exists(fig_dir):
            os.system('mkdir -p ' + fig_dir)
        figname = fig_dir + thre_scene + '_' + soil_T_thre + \
                '_US_landscape_diff_' +  month  + \
                '_overpass_sat_ColumnAmountNO2Trop_' + \
                'mod_NO2Trop' + ak + '_' + time_range[0:9] + '.png'
        plt.savefig(figname, format='png', dpi=300)




