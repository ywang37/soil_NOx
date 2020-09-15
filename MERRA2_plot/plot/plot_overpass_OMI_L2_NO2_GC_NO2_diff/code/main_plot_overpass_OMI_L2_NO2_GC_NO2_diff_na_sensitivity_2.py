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
from mylib.plot_utility import plot_comparision

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
import sn_p
from sn_plot import plot_compare_4_to_1_new, layout_2

#######################
# Start user parameters
#

data_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/resample/data/na_sensitivity_2/monthly/'

fig_dir = '../figure/na_sensitivity_2/'


year_list = ['2009', '2012']

month_list = ['06', '07', '08']


vmin = 0.0
vmax = 3
units = r'[10$^{15}$ molec cm$^{-2}$]'

diff_vmin = -1.5
diff_vmax = -diff_vmin
cb_diff_ticks = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0,1.5]

region_name = 'US'
region_limit = sn_p.region_dict[region_name]
xticks = (-120, -100, -80)
yticks = (30, 40, 50)

scene_tup = ('soil_T_obs', 'ori',
        'ori_all_scale_0.50',
        'soil_T_obs_all_scale_0.30',
        'soil_T_obs_all_scale_0.40',
        'soil_T_obs_all_scale_0.50')


title_dict = {}
title_dict['ori'] = 'Control'
title_dict['ori_all_scale_0.50'] = 'Control * 0.5'
title_dict['soil_T_obs'] = r'T$_{soil}$_S$_{new}$'
title_dict['soil_T_obs_all_scale_0.30'] = r'T$_{soil}$_S$_{new}$ * 0.3'
title_dict['soil_T_obs_all_scale_0.40'] = r'T$_{soil}$_S$_{new}$ * 0.4'
title_dict['soil_T_obs_all_scale_0.50'] = r'T$_{soil}$_S$_{new}$ * 0.5'

scale = 1e-15

lw = 0.5

verbose = True

#
# End user parameters
#####################

month_day_dict = {}
month_day_dict['06'] = '30'
month_day_dict['07'] = '31'
month_day_dict['08'] = '31'

#
for year in year_list:

    for month in month_list:

        time_range = year + '-' + month + '-01_' + year+ '-' + month + '-' + \
                month_day_dict[month]

        sat_varname = 'sat_ColumnAmountNO2Trop'
        varname = ['Latitude_e', 'Longitude_e', sat_varname]
        mod_varname_list = []
        for scene in scene_tup:

            mod_varname = 'mod_NO2Trop_AK_tp_sat_' + scene
            varname.append(mod_varname)
            mod_varname_list.append(mod_varname)

        infile = data_dir + 'model_satellite_' + time_range + '.nc'
 
        tmp_data = read_nc(infile, varname, verbose=verbose)

        data_dict = {}
        #data_dict['Latitude_e'] = tmp_data['Latitude_e']
        #data_dict['Longitude_e'] = tmp_data['Longitude_e']
        lat_e = tmp_data['Latitude_e']
        lon_e = tmp_data['Longitude_e']
        data_dict['OMI'] = tmp_data[sat_varname] * scale
        data_dict['ori'] = tmp_data[mod_varname_list[0]] * scale
        for i in range(len(scene_tup)):

            scene = scene_tup[i]
            mod_varname = mod_varname_list[i]
            data_dict[scene] = tmp_data[mod_varname] * scale

            var1 = data_dict['OMI'] 
            var2 = data_dict[scene]

            title = {}
            title['var1'] = 'OMI'
            title['var2'] = title_dict[scene]
            title['diff'] = title_dict[scene] + '- OMI'
            plot_comparision(var1, var2, lat_e, lon_e,
                    vmin=vmin, vmax=vmax,
                    diff_vmin=diff_vmin, diff_vmax=diff_vmax,
                    region_limit=region_limit,
                    xticks=xticks,
                    yticks=yticks,
                    mask_ocean=True,
                    cb_diff_ticks=cb_diff_ticks,
                    title_dict=title,
                    units=units)

            if not os.path.exists(fig_dir):
                os.system('mkdir -p ' + fig_dir)
            figname = fig_dir + 'US_diff_overpass_sat_ColumnAmountNO2Trop_' + \
                    'mod_NO2Trop_AK_' + time_range + '_' + \
                    scene + '.png'
            plt.savefig(figname, format='png', dpi=300)




