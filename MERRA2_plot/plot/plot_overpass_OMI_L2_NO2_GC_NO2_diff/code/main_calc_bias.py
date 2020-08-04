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
from mylib.colormap.emphasize_small_map import emphasize_small_map
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.grid_utility import region_ave_sum
from mylib.io import read_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
import sn_p
from sn_plot import plot_compare_4_to_1_new, layout_2

#######################
# Start user parameters
#

data_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/resample/data/multimonth_mean/'

fig_dir = '../figure/'

area_file = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_runs/GC_ori/\
runs/merra2_2x25_tropchem_spinup_20050601/OutputDir/\
HEMCO_diagnostics.200505160000.nc'

month_list = ['multimonth_mean', 'month_06', 'month_07', 'month_08']


cmap1 = deepcopy(emphasize_small_map)
vmin1 = 0.0
vmax1 = 8
units1 = r'[10$^{15}$ molec cm$^{-2}$]'

cmap_diff = deepcopy(gbcwpry_map)
vmin_diff = -2
vmax_diff = -vmin_diff
units_diff = r'[10$^{15}$ molec cm$^{-2}$]'


#ak_list = ['_AK', '']
ak_list = ['_AK']

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

title_dict = {}
title_dict['OMI'] = 'OMI'
title_dict['ori'] = 'Control'
title_dict['diff_ori'] = 'Control - OMI'
title_dict['diff_soil_T_ori'] = r'T$_{soil}$_S$_{old}$ - OMI'
title_dict['diff_surf_T_obs'] = r'T$_{air}$_S$_{new}$ - OMI'
title_dict['diff_soil_T_obs'] = r'T$_{soil}$_S$_{new}$ - OMI'

region_name_list = ['Central_US']

#scale = 1e-15

lw = 0.5

verbose = True

#
# End user parameters
#####################

# area
area_data = read_nc(area_file, ['AREA'], verbose=verbose)
area = area_data['AREA']


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
for month in month_list:

    sat_varname = 'sat_ColumnAmountNO2Trop_'+month+'_yearly_mean'
    varname = ['Latitude', 'Longitude', sat_varname]
    mod_varname_list = []
    for scene in scene_tup:

        mod_varname = 'mod_NO2Trop_AK_tp_sat_' + scene \
                + '_' + month + '_yearly_mean'
        varname.append(mod_varname)
        mod_varname_list.append(mod_varname)

    infile = data_dir + 'model_satellite_' + time_range + '.nc'
 
    tmp_data = read_nc(infile, varname, verbose=verbose)

    # latitude and longitude
    lat = tmp_data['Latitude']
    lon = tmp_data['Longitude']

    # loop region
    for region_name in region_name_list:

        print('month: ' + month + ', region: ' + region_name)

        region_limit = sn_p.region_dict[region_name]

        # get satellite data
        sat_data = deepcopy( tmp_data[sat_varname] )
        sat_region_dict = region_ave_sum(sat_data, in_weight=area, 
                lat=lat, lon=lon, region_limit=region_limit)
        sat_region_mean = sat_region_dict['mean']

        mod_region_mean_list = []
        for i in range(len(scene_tup)):
            scene = scene_tup[i]
            mod_varname = mod_varname_list[i]

            # get model data
            mod_data = deepcopy( tmp_data[mod_varname] )
            mod_region_dict = region_ave_sum(mod_data, in_weight=area,
                    lat=lat, lon=lon, region_limit=region_limit)
            mod_region_mean = mod_region_dict['mean']
            mod_region_mean_list.append(mod_region_mean)

        # print
        print('Source      Value      Bias      NMB')
        print('OMI        {:.2e}'.format(sat_region_mean))
        for i in range(len(scene_tup)):
            scene = scene_tup[i]
            mod_region_mean = mod_region_mean_list[i]
            print(scene.ljust(10) + ' {:.2e}  {:.2e}  {:.1f}%'.format(
                mod_region_mean, mod_region_mean-sat_region_mean,
                (mod_region_mean-sat_region_mean)/sat_region_mean*100.0))






