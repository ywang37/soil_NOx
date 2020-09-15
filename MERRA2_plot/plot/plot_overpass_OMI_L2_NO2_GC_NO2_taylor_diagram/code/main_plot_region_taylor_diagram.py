"""
Created on June 30, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import os
import sys, getopt

sys.path.append('/Dedicated/jwang-data/ywang/opt/anaconda3/lib/python3.7/site-packages')
from mylib.cartopy_plot import cartopy_plot
from mylib.colormap.colormap_utility import truncate_colormap
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.colormap.emphasize_small_map import emphasize_small_map
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.grid_utility import region_ave_sum
from mylib.grid_utility import get_center_index_latlon
from mylib.io import read_nc
from mylib.modified_taylor_diagrams import ModTaylorDiagram
from mylib.trend_analysis import trend_analysis

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
import sn_p
from sn_plot import plot_time_series

#######################
# Start user parameters
#

start_year = 2005
end_year   = 2019

xx = np.array(range(start_year,end_year+1))

time_range = str(start_year) + '-' + str(end_year) + '_06-08'

NO2_VCD_file = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/resample/data/multimonth_mean/\
model_satellite_' + time_range + '.nc'

#area_file = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_runs/GC_ori/\
#runs/merra2_2x25_tropchem_spinup_20050601/OutputDir/\
#HEMCO_diagnostics.200505160000.nc'

roor_fig_dir = '../figure/'

month_list = ['multimonth_mean', 'month_06', 'month_07', 'month_08']

verbose = True

ak_list = ['_AK', '']

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

model_label_dict = {}
model_label_dict['ori'] = 'Control'
model_label_dict['soil_T_ori'] = r'T$_{soil}$_S$_{old}$'
model_label_dict['surf_T_obs'] = r'T$_{air}$_S$_{new}$'
model_label_dict['soil_T_obs'] = r'T$_{soil}$_S$_{new}$'


region_name_list = ['Central_US', 'Eastern_US', 'US', 'CA_NV']

# whether or not plot region selection
region_flag_dict = {}
region_flag_dict['Central_US'] = True
region_flag_dict['Eastern_US'] = True
region_flag_dict['US'] = True
region_flag_dict['CA_NV'] = True

xlabel = 'Year'

ylabel = r'NO$_2$ VCD [molec cm$^{-2}$]'

ratio_ylabel = 'Percentage relative to ' + str(start_year) + ' [%]'

percentage_scale = 100.0

ylim = [0.0, 3.0e15]

ratio_ylim = [60.0, 140.0]

US_flag_file = '/Dedicated/jwang-data/ywang/opt/anaconda3/lib/python3.7/\
site-packages/mylib/old/pixel_in_contiguous_US/US_flag_2x25.nc'

#
# End user parameters
#####################

## area
#tmp_data = read_nc(area_file, ['AREA'], verbose=verbose)
#area = tmp_data['AREA']


# US flag mask
US_flag_mask = read_nc(US_flag_file, ['flag'], verbose=verbose)
US_flag_mask = US_flag_mask['flag']
US_flag_mask = US_flag_mask < 0.5

#
for month in month_list:

    # NO2 VCD
    sat_NO2_name = 'sat_ColumnAmountNO2Trop_' + month + \
            '_yearly_mean'
    for ak in ak_list:

        overpass_varname_list = ['Latitude_e', 'Longitude_e', 
                'Latitude', 'Longitude', sat_NO2_name]
        NO2_VCD_name_list = [sat_NO2_name]
        model_NO2_VCD_name_list = []
        label_dict = {}
        label_dict[sat_NO2_name] = 'OMI'
        for scene in scene_tup:

            mod_NO2_varname = 'mod_NO2Trop' + ak + '_tp_sat_' + scene \
                    + '_' + month + '_yearly_mean'
            overpass_varname_list.append(mod_NO2_varname)
            NO2_VCD_name_list.append(mod_NO2_varname)
            model_NO2_VCD_name_list.append(mod_NO2_varname)
            label_dict[mod_NO2_varname] = \
                    model_label_dict.get(scene)
 
        NO2_data_dict = \
                read_nc(NO2_VCD_file, overpass_varname_list, verbose=verbose)

        # mask outside US
        for NO2_VCD_name in NO2_VCD_name_list:
            NO2_data_dict[NO2_VCD_name][US_flag_mask] = np.nan

        # lat and lon
        lat = NO2_data_dict['Latitude']
        lon = NO2_data_dict['Longitude']
        lat_e = NO2_data_dict['Latitude_e']
        lon_e = NO2_data_dict['Longitude_e']

        # loop region
        for region_name in region_name_list:

            print('month: ' + month + ', AK: ' + ak + ', region: ' + 
                    region_name)

            region_limit = sn_p.region_dict[region_name]
            i1, i2, j1, j2 = \
                    get_center_index_latlon(lat_e[:,0], lon_e[0,:],
                            region_limit)
            i1 = i1 + 1
            j1 = j1 + 1
            r_lat = lat[i1:i2+1,j1:j2+1]
            r_lon = lon[i1:i2+1,j1:j2+1]
            r_lat_e = lat_e[i1:i2+2,j1:j2+2]
            r_lon_e = lon_e[i1:i2+2,j1:j2+2]
            print(r_lat.shape)

            fig_dir = roor_fig_dir + region_name + '/'
            if not os.path.exists(fig_dir):
                os.system('mkdir -p ' + fig_dir)

            # NO2 VCD (OMI and GC scenes)
            region_NO2_VCD_dict = {}
            for NO2_VCD_name in NO2_VCD_name_list:

                # get data in the region
                NO2_VCD = NO2_data_dict[NO2_VCD_name][i1:i2+1,j1:j2+1]
                region_NO2_VCD_dict[NO2_VCD_name] = NO2_VCD

                    # whether or not plot region selection                   
                if region_flag_dict.get(region_name, False):

                        region_flag_dict[region_name] = False

                        #cartopy_plot(r_lon_e, r_lat_e, NO2_VCD,
                        #        countries=True, states=True)
                        #plt.savefig(fig_dir + region_name + '_data.png', 
                        #        format='png',dpi=300)

            # convert data to 1-D
            flag = np.full(r_lat.shape, True)
            for NO2_VCD_name in NO2_VCD_name_list:
                NO2_VCD = region_NO2_VCD_dict[NO2_VCD_name]
                flag = np.logical_and(flag, np.logical_not(np.isnan(NO2_VCD)))
            print(flag)
            for NO2_VCD_name in NO2_VCD_name_list:
                region_NO2_VCD_dict[NO2_VCD_name] = \
                        np.array(region_NO2_VCD_dict[NO2_VCD_name][flag]) /1e15
                print(region_NO2_VCD_dict[NO2_VCD_name].shape)

            # begin plot
            fig = plt.figure(figsize=(7,7))
            plt.rcParams.update({'font.size': 12})

            mtd = ModTaylorDiagram(fig=fig, 
                    title_expected=label_dict[sat_NO2_name])

            sat_NO2_VCD = region_NO2_VCD_dict[sat_NO2_name]
            stringID_labels = []
            for i in range(len(model_NO2_VCD_name_list)):
                NO2_VCD_name = model_NO2_VCD_name_list[i]
                NO2_VCD = region_NO2_VCD_dict[NO2_VCD_name]
                mtd.add_prediction(sat_NO2_VCD, NO2_VCD, '', str(i+1), 'o')
                print(i, sat_NO2_VCD.shape, NO2_VCD.shape)
                print(type(sat_NO2_VCD), type(NO2_VCD))
                print(sat_NO2_VCD)
                print(NO2_VCD)
                label = (r'$'+str(i+1)+'$', label_dict[NO2_VCD_name])
                stringID_labels.append(label)

            mtd.plot()

            mtd.stringID_legend( stringID_labels, loc='upper left' )


            figname = fig_dir + 'taylor_diagram_' + \
                    region_name + '_' +  month  + \
                    '_NO2_VCD' + ak + '_' + time_range[0:9] + '.png'
            plt.savefig(figname, format='png', dpi=300)

