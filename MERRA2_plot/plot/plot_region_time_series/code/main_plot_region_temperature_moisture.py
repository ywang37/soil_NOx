"""
Created on June 30, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import os
import sys, getopt

from mylib.cartopy_plot import cartopy_plot
from mylib.colormap.colormap_utility import truncate_colormap
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.colormap.emphasize_small_map import emphasize_small_map
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.grid_utility import region_ave_sum
from mylib.io import read_nc
from mylib.time_series.plot_time_series import plot_time_series_twinx
from mylib.trend_analysis import trend_analysis

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
import sn_p

#######################
# Start user parameters
#

start_year = 2005
end_year   = 2019

xx = np.array(range(start_year,end_year+1))

time_range = str(start_year) + '-' + str(end_year) + '_06-08'

infile = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/resample/data/multimonth_mean/\
model_satellite_' + time_range + '.nc'

area_file = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_runs/GC_ori/\
runs/merra2_2x25_tropchem_spinup_20050601/OutputDir/\
HEMCO_diagnostics.200505160000.nc'

roor_fig_dir = '../figure/temp_moist/'

month_list = ['multimonth_mean', 'month_06', 'month_07', 'month_08']

verbose = True

region_name_list = ['Central_US', 'Eastern_US', 'US', 'CA_NV']

# whether or not plot region selection
region_flag_dict = {}
region_flag_dict['Central_US'] = True
region_flag_dict['Eastern_US'] = True
region_flag_dict['US'] = True
region_flag_dict['CA_NV'] = True

temp_varn_list = ['mod_Met_TS', 'mod_Met_TSOIL1']
temp_varn_dict = {}
temp_varn_dict['mod_Met_TS'] = '2-m'
temp_varn_dict['mod_Met_TSOIL1'] = 'Soil'

xlabel = 'Year'

l_ylim = [20.0, 40.0]
r_ylim = [0.0, 1.0]



US_flag_file = '/Dedicated/jwang-data/ywang/opt/anaconda3/lib/python3.7/\
site-packages/mylib/old/pixel_in_contiguous_US/US_flag_2x25.nc'

#
# End user parameters
#####################

# area
tmp_data = read_nc(area_file, ['AREA'], verbose=verbose)
area = tmp_data['AREA']


# US flag mask
US_flag_mask = read_nc(US_flag_file, ['flag'], verbose=verbose)
US_flag_mask = US_flag_mask['flag']
US_flag_mask = US_flag_mask < 0.5

coor_list = ['Latitude_e', 'Longitude_e', \
        'Latitude', 'Longitude']

var_prefix_list = ['mod_Met_GWETTOP', 'mod_Met_TS', \
        'mod_Met_TSOIL1']

#
for month in month_list:

    varname_list = []
    for i in range(len(var_prefix_list)):
        varname_list.append(var_prefix_list[i] + '_' + month)
    varname_list = varname_list + coor_list
 
    all_data_dict = \
            read_nc(infile, varname_list, verbose=verbose)

    # lat and lon
    lat = all_data_dict['Latitude']
    lon = all_data_dict['Longitude']
    lat_e = all_data_dict['Latitude_e']
    lon_e = all_data_dict['Longitude_e']

    # loop region
    for region_name in region_name_list:

        print('month: ' + month + ', region: ' + region_name)

        region_limit = sn_p.region_dict[region_name]

        fig_dir = roor_fig_dir + region_name + '/'
        if not os.path.exists(fig_dir):
            os.system('mkdir -p ' + fig_dir)

        l_data_dict = {}
        r_data_dict = {}
        for temp_varn in temp_varn_list:

            # temperature
            temp = all_data_dict[temp_varn + '_' + month] - 273.15
            dim = temp.shape
            l_data_dict[temp_varn] = \
                    np.zeros((dim[0],), dtype=np.float)

            # moisture
            moist = all_data_dict['mod_Met_GWETTOP_' + month]
            r_data_dict['mod_Met_GWETTOP'] = \
                    np.zeros((dim[0],), dtype=np.float)

            for i in range(dim[0]):

                # temperature
                region_temp = region_ave_sum(temp[i,:,:], in_weight=area,
                        flag_mask=US_flag_mask,
                        lat=lat, lon=lon, region_limit=region_limit)
                l_data_dict[temp_varn][i] = region_temp['mean']

                # moisture
                region_moist = region_ave_sum(moist[i,:,:], in_weight=area,
                        flag_mask=US_flag_mask,
                        lat=lat, lon=lon, region_limit=region_limit)
                r_data_dict['mod_Met_GWETTOP'][i] = region_moist['mean']

                # whether or not plot region selection                   
                if region_flag_dict.get(region_name, False):

                    region_flag_dict[region_name] = False

                    cartopy_plot(lon_e, lat_e, region_temp['data'],
                            countries=True, states=True)
                    plt.savefig(fig_dir + region_name + '_data.png', 
                            format='png',dpi=300)
                    cartopy_plot(lon_e, lat_e, region_temp['weight'],
                            countries=True, states=True)
                    plt.savefig(fig_dir + region_name + '_weight.png', 
                            format='png',dpi=300)
                    cartopy_plot(lon_e, lat_e, 
                            region_temp['final_flag'].astype(float),
                            countries=True, states=True)
                    plt.savefig(fig_dir + region_name + '_final_flag.png', 
                            format='png',dpi=300)

            # temp trends
            yy_temp = deepcopy(l_data_dict[temp_varn])
            ta_temp = trend_analysis()
            ta_temp.analysis_yearly(yy_temp)
            td_label_temp = '({:.3f}'.format(ta_temp.popt[1]) + \
                    u'\u00B1' + '{:.3f}'.format(ta_temp.trend_std) + \
                    u' \u00B0C' + r' a$^{-1}$)'
            print(td_label_temp)
            l_label_dict = {}
            l_label_dict[temp_varn] = temp_varn_dict[temp_varn] + \
                    ' temperature '+ td_label_temp

            # moisture trends
            yy_moist = deepcopy(r_data_dict['mod_Met_GWETTOP'])
            ta_moist = trend_analysis()
            ta_moist.analysis_yearly(yy_moist)
            td_label_moist = '({:.5f}'.format(ta_moist.popt[1]) + \
                    u'\u00B1' + '{:.5f}'.format(ta_moist.trend_std) + \
                    r' a$^{-1}$)'
            print(td_label_moist)
            r_label_dict = {}
            r_label_dict['mod_Met_GWETTOP'] = \
                    'Soil moisture ' + td_label_moist

            # begin plot
            l_ylabel = temp_varn_dict[temp_varn] + \
                    ' temperature ' + u'[\u00B0C]'
            r_ylabel = 'Soil moisture [unitless]'
            plot_time_series_twinx(l_data_dict, r_data_dict,
                    [temp_varn], ['mod_Met_GWETTOP'], xx,
                    l_label_dict=l_label_dict, r_label_dict=r_label_dict,
                    l_ylim=l_ylim, r_ylim=r_ylim,
                    xlabel=xlabel,
                    l_ylabel=l_ylabel, r_ylabel=r_ylabel)
            figname = fig_dir + region_name + '_' +  month  + \
                     '_' + temp_varn_dict[temp_varn] + '_temp_moist_' + \
                     time_range[0:9] + '.png'
            plt.savefig(figname, format='png', dpi=300)
