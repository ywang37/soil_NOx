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
from mylib.trend_analysis import trend_analysis

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_io import read_NOx_emissions
import sn_p
from sn_plot import plot_time_series

#######################
# Start user parameters
#

start_year = 2005
end_year   = 2019

xx = np.array(range(start_year,end_year+1))

year_range = str(start_year) + '-' + str(end_year)

time_range = year_range + '_06-08'

NO2_VCD_file = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/resample/data/multimonth_mean/\
model_satellite_' + time_range + '.nc'

area_file = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_runs/GC_ori/\
runs/old/merra2_2x25_tropchem_spinup_20050601/OutputDir/\
HEMCO_diagnostics.200505160000.nc'

emissions_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/emissions/data/'

root_fig_dir = '../figure_soil_contribution/'

month_list = ['multimonth_mean', 'month_06', 'month_07', 'month_08']

verbose = True

ak = '_AK'

scene_tup = ['ori', 'soil_T_obs']

model_label_dict = {}
model_label_dict['ori'] = 'Control'
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

NO_soil_name = 'EmisNO_Soil_yearly'
NO_total_name = 'EmisNO_Total_yearly'
NOx_total_name = 'EmisNOx_Total_yearly'
NO2_anthro_name = 'EmisNO2_Anthro_yearly'

NOx_emissions_name_list = []
for scene in scene_tup:
    NOx_emissions_name_list.append(NO_soil_name+'_'+scene)
    NOx_emissions_name_list.append(NOx_total_name+'_'+scene)


#
# End user parameters
#####################

#
month_list = ['multimonth_mean', 'month_06', 'month_07', 'month_08']
month_range_dict = {}
month_range_dict['multimonth_mean'] = '06-08'
month_range_dict['month_06']        = '06-06'
month_range_dict['month_07']        = '07-07'
month_range_dict['month_08']        = '08-08'

# area
tmp_data = read_nc(area_file, ['AREA'], verbose=verbose)
area = tmp_data['AREA']


# US flag mask
US_flag_mask = read_nc(US_flag_file, ['flag'], verbose=verbose)
US_flag_mask = US_flag_mask['flag']
US_flag_mask = US_flag_mask < 0.5

#
for month in month_list:

    month_range = month_range_dict[month]

    # NO2 VCD
    sat_NO2_name = 'sat_ColumnAmountNO2Trop_'+month

    overpass_varname_list = ['Latitude_e', 'Longitude_e', 
            'Latitude', 'Longitude', sat_NO2_name]
    NO2_VCD_name_list = [sat_NO2_name]
    all_NO2_VCD_name_list = [sat_NO2_name]
    model_NO2_VCD_name_list = []
    label_dict = {}
    label_dict[sat_NO2_name] = 'OMI'
    for scene in scene_tup:

        mod_NO2_varname = 'mod_NO2Trop_AK_tp_sat_' + scene \
                    + '_' + month
        overpass_varname_list.append(mod_NO2_varname)
        NO2_VCD_name_list.append(mod_NO2_varname)
        all_NO2_VCD_name_list.append(mod_NO2_varname)
        all_NO2_VCD_name_list.append(mod_NO2_varname+'_no_soil')
        model_NO2_VCD_name_list.append(mod_NO2_varname)
        label_dict[mod_NO2_varname] = \
               model_label_dict.get(scene)
        label_dict[mod_NO2_varname+'_no_soil'] = \
                model_label_dict.get(scene) + '(no_soil)'
 
    NO2_data_dict = \
            read_nc(NO2_VCD_file, overpass_varname_list, verbose=verbose)

    # lat and lon
    lat = NO2_data_dict['Latitude']
    lon = NO2_data_dict['Longitude']
    lat_e = NO2_data_dict['Latitude_e']
    lon_e = NO2_data_dict['Longitude_e']

    # NOx emissions
    NOx_emi_data_dict = {}
    for scene in scene_tup:
        emissions_file = emissions_dir + scene + '/' + \
                'HEMCO_diagnostics.' + time_range + '.nc'

        tmp_NOx_emi_data_dict = read_NOx_emissions(emissions_file,
                verbose=verbose, suffix='_yearly')

        # total
        NOx_emi_data_dict[NOx_total_name+'_'+scene] = \
                tmp_NOx_emi_data_dict[NO_total_name] + \
                tmp_NOx_emi_data_dict[NO2_anthro_name]

        # soil
        NOx_emi_data_dict[NO_soil_name+'_'+scene] = \
                tmp_NOx_emi_data_dict[NO_soil_name]

    # loop region
    for region_name in region_name_list:

        print('month: ' + month + ', region: ' + 
                region_name)

        region_limit = sn_p.region_dict[region_name]

        fig_dir = root_fig_dir + region_name + '/'
        if not os.path.exists(fig_dir):
            os.system('mkdir -p ' + fig_dir)

        # NO2 VCD (OMI and GC scenes)
        region_NO2_VCD_dict = {}
        for NO2_VCD_name in NO2_VCD_name_list:

            NO2_VCD = NO2_data_dict[NO2_VCD_name]
            dim = NO2_VCD.shape

            region_NO2_VCD_dict[NO2_VCD_name] = \
                    np.zeros((dim[0],), dtype=np.float)

            # loop time
            for i in range(dim[0]):

                region_NO2_VCD = \
                        region_ave_sum(NO2_VCD[i,:,:], in_weight=area,
                        flag_mask=US_flag_mask,
                        lat=lat, lon=lon, region_limit=region_limit)

                # save data
                region_NO2_VCD_dict[NO2_VCD_name][i] = \
                        region_NO2_VCD['mean']

                # whether or not plot region selection                   
                if region_flag_dict.get(region_name, False):

                    region_flag_dict[region_name] = False

                    cartopy_plot(lon_e, lat_e, region_NO2_VCD['data'],
                            countries=True, states=True)
                    plt.savefig(fig_dir + region_name + '_data.png', 
                            format='png',dpi=300)
                    cartopy_plot(lon_e, lat_e, region_NO2_VCD['weight'],
                            countries=True, states=True)
                    plt.savefig(fig_dir + region_name + '_weight.png', 
                            format='png',dpi=300)
                    cartopy_plot(lon_e, lat_e, 
                            region_NO2_VCD['final_flag'].astype(float),
                            countries=True, states=True)
                    plt.savefig(fig_dir + region_name + '_final_flag.png', 
                            format='png',dpi=300)

        # NOx emissions
        region_NOx_emissions_dict = {}
        for NOx_emissions_name in NOx_emissions_name_list:

            NOx_emissions = NOx_emi_data_dict[NOx_emissions_name]
            dim = NOx_emissions.shape

            region_NOx_emissions_dict[NOx_emissions_name] = \
                    np.zeros((dim[0],), dtype=np.float)

            # loop time
            for i in range(dim[0]):

                region_NOx_emissions = \
                        region_ave_sum(NOx_emissions[i,:,:],
                                in_weight=area,
                                flag_mask=US_flag_mask,
                                lat=lat, lon=lon, region_limit=region_limit)

                # save data
                region_NOx_emissions_dict[NOx_emissions_name][i] = \
                        region_NOx_emissions['sum']

        # soil ratio
        soil_ratio_dict = {}
        for scene in scene_tup:

            soil_ratio = np.zeros((dim[0],), dtype=np.float)
            for i in range(dim[0]):
                soil_ratio[i] = \
                    region_NOx_emissions_dict[NO_soil_name+'_'+scene][i] / \
                    region_NOx_emissions_dict[NOx_total_name+'_'+scene][i]

            soil_ratio_dict['soil_ratio_'+scene] = soil_ratio


        for NO2_VCD_name in model_NO2_VCD_name_list:

            region_NO2_VCD_dict[NO2_VCD_name+'_no_soil'] = \
                    region_NO2_VCD_dict[NO2_VCD_name] * \
                    (1.0 - soil_ratio)

        # trend without soil
        for scene in scene_tup:
            mod_NO2_varname = 'mod_NO2Trop_AK_tp_sat_' + scene \
                    + '_' + month
            region_NO2_VCD_dict[NO2_VCD_name+'_no_soil'] = \
                    region_NO2_VCD_dict[NO2_VCD_name] * \
                    (1.0 - soil_ratio_dict['soil_ratio_'+scene])

        # begin plot

        # original value
        plot_time_series(region_NO2_VCD_dict, all_NO2_VCD_name_list,
                xx, label_dict=label_dict,
                xlabel=xlabel, ylabel=ylabel,
                ylim=ylim)
        figname = fig_dir + region_name + '_' +  month  + \
                 '_NO2_VCD_AK_' + time_range[0:9] + '.png'
        plt.savefig(figname, format='png', dpi=300)

        # percentage trend
        pert_label_dict = {}
        for NO2_VCD_name in all_NO2_VCD_name_list:

            # get data
            yy = deepcopy(region_NO2_VCD_dict[NO2_VCD_name])
            yy = yy / yy[0] * percentage_scale
                
            # trend analysis
            ta = trend_analysis()
            ta.analysis_yearly(yy)
            print(ta.popt[1], ta.trend_std)
            trend_label = '({:.2f}'.format(ta.popt[1]) + \
                    u'\u00B1' + '{:.2f}'.format(ta.trend_std) + \
                    r' % a$^{-1}$)'
            pert_label_dict[NO2_VCD_name] = \
                    label_dict[NO2_VCD_name] + ' ' + \
                    trend_label

        # percentage
        plot_time_series(region_NO2_VCD_dict, all_NO2_VCD_name_list,
                xx, label_dict=pert_label_dict, ratio=True,
                scale=percentage_scale,
                xlabel=xlabel, ylabel=ratio_ylabel,
                ylim=ratio_ylim)
        figname = fig_dir + 'ratio_' + region_name + '_' +  month  + \
                '_NO2_VCD_AK_' + time_range[0:9] + '.png'
        plt.savefig(figname, format='png', dpi=300)






