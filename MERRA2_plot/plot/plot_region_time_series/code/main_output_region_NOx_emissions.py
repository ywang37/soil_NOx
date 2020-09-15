"""
Created on June 30, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import sys, getopt

from mylib.cartopy_plot import cartopy_plot
from mylib.constants import day_to_second
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
from sn_plot import plot_stackplot

#######################
# Start user parameters
#

start_year = 2012
end_year   = 2012

#xx = np.array(range(start_year,end_year+1))

year_range = str(start_year) + '-' + str(end_year)

#month_range_list = ['06-06', '07-07', '08-08', '06-08']
month_range_list = ['06-08']

emissions_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/emissions/data/'

area_file = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_runs/GC_ori/\
runs/merra2_2x25_tropchem_spinup_20050601/OutputDir/\
HEMCO_diagnostics.200505160000.nc'

#root_fig_dir = '../figure/NOx_emissions/'

verbose = True

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs', \
        'soil_T_obs_fit_18_degree',
        'soil_T_obs_fit_18_degree_scale_0.60',
        'soil_T_obs_fit_18_degree_scale_0.70',
        'soil_T_obs_fit_18_degree_scale_0.80',
        'soil_T_obs_fit_OMI_NO2',
        'soil_T_obs_scale_0.60',
        'soil_T_obs_scale_0.65',
        'soil_T_obs_scale_0.70',
        'soil_T_obs_scale_0.75',
        'soil_T_obs_scale_0.80'
        ]

NO_sectors = ['EmisNO_BioBurn', 'EmisNO_Soil', \
        'EmisNO_Anthro', 'EmisNO_Lightning']

NO2_setctors = ['EmisNO2_Anthro']

NOx_emissions_name_list = ['EmisNOx_Anthro', 
        'EmisNO_Soil', 'EmisNO_Lightning',
        'EmisNO_BioBurn']
NOx_emissions_total_name = 'EmisNOx_Total'

#label_dict = {}
#label_dict['EmisNOx_Anthro_yearly'] = 'Anthropogenic'
#label_dict['EmisNO_Soil_yearly'] = 'Soil'
#label_dict['EmisNO_Lightning_yearly'] = 'Lightning'
#label_dict['EmisNO_BioBurn_yearly'] = 'Fires'

#region_name_list = ['Central_US', 'Eastern_US', 'US', 'CA_NV']
region_name_list = ['Central_US']

# whether or not plot region selection
region_flag_dict = {}
region_flag_dict['Central_US'] = False
region_flag_dict['Eastern_US'] = False
region_flag_dict['US'] = False
region_flag_dict['CA_NV'] = False

#xlabel = 'Year'

#ylabel = r'NO$_x$ emissions [Gg N]'

#xminorticks = np.array(range(start_year, end_year+1))

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

coor_varns = ['Latitude', 'Latitude_e', 'Longitude', 'Longitude_e']

#
for month_range in month_range_list:

    data_dict = {}
    for region_name in region_name_list:
        data_dict[region_name] = {}
        data_dict[region_name]['scene'] = []
        for scene in scene_tup:
            data_dict[region_name]['scene'].append(scene)
        for  NOx_emissions_name in \
                NOx_emissions_name_list:
            data_dict[region_name][NOx_emissions_name] = []


    for scene in scene_tup:

        time_range = year_range + '_' + month_range

        emissions_file = emissions_dir + scene + '/' + \
                'HEMCO_diagnostics.' + time_range + '.nc'
 
        NOx_data_dict = read_NOx_emissions(emissions_file, 
                coor_varns=coor_varns, verbose=verbose,
                suffix='')

        # lat and lon
        lat = NOx_data_dict['Latitude']
        lon = NOx_data_dict['Longitude']
        lat_e = NOx_data_dict['Latitude_e']
        lon_e = NOx_data_dict['Longitude_e']

        NOx_data_dict['EmisNOx_Anthro'] = \
                NOx_data_dict['EmisNO_Anthro'] + \
                NOx_data_dict['EmisNO2_Anthro']

        # loop region
        for region_name in region_name_list:

            print('scene: ' + scene + ', month: ' + month_range + 
                    ', region: ' + region_name)

            region_limit = sn_p.region_dict[region_name]

            # NOx emissions
            for NOx_emissions_name in NOx_emissions_name_list:

                NOx_emissions = NOx_data_dict[NOx_emissions_name]


                region_NOx_emissions = \
                        region_ave_sum(NOx_emissions, 
                        in_weight=area,
                        flag_mask=US_flag_mask,
                        lat=lat, lon=lon, region_limit=region_limit)

                tmp = region_NOx_emissions['sum']

                # conversion
                # ng N/s => Gg N
                tmp = tmp * day_to_second * \
                        sn_p.days_dict[month_range] / 1e18

                # save data
                data_dict[region_name][NOx_emissions_name].append(tmp)


    # DataFrame
    for region_name in region_name_list:
        print('-------------------- ' + region_name  + ' --------------------')
        data_dict[region_name] = pd.DataFrame(data_dict[region_name])
        print(data_dict[region_name])



#            figname = fig_dir + region_name + '_' + scene + \
#                    '_NOx_emissions' + '_' + time_range + '.png'
#            plt.savefig(figname, format='png', dpi=300)

