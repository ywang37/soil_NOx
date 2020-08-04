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

start_year = 2005
end_year   = 2019

xx = np.array(range(start_year,end_year+1))

year_range = str(start_year) + '-' + str(end_year)

month_range_list = ['06-06', '07-07', '08-08', '06-08']

emissions_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/emissions/data/'

area_file = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_runs/GC_ori/\
runs/merra2_2x25_tropchem_spinup_20050601/OutputDir/\
HEMCO_diagnostics.200505160000.nc'

root_fig_dir = '../figure/NOx_emissions/'

verbose = True

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

NO_sectors = ['EmisNO_BioBurn_yearly', 'EmisNO_Soil_yearly', \
        'EmisNO_Anthro_yearly', 'EmisNO_Lightning_yearly']

NO2_setctors = ['EmisNO2_Anthro_yearly']

NOx_emissions_name_list = ['EmisNOx_Anthro_yearly', 
        'EmisNO_Soil_yearly', 'EmisNO_Lightning_yearly',
        'EmisNO_BioBurn_yearly']

label_dict = {}
label_dict['EmisNOx_Anthro_yearly'] = 'Anthropogenic'
label_dict['EmisNO_Soil_yearly'] = 'Soil'
label_dict['EmisNO_Lightning_yearly'] = 'Lightning'
label_dict['EmisNO_BioBurn_yearly'] = 'Fires'

region_name_list = ['Central_US', 'Eastern_US', 'US', 'CA_NV']

# whether or not plot region selection
region_flag_dict = {}
region_flag_dict['Central_US'] = True
region_flag_dict['Eastern_US'] = True
region_flag_dict['US'] = True
region_flag_dict['CA_NV'] = True

xlabel = 'Year'

ylabel = r'NO$_x$ emissions [Gg N]'

xminorticks = np.array(range(start_year, end_year+1))

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
for scene in scene_tup:
    
    for month_range in month_range_list:

        time_range = year_range + '_' + month_range

        emissions_file = emissions_dir + scene + '/' + \
                'HEMCO_diagnostics.' + time_range + '.nc'

 
        NOx_data_dict = read_NOx_emissions(emissions_file, 
                coor_varns=coor_varns, verbose=verbose,
                suffix='_yearly')

        # lat and lon
        lat = NOx_data_dict['Latitude']
        lon = NOx_data_dict['Longitude']
        lat_e = NOx_data_dict['Latitude_e']
        lon_e = NOx_data_dict['Longitude_e']

        NOx_data_dict['EmisNOx_Anthro_yearly'] = \
                NOx_data_dict['EmisNO_Anthro_yearly'] + \
                NOx_data_dict['EmisNO2_Anthro_yearly']

        # loop region
        for region_name in region_name_list:

            print('scene: ' + scene + ', month: ' + month_range + 
                    ', region: ' + region_name)

            region_limit = sn_p.region_dict[region_name]

            fig_dir = root_fig_dir + region_name + '/'
            if not os.path.exists(fig_dir):
                os.system('mkdir -p ' + fig_dir)

            # NOx emissions
            region_NOx_emissions_dict = {}
            for NOx_emissions_name in NOx_emissions_name_list:

                NOx_emissions = NOx_data_dict[NOx_emissions_name]
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

                    # conversion
                    # ng N/s => Gg N
                    region_NOx_emissions_dict[NOx_emissions_name][i] = \
                            region_NOx_emissions_dict[NOx_emissions_name][i] \
                            * day_to_second * \
                            sn_p.days_dict[month_range] / 1e18

                    # whether or not plot region selection                   
                    if region_flag_dict.get(region_name, False):

                        region_flag_dict[region_name] = False

                        cartopy_plot(lon_e, lat_e, 
                                region_NOx_emissions['data'],
                                countries=True, states=True)
                        plt.savefig(fig_dir + region_name + '_data.png', 
                                format='png',dpi=300)
                        cartopy_plot(lon_e, lat_e, 
                                region_NOx_emissions['weight'],
                                countries=True, states=True)
                        plt.savefig(fig_dir + region_name + '_weight.png', 
                                format='png',dpi=300)
                        cartopy_plot(lon_e, lat_e, 
                                region_NOx_emissions['final_flag'].astype(float),
                                countries=True, states=True)
                        plt.savefig(fig_dir + region_name + '_final_flag.png', 
                                format='png',dpi=300)

            # begin plot

            # original value
            plot_stackplot(region_NOx_emissions_dict, 
                    NOx_emissions_name_list,
                    xx, label_dict=label_dict,
                    xlabel=xlabel, ylabel=ylabel,
                    xminorticks=xminorticks)
            figname = fig_dir + region_name + '_' + scene + \
                    '_NOx_emissions' + '_' + time_range + '.png'
            plt.savefig(figname, format='png', dpi=300)

