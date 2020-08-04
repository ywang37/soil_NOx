"""
Created on July 13, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

from mylib.grid_utility import generate_grid_gc_2x25
from mylib.grid_utility import region_ave_sum, fill_region_color
from mylib.io import read_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
import sn_p

#######################
# Start user parameters
#

US_flag_file = '/Dedicated/jwang-data/ywang/opt/anaconda3/lib/python3.7/\
site-packages/mylib/old/pixel_in_contiguous_US/US_flag_2x25.nc'

region_name_list = ['US', 'CA_NV', 'Central_US', 'Eastern_US']

color_list = ['darkgrey', 'C0', 'C1', 'C2']

verbose=True

xtick = (-120, -100, -80)
ytick = (30, 40, 50)

#
# End user parameters
#####################

# US flag mask
US_flag_mask = read_nc(US_flag_file, ['flag'], verbose=verbose)
US_flag_mask = US_flag_mask['flag']
US_flag_mask = US_flag_mask < 0.5


lat_e, lon_e, lat, lon = generate_grid_gc_2x25()
lon_e, lat_e = np.meshgrid(lon_e, lat_e)
lon, lat = np.meshgrid(lon, lat)
tmp0 = np.zeros_like(lon)

# region
region_flag_list = []
for i in range(len(region_name_list)):

    region_name = region_name_list[i]

    region_limit = sn_p.region_dict[region_name]

    tmp1 = region_ave_sum(tmp0, flag_mask=US_flag_mask,
            lat=lat, lon=lon, region_limit=region_limit)

    region_flag = np.logical_not(tmp1['final_flag'])
    region_flag_list.append(region_flag)

# plot
fill_region_color(lon_e, lat_e, region_flag_list, color_list,
        xtick=xtick, ytick=ytick,
        region_limit=sn_p.region_dict['US'])
plt.savefig('../figure/region_color_shaded.png', format='png', dpi=300)


