"""
Created on September 3, 2020

@author: Yi Wang
"""

import numpy as np

from mylib.borders.select_pixel import pixel_in_US_state
from mylib.grid_utility import generate_grid_2
from mylib.io import write_nc

# define grid
lat_min =    9.75
lat_max =   70.25
lat_res =    0.5
lon_min = -140.3125
lon_max =  -39.6875
lon_res =    0.625

# generate grid
lat_e, lon_e, lat_c, lon_c = generate_grid_2(lat_res, lat_min, lat_max,
        lon_res, lon_min, lon_max)
lon_c, lat_c = np.meshgrid(lon_c, lat_c)
lon_e, lat_e = np.meshgrid(lon_e, lat_e)
pixel = np.array([lon_c.flatten(), lat_c.flatten()]).T

# find pixels in CA
result = pixel_in_US_state(pixel, 'California')
flag = result['flag']

# output
data_dict = {}
data_dict['Latitude'] = lat_c
data_dict['Latitude_e'] = lat_e
data_dict['Longitude'] = lon_c
data_dict['Longitude_e'] = lon_e
data_dict['flag'] = flag
write_nc('../data/CA_flag_05x0625.nc', data_dict)
