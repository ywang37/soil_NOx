"""
Created on March 29, 2020

@author: Yi Wang
"""

import numpy as np

def generate_start_end(month):
    """
    month is 'YYYYMM'

    get 'YYYY-MM-01_YYYY-MM-end'
    """

    days_dict = {
            '06' : '30',
            '07' : '31',
            '08' : '31'
            }

    YYYY = month[0:4]
    MM   = month[4:6]

    result = YYYY + '-' + MM + '-01_' + YYYY + '-' + MM + '-' + days_dict[MM]

    return result

def generate_latlon_label(lat, lon):
    """
    """

    # latitide
    if (lat >= 0.0):
        lat_dir = 'N'
    else:
        lat_dir = 'S'
    lat_label = str(np.absolute(lat)) + '\u00b0' + lat_dir

    # longitude
    if (lon >= 0.0):
        lon_dir = 'E'
    else:
        lon_dir = 'W'
    lon_label = str(np.absolute(lon)) + '\u00b0' + lon_dir

    return lat_label, lon_label
