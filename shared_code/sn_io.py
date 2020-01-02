"""
Created on January 2, 2020

@author: Yi Wang
"""

import numpy as np

from mylib.grid_utility import generate_grid_gc_2x25
from mylib.io import read_nc

def read_nc_emissions_multifiles(root_dir, scene_tup, month,
        varname, gc_run='geosfp_2x25_tropchem', res='2x25', 
        verbose=True):
    """
    """

    out_dict = {}
    emi_dict = {}
    for i in range(len(scene_tup)):

        scene = scene_tup[i]

        filename = root_dir + 'GEOS-Chem_' + scene + '/runs/' + \
                gc_run + '_' + month + '/' + 'HEMCO_diagnostics.' + \
                month + '010000.nc'

        data_dict = read_nc(filename, [varname], verbose=verbose, squeeze=True)
        emi_dict[scene] = data_dict[varname]

    if res == '2x25':
        lat_e, lon_e, lat_c, lon_c = generate_grid_gc_2x25()
    lon_e, lat_e = np.meshgrid(lon_e, lat_e)


    out_dict['scene_tup'] = scene_tup
    out_dict['emi_dict'] = emi_dict
    out_dict['lat_e'] = lat_e
    out_dict['lon_e'] = lon_e

    return out_dict
