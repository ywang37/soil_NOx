"""
Created on January 2, 2020

@author: Yi Wang
"""

from netCDF4 import Dataset
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
        if emi_dict[scene].ndim == 3:
            emi_dict[scene] = np.sum(emi_dict[scene], axis=0)

    if res == '2x25':
        lat_e, lon_e, lat_c, lon_c = generate_grid_gc_2x25()
    lon_e, lat_e = np.meshgrid(lon_e, lat_e)


    out_dict['scene_tup'] = scene_tup
    out_dict['emi_dict'] = emi_dict
    out_dict['lat_e'] = lat_e
    out_dict['lon_e'] = lon_e

    return out_dict


def get_nc_NO_emi_ratio_multifiles(root_dir, scene_tup, month,
        varname,
        gc_run='geosfp_2x25_tropchem', res='2x25', 
        verbose=True):
    """
    The ratio *varname* NO emisions to total NO emissions
    """

    varname_all = [varname, 'EmisNO_Total']

    out_dict = {}
    emi_ratio_dict = {}
    for i in range(len(scene_tup)):

        scene = scene_tup[i]

        filename = root_dir + 'GEOS-Chem_' + scene + '/runs/' + \
                gc_run + '_' + month + '/' + 'HEMCO_diagnostics.' + \
                month + '010000.nc'

        data_dict = read_nc(filename, varname_all, 
                verbose=verbose, squeeze=True)

        # soil
        emi_soil = data_dict[varname_all[0]]

        # total
        emi_total = np.sum(data_dict[varname_all[1]], axis=0)

        # ratio
        ratio = np.full_like(emi_total, np.nan)
        flag = emi_total > 1e-6
        ratio[flag] = emi_soil[flag] / emi_total[flag]
        emi_ratio_dict[scene] = ratio


    if res == '2x25':
        lat_e, lon_e, lat_c, lon_c = generate_grid_gc_2x25()
    lon_e, lat_e = np.meshgrid(lon_e, lat_e)


    out_dict['scene_tup'] = scene_tup
    out_dict['emi_ratio_dict'] = emi_ratio_dict
    out_dict['lat_e'] = lat_e
    out_dict['lon_e'] = lon_e

    return out_dict


def save_ave(filename, data_dict, verbose=True):
    """ save avarage data
    """

    coord_name_list = [
            'Latitude',
            'Latitude_e',
            'Longitude',
            'Longitude_e'
            ]

    if verbose:
        print(' - save_ave: output ' + filename)

    # open file
    nc_f = Dataset(filename, 'w')

    # grid, _e means edge
    Latitude    = data_dict['Latitude']
    Longitude   = data_dict['Longitude']
    Latitude_e  = data_dict['Latitude_e']
    Longitude_e = data_dict['Longitude_e']

    # Dimensions of a netCDF file
    dim_lat = nc_f.createDimension('Latitude',  Latitude.shape[0])
    dim_lon = nc_f.createDimension('Longitude', Latitude.shape[1])
    dim_lat_e = nc_f.createDimension('Latitude_e',  Latitude_e.shape[0])
    dim_lon_e = nc_f.createDimension('Longitude_e', Latitude_e.shape[1])

    # create variables in a netCDF file

    # lat and lon
    Latitude_v = nc_f.createVariable('Latitude', 'f4', 
            ('Latitude', 'Longitude'))
    Longitude_v = nc_f.createVariable('Longitude', 'f4',
            ('Latitude', 'Longitude'))
    Latitude_e_v = nc_f.createVariable('Latitude_e', 'f4',
            ('Latitude_e', 'Longitude_e'))
    Longitude_e_v = nc_f.createVariable('Longitude_e', 'f4',
            ('Latitude_e', 'Longitude_e'))

    # variables
    nc_var_dict = {}
    for varname in data_dict:
        if not (varname in coord_name_list):
            if varname is not 'count':
                nc_var = nc_f.createVariable(varname, 'f4',
                        ('Latitude', 'Longitude'))
            else:
                nc_var = nc_f.createVariable(varname, int,
                        ('Latitude', 'Longitude'))
            nc_var_dict[varname] = nc_var

    # write variables

    # lat and lon
    Latitude_v[:]    = Latitude
    Longitude_v[:]   = Longitude
    Latitude_e_v[:]  = Latitude_e
    Longitude_e_v[:] = Longitude_e

    for varname in nc_var_dict:
        nc_var_dict[varname][:] = data_dict[varname]

    # close file
    nc_f.close()

