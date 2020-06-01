"""
Created on January 2, 2020

@author: Yi Wang
"""

from calendar import monthrange
from copy import deepcopy
import datetime
from netCDF4 import Dataset
import numpy as np
import os

from mylib.constants import molec_to_kgN
from mylib.grid_utility import generate_grid_gc_2x25, get_center_index
from mylib.grid_utility import get_center_index_latlon
from mylib.io import read_nc


#
#------------------------------------------------------------------------------
#
def read_nc_emissions_multifiles(root_dir, scene_tup, month,
        varname, gc_run='geosfp_2x25_tropchem', res='2x25',
        scene_prefix='GEOS-Chem_',
        subdir = '',
        path=1,
        verbose=True):
    """
    """

    out_dict = {}
    emi_dict = {}
    for i in range(len(scene_tup)):

        scene = scene_tup[i]

        if path == 1:
            filename = root_dir + scene_prefix + scene + '/runs/' + \
                    gc_run + '_' + month + '/' + subdir + 'HEMCO_diagnostics.' + \
                    month + '010000.nc'
        elif path == 2:
            filename = root_dir + scene + '/HEMCO_diagnostics.' + month  + '.nc'

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
#
#------------------------------------------------------------------------------
#
def read_NO_emissions(filename, varns=[], verbose=True,
        amount=True, days=None):
    """
    (ywang, 04/07/2020)
    Parameters
    ----------
    amount : logical
        If True, convert emission rates to emission amount
    """

    if len(varns)  == 0:
        varnames = ['EmisNO_Soil', 'EmisNO_Total']
    else:
        varnames = deepcopy(varns)
    varnames.append('AREA')
    varnames = list(set(varnames))

    if verbose:
        print(' - read_NO_emissions: reading ' + filename)

    # read data
    out_data = read_nc(filename, varnames, squeeze=True)
    for var in out_data:
        if out_data[var].ndim == 3:
            out_data[var] = np.sum(out_data[var], axis=0)

    # convert emission rates to emission amount
    if amount:

        # assume calcualte emissions for one month
        # and get month from filename
        if days is None:
            time_c = filename.split('/').split('.')[1]
            year = int(time_c[0:4])
            month = int(time_c[4:6])
            days = monthrange(year, month)[1]

        day_sceonds = 24.0 * 3600.0

        # molec/cm2/s => kg N
        AREA = out_data['AREA']
        for var in out_data:
            if var != 'AREA':
                out_data[var+'_amount'] = out_data[var] * day_sceonds * AREA \
                        * molec_to_kgN

    return out_data
#
#------------------------------------------------------------------------------
#
def get_nc_NO_emi_ratio_multifiles(root_dir, scene_tup, month,
        varname,
        gc_run='geosfp_2x25_tropchem', res='2x25', 
        region_limit=None,
        verbose=True):
    """
    The ratio *varname* NO emisions to total NO emissions
    """

    varname_all = [varname, 'EmisNO_Total']

    if res == '2x25':
        lat_e, lon_e, lat_c, lon_c = generate_grid_gc_2x25()

    # region_limit
    if region_limit is not None:
        i1, i2, j1, j2 = \
                get_center_index_latlon(lat_e, lon_e, region_limit)


    out_dict = {}
    emi_ratio_dict = {}
    emi_varname_dict = {}
    emi_total_dict = {}
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

        if region_limit is not None:
            emi_soil  = emi_soil[i1:i2+1,j1:j2+1]
            emi_total = emi_total[i1:i2+1,j1:j2+1]

        # ratio
        ratio = np.full_like(emi_total, np.nan)
        flag = emi_total > 1e-6
        ratio[flag] = emi_soil[flag] / emi_total[flag]
        emi_ratio_dict[scene] = ratio

        emi_varname_dict[scene] = emi_soil

        emi_total_dict[scene] = emi_total


    lon_e, lat_e = np.meshgrid(lon_e, lat_e)
    if region_limit is not None:
        lat_e = lat_e[i1:i2+2,j1:j2+2]
        lon_e = lon_e[i1:i2+2,j1:j2+2]

    out_dict['scene_tup'] = scene_tup
    out_dict['emi_ratio_dict'] = emi_ratio_dict
    out_dict[varname+'_dict'] = emi_varname_dict
    out_dict['emi_total_dict'] = emi_total_dict
    out_dict['lat_e'] = lat_e
    out_dict['lon_e'] = lon_e

    return out_dict
#
#------------------------------------------------------------------------------
#
def read_VCD_multiscenes(filename, scene_tup, mod_scene_var_list=[], 
        sat_var_list=[], mod_var_list=[], latlon=True, verbose=True):
    """
    """

    # all variable names
    varname_list = []

    # mod_scene_var_list
    for i in range(len(scene_tup)):

        scene = scene_tup[i]

        for j in range(len(mod_scene_var_list)):

            mod_scene_var = mod_scene_var_list[j]

            varname_list.append(mod_scene_var + scene)

    # sat_var_list
    varname_list.extend(sat_var_list)

    # mod_var_list
    varname_list.extend(mod_var_list)

    # latlon
    if latlon:
        varname_list.extend(['Latitude', 'Longitude', 
            'Latitude_e', 'Longitude_e'])

    # get all varilables
    out_dict = read_nc(filename, varname_list, verbose=verbose)

    return out_dict
#
#------------------------------------------------------------------------------
#
def read_VCD_multiscenes_days(file_str, data_dir, startDate, endDate,
        scene_tup, mod_scene_var_list=[],
        sat_var_list=[], mod_var_list=[], latlon=True,
        sat_cv_name='sat_ColumnAmountNO2Trop',
        sat_cv_flag=True, sat_cv_thre=None, sat_ave_thre=None,
        region_limit=None, verbose=True):
    """
    (ywang, 03/30/20)
    """

    # all variable names
    varname_list = []

    # mod_scene_var_list
    for i in range(len(scene_tup)):

        scene = scene_tup[i]

        for j in range(len(mod_scene_var_list)):

            mod_scene_var = mod_scene_var_list[j]

            varname_list.append(mod_scene_var + scene)

    # sat_var_list
    varname_list.extend(sat_var_list)

    # mod_var_list
    varname_list.extend(mod_var_list)

    # data_dir
    if data_dir[-1] != '/':
        data_dir = data_dir + '/'

    # Date
    currDate   = startDate
    currDate_D = datetime.datetime.strptime(currDate, '%Y-%m-%d')
    endDate_D  = datetime.datetime.strptime(endDate,  '%Y-%m-%d')

    first_flag = True
    out_dict = {}
    for varname in varname_list:
        out_dict[varname] = []
    while currDate_D <= endDate_D:
        
        # current date
        currDate = str(currDate_D)[0:10]
        print(''.join(np.full((79,), '-')))
        print('processing ' + currDate)

        filename = data_dir + file_str.replace('YYYY-MM-DD', currDate)

        if os.path.exists(filename):
            one_dict = read_VCD_multiscenes(filename, scene_tup, 
                    mod_scene_var_list=mod_scene_var_list,
                    sat_var_list=sat_var_list, 
                    mod_var_list=mod_var_list, latlon=first_flag, 
                    verbose=True)

            if first_flag:
                out_dict['Latitude']    = one_dict['Latitude']
                out_dict['Latitude_e']  = one_dict['Latitude_e']
                out_dict['Longitude']   = one_dict['Longitude']
                out_dict['Longitude_e'] = one_dict['Longitude_e']
                first_flag = False

            for varname in varname_list:
                out_dict[varname].append(one_dict[varname])

        else:
            exit()

        # go to next day
        currDate_D = currDate_D + datetime.timedelta(days=1)

    for varname in varname_list:
        out_dict[varname] = np.array(out_dict[varname])
        print(varname, out_dict[varname].shape)

    if region_limit is not None:

        lat_e_1D = out_dict['Latitude_e'][:,0]
        lat_min = region_limit[0]
        lat_max = region_limit[2]
        i1 = get_center_index(lat_e_1D, lat_min)
        i2 = get_center_index(lat_e_1D, lat_max)

        lon_e_1D = out_dict['Longitude_e'][0,:]
        lon_min = region_limit[1]
        lon_max = region_limit[3]
        j1 = get_center_index(lon_e_1D, lon_min)
        j2 = get_center_index(lon_e_1D, lon_max)

        for varname in varname_list:
            out_dict[varname] = out_dict[varname][:,i1:i2+1,j1:j2+1]

        out_dict['Latitude']    = out_dict['Latitude'][i1:i2+1,j1:j2+1]
        out_dict['Latitude_e']  = out_dict['Latitude_e'][i1:i2+2,j1:j2+2]
        out_dict['Longitude']   = out_dict['Longitude'][i1:i2+1,j1:j2+1]
        out_dict['Longitude_e'] = out_dict['Longitude_e'][i1:i2+2,j1:j2+2]

    # coefficient of variance
    if (sat_cv_thre is not None) or (sat_ave_thre is not None):
        sat_cv_flag = True
    if sat_cv_flag:
        out_dict[sat_cv_name + '_ave'] = \
                np.nanmean(out_dict[sat_cv_name], axis=0)
        out_dict[sat_cv_name + '_std'] = \
                np.nanstd(out_dict[sat_cv_name], axis=0)
        out_dict[sat_cv_name + '_cv'] = out_dict[sat_cv_name + '_std'] / \
                out_dict[sat_cv_name + '_ave']
        flag = np.full_like(out_dict[sat_cv_name + '_cv'], True)

    # sat_cv_thre
    if sat_cv_thre is not None:
        flag = np.logical_and(flag, \
                out_dict[sat_cv_name + '_cv'] > sat_cv_thre)

    # sat_ave_thre
    if sat_ave_thre is not None:
        flag = np.logical_and(flag, \
                out_dict[sat_cv_name + '_ave'] > sat_ave_thre)

    if (sat_cv_thre is not None) or (sat_ave_thre is not None):
        for varname in varname_list:
            out_dict[varname][:,np.logical_not(flag)] = np.nan




    return out_dict
#
#------------------------------------------------------------------------------
#
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

