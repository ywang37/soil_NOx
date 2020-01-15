"""
Created on January 13, 2020

@author: Yi Wang
"""

import datetime
import glob
import numpy as np

from mylib.conversion import vmr_to_molec_cm2
from mylib.gc_io.read_nd49 import read_nd49_resample
from mylib.pro_omi_no2_l2.io_omi_no2_l2 import read_OMI_NO2_L2
from mylib.pro_omi_no2_l2.pro_omi_no2_l2 import QC_OMI_NO2_L2
from mylib.pro_satellite.sat_model_sample import sat_model_sample
from mylib.pro_satellite.sat_model_sample import save_sat_model_sample

#######################
# Start user parameters
#


gc_dir = '/Dedicated/jwang-data/ywang/soil_NOx/GEOS-Chem_ori/\
runs/geosfp_2x25_tropchem_201806/ND49/'

sat_dir = '/Dedicated/jwang-data/shared_satData/OMI_NO2_L2/2018/06/'

out_dir = '../data/'

startDate = '2018-06-28'
endDate   = '2018-06-28'




mod_varname_list = ['IJ-AVG-$_NO2', 'TIME-SER_AIRDEN', \
        'BXHGHT-$_BXHEIGHT', 'PEDGE-$_PSURF']

# mod_coord_dict
mod_coord_dict = {}
mod_coord_dict['coord_format']  = 'ses'
mod_coord_dict['mod_lat_start'] =  -91.0
mod_coord_dict['mod_lat_end']   =   91.0
mod_coord_dict['mod_lat_step']  =    2.0
mod_coord_dict['mod_lon_start'] = -181.25
mod_coord_dict['mod_lon_end']   =  178.75
mod_coord_dict['mod_lon_step']  =    2.5

sat_varname_list = [\
        '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/TropopausePressure'
        ]


verbose = True

#
# End user parameters
#####################

# Date
currDate   = startDate
currDate_D = datetime.datetime.strptime(currDate, '%Y-%m-%d')
endDate_D  = datetime.datetime.strptime(endDate,  '%Y-%m-%d')


while currDate_D <= endDate_D:

    # current date
    currDate = str(currDate_D)[0:10]
    print(''.join(np.full((79,), '-')))
    print('processing ' + currDate)

    # A satelite file may span two days. Thus, we need to
    # find model data from previous day and next day
    preDate_D  = currDate_D + datetime.timedelta(days=-1)
    nextDate_D = currDate_D + datetime.timedelta(days=1)
    preDate  = str(preDate_D)
    nextDate = str(nextDate_D)

    # all model files
    p_date = preDate[0:4]  + preDate[5:7]  + preDate[8:10]
    c_date = currDate[0:4] + currDate[5:7] + currDate[8:10]
    n_date = nextDate[0:4] + nextDate[5:7] + nextDate[8:10]
    pre_files  = gc_dir + 'ts' + p_date + '.bpch'
    curr_files = gc_dir + 'ts' + c_date + '.bpch'
    next_files = gc_dir + 'ts' + n_date + '.bpch'

    # read model data
    model_files = [pre_files, curr_files, next_files]
    model_data = read_nd49_resample(model_files, mod_varname_list)

    print(model_data['IJ-AVG-$_NO2'].shape)

    for i in range(len(model_data['Time'])):
        print(i, model_data['Time'][i])

    # find satellite files
    sat_wildcard = sat_dir + 'OMI-Aura_L2-OMNO2_' + currDate[0:4] \
            + 'm' + currDate[5:7] + currDate[8:10] + '*_v003-*.he5'
    all_sat_files = glob.glob(sat_wildcard)
    all_sat_files.sort()
    print('All satellite files on ' + currDate)
    for i in range(len(all_sat_files)):
        print('  ' + all_sat_files[i])

    # process satellite files
    print('Process satellites on ' + currDate)
    #for i in range(len(all_sat_files)):
    for i in [3]:

        # read satellite file
        sat_file = all_sat_files[i]
        print('  reading ' + sat_file)
        sat_data = read_OMI_NO2_L2(sat_file, varnames=sat_varname_list,
                verbose=verbose)

        # quality control
        sat_NO2_trop = sat_data['ColumnAmountNO2Trop']
        sza          = sat_data['SolarZenithAngle']
        vza          = sat_data['ViewingZenithAngle']
        CF           = sat_data['CloudFraction']
        XTtrackQ     = sat_data['XTrackQualityFlags']
        vcdQ         = sat_data['VcdQualityFlags']
        TerrRef      = sat_data['TerrainReflectivity']
        sat_flag = QC_OMI_NO2_L2(sat_NO2_trop, sza, vza, CF, 
                XTtrackQ, vcdQ, TerrRef)

        # unit conversion (ppbv => molec/cm2)
        # at every layer
        mod_NO2 = vmr_to_molec_cm2(model_data['IJ-AVG-$_NO2'],
                model_data['TIME-SER_AIRDEN'], model_data['BXHGHT-$_BXHEIGHT'])


        # Sample model results according satellite observations
        # and regrid satellite observations to model grids.
        mod_TAI93 = model_data['TAI93']
        mod_var_dict = {}
        mod_var_dict['NO2'] = mod_NO2
        mod_var_dict['PEdge_Bot'] = model_data['PEDGE-$_PSURF']
        sat_lat = sat_data['Latitude']
        sat_lon = sat_data['Longitude']
        sat_TAI93 = np.tile(sat_data['Time'], sat_lat.shape[1])
        sat_TAI93 = sat_TAI93.reshape(sat_lat.shape[::-1])
        sat_TAI93 = sat_TAI93.T
        print(sat_TAI93.shape)
        print(sat_TAI93[:,0])
        print(sat_TAI93[1,:])
        sat_obs_dict = {}
        sat_obs_dict['ColumnAmountNO2Trop'] = sat_data['ColumnAmountNO2Trop']
        sat_obs_dict['TropopausePressure'] = sat_data['TropopausePressure']
        sat_obs_dict['AmfTrop'] = sat_data['AmfTrop']
        sat_obs_dict['ScatteringWeight'] = sat_data['ScatteringWeight']
        sat_mod_dict = \
                sat_model_sample(mod_coord_dict, mod_TAI93, mod_var_dict,
                sat_lat, sat_lon, sat_TAI93, sat_obs_dict,
                sat_flag=sat_flag)

        # save data
        lon, lat = np.meshgrid(model_data['longitude'], model_data['latitude'])
        sat_mod_dict['Latitude']  = lat
        sat_mod_dict['Longitude'] = lon
        lon_e, lat_e = \
                np.meshgrid(model_data['longitude_e'], 
                        model_data['latitude_e'])
        sat_mod_dict['Latitude_e']   = lat_e
        sat_mod_dict['Longitude_e']  = lon_e
        out_file = out_dir + 'model_satellite_' + \
                sat_file.split('/')[-1][18:32] + '.nc'
        save_sat_model_sample(out_file, sat_mod_dict)


    # go to next day
    currDate_D = currDate_D + datetime.timedelta(days=1)

print('Done')
