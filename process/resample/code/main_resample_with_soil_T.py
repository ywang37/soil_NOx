"""
Created on January 13, 2020

@author: Yi Wang
"""

import datetime
import glob
import numpy as np

from mylib.amf.amf import AMF_trop
from mylib.conversion import vmr_to_molec_cm2
from mylib.gc_io.met import get_geosfp_hourly_A1_3days
from mylib.gc_io.read_nd49 import read_nd49_resample
from mylib.pro_omi_no2_l2.io_omi_no2_l2 import read_OMI_NO2_L2
from mylib.pro_omi_no2_l2.pro_omi_no2_l2 import QC_OMI_NO2_L2
from mylib.pro_satellite.sat_model_sample import sat_model_sample
from mylib.pro_satellite.sat_model_sample import save_sat_model_sample

#######################
# Start user parameters
#

startDate = '2018-08-31'
endDate   = '2018-08-31'

yyyy = startDate[0:4]
mm   = startDate[5:7]
yyyymm = yyyy + mm

scene_list = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

gc_root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/'

gc_rel_path = 'runs/geosfp_2x25_tropchem_' + yyyymm  + '/ND49/'

sat_dir = '/Dedicated/jwang-data/shared_satData/OMI_NO2_L2/' + \
        yyyy + '/' + mm + '/'

met_root_dir = '/Dedicated/jwang-data/GCDATA/GEOS_2x2.5/GEOS_FP_soil_T/'
get_geosfp_A1 = get_geosfp_hourly_A1_3days

out_dir = '../data/granule/'


# species names
mod_spename = 'IJ-AVG-$_NO2'

# all model variable names
mod_varname_list = [mod_spename, 'TIME-SER_AIRDEN', \
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

# soil temperature
soil_T_flag = True

# calculate air mass factor
amf_flag = True

verbose = True

flag_2D = True
flag_1D = True

# use tropopause pressure from satellite data
tp_sat_flag = True

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
    # all dates
    p_date = preDate[0:4]  + preDate[5:7]  + preDate[8:10]
    c_date = currDate[0:4] + currDate[5:7] + currDate[8:10]
    n_date = nextDate[0:4] + nextDate[5:7] + nextDate[8:10]

    mod_var_dict = {}

    # soil temperature
    if soil_T_flag:
        no_pre  = False
        no_next = False
        if preDate[5:7] != currDate[5:7]:
            no_pre = True
        if nextDate[5:7] != currDate[5:7]:
            no_next = True
        soil_T = get_geosfp_A1(met_root_dir, c_date, ['TSOIL1'],
                no_pre=no_pre, no_next=no_next, verbose=verbose)
        mod_var_dict['TSOIL1'] = soil_T['TSOIL1']

    # get filenames
    for scene in scene_list:

        print(' get GC data from ' + scene)

        # data directory
        gc_dir = gc_root_dir + 'GEOS-Chem_' + scene + '/' + gc_rel_path

        # filenames for previous, current, and next day.
        pre_files  = gc_dir + 'ts' + p_date + '.bpch'
        curr_files = gc_dir + 'ts' + c_date + '.bpch'
        next_files = gc_dir + 'ts' + n_date + '.bpch'

        # read model data
        model_files = [pre_files, curr_files, next_files]
        if (scene == scene_list[0]):
            # read all model variables
            model_data = \
                    read_nd49_resample(model_files, mod_varname_list)
            # add scene suffix to species name
            model_data[mod_spename + '_' + scene] = \
                    model_data.pop(mod_spename)
        else:
            # only read a species
            model_species = \
                    read_nd49_resample(model_files, [mod_spename])
            # add species of a centain scene to *model_data*
            model_data[mod_spename + '_' + scene] = \
                    model_species.pop(mod_spename)
                    
    # unit conversion (ppbv => molec/cm2)
    # at every layer
    # and prepare model data for resmapling
    for scene in scene_list:
        mod_NO2 = vmr_to_molec_cm2(model_data[mod_spename + '_' + scene],
                model_data['TIME-SER_AIRDEN'], model_data['BXHGHT-$_BXHEIGHT'])
        mod_var_dict['NO2_' + scene] = mod_NO2

    # preprae additional data for resampling
    mod_TAI93 = model_data['TAI93']
    mod_var_dict['PEdge_Bot'] = model_data['PEDGE-$_PSURF']


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
    for i in range(len(all_sat_files)):
    #for i in [14]:

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

        # prepare satellite data for resmapling
        sat_lat = sat_data['Latitude']
        sat_lon = sat_data['Longitude']
        sat_TAI93 = np.tile(sat_data['Time'], sat_lat.shape[1])
        sat_TAI93 = sat_TAI93.reshape(sat_lat.shape[::-1])
        sat_TAI93 = sat_TAI93.T
        sat_obs_dict = {}
        sat_obs_dict['ColumnAmountNO2Trop'] = sat_data['ColumnAmountNO2Trop']
        sat_obs_dict['TropopausePressure'] = sat_data['TropopausePressure']
        sat_obs_dict['AmfTrop'] = sat_data['AmfTrop']
        sat_obs_dict['ScatteringWeight'] = sat_data['ScatteringWeight']

        # Sample model results according satellite observations
        # and regrid satellite observations to model grids.
        sat_mod_dict = \
                sat_model_sample(mod_coord_dict, mod_TAI93, mod_var_dict,
                sat_lat, sat_lon, sat_TAI93, sat_obs_dict,
                sat_flag=sat_flag)

        # skip if there is no satellite data that that is
        # in the region
        if sat_mod_dict['valid_sat'] == 0:
            print('There is no satellite data that that is in the region.')
            continue

        if amf_flag:

            # scattering weight pressure
            SW_AK_press = sat_data['ScatteringWtPressure']

            if flag_2D:

                # prepare data for AMF_trop function
                mod_grid_dict = sat_mod_dict['mod_grid_dict']
                sat_grid_dict = sat_mod_dict['sat_grid_dict']
                PEdge_Bot_2D = mod_grid_dict['PEdge_Bot']
                SW_AK_2D     = sat_grid_dict['ScatteringWeight']
                flag         = (sat_mod_dict['count'] > 0)

                # use tropopause pressure from satellite
                if tp_sat_flag:

                    # prepare data for AMF_trop function
                    ind_l_2D        = None
                    P_tropopause_2D = sat_grid_dict['TropopausePressure']

                    # loop different scenes
                    for scene in scene_list:
                        
                        # prepare species from AMF_trop function
                        layer_val_2D = mod_grid_dict['NO2_' + scene]

                        print('Calculate AMF at grids for ' + scene + 
                                ', use tropopause from satellite ')

                        # calculate AMF
                        amf_trop_data_2D_satp = AMF_trop(layer_val_2D,
                                PEdge_Bot_2D, SW_AK_2D, SW_AK_press,
                                ind_l_arr=ind_l_2D,
                                P_tropopause_arr=P_tropopause_2D,
                                var='AMF', flag=flag)

                        # add data to sat_mod_dict through
                        # mod_grid_dict and sat_grid_dict
                        mod_grid_dict['AmfTrop_tp_sat_' + scene] = \
                                amf_trop_data_2D_satp['AMF']
                        mod_grid_dict['NO2Trop_tp_sat_' + scene] = \
                                amf_trop_data_2D_satp['VCD']

                        # convert VCD to VCD_AK and save to
                        # mod_grid_dict
                        mod_grid_dict['NO2Trop_AK_tp_sat_' + scene] = \
                                mod_grid_dict['NO2Trop_tp_sat_' + scene] \
                                * mod_grid_dict['AmfTrop_tp_sat_' + scene] \
                                / sat_grid_dict['AmfTrop']


            if flag_1D:

                # prepare data for AMF_trop function
                mod_1D_dict = sat_mod_dict['mod_1D_dict']
                sat_1D_dict = sat_mod_dict['sat_1D_dict']
                PEdge_Bot_1D = mod_1D_dict['PEdge_Bot']
                SW_AK_1D     = sat_1D_dict['ScatteringWeight']
                flag         = None

                # use tropopause pressure from satellite
                if tp_sat_flag:

                    # prepare data for AMF_trop function
                    ind_l_1D        = None
                    P_tropopause_1D = sat_1D_dict['TropopausePressure']

                    # loop different scenes
                    for scene in scene_list:

                        layer_val_1D = mod_1D_dict['NO2_' + scene]

                        print('Calculate AMF at stations for ' + scene + 
                                ', use tropopause from satellite ')

                        # calculate AMF
                        amf_trop_data_1D_satp = AMF_trop(layer_val_1D,
                                PEdge_Bot_1D, SW_AK_1D, SW_AK_press,
                                ind_l_arr=ind_l_1D,
                                P_tropopause_arr=P_tropopause_1D,
                                var='AMF', flag=flag)

                        # add data to sat_mod_dict through
                        # mod_1D_dict and sat_1D_dict
                        mod_1D_dict['AmfTrop_tp_sat_' + scene] = \
                                amf_trop_data_1D_satp['AMF']
                        mod_1D_dict['NO2Trop_tp_sat_' + scene] = \
                                amf_trop_data_1D_satp['VCD']

                        # convert VCD to VCD_AK and save to
                        # mod_1D_dict
                        mod_1D_dict['NO2Trop_AK_tp_sat_' + scene] = \
                                mod_1D_dict['NO2Trop_tp_sat_' + scene] \
                                * mod_1D_dict['AmfTrop_tp_sat_' + scene] \
                                / sat_1D_dict['AmfTrop']


        # save data
        if ( sat_mod_dict['valid_sat'] > 0 ):
            lon, lat = np.meshgrid(model_data['longitude'], 
                    model_data['latitude'])
            sat_mod_dict['Latitude']  = lat
            sat_mod_dict['Longitude'] = lon
            lon_e, lat_e = \
                    np.meshgrid(model_data['longitude_e'], 
                            model_data['latitude_e'])
            sat_mod_dict['Latitude_e']   = lat_e
            sat_mod_dict['Longitude_e']  = lon_e
            sat_mod_dict['sat_ScatteringWtPressure'] = \
                    sat_data['ScatteringWtPressure']
            out_file = out_dir + 'model_satellite_' + \
                    sat_file.split('/')[-1][18:32] + '.nc'
            save_sat_model_sample(out_file, sat_mod_dict)


    # go to next day
    currDate_D = currDate_D + datetime.timedelta(days=1)

print('Done')
