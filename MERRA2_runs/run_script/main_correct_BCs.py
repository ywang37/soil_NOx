"""
Created on August 26, 2020

@author: Yi Wang
"""

import datetime
import glob
import os

from mylib.model.gc_data import correct_BC

#######################
# Start user parameters
#

#year = 2005
#month = 8
#yyyymm_list = []
#for year in range(2006, 2014):
#    yyyymm_list.append(str(year) + '05')

#yyyymm_list = ['200905', '200906', '200907', 
#        '201205', '201206', '201207']
#yyyymm_list = ['200908', '201208']
yyyymm_list = ['200908', '201208']

# spinup month and day
spinup_month = '05'
spinup_day   = '16'

#scene_list = ['GC_ori', 'GC_soil_T_obs', 'GC_soil_T_ori', 'GC_surf_T_obs']
#scene_list = ['GC_ori', 'GC_soil_T_obs']
#scene_list = ['GC_ori_all_scale_0.50', 'GC_soil_T_obs_all_scale_0.50']
scene_list = ['GC_ori_all_scale_0.50', 'GC_soil_T_obs_all_scale_0.50',
        'GC_ori', 'GC_soil_T_obs']

root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_runs/'


#
# End user parameters
#####################

for yyyymm in yyyymm_list:

    yyyy = yyyymm[0:4]
    mm   = yyyymm[4:6]

    for scene in scene_list:

        base_path = root_dir + scene + '/' + 'runs/merra2_2x25_tropchem'

        if (mm == spinup_month):
            next_mm = str(int(mm) + 1).zfill(2)
            curr_file = base_path + '_spinup_' + yyyy + next_mm + '01/' + \
                    'OutputDir/' + \
                    'GEOSChem.BoundaryConditions.' + yyyy + mm + spinup_day + \
                    '_0000z.nc4'
            pre_file = None
        else:
            currDate_D = datetime.datetime.strptime(yyyymm + '01', '%Y%m%d')
            preDate_D  = currDate_D + datetime.timedelta(days=-1)
            preDate    = str(preDate_D)
            if (preDate[5:7] == spinup_month):
                pre_dir = 'spinup_' + yyyymm + '01'
            else:
                pre_dir = preDate[0:4] + preDate[5:7]
            curr_file = base_path + '_' + yyyymm + '/' + 'OutputDir/' + \
                    'GEOSChem.BoundaryConditions.' + yyyymm + '01_0000z.nc4'
            pre_file = base_path + '_' + pre_dir + '/' + 'OutputDir/' + \
                    'GEOSChem.BoundaryConditions.' + yyyymm + '01_0000z.nc4'

        # correct boundary condition files
        correct_BC(curr_file, pre_file)

