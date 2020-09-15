"""
Created on January 1, 2020

@author: Yi Wang
"""

import sys

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_process import daily_to_monthly_emissions

#######################
# Start user parameters
#

in_root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_runs/'

out_root_dir = '../data/'

scene_tup = [ \
        'soil_T_obs_fit_18_degree',
        'soil_T_obs_fit_18_degree_scale_0.60',
        'soil_T_obs_fit_18_degree_scale_0.70',
        'soil_T_obs_fit_18_degree_scale_0.80',
        'soil_T_obs_fit_OMI_NO2',
        'soil_T_obs_scale_0.60',
        'soil_T_obs_scale_0.65',
        'soil_T_obs_scale_0.70',
        'soil_T_obs_scale_0.75',
        'soil_T_obs_scale_0.80',
        ]

year_list = [2009, 2012]

month_list = [6, 7, 8]

#year_list = list(range(2005,2020))

#
# End user parameters
#####################


for year in year_list:
    for month in month_list:
        for scene in scene_tup:

            inDir = in_root_dir + 'GC_' + scene + '/runs/'

            outDir = out_root_dir + scene + '/'

            daily_to_monthly_emissions(inDir, outDir, year, month)
