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

out_root_dir = '../data/na_05x0625/'

scene_tup = ('soil_T_obs', 'ori',
        'ori_all_scale_0.50',
        'soil_T_obs_all_scale_0.30',
        'soil_T_obs_all_scale_0.40',
        'soil_T_obs_all_scale_0.50')

year_list = [2009, 2012]

month_list = [6, 7, 8]

#year_list = list(range(2005,2020))

res = 'na_05x0625'

#
# End user parameters
#####################


for year in year_list:
    for month in month_list:
        for scene in scene_tup:

            inDir = in_root_dir + 'GC_' + scene + '/runs/'

            outDir = out_root_dir + scene + '/'

            daily_to_monthly_emissions(inDir, outDir, year, month,
                    res=res)
