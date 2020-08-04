"""
Created on January 1, 2020

@author: Yi Wang
"""

import sys, getopt

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_process import average_of_multi_month_emissions

#######################
# Start user parameters
#

in_root_dir = '../data/'

out_root_dir = '../data/'

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

start_year = 2005
end_year = 2019

start_month = 6
end_month = 8


#
# End user parameters
#####################

argv = sys.argv[1:]

opts = 'sy:ey:sm:em'
long_opts = ['start_year=', 'end_year=', 'start_month=', 'end_month=']

try:
    opts, args = getopt.getopt(argv, opts, long_opts)
except getopt.GetoptError:
    print('python main_calc_ave_emi.py -sy <start_year> ' +
            '-ey <end_year> -sm <start_month> -em <end_month>')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-sy', '--start_year'):
        start_year = int(arg)
    if opt in ('-ey', '--end_year'):
        end_year = int(arg)
    if opt in ('-sm', '--start_month'):
        start_month = int(arg)
    if opt in ('-em', '--end_month'):
        end_month = int(arg)



for scene in scene_tup:

    inDir = in_root_dir + scene + '/'

    outDir = out_root_dir + scene + '/'

    average_of_multi_month_emissions(inDir, outDir, 
            start_year, end_year,
            start_month, end_month, res='2x2.5',
            month_day_flag=True)
