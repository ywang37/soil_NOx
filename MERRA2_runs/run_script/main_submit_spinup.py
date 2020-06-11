"""
Created on May 27, 2019

@author: Yi Wang
"""

import glob
import os

#######################
# Start user parameters
#

year = 2009

scene_list = ['GC_ori', 'GC_soil_T_obs', 'GC_soil_T_ori', 'GC_surf_T_obs']
#scene_list = ['GC_ori']

root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_runs/'

queue_dict = {}
queue_dict['GC_ori']        = 'CGRER'
queue_dict['GC_soil_T_obs'] = 'CGRER'
queue_dict['GC_soil_T_ori'] = 'CGRER'
queue_dict['GC_surf_T_obs'] = 'CGRER'
#queue_dict['GC_ori']        = 'INFORMATICS'
#queue_dict['GC_soil_T_obs'] = 'INFORMATICS'
#queue_dict['GC_soil_T_ori'] = 'INFORMATICS'
#queue_dict['GC_surf_T_obs'] = 'INFORMATICS'


#
# End user parameters
#####################

year_c = str(year)

mmdd = '0601'

for scene in scene_list:

    ori_dir = root_dir + scene + '/' + 'runs/merra2_2x25_tropchem_spinup'

    new_dir = ori_dir + '_' + year_c + mmdd

    # copy directroy
    if not os.path.isdir(new_dir):
        os.system('cp -r ' + ori_dir + ' ' + new_dir)

    ##################
    # edit input.geos
    ##################
    input_file = new_dir + '/input.geos'

    line3_ori = '"Start YYYYMMDD, hhmmss  : 20180516 000000"'
    line3_new = '"Start YYYYMMDD, hhmmss  : ' + year_c + '0516 000000"'
    line3_cmd = 'sed -i 3s/' + line3_ori + '/' + line3_new + '/ ' + input_file
    os.system(line3_cmd)

    line4_ori = '"End   YYYYMMDD, hhmmss  : 20180601 000000"'
    line4_new = '"End   YYYYMMDD, hhmmss  : ' + year_c + '0601 000000"'
    line4_cmd = 'sed -i 4s/' + line4_ori + '/' + line4_new + '/ ' + input_file
    os.system(line4_cmd)

    ###########
    # edit run
    ###########
    run_file = new_dir + '/run'

    f_run = open(run_file, 'r')
    run_lines = f_run.readlines()
    f_run.close()

    # queue
    queue = queue_dict[scene]
    run_lines[3] = '#$ -q ' + queue + '\n'

    # working directory
    wkdir = 'cd ' + new_dir + '\n'
    run_lines[13] = wkdir

    f_new_run = open(run_file, 'w+')
    for line in run_lines:
        f_new_run.write(line)
    f_new_run.close()

    ###############
    # restart file
    ###############
    os.chdir(new_dir)
    os.system('mv GEOSChem.Restart.20180516_0000z.nc4 ' + \
            'GEOSChem.Restart.' + year_c + '0516_0000z.nc4')

    #############
    # submit job
    #############
    os.system('qsub run')



