"""
Created on March 28, 2020

@author: Yi Wang
"""

import os
import sys

#######################
# Start user parameters
#

# usage:
# python main_submit_daily_regrid.py start_year end_year

start_year = int(sys.argv[1])
end_year   = int(sys.argv[2])

#
# End user parameters
#####################



def replace_run_sh(in_file, out_file, year):
    """
    """

    # in
    fin = open(in_file, 'r')
    data = fin.readlines()
    fin.close()

    # out
    startDate = year + '-01-01' 
    endDate   = year + '-12-31'
    data[1]  = data[1].replace('regrid', 'y' + year + 'regrid')
    data[-1] = data[-1].replace('startDate', startDate)
    data[-1] = data[-1].replace('endDate',   endDate  )
    fout = open(out_file, 'w')
    for line in data:
        fout.write(line)
    fout.close()

cwd_dir = os.getcwd()
for year in range(start_year, end_year+1):

    year_c = str(year)
    print('proesss ' + year_c)

    # copy run
    os.chdir(cwd_dir)
    run_dir = cwd_dir + '/tmp/' + year_c
    if os.path.isdir(run_dir):
        os.system('rm -rf ' + run_dir)
    os.system('mkdir ' + run_dir)
    run_file_ori = './run_daily_regrid.sh'
    run_file = run_dir + \
            '/run_daily_regrid_' + year_c  + '.sh'
    replace_run_sh(run_file_ori, run_file, year_c)

    # submit
    os.chdir(run_dir)
    os.system('qsub ' + run_file)
