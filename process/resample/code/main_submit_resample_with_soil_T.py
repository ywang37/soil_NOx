"""
Created on March 28, 2020

@author: Yi Wang
"""

import os

#######################
# Start user parameters
#

month_list = ['06', '07', '08']

year_list = ['2016', '2017', '2018', '2019']

day_dict = {
        '06' : '30',
        '07' : '31',
        '08' : '31'
        }


#
# End user parameters
#####################



def replace_run_sh(in_file, out_file,  year, month):
    """
    """

    # in
    fin = open(in_file, 'r')
    data = fin.readlines()
    fin.close()

    # out
    startDate = year + '-' + month + '-01' 
    endDate   = year + '-' + month + '-' + day_dict[month]
    data[1]  = data[1].replace('soil_NOx', 'm' + year[2:4] + month + \
            '_soil_NOx')
    data[-1] = data[-1].replace('startDate', startDate)
    data[-1] = data[-1].replace('endDate',   endDate  )
    fout = open(out_file, 'w')
    for line in data:
        fout.write(line)
    fout.close()

cwd_dir = os.getcwd()
for year in year_list:
    for month in month_list:

        YYYYMM = year + month

        print('proesss ' + YYYYMM)

        # copy run
        os.chdir(cwd_dir)
        run_dir = cwd_dir + '/tmp/' + YYYYMM
        if os.path.isdir(run_dir):
            os.system('rm -rf ' + run_dir)
        os.system('mkdir ' + run_dir)
        run_file_ori = './run_resample_with_soil_T.sh'
        run_file = run_dir + \
                '/run_resample_with_soil_T_' + YYYYMM  + '.sh'
        replace_run_sh(run_file_ori, run_file, year, month)

        # submit
        os.chdir(run_dir)
        os.system('qsub ' + run_file)
