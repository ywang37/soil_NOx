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


cwd_dir = os.getcwd()
for year in year_list:
    for month in month_list:

        YYYYMM = year + month

        print('proesss ' + YYYYMM)

        # submit
        startDate = year + '-' + month + '-01'
        endDate   = year + '-' + month + '-' + day_dict[month]
        cmd = 'python main_plot_NO2_VCD.py ' + \
                'model_satellite_' + startDate  +'_' + endDate + ' ' + \
                YYYYMM
        print(cmd)
        os.system(cmd)
