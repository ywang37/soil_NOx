"""
Created on March 28, 2020

@author: Yi Wang
"""

import os

#######################
# Start user parameters
#

month_list = ['06', '07', '08']

#year_list = ['2016', '2017', '2018', '2019']
#year_list = ['2010', '2012']
year_list = ['2009', '2012']

#year_list = []
#for year in range(2005, 2020):
#    year_list.append(str(year))

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
        cmd = 'python main_granule_to_daily_sensitivity_2.py ' + \
                startDate + ' ' + endDate
        print(cmd)
        os.system(cmd)
