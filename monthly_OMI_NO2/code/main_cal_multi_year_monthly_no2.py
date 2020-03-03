"""
Created on September 16, 2019

@author: Yi Wang
"""

import datetime

from mylib.pro_omi_no2_l3.io_omi_no2_l3 import read_month_OMI_NO2_L3
from mylib.pro_omi_no2_l3.pro_omi_no2_l3 import cal_multi_year_month_OMI_NO2_L3

#######################
# Start user parameters
#

start_year = 2005
end_year  = 2019

in_dir = '../data/monthly_mean/'

out_dir = '../data/multi_year_monthly_mean/'

month_list = [
        '01', '02', '03', '04', '05', '06',
        '07', '08', '09', '10', '11', '12'
        ]

verbose = True

#
# End user parameters
#####################

for month in month_list:

    print('------------------------------------------------------')
    print('processing ' + month)

    # Calculate multi-year mean of monthly OMI L3 NO2 and output.
    NO2 = cal_multi_year_month_OMI_NO2_L3(month, start_year, end_year, 
            in_dir, out_dir, verbose=verbose)
