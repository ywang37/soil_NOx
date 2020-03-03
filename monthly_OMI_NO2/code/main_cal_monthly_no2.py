"""
Created on September 16, 2019

@author: Yi Wang
"""

import datetime

from mylib.pro_omi_no2_l3.io_omi_no2_l3 import output_month_OMI_NO2_L3
from mylib.pro_omi_no2_l3.pro_omi_no2_l3 import cal_month_OMI_NO2_L3

#######################
# Start user parameters
#

start_month = '2019-01'
end_month   = '2019-12'
curr_month  = start_month

in_root_dir = '/Dedicated/jwang-data/ywang/OMI_NO2/level3/'

out_dir = '../data/monthly_mean/'

verbose = True

#
# End user parameters
#####################

curr_month_d = datetime.datetime.strptime(curr_month+'-01', '%Y-%m-%d')
end_month_d  = datetime.datetime.strptime(end_month+'-01',  '%Y-%m-%d')

while curr_month_d <= end_month_d:

    curr_month = str(curr_month_d)[0:7]
    yyyy = curr_month[0:4]
    mm = curr_month[5:7]
    yyyymm = yyyy + mm
    print('------------------------------------------------------')
    print('processing ' + yyyymm)

    # calculate monthly mean
    in_dir = in_root_dir + yyyy + '_ori/'
    NO2 = cal_month_OMI_NO2_L3(yyyymm, in_dir, verbose=verbose)

    # output monthly mean
    output_month_OMI_NO2_L3(out_dir, curr_month, NO2, verbose=verbose)

    # go to next month
    curr_month_d = \
            datetime.datetime(curr_month_d.year+(curr_month_d.month//12), 
                    ((curr_month_d.month%12)+1), 1)




