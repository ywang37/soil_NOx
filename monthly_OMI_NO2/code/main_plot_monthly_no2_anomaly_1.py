"""
Created on September 17, 2019

@author: Yi Wang
"""

import datetime

from mylib.pro_omi_no2_l3.plot_omi_no2_l3 import plot_NO2_combine

#######################
# Start user parameters
#

start_year = '2005'
end_year  = '2019'

mean_dir = '../data/multi_year_monthly_mean/'

month_dir = '../data/monthly_mean/'

fig_dir = '../figure/anomaly_1/'

year = '2012'
#year = '2006'

month_list = [
        '01', '02', '03', '04', '05', '06',
        '07', '08', '09', '10', '11', '12'
        ]
#month_list = ['01']

region_limit = [25.0, -126.0, 50.0, -66.0]

vmax_mean_dict = {
        '01' : 0.4,
        '02' : 0.4,
        '03' : 0.4,
        '04' : 0.3,
        '05' : 0.2,
        '06' : 0.2,
        '07' : 0.2,
        '08' : 0.2,
        '09' : 0.2,
        '10' : 0.3,
        '11' : 0.3,
        '12' : 0.4,
        }

vmax_ano_dict = {
        '01' : 0.1,
        '02' : 0.1,
        '03' : 0.1,
        '04' : 0.05,
        '05' : 0.05,
        '06' : 0.05,
        '07' : 0.05,
        '08' : 0.05,
        '09' : 0.05,
        '10' : 0.1,
        '11' : 0.1,
        '12' : 0.1,
        }

verbose = True

#
# End user parameters
#####################

for month in month_list:

    print('------------------------------------------------------')
    print('processing ' + month)

    vmax_mean = vmax_mean_dict[month]
    vmax_ano = vmax_ano_dict[month]
    vmin_ano = -vmax_ano

    # Calculate multi-year mean of monthly OMI L3 NO2 and output.
    plot_NO2_combine(mean_dir, month_dir, 
            start_year, end_year, year, month, verbose=verbose,
            dir_fig=fig_dir, region_limit=region_limit,
            vmax_mean=vmax_mean,
            vmin_ano=vmin_ano, vmax_ano=vmax_ano)


