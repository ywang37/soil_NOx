"""
Created on May 21, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import sys, getopt

from mylib.correlation import plot_pearsonr_map

#######################
# Start user parameters
#

data_dir = '../data/'

fig_dir = '../figure/'


r_vmin = -0.8
r_vmax = 0.8

varname_list = ['multimonth_mean', 'monthly',\
        'month_01', 'month_02', 'month_12']

#
# End user parameters
#####################

time_range = '2005-2019_winter'

#
argv = sys.argv[1:]
opts = 'tr'
long_opts = ['time_range=']
try:
    opts, args = getopt.getopt(argv, opts, long_opts)
except getopt.GetoptError:
    print('python main_plot_soil_T_OMI_L3_NO2_corr.py -tr <time_range> ')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-tr', '--time_range'):
        time_range = arg
#

name = 'corr_GWETTOP_NO2_Trop_CS_' + time_range

corr_file = data_dir + name + '.nc'

for varname in varname_list:
    plot_pearsonr_map(corr_file, fig_dir, varname,
            name=name, r_vmin=r_vmin, r_vmax=r_vmax)
