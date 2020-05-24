"""
Created on May 21, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import sys, getopt

from mylib.trend_analysis import plot_trend_map

#######################
# Start user parameters
#

data_dir = '../data/'

fig_dir = '../figure/'


mean_vmin = -30.0
mean_vmax = 30.0
mean_units = u'2-m temperature [\u00B0C]'

trend_vmax = 0.2
trend_vmin = -trend_vmax
trend_units = u'2-m temperature trend [\u00B0C '  r'a$^{-1}$]'

sigma_units = \
        u'Standard deviation of 2-m temperature trend [\u00B0C '  r'a$^{-1}$]'

#
# End user parameters
#####################

time_range = '2005-2019_06-08'

#
argv = sys.argv[1:]
opts = 'tr'
long_opts = ['time_range=']
try:
    opts, args = getopt.getopt(argv, opts, long_opts)
except getopt.GetoptError:
    print('python main_plot_soil_T_trends.py -tr <time_range> ')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-tr', '--time_range'):
        time_range = arg
#

name = 'T2M_' + time_range

trend_file = data_dir + name + '.nc'


plot_trend_map(trend_file, fig_dir,
        name=name,
        mean_vmin=mean_vmin, mean_vmax=mean_vmax,
        mean_units=mean_units,
        trend_vmin=trend_vmin, trend_vmax=trend_vmax,
        trend_units=trend_units,
        sigma_units=sigma_units)
