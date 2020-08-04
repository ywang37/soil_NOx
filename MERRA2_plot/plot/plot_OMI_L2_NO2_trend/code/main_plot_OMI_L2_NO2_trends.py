"""
Created on May 21, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import sys, getopt

from mylib.colormap.emphasize_small_map import emphasize_small_map
from mylib.trend_analysis import plot_trend_map

#######################
# Start user parameters
#

# python main_plot_OMI_L2_NO2_trends.py --time_range=2005-2019_06-08

data_dir = '../data/'

fig_dir = '../figure/'

mean_cmap = deepcopy(emphasize_small_map)
mean_vmin = 0
mean_vmax = 8
mean_units = r'OMI L2 NO$_2$ [10$^{15}$  molec cm$^{-2}$]'

trend_vmax = 0.1
trend_vmin = -trend_vmax
trend_units = r'OMI L2 NO$_2$ trend [10$^{15}$ molec cm$^{-2}$ a$^{-1}$]'

sigma_units = r'Standard deviation of OMI L2 NO$_2$ trend ' + \
        r'[10$^{15}$ molec cm$^{-2}$ a$^{-1}$]'

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
    print('python main_plot_OMI_L2_NO2_trends.py -tr <time_range> ')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-tr', '--time_range'):
        time_range = arg
#

name = 'OMI_L2_NO2_' + time_range

trend_file = data_dir + name + '.nc'


plot_trend_map(trend_file, fig_dir,
        name=name,
        mean_vmin=mean_vmin, mean_vmax=mean_vmax,
        mean_units=mean_units, mean_cmap=mean_cmap,
        trend_vmin=trend_vmin, trend_vmax=trend_vmax,
        trend_units=trend_units,
        sigma_units=sigma_units)
