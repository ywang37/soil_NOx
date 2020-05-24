"""
Created on May 21, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import sys, getopt

from mylib.trend_analysis import plot_trend_map
from mylib.io import read_nc

#######################
# Start user parameters
#

data_dir = '../data/'

fig_dir = '../figure/'


mean_vmin = 0
mean_vmax = 1.0
mean_units = 'Soil moisture [unitless]'

trend_vmax = 0.02
trend_vmin = -trend_vmax
trend_units = r'Soil moisture trend [uniteless a$^{-1}$]'

sigma_units = \
        r'Standard deviation of soil moisture trend [unitless a$^{-1}$]'

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
    print('python main_plot_soil_moisture_trends.py -tr <time_range> ')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-tr', '--time_range'):
        time_range = arg
#

name = 'GWETTOP_' + time_range

trend_file = data_dir + name + '.nc'

# mask
tmp_data = read_nc(trend_file, ['mean'], verbose=True)
mean_mask = tmp_data['mean'] == 1.0


plot_trend_map(trend_file, fig_dir,
        name=name,
        mean_mask=mean_mask,
        mean_vmin=mean_vmin, mean_vmax=mean_vmax,
        mean_units=mean_units,
        trend_vmin=trend_vmin, trend_vmax=trend_vmax,
        trend_units=trend_units,
        sigma_units=sigma_units)
