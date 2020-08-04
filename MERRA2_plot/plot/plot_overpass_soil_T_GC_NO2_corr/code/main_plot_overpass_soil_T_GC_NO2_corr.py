"""
Created on May 21, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import sys, getopt

from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.correlation import plot_pearsonr_map

#######################
# Start user parameters
#

data_dir = '../data/'

fig_dir = '../figure/'


r_vmin = -0.8
r_vmax = 0.8

varname_list = ['multimonth_mean', 'monthly',\
        'month_06', 'month_07', 'month_08']

r_signi_cmap = deepcopy(gbcwpry_map)

varn_list = ['mod_NO2Trop_AK_tp_sat_ori',
             'mod_NO2Trop_AK_tp_sat_soil_T_ori',
             'mod_NO2Trop_AK_tp_sat_surf_T_obs',
             'mod_NO2Trop_AK_tp_sat_soil_T_obs',
             'mod_NO2Trop_tp_sat_ori',
             'mod_NO2Trop_tp_sat_soil_T_ori',
             'mod_NO2Trop_tp_sat_surf_T_obs',
             'mod_NO2Trop_tp_sat_soil_T_obs']

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
    print('python main_plot_soil_T_OMI_L3_NO2_corr.py -tr <time_range> ')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-tr', '--time_range'):
        time_range = arg
#

for varn in varn_list:

    name = 'corr_overpass_mod_Met_TSOIL1_' + varn + '_' + time_range

    corr_file = data_dir + name + '.nc'

    for varname in varname_list:
        plot_pearsonr_map(corr_file, fig_dir, varname,
                name=name, r_vmin=r_vmin, r_vmax=r_vmax,
                r_signi_cmap=r_signi_cmap)
