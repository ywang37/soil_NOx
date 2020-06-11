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

NH_label_dict={}
NH_label_dict['multimonth_mean'] = 'NH summer (JJA)'
NH_label_dict['monthly'] = 'NH summer (JJA)'
NH_label_dict['month_06'] = 'NH June'
NH_label_dict['month_07'] = 'NH July'
NH_label_dict['month_08'] = 'NH August'

SH_label_dict={}
SH_label_dict['multimonth_mean'] = 'SH summer (DJF)'
SH_label_dict['monthly'] = 'SH summer (DJF)'
SH_label_dict['month_06'] = 'SH December'
SH_label_dict['month_07'] = 'SH January'
SH_label_dict['month_08'] = 'SH Feburary'

r_cmap = deepcopy(gbcwpry_map)
#r_cmap = plt.get_cmap('seismic')

#
# End user parameters
#####################

time_range = '2005-2019_combine_NH_SH_summer'

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

name = 'corr_TSOIL1_NO2_Trop_CS_' + time_range

corr_file = data_dir + name + '.nc'

for varname in varname_list:
    NH_label = NH_label_dict[varname]
    SH_label = SH_label_dict[varname]
    plot_pearsonr_map(corr_file, fig_dir, varname,
            name=name, r_vmin=r_vmin, r_vmax=r_vmax,
            r_cmap=deepcopy(r_cmap),
            equator=True, NH_label=NH_label,
            SH_label=SH_label)
