"""
Created on January 20, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

from mylib.io import read_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_io import read_VCD_multiscenes
from sn_plot import plot_two_month_diff
from sn_utility import generate_start_end

#######################
# Start user parameters
#

# usage:
# python main_plot_NO2_VCD.py month1 month2

month1 = sys.argv[1]
month2 = sys.argv[2]

file_str = 'model_satellite'

NO2_VCD_dir = '/Dedicated/jwang-data/ywang/soil_NOx/process/\
resample/data/monthly/'

filename1 = NO2_VCD_dir + file_str + '_' + \
        generate_start_end(month1) + '.nc'
filename2 = NO2_VCD_dir + file_str + '_' + \
        generate_start_end(month2) + '.nc'

fig_dir = '../figure/'

sat_var_list = ['sat_ColumnAmountNO2Trop']

mod_var_list = ['mod_TSOIL1']

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

AK = '_AK'
#AK = ''

region_limit = [40, -104, 50, -91]

mod_varname_str = 'mod_NO2Trop' + AK + '_tp_sat_'

mod_scene_var_list = [mod_varname_str]

#mod_sim_str = 'GC' + AK + '_'

#sat_sim_name = 'OMI'

xticks=np.arange(-180.0, 180.1, 5.0)
yticks=np.arange(-90.0, 90.1, 5.0)

min_VCD = 0.0
max_VCD = 4.0e15
label_VCD = r'[molec $cm_{-2}$]'

max_diff = 4e15
min_diff = -max_diff
label_diff = r'[molec $cm_{-2}$]'


lw=0.5

#
# End user parameters
#####################

# read_data
mon1_dict = read_VCD_multiscenes(filename1, scene_tup, 
        mod_scene_var_list=mod_scene_var_list, sat_var_list=sat_var_list,
        mod_var_list=mod_var_list)
mon2_dict = read_VCD_multiscenes(filename2, scene_tup, 
         mod_scene_var_list=mod_scene_var_list, sat_var_list=sat_var_list,
         mod_var_list=mod_var_list)

# plot
plot_two_month_diff(mon1_dict, mon2_dict, sat_var_list[0],
        mod_varname_str, scene_tup, region_limit=region_limit,
        vmin=min_VCD, vmax=max_VCD,
        diff_min=min_diff, diff_max=max_diff,
        xticks=xticks, yticks=yticks)

figname = fig_dir + file_str + AK + '_' + month1 + '_' + month2  + '.png'
plt.savefig(figname, format='png', dpi=300)
