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
from sn_plot import plot_NO2_T
from sn_utility import generate_start_end

#######################
# Start user parameters
#

# usage:
# python main_plot_NO2_VCD.py month

month = sys.argv[1]

file_str = 'model_satellite'

NO2_VCD_dir = '/Dedicated/jwang-data/ywang/soil_NOx/process/\
resample/data/monthly/'

filename = NO2_VCD_dir + file_str + '_' + \
        generate_start_end(month) + '.nc'

fig_dir = '../figure/'

sat_var_list = ['sat_ColumnAmountNO2Trop']

mod_var_list = ['mod_TSOIL1']

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

T_name = 'mod_TSOIL1'

AK = '_AK'
#AK = ''

region_limit = [20, 65, 40, 95]

mod_varname_str = 'mod_NO2Trop' + AK + '_tp_sat_'

mod_scene_var_list = [mod_varname_str]

#mod_sim_str = 'GC' + AK + '_'

#sat_sim_name = 'OMI'

xticks=np.arange(-180.0, 180.1, 10.0)
yticks=np.arange(-90.0, 90.1, 5.0)

min_VCD = 0.0
max_VCD = 1.0e16
label_VCD = r'[molec $cm_{-2}$]'

max_diff = 4e15
min_diff = -max_diff
label_diff = r'[molec $cm_{-2}$]'


lw=0.5

#
# End user parameters
#####################

# read_data
data_dict = read_VCD_multiscenes(filename, scene_tup, 
        mod_scene_var_list=mod_scene_var_list, sat_var_list=sat_var_list,
        mod_var_list=mod_var_list)

# plot
plot_NO2_T(data_dict, sat_var_list[0],
        mod_varname_str, scene_tup, T_name=T_name)

figname = fig_dir + 'India_NO2_T_' + \
        file_str + AK + '_' + month + '.png'
plt.savefig(figname, format='png', dpi=300)
