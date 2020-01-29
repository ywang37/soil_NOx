"""
Created on January 20, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import sys

from mylib.io import read_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_io import get_nc_NO_emi_ratio_multifiles
from sn_plot import plot_compare_4_to_1

#######################
# Start user parameters
#

file_str = 'model_satellite_2018-06-01_2018-06-30'

NO2_VCD_dir = '/Dedicated/jwang-data/ywang/soil_NOx/process/\
resample/data/monthly/'

emi_root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/'

month = '201806'

emi_ratio_name = 'EmisNO_Soil'

fig_dir = '../figure/'

sat_varname = 'sat_ColumnAmountNO2Trop'

sat_cv_varname = 'sat_ColumnAmountNO2Trop_cv'

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

AK = '_AK'
#AK = ''

mod_varname_str = 'mod_NO2Trop' + AK + '_tp_sat_'

mod_sim_str = 'GC' + AK + '_'

sat_sim_name = 'OMI'

min_VCD = 0.0
max_VCD = 1.0e16
label_VCD = r'[molec $cm_{-2}$]'

max_diff = 3e15
min_diff = -max_diff
label_diff = r'[molec $cm_{-2}$]'

max_ratio = 2
min_ratio = -max_ratio
label = ''

lw=0.5

#
# End user parameters
#####################

# emission ratio
emi_ratio = get_nc_NO_emi_ratio_multifiles(emi_root_dir,
        scene_tup, month, emi_ratio_name)
emi_ratio = emi_ratio['emi_ratio_dict']
emi_ratio = emi_ratio[scene_tup[0]]

# all variable to be read
mod_varname_list = []
for i in range(len(scene_tup)):
    mod_varname_list.append( mod_varname_str + scene_tup[i] )
all_varname_list = mod_varname_list + \
        ['Latitude', 'Longitude', 'Latitude_e', 'Longitude_e']
all_varname_list.append(sat_varname)
all_varname_list.append(sat_cv_varname)

mod_sim_name_dict = {}
for i in range(len(mod_varname_list)):
    mod_varname = mod_varname_list[i]
    mod_sim_name_dict[mod_varname] = mod_sim_str + scene_tup[i]

# read_data
NO2_VCD_file = NO2_VCD_dir + file_str + '.nc'
data_dict = read_nc(NO2_VCD_file, all_varname_list)

# sat_cv
sat_cv = data_dict[sat_cv_varname]

plot_compare_4_to_1(data_dict, sat_varname, mod_varname_list,
        sat_sim_name, mod_sim_name_dict, lw=lw,
        sat_cv=sat_cv,
        soil_ratio=emi_ratio,
        min_VCD=min_VCD, max_VCD=max_VCD,
        min_diff=min_diff, max_diff=max_diff,
        min_ratio=min_ratio, max_ratio=max_ratio)

figname = fig_dir + file_str + AK + '.png'
plt.savefig(figname, format='png', dpi=300)
