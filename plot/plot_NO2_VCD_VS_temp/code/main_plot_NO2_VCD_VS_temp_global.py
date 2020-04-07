"""
Created on January 20, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

from mylib.io import read_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_io import read_VCD_multiscenes_days, get_nc_NO_emi_ratio_multifiles
import sn_p
from sn_plot import plot_NO2_VS_T, plot_NO2_CV_emi_ratio
from sn_utility import generate_start_end

#######################
# Start user parameters
#

# usage:
# python main_plot_NO2_VCD_VS_temp.py startDate endDate emi_month region \
#        ratio_thre

startDate = sys.argv[1]
endDate = sys.argv[2]

emi_month = sys.argv[3]

region_name = sys.argv[4]

emi_ratio_thre = float(sys.argv[5])

emi_value_thre = 1e10

file_str0 = 'model_satellite'
file_str = file_str0 + '_YYYY-MM-DD.nc'

NO2_VCD_dir = '/Dedicated/jwang-data/ywang/soil_NOx/process/\
resample/data/daily/'

emi_root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/'

fig_dir = '../figure/'

sat_var_list = ['sat_ColumnAmountNO2Trop']

mod_var_list = ['mod_TSOIL1']

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

T_name = 'mod_TSOIL1'

emi_ratio_name = 'EmisNO_Soil'

AK = '_AK'
#AK = ''


region_limit = sn_p.region_dict[region_name]

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

xticks=np.arange(-180.0, 180.1, 60)
yticks=np.arange(-90.0, 90.1, 30)

lw=0.5

#
# End user parameters
#####################

# read NO2 VCD data
data_dict = read_VCD_multiscenes_days(file_str, NO2_VCD_dir, 
        startDate, endDate, scene_tup, 
        mod_scene_var_list=mod_scene_var_list, sat_var_list=sat_var_list,
        mod_var_list=mod_var_list,
        region_limit=region_limit,
        )

# read emisssions data
emi_dict = get_nc_NO_emi_ratio_multifiles(emi_root_dir,
        scene_tup, emi_month, emi_ratio_name,region_limit=region_limit)
emi_ratio = emi_dict['emi_ratio_dict']
emi_ratio = emi_ratio[scene_tup[0]]
emi_ratio_dict = {}
emi_ratio_dict['emi_ratio'] = emi_ratio
emi_ratio_dict['lat_e'] = emi_dict['lat_e']
emi_ratio_dict['lon_e'] = emi_dict['lon_e']
emi_ratio_dict['title'] = emi_month + ' soil NO emi ratio'
emi_soil = emi_dict[emi_ratio_name+'_dict']
emi_soil = emi_soil[scene_tup[0]]
emi_value_flag = emi_soil > emi_value_thre
emi_ratio_dict['emi_ratio'][np.logical_not(emi_value_flag)] = np.nan

# plot_NO2_CV
plot_NO2_CV_emi_ratio(data_dict, emi_dict=emi_ratio_dict)
figname = fig_dir + region_name  +'_NO2_CV_' + \
        file_str0 + AK + '_' + startDate + '_' + endDate + \
        '_ratio_' + str(emi_ratio_thre) + '.png'
plt.savefig(figname, format='png', dpi=300)

# plot_NO2_T_series
plot_NO2_VS_T(data_dict, sat_var_list[0],
        mod_varname_str, scene_tup, T_name=T_name,
        emi_ratio=emi_ratio_dict['emi_ratio'],
        emi_ratio_thre=emi_ratio_thre,
        xticks=xticks, yticks=yticks,
        )
figname = fig_dir + region_name + '_NO2_VS_T_' + \
        file_str0 + AK + '_' + startDate + '_' + endDate + \
        '_ratio_' + str(emi_ratio_thre) + '.png'
plt.savefig(figname, format='png', dpi=300)

