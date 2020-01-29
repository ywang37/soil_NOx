"""
Created on January 1, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import sys

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_plot import plot_panel_variables, layout_2

#######################
# Start user parameters
#

root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/'

fig_dir = '../figure/'

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

month = '201806'

varname_Soil = 'EmisNO_Soil'

varname_Total = 'EmisNO_Total'

tot_vmin = 0.0
tot_vmax = 8e11
tot_units = r'[10$^{11}$ molec cm$^{-2}$ s$^{-1}$]'
#tot_cb_ticks = [0, 2e11, 4e11, 6e11, 8e11]
#tot_cb_ticklabels = ['0', '2', '4', '6', '8']
tot_cb_ticks = None
tot_cb_ticklabels = None

lw=0.5


#
# End user parameters
#####################


# read emissions


# emissions
plot_panel_variables(root_dir, scene_tup, month,
        varname=varname_Total, 
        vmin=tot_vmin, vmax=tot_vmax, units=tot_units,
        lw=lw,
        cb_ticks=tot_cb_ticks, cb_ticklabels=tot_cb_ticklabels)
figname = fig_dir + varname_Total + '_' + month + '.png'
plt.savefig(figname, format='png', dpi=300)

