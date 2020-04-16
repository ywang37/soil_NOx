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

varname = 'EmisNO_Soil'

vmin = 0.0
vmax = 2e11
units = r'[10$^{11}$ molec cm$^{-2}$ s$^{-1}$]'
cb_ticks = [0, 0.5e11, 1e11, 1.5e11, 2e11]
cb_ticklabels = ['0.0', '0.5', '1.0', '1.5', '2.0']


vmax_diff = 1e11
vmin_diff = -vmax_diff
units_diff = r'[10$^{11}$ molec cm$^{-2}$ s$^{-1}$]'
cb_ticks_diff = [-1.0e11, -0.5e11, 0, 0.5e11, 1e11]
cb_ticklabels_diff = ['-1.0', '-0.5', '0.0', '0.5', '1.0']

lw=0.5

diff_title_list = \
        [
        'Control', \
        'Soil_T_old - Control', \
        'Air_T_new - Control', \
        'Soil_T_new - Control', \
        ]

#
# End user parameters
#####################


# emissions
plot_panel_variables(root_dir, scene_tup, month,
        varname=varname, vmin=vmin, vmax=vmax, units=units,
        lw=lw,
        cb_ticks=cb_ticks, cb_ticklabels=cb_ticklabels)
figname = fig_dir + varname + '_' + month + '.png'
plt.savefig(figname, format='png', dpi=300)

# emission differences
plot_panel_variables(root_dir, scene_tup, month,
        varname=varname,  vmin=vmin, vmax=vmax, units=units, 
        vmin_diff=vmin_diff, vmax_diff=vmax_diff, units_diff=units_diff,
        diff=True, layout=layout_2,
        lw=lw,
        title_list=diff_title_list,
        cb_ticks=cb_ticks, cb_ticklabels=cb_ticklabels,
        cb_ticks_diff=cb_ticks_diff, cb_ticklabels_diff=cb_ticklabels_diff,
        )
figname = fig_dir + varname + '_diff_' + month + '.png'
plt.savefig(figname, format='png', dpi=300)
