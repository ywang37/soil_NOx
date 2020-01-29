"""
Created on January 1, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import sys

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_io import get_nc_NO_emi_ratio_multifiles
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
vmax = 0.5
units = ''
cb_ticks = None
cb_ticklabels = None

lw=0.5


#
# End user parameters
#####################

title_list = []
for secene in scene_tup:
    title_list.append(secene + '_NO_emi_ratio')


# read emissions


# emissions
plot_panel_variables(root_dir, scene_tup, month,
        varname=varname, read_func_varname='emi_ratio_dict',
        read_func=get_nc_NO_emi_ratio_multifiles,
        vmin=vmin, vmax=vmax, units=units,
        lw=lw, title_list=title_list,
        cb_ticks=cb_ticks, cb_ticklabels=cb_ticklabels)
figname = fig_dir + varname + '_ratio_' + month + '.png'
plt.savefig(figname, format='png', dpi=300)

