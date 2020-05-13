"""
Created on January 1, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import sys

from mylib.constants import kg_NO_to_ng_N

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_plot import plot_panel_variables, layout_2

#######################
# Start user parameters
#

root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_runs/'

fig_dir = '../figure/'

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

month = '201806'

varname = 'EmisNO_Soil'

gc_run = 'merra2_2x25_tropchem'

subdir = 'OutputDir/'

scene_prefix='GC_'

scale = kg_NO_to_ng_N

vmin = 0.0
vmax = 20.0
units = r'[ng N m$^{-2}$ s$^{-1}$]'
cb_ticks = None
cb_ticklabels = None


vmax_diff = 10.0
vmin_diff = -vmax_diff
units_diff = r'[ng N m$^{-2}$ s$^{-1}$]'
cb_ticks_diff = None
cb_ticklabels_diff = None

#vmin = None
#vmax = None
#units = ''
#cb_ticks = None
#cb_ticklabels = None


#vmax_diff = None
#vmin_diff = None
#units_diff = ''
#cb_ticks_diff = None
#cb_ticklabels_diff = None

lw=0.5

title_list = \
        [
        'Control', \
        'Soil_T_old', \
        'Air_T_new', \
        'Soil_T_new', \
        ]

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
        gc_run=gc_run,
        subdir=subdir,
        scene_prefix=scene_prefix,
        scale=scale,
        lw=lw,
        title_list=title_list,
        cb_ticks=cb_ticks, cb_ticklabels=cb_ticklabels)
figname = fig_dir + varname + '_' + month + '.png'
plt.savefig(figname, format='png', dpi=300)

# emission differences
plot_panel_variables(root_dir, scene_tup, month,
        varname=varname,  vmin=vmin, vmax=vmax, units=units, 
        vmin_diff=vmin_diff, vmax_diff=vmax_diff, units_diff=units_diff,
        gc_run=gc_run,
        subdir=subdir,
        scene_prefix=scene_prefix,
        scale=scale,
        diff=True, layout=layout_2,
        lw=lw,
        title_list=diff_title_list,
        cb_ticks=cb_ticks, cb_ticklabels=cb_ticklabels,
        cb_ticks_diff=cb_ticks_diff, cb_ticklabels_diff=cb_ticklabels_diff,
        )
figname = fig_dir + varname + '_diff_' + month + '.png'
plt.savefig(figname, format='png', dpi=300)
