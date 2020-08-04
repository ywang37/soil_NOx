"""
Created on January 1, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyval

from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.layout import multiFigure, h_2_ax, panel_tick_label
from mylib.pro_satellite.pro_satellite import calculate_pixel_edge2

# Hudman's parameters
e_coef = 0.103

# Olikawa's parameters
a = -0.009
b = 0.837
c = -22.527
d = 196.149

temp_del = 0.5

# temperature for Hudman's
temp1 = np.arange(0, 40.01, temp_del)
H_temp = np.array(temp1)
H_temp[H_temp>=30.0] = 30.0
H_temp_func = np.exp(e_coef * H_temp)

# temperature for Olikawa's
temp2 = np.arange(20.0, 40.01, temp_del)
O_temp = np.array(temp2)
O_temp_func = polyval(O_temp, [d,c,b,a])
ind = np.where(H_temp < 20.0)
O_temp_func = np.hstack([H_temp_func[ind], O_temp_func])

# soil moisture
arid_a, arid_b = 8.24, 12.5
nonarid_a, nonarid_b = 5.5, 5.55
moist_del = 0.025
moist = np.array(np.arange(0, 1.001, moist_del))
arid_func    =    arid_a * moist * np.exp(   -arid_b * moist * moist)
nonarid_func = nonarid_a * moist * np.exp(-nonarid_b * moist * moist)


# temp_moisture
temp_moist_name_list = ['old_arid', 'old_nonarid', 'new_arid', 'new_nonarid']
temp_moist_dict = {}
temp_moist_dict['old_arid']    = np.outer(arid_func,    H_temp_func)
temp_moist_dict['old_nonarid'] = np.outer(nonarid_func, H_temp_func)
temp_moist_dict['new_arid']    = np.outer(arid_func,    O_temp_func)
temp_moist_dict['new_nonarid'] = np.outer(nonarid_func, O_temp_func)

# title
title_dict = {}
title_dict['old_arid'] = 'Arid soil, old T scheme'
title_dict['old_nonarid'] = 'Nonarid soil, old T scheme'
title_dict['new_arid'] = 'Arid soil, new T scheme'
title_dict['new_nonarid'] = 'Nonarid soil, new T scheme'

xx, yy = np.meshgrid(temp1, moist*100)
xlim = [xx[0,0], xx[0,-1]]
ylim = [yy[0,0], yy[-1,0]]
yy, xx = calculate_pixel_edge2(yy, xx)

# plot
vmin = 0.0
vmax = 60.0
p_dict = multiFigure(2, 2, bottom=0.2, top=0.92, hspace=0.2, wspace=0.05)
ax_list = p_dict['axes']
mesh_list = []
for i in range(len(temp_moist_name_list)):

    ax = ax_list[i]

    temp_moist_name = temp_moist_name_list[i]

    ax.set_title(title_dict[temp_moist_name])

    ax.set_xticks([])
    ax.set_yticks([])

    temp_moist = temp_moist_dict[temp_moist_name]

    mesh = ax.pcolormesh(xx, yy, temp_moist, cmap=deepcopy(WhGrYlRd_map),
            vmin=vmin, vmax=vmax)
    mesh_list.append(mesh)


panel_tick_label(ax_list, 2, xticks=np.arange(0, 40.1, 5),
        yticks=np.arange(0, 100.1, 20), xlabel=u'Soil temperature [\u00B0C]',
        ylabel='Soil water-filled pore space [%]')

for ax in ax_list:
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

cax = h_2_ax(p_dict['fig'], ax_list[2], ax_list[3], y_off=-0.11)
cbar = plt.colorbar(mesh_list[3], cax=cax, orientation='horizontal')
cbar.set_label(r'Soil NO$_x$ emissions dependence on ' + \
        'soil moisture and temperature')

plt.savefig('../figure/temp_func_moist_func.png', format='png', dpi=300)



