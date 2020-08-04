"""
Created on July 15, 2020

@author: Yi Wang
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

from mylib.cartopy_plot import cartopy_plot
from mylib.io import read_nc
from mylib.layout import h_1_ax

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
import sn_p

filename = '/Dedicated/jwang-data/GCDATA/ExtData/\
HEMCO/SOILNOX/v2014-07/soilNOx.climate.generic.05x05.nc'

# read data
varnames = ['lon', 'lat', 'NON_ARID', 'ARID']
data_dict = read_nc(filename, varnames, verbose=True, squeeze=True)
lon = data_dict['lon']
lat = data_dict['lat']
lon, lat = np.meshgrid(lon, lat)
dat = np.full_like(lon, np.nan)
dat[data_dict['NON_ARID'] > 0.5] = 0.5
dat[data_dict['ARID'] > 0.5]     = 1.5

# plot
xtick = (-120, -100, -80)
ytick = (30, 40, 50)
pout = cartopy_plot(lon, lat, dat, vmin=0.0, vmax=2.0,
        cmap=matplotlib.colors.ListedColormap(['C0', 'C1']),
        bad_c='white', cbar=False,
        xtick=xtick, ytick=ytick,
        region_limit=sn_p.region_dict['US'],
        countries=True, states=True)
plt.subplots_adjust(left=0.25, right=0.75)
cax = h_1_ax(pout['fig'], pout['ax'], y_off=-0.09)
cb = plt.colorbar(pout['mesh'], cax=cax, orientation='horizontal',
        ticks=[0.5, 1.5])
cb.ax.set_xticklabels(['Non-arid', 'Arid'])

plt.savefig('../figure/arid_region.png', format='png', dpi=300)
