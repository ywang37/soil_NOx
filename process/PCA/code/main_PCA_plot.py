from copy import deepcopy
import datetime
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
from sklearn.decomposition import PCA

from mylib.io import read_nc
from mylib.pca.plot import plot_2D_components_and_coeffs
from mylib.pca.plot import plot_ave
from mylib.pca.plot import plot_explained_variance_ratio


#######################
# Start user parameters
#


dataRootDir = '../data/PCA/'
figRootDir = '../figure/'

start_year = 2005
end_year   = 2019

scene = 'summer'
#scene = 'all_season'

land_ocean = 'land'

deseasonal = True

soil_emi_flag = True
soil_emi_ratio_thre = 0.5

res = '2x25'

n_components = 6

verbose = True

#
# End user parameters
#####################

# deseaonal
deseasonal_c = ''
if deseasonal:
    deseasonal_c = '_deseasonal'

# soil_emi_flag
soil_emi_c = ''
if soil_emi_flag:
    soil_emi_c = '_soil_emi{}'.format(soil_emi_ratio_thre)

if scene == 'all_season':
    xticks = np.array(range(0, (end_year-start_year+1)*12, 12))
    xticklablels = []
    for year in range(start_year, end_year+1):
        xticklablels.append(str(year) + '01')
elif scene == 'summer':
    xticks = np.array(range(0, (end_year-start_year+1)*3, 3))
    xticklablels = []
    for year in range(start_year, end_year+1):
        xticklablels.append(str(year) + '06')

time_ticks = (xticks, xticklablels)

title = 'PCA_' + str(start_year) + '-' + str(end_year) \
        + '_' + scene + '_' + land_ocean + '_' + res \
        + deseasonal_c + soil_emi_c

fig_dir = figRootDir + title + '/'
if not os.path.isdir(fig_dir):
    os.system('mkdir -p ' + fig_dir)

# read data
in_file = dataRootDir + title + '.nc'
varnames = ['Latitude_e', 'Longitude_e', 'Latitude', 'Longitude',
        'ori_data', 'ave', 'components', 'coeffs', 
        'explained_variance_ratio', 'explained_variance']
data_dict = read_nc(in_file, varnames, verbose=True)
lat_e = data_dict['Latitude_e']
lon_e = data_dict['Longitude_e']
components = data_dict['components']
coeffs = data_dict['coeffs']
explained_variance_ratio = data_dict['explained_variance_ratio']
ave = data_dict['ave']

# plot components and coeff
plot_2D_components_and_coeffs(components, lat_e, lon_e, coeffs,
        explained_variance_ratio,
        fig_dir, name=title, n_components=n_components,
        time_ticks=time_ticks)

# plot_explained_variance_ratio
plot_explained_variance_ratio(explained_variance_ratio[:n_components])
figname = fig_dir + title + '_contributions.png'
plt.savefig(figname, format='png', dpi=300)

# plot average
plot_ave(ave, lat_e, lon_e, units=r'[molec cm$^{-2}$]')
figname = fig_dir + title + '_ave.png'
plt.savefig(figname, format='png', dpi=300)
