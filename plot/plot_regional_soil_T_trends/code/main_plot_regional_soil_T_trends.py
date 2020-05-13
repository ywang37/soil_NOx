"""
Created on April 22, 2020

@author: Yi Wang
"""

from calendar import monthrange
import matplotlib.pyplot as plt
import numpy as np
import sys, getopt

from mylib.grid_utility import generate_grid_gc_2x25, get_center_index_latlon
from mylib.io import read_nc
from mylib.pca.plot import plot_2D_components_and_coeffs

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_io import read_NO_emissions
import sn_p
from sn_plot import plot_ave_series

#######################
# Start user parameters
#


pca_file = '/Dedicated/jwang-data/ywang/soil_NOx/process/PCA/data/PCA/\
PCA_2005-2019_summer_land_2x25_deseasonal_soil_emi0.3.nc'

root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/process/soil_T/data/monthly/'

area_file = '/Dedicated/jwang-data/ywang/soil_NOx/GEOS-Chem_ori/runs/' \
        + 'geosfp_2x25_tropchem_201406/HEMCO_diagnostics.201406010000.nc'

fig_dir = '../figure/'

region_list = ['Kazakhstan', 'India', 'Sahel', 'US']

units = u'Soil temperature [\u00B0C]'

verbose = True

pca_st_yr = 2005
pca_ed_yr = 2019
pca_xticks = np.array(range(0, (pca_ed_yr-pca_st_yr+1)*3, 3))
pca_xticklablels = []
for pca_year in range(pca_st_yr, pca_ed_yr+1):
    pca_xticklablels.append(str(pca_year) + '06')
pca_time_ticks = (pca_xticks, pca_xticklablels)

res = '2x25'

pca_title = 'PCA'

m_vmin = None

#
# End user parameters
#####################

# agruements
start_year = 2014
end_year = 2019

start_month = 6
end_month = 8

n_mon = end_month - start_month + 1
xticks = np.array(range(0, (end_year-start_year+1)*n_mon, n_mon))
xticklablels = []
for year in range(start_year, end_year+1):
    xticklablels.append(str(year) + str(start_month).zfill(2))
    time_ticks = (xticks, xticklablels)

argv = sys.argv[1:]

opts = 'sy:ey:sm:em'
long_opts = ['start_year=', 'end_year=', 'start_month=', 'end_month=']

try:
    opts, args = getopt.getopt(argv, opts, long_opts)
except getopt.GetoptError:
    print('python main_calc_ave_emi.py -sy <start_year> ' + 
            '-ey <end_year> -sm <start_month> -em <end_month>')
    sys.exit(2)
for opt, arg in opts:
    if opt in ('-sy', '--start_year'):
        start_year = int(arg)
    if opt in ('-ey', '--end_year'):
        end_year = int(arg)
    if opt in ('-sm', '--start_month'):
        start_month = int(arg)
    if opt in ('-em', '--end_month'):
        end_month = int(arg)

# year and month
sy_c = str(start_year)
ey_c = str(end_year)
sm_c = str(start_month).zfill(2)
em_c = str(end_month).zfill(2)

# variable name
varn = 'TSOIL1'

# loop year and month
T_mon = []
coor_flag = True
for iyr in range(start_year, end_year+1):
    for imo in range(start_month, end_month+1):

        yyyy   = str(iyr)
        mm     = str(imo).zfill(2)
        yyyymm = yyyy + mm

        print('---- process ' + yyyymm  + ' ----')

        T_file = root_dir + yyyy + '/GEOSFP.' + yyyymm + '.' + res + '.nc'

        # read data
        indata = read_nc(T_file, varnames=[varn])
        T_mon.append(indata[varn] - 273.15)

# monthly soil temperature
T_mon = np.array(T_mon)

# PCA data
pca_varns = ['Latitude_e', 'Longitude_e', 'Latitude', 'Longitude',
        'components', 'coeffs', 'explained_variance_ratio']
pca_dict = read_nc(pca_file, pca_varns, verbose)
lat_e = pca_dict['Latitude_e']
lon_e = pca_dict['Longitude_e']
components = pca_dict['components']
coeffs = pca_dict['coeffs']
explained_variance_ratio = pca_dict['explained_variance_ratio']
comp1 = components[:,:,0]

# area
area_data = read_nc(area_file, varnames=['AREA'])
area = area_data['AREA']

# plot components and coeff
plot_2D_components_and_coeffs(components, lat_e, lon_e, coeffs,
        explained_variance_ratio,
        fig_dir, name=pca_title, n_components=1,
        time_ticks=pca_time_ticks)

# only keep data over regions with soil NOx emissions
# account for at least 30%
flag_nan = np.isnan(comp1)
for im in range(T_mon.shape[0]):
    T_mon[im,flag_nan] = np.nan
area[flag_nan] = np.nan

# loop regions
for region in region_list:

    # get region limit
    region_limit = sn_p.region_dict[region]

    # get region indexes
    i1, i2, j1, j2 = \
            get_center_index_latlon(lat_e[:,0], lon_e[0,:], region_limit)

    # region data
    r_area = area[i1:i2+1,j1:j2+1]
    r_T_mon = T_mon[:,i1:i2+1,j1:j2+1]
    r_lat_e = lat_e[i1:i2+2,j1:j2+2]
    r_lon_e = lon_e[i1:i2+2,j1:j2+2]

    plot_ave_series(r_T_mon, r_area, r_lat_e, r_lon_e,
            units=units, m_vmin=m_vmin, time_ticks=time_ticks)

    #
    sub_title = 'Soil_T_' + str(start_year) + '-' + str(end_year) \
            + '_' + region + '.png'
    figname = fig_dir + sub_title
    plt.savefig(figname, format='png', dpi=300)







