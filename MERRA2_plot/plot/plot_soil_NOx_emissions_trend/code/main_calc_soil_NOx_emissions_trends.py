"""
Created on April 22, 2020

@author: Yi Wang
"""

from calendar import monthrange
import matplotlib.pyplot as plt
import numpy as np
import sys, getopt
from tqdm import tqdm

from mylib.constants import kg_NO_to_ng_N
from mylib.grid_utility import generate_grid_gc_2x25, get_center_index_latlon
from mylib.io import read_nc, write_nc
from mylib.pca.plot import plot_2D_components_and_coeffs
from mylib.trend_analysis import trend_analysis

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_io import read_NO_emissions
import sn_p
from sn_plot import plot_ave_series

#######################
# Start user parameters
#

in_data_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/process/\
emissions/data/'

out_data_dir = '../data/'

verbose = True

varn = 'EmisNO_Soil'

scene_tup = ['ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

#
# End user parameters
#####################

# agruements
start_year = 2005
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

month_days_dict = {}
month_days_dict['06'] = '30'
month_days_dict['07'] = '31'
month_days_dict['08'] = '31'

# loop scene
for scene in scene_tup:

    # loop year and month
    T_year = []
    coor_flag = True
    coor_var = ['Latitude', 'Longitude', 'Latitude_e', 'Longitude_e']
    for iyr in range(start_year, end_year+1):

        T_mon = []
        for imo in range(start_month, end_month+1):

            yyyy   = str(iyr)
            mm     = str(imo).zfill(2)
            yyyymm = yyyy + mm

            print('---- process ' + yyyymm  + ' ----')

            T_file = in_data_dir + scene + '/' + \
                    'HEMCO_diagnostics.' + yyyymm + '.nc'

            # read data
            if coor_flag:
                coor_flag = False
                indata = read_nc(T_file, varnames=[varn]+coor_var, 
                        squeeze=True)
                lat   = indata['Latitude']
                lon   = indata['Longitude']
                lat_e = indata['Latitude_e']
                lon_e = indata['Longitude_e']
            else:
                indata = read_nc(T_file, varnames=[varn], squeeze=True)
            T_mon.append(indata[varn] * kg_NO_to_ng_N)

        T_mon = np.array(T_mon)
        T_mon_mean = np.mean(T_mon, axis=0)

        T_year.append(T_mon_mean)

    # yearly soil temperature
    T_year = np.array(T_year)

    # mean
    T_mean = np.mean(T_year, axis=0)

    # loop grids to calculate trends
    trend = np.full_like(lat, np.nan)
    trend_sigma = np.full_like(lat, np.nan)
    for i in tqdm(range(lat.shape[0])):
        for j in range(lat.shape[1]):

            # data in one grid
            y = T_year[:,i,j]

            # trend analysis
            if not np.any(np.isnan(y)):
                ta = trend_analysis()
                ta.analysis_yearly(y)
                trend[i,j] = ta.popt[1]
                trend_sigma[i,j] = ta.trend_std

    # prepare coordinate data
    coor_dict = {}
    coor_dict['Latitude']    = lat
    coor_dict['Longitude']   = lon
    coor_dict['Latitude_e']  = lat_e
    coor_dict['Longitude_e'] = lon_e

    # prepare output data
    out_dict = {}
    out_dict['mean'] = T_mean
    out_dict['trend'] = trend
    out_dict['trend_sigma'] = trend_sigma

    # add coordinate
    out_dict.update(coor_dict)

    # units
    units_dict = {}
    units_dict['mean'] = 'ng N/m2/s'
    units_dict['trend'] = 'ng N/m2/s/a'
    units_dict['trend_sigma'] = 'ng N/m2/s/a'

    # output
    outfile = out_data_dir + scene + '_' + varn + '_' + \
            sy_c + '-' + ey_c + '_' + \
            sm_c + '-' + em_c + '.nc'
    write_nc(outfile, out_dict, units_dict=units_dict, verbose=True)
