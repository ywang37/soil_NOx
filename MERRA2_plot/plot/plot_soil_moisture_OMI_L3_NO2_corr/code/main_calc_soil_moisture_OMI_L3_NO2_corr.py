"""
Created on April 22, 2020

@author: Yi Wang
"""

from calendar import monthrange
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import sys, getopt
from tqdm import tqdm

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

root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/soil_T/data/monthly/'

NO2_root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/process/\
PCA_OMI_L3_NO2/data/monthly/'

data_dir = '../data/'

verbose = True

res = '2x25'

# variable name
varn = 'GWETTOP'

# NO2 variable name
NO2_varn = 'NO2_Trop_CS'

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


# loop year and month
T_year = []
T_monthly = []
NO2_year = []
NO2_monthly = []
T_month_dict = {}
NO2_month_dict = {}
for imo in range(start_month, end_month+1):
    mm = str(imo).zfill(2)
    T_month_dict[mm] = []
    NO2_month_dict[mm] = []
coor_flag = True
coor_var = ['Latitude', 'Longitude', 'Latitude_e', 'Longitude_e']
for iyr in range(start_year, end_year+1):

    T_mon = []
    NO2_mon = []
    for imo in range(start_month, end_month+1):

        yyyy   = str(iyr)
        mm     = str(imo).zfill(2)
        yyyymm = yyyy + mm

        print('---- process ' + yyyymm  + ' ----')

        # read temperature data
        T_file = root_dir + yyyy + '/MERRA2.' + yyyymm + '.' + res + '.nc'
        if coor_flag:
            coor_flag = False
            indata = read_nc(T_file, varnames=[varn]+coor_var)
            lat   = indata['Latitude']
            lon   = indata['Longitude']
            lat_e = indata['Latitude_e']
            lon_e = indata['Longitude_e']
        else:
            indata = read_nc(T_file, varnames=[varn])
        T_tmp = indata[varn] - 273.15
        T_mon.append(T_tmp)
        T_monthly.append(T_tmp)
        T_month_dict[mm].append(T_tmp)

        # read NO2 data
        NO2_file = NO2_root_dir + 'OMI_NO2_' + yyyy + '-' + mm + '_' + res \
                + '.nc'
        NO2_indata = read_nc(NO2_file, varnames=[NO2_varn])
        # molec/cm2 => 1e15 molec/cm2
        NO2_tmp = NO2_indata[NO2_varn] / 1e15
        NO2_mon.append(NO2_tmp) 
        NO2_monthly.append(NO2_tmp)
        NO2_month_dict[mm].append(NO2_tmp)


    # tempearture data
    T_mon = np.array(T_mon)
    T_mon_mean = np.mean(T_mon, axis=0)
    T_year.append(T_mon_mean)

    # NO2 data
    NO2_mon = np.array(NO2_mon)
    NO2_mon_mean = np.mean(NO2_mon, axis=0)
    NO2_year.append(NO2_mon_mean)

# yearly
T_year = np.array(T_year)
NO2_year = np.array(NO2_year)

# mean
T_mean = np.mean(T_year, axis=0)
NO2_mean = np.mean(NO2_year, axis=0)

# monthly
T_monthly = np.array(T_monthly)
NO2_monthly = np.array(NO2_monthly)

# every month
for imo in range(start_month, end_month+1):
    mm = str(imo).zfill(2)
    T_month_dict[mm] = np.array(T_month_dict[mm])
    NO2_month_dict[mm] = np.array(NO2_month_dict[mm])


# loop grids to calculate linear correlation
month_mean_r = np.full_like(lat, np.nan)
month_mean_p = np.full_like(lat, np.nan)
monthly_r = np.full_like(lat, np.nan)
monthly_p = np.full_like(lat, np.nan)
month_r_dict = {}
month_p_dict = {}
for imo in range(start_month, end_month+1):
    mm = str(imo).zfill(2)
    month_r_dict[mm] = np.full_like(lat, np.nan)
    month_p_dict[mm] = np.full_like(lat, np.nan)
for i in tqdm(range(lat.shape[0])):
    for j in range(lat.shape[1]):


        # multimonth mean
        x_yr = NO2_year[:,i,j]
        y_yr = T_year[:,i,j]
        if (not np.any(np.isnan(x_yr))) and (not np.any(np.isnan(y_yr))):
            corr_yr = pearsonr(x_yr, y_yr)
            month_mean_r[i,j] = corr_yr[0]
            month_mean_p[i,j] = corr_yr[1]

        # monhtly
        x_mo = NO2_monthly[:,i,j]
        y_mo = T_monthly[:,i,j]
        if (not np.any(np.isnan(x_mo))) and (not np.any(np.isnan(y_mo))):
            corr_mo = pearsonr(x_mo, y_mo)
            monthly_r[i,j] = corr_mo[0]
            monthly_p[i,j] = corr_mo[1]

        # every month
        for imo in range(start_month, end_month+1):
            mm = str(imo).zfill(2)
            x_ev_mo = NO2_month_dict[mm][:,i,j]
            y_ev_mo = T_month_dict[mm][:,i,j]
            if (not np.any(np.isnan(x_ev_mo))) and \
                    (not np.any(np.isnan(y_ev_mo))):
                corr_ev_mo = pearsonr(x_ev_mo, y_ev_mo)
                month_r_dict[mm][i,j] = corr_ev_mo[0]
                month_p_dict[mm][i,j] = corr_ev_mo[1]


# prepare coordinate data
coor_dict = {}
coor_dict['Latitude']    = lat
coor_dict['Longitude']   = lon
coor_dict['Latitude_e']  = lat_e
coor_dict['Longitude_e'] = lon_e

# prepare output data
out_dict = {}
out_dict['mean_T'] = T_mean
out_dict['mean_NO2'] = NO2_mean
out_dict['r_multimonth_mean'] = month_mean_r
out_dict['p_multimonth_mean'] = month_mean_p
out_dict['r_monthly'] = monthly_r
out_dict['p_monthly'] = monthly_p
for imo in range(start_month, end_month+1):
    mm = str(imo).zfill(2)
    out_dict['r_month_'+mm] = month_r_dict[mm]
    out_dict['p_month_'+mm] = month_p_dict[mm]

# add coordinate
out_dict.update(coor_dict)

# units
units_dict = {}
units_dict['mean_T'] = u'\u00B0C'
units_dict['mean_NO2'] = '1e15 molec/cm2'
units_dict['r_multimonth_mean'] = 'unitless'
units_dict['p_multimonth_mean'] = 'unitless'
units_dict['r_monthly'] = 'unitless'
units_dict['p_monthly'] = 'unitless'
for imo in range(start_month, end_month+1):
    mm = str(imo).zfill(2)
    units_dict['r_month_'+mm] = 'unitless'
    units_dict['p_month_'+mm] = 'unitless'

# output
outfile = data_dir + 'corr_' + varn + '_' + NO2_varn + '_' + \
        sy_c + '-' + ey_c + '_' + sm_c + '-' + em_c + '.nc'
write_nc(outfile, out_dict, units_dict=units_dict, verbose=True)
