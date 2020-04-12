"""
Created on April 8, 2020

@author: Yi Wang
"""

from calendar import monthrange
import numpy as np
import sys, getopt

from mylib.grid_utility import generate_grid_gc_2x25
from mylib.io import write_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_io import read_NO_emissions


#######################
# Start user parameters
#

root_dir = '/Dedicated/jwang-data/ywang/soil_NOx/GEOS-Chem_ori/runs/'

out_dir = '../data/'

#
# End user parameters
#####################

# agruements
start_year = 2014
end_year = 2019

start_month = 6
end_month = 8

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

# emission category
emi_caty = ['EmisNO_Soil', 'EmisNO_Total']
emi_mon = {}
for emi_c in emi_caty:
    emi_mon[emi_c] = []


# loop year and month
days_mon  = []
for iyr in range(start_year, end_year+1):
    for imo in range(start_month, end_month+1):

        yyyy   = str(iyr)
        mm     = str(imo).zfill(2)
        yyyymm = yyyy + mm

        print('---- process ' + yyyymm  + ' ----')

        emi_file = root_dir + 'geosfp_2x25_tropchem_' + yyyymm \
                + '/HEMCO_diagnostics.' + yyyymm + '010000.nc'

        # read data
        NO_emi = read_NO_emissions(emi_file, varns=emi_caty, amount=False)
        for emi_c in emi_caty:
            emi_mon[emi_c].append(NO_emi[emi_c])

        # days in a month
        days = monthrange(iyr, imo)[1]
        days_mon.append(days)

# list to array
for emi_c in emi_caty:
    emi_mon[emi_c] = np.array(emi_mon[emi_c])
days_mon = np.array(days_mon)

# average
emi_ave = {}
for emi_c in emi_caty:
    tmp = []
    for i in range(len(days_mon)):
        tmp.append( emi_mon[emi_c][i,:,:] * days_mon[i] )
    tmp = np.array(tmp)
    emi_ave[emi_c] = np.sum(tmp, axis=0) / np.sum(days_mon)

# ratio
emi_ratio = {}
for emi_c in emi_caty:
    emi_ratio[emi_c + '_ratio'] = \
            emi_ave[emi_c] / emi_ave['EmisNO_Total']

# Latitude and Longitude
lat_e, lon_e, lat_c, lon_c = generate_grid_gc_2x25()
lon_e, lat_e = np.meshgrid(lon_e, lat_e)
lon_c, lat_c = np.meshgrid(lon_c, lat_c)

# output 
data_dict = {}
data_dict['Latitude'] = lat_c
data_dict['Longitude'] = lon_c
data_dict['Latitude_e'] = lat_e
data_dict['Longitude_e'] = lon_e
data_dict.update(emi_ave)
data_dict.update(emi_ratio)

units_dict = {}
for emi_c in emi_caty:
    units_dict[emi_c] = 'molec/cm2'
        
out_file = out_dir + 'NO_emi_ave_' + sy_c + '-' + ey_c \
        + '_' + sm_c + '-' + em_c + '.nc'
write_nc(out_file, data_dict, units_dict=units_dict)








