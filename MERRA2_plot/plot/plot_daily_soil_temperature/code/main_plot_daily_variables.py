"""
Created on July 27, 2020

@author: Yi Wang
"""

from copy import deepcopy
import datetime
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

from mylib.grid_utility import generate_grid_gc_2x25
from mylib.io import read_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_process import get_grid_data
from sn_utility import generate_latlon_label

#######################
# Start user parameters
#

startDate = '2019-06-01'
endDate   = '2019-08-31'

root_data_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/soil_T/data/daily/'

fig_dir = '../figure/'

varname_list = ['GWETTOP', 'T2M', 'TSOIL1']

T_list = ['T2M', 'TSOIL1']

nlat = 91
nlon = 144

res = '2x25'

grid_list = [ (34, -100), (36, -100), 
        (38, -100), (40, -100), (42, -100),
        (44, -100), (46, -100) ]

mo_dict = {}
mo_dict['06'] = 'JUN'
mo_dict['07'] = 'JUL'
mo_dict['08'] = 'AUG'

ylabel_dict = {}
ylabel_dict['GWETTOP'] = 'Soil water-filled pore space [%]'
ylabel_dict['T2M'] = u'2-m temperature [\u00B0C]'
ylabel_dict['TSOIL1'] = u'Soil temperature [\u00B0C]'

#
# End user parameters
#####################

label_list = []
for grid in grid_list:
    lat_label, lon_label = generate_latlon_label(grid[0], grid[1])
    label_list.append(lat_label + ' ' +  lon_label)

nan_arr = np.full((nlat, nlon), np.nan)

# Date
currDate   = startDate
currDate_D = datetime.datetime.strptime(currDate, '%Y-%m-%d')
endDate_D  = datetime.datetime.strptime(endDate,  '%Y-%m-%d')

indata_dict = {}
for varn in varname_list:
    indata_dict[varn] = []

xx = []
xticks = []
xticklabels = []
day = 0
while currDate_D <= endDate_D:


    # current date
    currDate = str(currDate_D)[0:10]
    yyyy = currDate[0:4]
    mm   = currDate[5:7]
    dd   = currDate[8:10]
    yyyymmdd = yyyy + mm + dd

    if ((dd == '01') or (dd == '16')):
        xticks.append( day )
        xticklabels.append( mo_dict[mm] + ' '  + str(int(dd)))

    xx.append(day)
    day += 1

    in_filename = root_data_dir + yyyy + '/' + mm + '/' + \
            'MERRA2.' + yyyymmdd  + '.2x25.nc'

    if os.path.exists(in_filename):

        indata = read_nc(in_filename, varname_list, verbose=True)

        for varn in varname_list:

            if varn in T_list:
                indata_dict[varn].append( indata.pop(varn) - 273.15)
            elif varn == 'GWETTOP':
                indata_dict[varn].append( indata.pop(varn) * 100.0)
            else:
                indata_dict[varn].append( indata.pop(varn) )
        
    else:

        print(in_filename + ' does not exists.')

        for varn in varname_list:
            indata_dict[varn].append( deepcopy(nan_arr) )

    # go to next day
    currDate_D = currDate_D + datetime.timedelta(days=1)

for varn in varname_list:
    indata_dict[varn] = np.array(indata_dict[varn])

# get data
if res == '2x25':
    lat_edge, lon_edge, lat_c, lon_c = generate_grid_gc_2x25()
selection_dict = \
        get_grid_data(indata_dict, varname_list, grid_list, lon_edge, lat_edge)

#for i in range(len(selection_dict['mod_Met_TSOIL1'])):
#    print(selection_dict['mod_Met_TSOIL1'][i])

# plot
for varn in varname_list:

    data = selection_dict[varn]

    fig = plt.figure()

    ax = fig.add_subplot(111)

    for i in range(len(data)):
        ax.plot(xx, data[i], 'o-', markersize=3, label=label_list[i])

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    ax.set_xlabel('Day')

    ylabel = ylabel_dict[varn]
    ax.set_ylabel(ylabel)

    plt.legend(loc='best')

    #plt.subplots_adjust()

    figname = fig_dir + 'daily_' + varn + '_' + startDate + '_' + \
            endDate + '.png'
    plt.savefig(figname, format='png', dpi=300)




