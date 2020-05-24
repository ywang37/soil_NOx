import datetime
import glob
import numpy as np
import os

from mylib.grid_utility import generate_grid_gc_2x25
from mylib.io import read_nc, write_nc


#######################
# Start user parameters
#

inRootDir  = '../data/daily/'
outRootDir = '../data/monthly/'

startYear = 2014
endYear   = 2019

month_list = ['06', '07', '08']

res = '2x25'

varnames = ['TSOIL1', 'T2M']

coor_names = ['Latitude', 'Longitude', 'Latitude_e', 'Longitude_e']

nan_thre = 200.0 # K

#
# End user parameters
#####################

# units_dict
units_dict = {}
for varn in varnames:
    units_dict[varn] = 'K'

coor_flag = True
coor_dict = {}
for year in range(startYear, endYear+1):

    year_c = str(year)

    # outDir
    outDir = outRootDir + year_c + '/' 
    if not os.path.isdir(outDir):
        os.system('mkdir -p ' + outDir)
    
    for month in month_list:

        print('--- ' + year_c + month + ' ---')

        indata_dict = {}
        for varn in varnames:
            indata_dict[varn] = []

        # inDir
        inDir = inRootDir + year_c + '/' + month + '/'

        # loop days
        startDate = year_c + month + '01'
        currDate_D = datetime.datetime.strptime(startDate, '%Y%m%d')
        while True:

            currDate = str(currDate_D)
            yyyy = currDate[0:4]
            mm   = currDate[5:7]
            dd   = currDate[8:10]
            yyyymmdd = yyyy + mm + dd
            if (mm != month):
                break

            # infile
            infile = inDir + 'GEOSFP.' + yyyymmdd + '.' + res + '.nc'

            # read data
            if coor_flag:
                indata = read_nc(infile, varnames+coor_names, verbose=True)
                coor_flag = False
                for coorn in coor_names:
                    coor_dict[coorn] = indata.pop(coorn)
            else:
                indata = read_nc(infile, varnames, verbose=True)
            for varn in varnames:
                indata_dict[varn].append(indata[varn])

            # go to next day
            currDate_D = currDate_D + datetime.timedelta(days=1)

        # calculate mean
        out_dict = {}
        for varn in varnames:
            out_dict[varn] = np.array(indata_dict[varn])
            out_dict[varn] = np.mean(out_dict[varn], axis=0)

        # add coordinate
        out_dict.update(coor_dict)

        # outfile
        outfile = outDir + 'GEOSFP.' + year_c + month + '.' + res + '.nc'
        write_nc(outfile, out_dict, units_dict=units_dict, verbose=True)




