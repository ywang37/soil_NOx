import datetime
import glob
import numpy as np
import os

from mylib.grid_utility import generate_grid_gc_2x25
from mylib.io import read_nc, write_nc


#######################
# Start user parameters
#

inRootDir  = '/Dedicated/jwang-data/GCDATA/GEOS_2x2.5/GEOS_FP_soil_T/'
outRootDir = '../data/daily/'

startYear = 2014
endYear   = 2019

month_list = ['06', '07', '08']

res = '2x25'

varnames = ['TSOIL1', 'T2M']

nan_thre = 200.0 # K

#
# End user parameters
#####################

if (res == '2x25'):
    lat_e, lon_e, lat_c, lon_c = generate_grid_gc_2x25()


lon_c, lat_c = np.meshgrid(lon_c, lat_c)
lon_e, lat_e = np.meshgrid(lon_e, lat_e)

# units_dict
units_dict = {}
for varn in varnames:
    units_dict[varn] = 'K'

for year in range(startYear, endYear+1):

    year_c = str(year)
    
    for month in month_list:

        print('--- ' + year_c + month + ' ---')

        # inDir
        inDir = inRootDir + year_c + '/' + month + '/'

        # outDir
        outDir = outRootDir + year_c + '/' + month + '/'
        if not os.path.isdir(outDir):
            os.system('mkdir -p ' + outDir)

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
            infile = inDir + '/GEOSFP.' + yyyymmdd \
                    + '.A1.' + res + '.nc'

            # read data
            indata = read_nc(infile, varnames, verbose=True)

            # nan
            for varn in varnames:
                flag = (indata[varn] <= nan_thre)
                indata[varn][flag] = np.nan

            # average
            outdata = {}
            for varn in varnames:
                outdata[varn] = np.mean(indata[varn], axis=0)


            # outfile
            outfile = outDir + '/GEOSFP.' + yyyymmdd \
                    + '.' + res + '.nc'

            # output
            outdata['Latitude']    = lat_c
            outdata['Longitude']   = lon_c
            outdata['Latitude_e']  = lat_e
            outdata['Longitude_e'] = lon_e
            write_nc(outfile, outdata, units_dict=units_dict, verbose=True)

            # go to next day
            currDate_D = currDate_D + datetime.timedelta(days=1)




