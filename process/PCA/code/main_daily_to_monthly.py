import datetime
import glob
import numpy as np

from mylib.io import read_nc, write_nc

#######################
# Start user parameters
#

min_pixels = 8

min_days = 3

inRootDir  = '../data/daily/'
outRootDir = '../data/monthly/'

startMonth = '2005-01'
endMonth   = '2019-12'

res = '2x25'

#
# End user parameters
#####################

varname_list = ['NO2_Trop_CS', 'sat_grid_count']

coord_name_list = ['Latitude', 'Longitude', 'Latitude_e', 'Longitude_e']

currMonth  = startMonth
currMonth_D = datetime.datetime.strptime(currMonth+'-01', '%Y-%m-%d')
endMonth_D  = datetime.datetime.strptime(endMonth+'-01',  '%Y-%m-%d')

while currMonth_D <= endMonth_D:

    currMonth = str(currMonth_D)[0:7]
    print('------------------------------------------------------')
    print('processing ' + currMonth)

    currMonth_D = datetime.datetime(currMonth_D.year+(currMonth_D.month//12), \
            ((currMonth_D.month%12)+1), 1)

    # get all filenames of current month
    all_files = glob.glob(inRootDir + currMonth[0:4] + 
            '/OMI_NO2_'+ currMonth + '-*' + '_' + res + '.nc')
    all_files.sort()

    out_dict = {}
    latlon_flag = True
    in_var_list = []
    in_count_list = []
    if len(all_files) > 0:

        # read all data
        for filename in all_files:

            if latlon_flag:
                latlon_flag = False
                in_data = read_nc(filename, varname_list + coord_name_list,
                        verbose=True)
                for varname in coord_name_list:
                    out_dict[varname] = in_data.pop(varname)
            else:
                in_data = read_nc(filename, varname_list, verbose=True)
            in_var = in_data['NO2_Trop_CS']
            in_count = in_data['sat_grid_count']
            in_var[in_count<min_pixels] = np.nan
            in_var_list.append(in_var)
            in_count_list.append(in_count)

        # days
        in_count_sum = np.array(in_count_list, dtype=int)
        in_count_sum = (in_count_sum >= min_pixels)
        in_count_sum = np.sum(in_count_sum, axis=0)
        out_dict['days'] = in_count_sum

        # monthly average
        in_var_all = np.array(in_var_list)
        in_var_ave = np.nanmean(in_var_all, axis=0)
        in_var_ave[in_count_sum<min_days] = np.nan
        out_dict['NO2_Trop_CS'] = in_var_ave

        # output
        out_file = outRootDir + 'OMI_NO2_' + currMonth + '_' + res + '.nc'
        write_nc(out_file, out_dict)

    else:
        print('No data in current month')

