from copy import deepcopy
import datetime
import glob
import numpy as np

from mylib.io import read_nc, write_nc

#######################
# Start user parameters
#


inRootDir  = '../data/monthly/'

start_year = 2005
end_year   = 2019

month_list = ['06', '07', '08']

res = '2x25'

# *ratio* data should be available
ratio = 0.8

#
# End user parameters
#####################

# threshold (at least *ratio* data are available)
threshold = int((end_year - start_year + 1) * len(month_list) * ratio)

varname_list = ['NO2_Trop_CS', 'days']

coord_name_list = ['Latitude', 'Longitude', 'Latitude_e', 'Longitude_e']

# generate month
yyyymm_list = []
for year in range(start_year, end_year+1):
    for month in month_list:
        yyyymm_list.append(str(year) + '-' + month)

latlon_flag = True
in_var_list = []
out_dict = {}
for currMonth in yyyymm_list:

    # read data
    infile = inRootDir + 'OMI_NO2_' + currMonth + '_' + res + '.nc'
    if latlon_flag:
        latlon_flag = False
        in_data = read_nc(infile, varname_list + coord_name_list,
                verbose=True)
        for varname in coord_name_list:
            out_dict[varname] = in_data.pop(varname)
    else:
        in_data = read_nc(infile, varname_list, verbose=True)
    in_var = in_data['NO2_Trop_CS']
    in_var_list.append(in_var)

# transpose
in_var_all = np.array(in_var_list)
in_var_all = np.transpose(in_var_all, (1,2,0))

# detect available grid
# (at least *threshold* months are available)
avail_month = np.logical_not( np.isnan(in_var_all) )
count_month = np.sum(avail_month, axis=2)
flag = (count_month >= threshold)

# calcaulte average for selected grid
ave = np.mean(in_var_all, axis=2)
ave[np.logical_not(flag)] = np.nan

# prepare data for PCA
ii_ind = []
jj_ind = []
PCA_data = []
for i in range(flag.shape[0]):
    for j in range(flag.shape[1]):

        if (flag[i,j]):

            # record index
            ii_ind.append(i)
            jj_ind.append(j)

            # store data
            tmp = deepcopy(in_var_all[i,j,:])
            tmp_nan = np.isnan(tmp)
            tmp[tmp_nan] = ave[i,j]
            PCA_data.append(tmp)
ii_ind = np.array(ii_ind) 
jj_ind = np.array(jj_ind) 
PCA_data = np.array(PCA_data)
print(ii_ind.shape)
print(jj_ind.shape)
print(PCA_data.shape)
print(np.sum(flag))




