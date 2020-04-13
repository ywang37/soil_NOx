from copy import deepcopy
import datetime
import glob
import numpy as np
from sklearn.decomposition import PCA

from mylib.io import read_nc
from mylib.land_ocean.land_ocean import is_land
from mylib.pca.io import write_2D_PCA_nc

#######################
# Start user parameters
#


inRootDir  = '../data/monthly/'
outRootDir = '../data/PCA/'

start_year = 2005
end_year   = 2019

scene = 'summer'
#scene = 'all_season'

land_ocean = 'land'

res = '2x25'

# *ratio* data should be available
ratio = 0.8


#
# End user parameters
#####################

if scene == 'summer':
    month_list = ['06', '07', '08']
elif scene == 'all_season':
    month_list = ['01', '02', '03', '04', '05', '06',
                  '07', '08', '09', '10', '11', '12']


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
ave = np.nanmean(in_var_all, axis=2)
ave[np.logical_not(flag)] = np.nan

# land_ocean
lon = deepcopy(out_dict['Longitude'])
lat = deepcopy(out_dict['Latitude'])
if land_ocean == 'land':
    land_ocean_flag = is_land(lon, lat)

# mask
ave[np.logical_not(land_ocean_flag)] = np.nan

# prepare data for PCA
ii_ind = []
jj_ind = []
PCA_data = []
for i in range(flag.shape[0]):
    for j in range(flag.shape[1]):

        if (flag[i,j] and land_ocean_flag[i,j]):

            # record index
            ii_ind.append(i)
            jj_ind.append(j)

            # store data
            tmp = in_var_all[i,j,:]
            tmp_nan = np.isnan(tmp)
            tmp[tmp_nan] = ave[i,j]
            PCA_data.append(tmp)
in_var_all[np.isnan(ave)] = np.nan
ii_ind = np.array(ii_ind) 
jj_ind = np.array(jj_ind) 
PCA_data = np.transpose(np.array(PCA_data))

# PCA
n_components=PCA_data.shape[0]-1
pca = PCA(n_components=n_components)
coeffs = pca.fit_transform(PCA_data)
#print(pca.explained_variance_ratio_)
#print(np.sum(pca.explained_variance_ratio_))
#print(pca.explained_variance_ / np.sum(pca.explained_variance_))
#print(coef.shape)

# map components back to grids
comp_ori = pca.components_
#print(comp_ori.shape)
components = np.full((*(flag.shape),n_components), np.nan)
#print(components.shape)
#print(components)
for ifs in range(comp_ori.shape[1]):
    
    # get index
    i = ii_ind[ifs]
    j = jj_ind[ifs]

    components[i,j,:] = comp_ori[:,ifs]

# prepare output
out_dict['ori_data'] = in_var_all
out_dict['ave'] = ave
out_dict['components'] = components
out_dict['coeffs'] = coeffs
out_dict['explained_variance_ratio'] = pca.explained_variance_ratio_
out_dict['explained_variance'] = pca.explained_variance_

# output
out_file = outRootDir + 'PCA_' + str(start_year) + '-' + str(end_year) \
        + '_' + scene + '_' + land_ocean + '_' + res + '.nc'
write_2D_PCA_nc(out_file, out_dict)






