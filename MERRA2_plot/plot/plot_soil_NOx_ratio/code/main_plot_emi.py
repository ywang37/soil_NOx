"""
Created on April 8, 2020

@author: Yi Wang
"""

#from calendar import monthrange
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import sys

from mylib.constants import kg_NO_to_ng_N
from mylib.io import read_nc, write_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_plot import plot_NOx_emi_ratio


#######################
# Start user parameters
#

in_data_dir = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/\
process/emissions/data/na_05x0625/soil_T_obs/'

out_file = '../data/soil_ratio_2009_and_2012.nc'

figname = '../figure/soil_ratio_2009_and_2012.png'

year_list = ['2009', '2012']

month_list = ['06', '07', '08']

#soil_thre = 1e10

verbose = True

#
# End user parameters
#####################

varnames = ['EmisNO_Total', 'EmisNO_Soil',
        'Latitude', 'Longitude', 'Latitude_e', 'Longitude_e']

emi_total = []
emi_soil = []
for year in year_list:
    for month in month_list:

        filename = in_data_dir + 'HEMCO_diagnostics.' + year + month + '.nc'

        tmp_data = read_nc(filename, varnames, verbose=verbose,
                squeeze=True)

        emi_total.append(np.sum(tmp_data['EmisNO_Total'], axis=0))

        emi_soil.append(tmp_data['EmisNO_Soil'])

emi_total = np.mean(np.array(emi_total), axis=0)
emi_soil  = np.mean(np.array(emi_soil),  axis=0)

soil_ratio = emi_soil / emi_total 


plot_dict = {}
plot_dict['Latitude']   = tmp_data['Latitude']
plot_dict['Longitude']  = tmp_data['Longitude']
plot_dict['Latitude_e']   = tmp_data['Latitude_e']
plot_dict['Longitude_e']  = tmp_data['Longitude_e']
plot_dict['EmisNO_Soil']  = emi_soil * kg_NO_to_ng_N
plot_dict['EmisNO_Total'] = emi_total * kg_NO_to_ng_N
plot_dict['EmisNO_Soil_ratio'] = deepcopy(soil_ratio)

# out data
units_dict = {}
units_dict['EmisNO_Soil'] = 'ng N m-2 s-1'
units_dict['EmisNO_Total'] = 'ng N m-2 s-1'
write_nc(out_file, plot_dict, units_dict=units_dict)

# plot
plot_NOx_emi_ratio(plot_dict, NO2_VCD_unit=r'[ng N m$^{-2}$ s$^{-1}$]')
plt.savefig(figname, format='png', dpi=300)

















