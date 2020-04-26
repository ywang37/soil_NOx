"""
Created on April 8, 2020

@author: Yi Wang
"""

#from calendar import monthrange
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import sys

from mylib.io import read_nc

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_plot import plot_NOx_emi_ratio


#######################
# Start user parameters
#

data_dir = '../data/'

fig_dir = '../figure/'

start_year = 2014
end_year = 2019

start_month = 6
end_month = 8

soil_thre = 1e10

verbose = True

#
# End user parameters
#####################

# year and month
sy_c = str(start_year)
ey_c = str(end_year)
sm_c = str(start_month).zfill(2)
em_c = str(end_month).zfill(2)

# read data
varnames = ['EmisNO_Soil', 'EmisNO_Soil_ratio', \
        'EmisNO_Total', 'Latitude_e', 'Longitude_e']
infile = data_dir + 'NO_emi_ave_' + sy_c + '-' + ey_c \
        + '_' + sm_c + '-' + em_c + '.nc'
data_dict = read_nc(infile, varnames, verbose=verbose)

# plot
flag = data_dict['EmisNO_Soil'] < soil_thre
data_dict['EmisNO_Soil_ratio'][flag] = np.nan
plot_NOx_emi_ratio(data_dict)
figname = fig_dir + 'NO_emi_ave_' + sy_c + '-' + ey_c \
        + '_' + sm_c + '-' + em_c + '.png'
plt.savefig(figname, format='png', dpi=300)






