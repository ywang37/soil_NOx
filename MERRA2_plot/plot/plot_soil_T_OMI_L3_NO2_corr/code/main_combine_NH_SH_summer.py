"""
Created on June 3, 2020

@author: Yi Wang
"""

import sys

sys.path.append('/Dedicated/jwang-data/ywang/soil_NOx/shared_code')
from sn_process import combine_NH_SH

#######################
# Start user parameters
#

name = 'corr_TSOIL1_NO2_Trop_CS_2005-2019'

data_dir = '../data/'

north_file = data_dir + name + '_06-08.nc'

south_file = data_dir + name + '_winter.nc'

combine_file = data_dir + name + '_combine_NH_SH_summer.nc'

#
# End user parameters
#####################


combine_NH_SH(north_file, south_file, combine_file)
