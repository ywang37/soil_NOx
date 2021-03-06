"""
Created on March 19, 2020

@author: Yi Wang
"""

from mylib.model.gc_data import link_for_soil_temp

#######################
# Start user parameters
#

res = 'GEOS_2x2.5'

model = 'GEOS_FP'

yyyymm_list = ['201405','201406', \
        '201407', '201408', '201409']
#yyyymm_list = ['201809']

#
# End user parameters
#####################

print('res is ' + res)
print('model is ' + model)


for yyyymm in yyyymm_list:

    year  = yyyymm[0:4]
    month = yyyymm[4:6]

    print('------ link ' + yyyymm + ' ------')
    link_for_soil_temp(res, model, year, month)
