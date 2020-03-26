"""
Created on March 19, 2020

@author: Yi Wang
"""

from mylib.model.gc_data import download_met

#######################
# Start user parameters
#

res = 'GEOS_2x2.5'

model = 'GEOS_FP'

#yyyymm_list = ['201809']
yyyymm_list = ['201605', '201606', '201607', '201608', '201609']

#
# End user parameters
#####################

print('res is ' + res)
print('model is ' + model)


for yyyymm in yyyymm_list:

    year  = yyyymm[0:4]
    month = yyyymm[4:6]

    print('Download ' + yyyymm)
    download_met(res, model, year, month)
