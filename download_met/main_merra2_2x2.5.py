"""
Created on March 19, 2020

@author: Yi Wang
"""

from mylib.model.gc_data import download_met

#######################
# Start user parameters
#

res = 'GEOS_2x2.5'

model = 'MERRA2'

#yyyymm_list = ['201809']
yyyymm_list = ['201905', '201906', '201907', '201908', '201909']

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
