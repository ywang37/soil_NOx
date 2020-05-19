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
#yyyymm_list = ['200506', '200507', '200508', '200509']
#yyyymm_list = ['200605', '200606', '200607', '200608', '200609']

yyyymm_list = []
month = ['05', '06', '07', '08', '09']
for year in range(2007, 2018):
    for mo in month:
        yyyymm_list.append(str(year) + mo)


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
