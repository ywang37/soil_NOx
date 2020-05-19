"""
Created on March 19, 2020

@author: Yi Wang
"""

from mylib.model.gc_data import link_for_soil_temp

#######################
# Start user parameters
#

res = 'GEOS_2x2.5'

model = 'MERRA2'

#yyyymm_list = ['201905','201906', \
#        '201907', '201908', '201909']
#yyyymm_list = ['201501']

yyyymm_list = []
month = ['05', '06', '07', '08', '09']
for year in range(2006, 2018):
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

    print('------ link ' + yyyymm + ' ------')
    link_for_soil_temp(res, model, year, month)
