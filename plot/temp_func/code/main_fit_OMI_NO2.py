"""
Created on January 1, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyval
import pandas as pd
from scipy.optimize import curve_fit

infile = '/Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/plot/\
plot_overpass_NO2_VCD_VS_soil_T/data/sat_model_NO2.csv'


scene = ['OMI', 'ori', 'soil_T_ori', 'surf_T_obs', 'soil_T_obs']

label_dict = {}
label_dict['OMI'] = 'OMI'
label_dict['ori'] = 'Control'
label_dict['soil_T_ori'] = r'T$_{soil}$_S$_{old}$'
label_dict['surf_T_obs'] = r'T$_{air}$_S$_{new}$'
label_dict['soil_T_obs'] = r'T$_{soil}$_S$_{new}$'

df = pd.read_csv(infile)


# soil temperature
soil_temp = df['soil_temp']


# scale data for plot
for i in range(len(scene)):
    df[scene[i]] = df[scene[i]] / 1e15

# fit OMI data
def fit_model(T, a, b, c, d):
    return a + b * T + c * T ** 2 + d * T ** 3

popt, pcov = curve_fit(fit_model, soil_temp, df[scene[0]])
omi_fit = fit_model(soil_temp, popt[0], popt[1], popt[2], popt[3])
print(popt)
print(soil_temp)
print(omi_fit)

# plot
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(soil_temp, df[scene[0]], 'o', c='C0', label=label_dict[scene[0]])
ax.plot(soil_temp, omi_fit, c='C0', 
        label=label_dict[scene[0]]+' fit')

for i in range(1 ,len(scene)):

    label = label_dict[scene[i]]

    if scene[i] == 'OMI':
        marker_line = 'o'
    else:
        marker_line = '-o'

    ax.plot(soil_temp, df[scene[i]], marker_line, label=label,
            c='C'+str(i))


ax.set_xlabel(u'Soil temperature [\u00B0C]')
ax.set_ylabel(r'NO$_2$ VCD [10$^{15}$ molec cm$^{-2}$]')
ax.set_ylim([0,3])
plt.legend()
plt.subplots_adjust(top=0.75, bottom=0.25, left=0.23, right=0.77)
plt.savefig('../figure/temp_func_fit_OMI_NO2.png', format='png', dpi=300)
