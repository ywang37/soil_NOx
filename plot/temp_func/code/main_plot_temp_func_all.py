"""
Created on January 1, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyval

# Hudman's parameters
e_coef = 0.103

# Olikawa's parameters
a = -0.009
b = 0.837
c = -22.527
d = 196.149

# temperature for Hudman's
temp1 = np.arange(0, 50.01, 0.1)
H_temp = np.array(temp1)
H_temp[H_temp>=30.0] = 30.0
H_temp_func = np.exp(e_coef * H_temp)

# temperature for Olikawa's
temp2 = np.arange(20.0, 50.01, 0.1)
O_temp = np.array(temp2)
O_temp[O_temp>=40.0] = 40.0
O_temp_func = polyval(O_temp, [d,c,b,a])

# plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(temp1, H_temp_func, label=r'GEOS-Chem (e$^{'+str(e_coef)+'T}$)')
ax.plot(temp2, O_temp_func, 
        label=
        'Fit_to_20\n('+str(a)+'T$^3$+'+str(b)+'T$^2$'+str(c)+'T+'+str(d)+')')

# Olikawa's parameters
a = -0.0061
b =  0.5708
c = -15.3634
d = 133.77

# temperature for Olikawa's
temp2 = np.arange(18.0, 50.01, 0.1)
O_temp = np.array(temp2)
O_temp[O_temp>=40.0] = 40.0
O_temp_func = polyval(O_temp, [d,c,b,a])
ax.plot(temp2, O_temp_func,
        label=
        'Fit_to_18\n('+str(a)+'T$^3$+'+str(b)+'T$^2$'+str(c)+'T+'+str(d)+')')

a = 0.001367
b = -0.11
c = 3.066
d = -20.43

# temperature for Olikawa's
temp2 = np.arange(20.0, 50.01, 0.1)
O_temp = np.array(temp2)
O_temp[O_temp>=40.0] = 40.0
O_temp_func = polyval(O_temp, [d,c,b,a])

ax.plot(temp2, O_temp_func,
        label=
        'Fit OMI NO$_2$\n('+str(a)+'T$^3$'+str(b)+'T$^2$+'+str(c)+'T'+str(d)+')')


ax.set_xlim([0, 50])
ax.set_ylim([0, 100])
ax.set_xlabel(u'Soil temperature, T [\u00B0C]')
ax.set_ylabel('Temperature function, f(T)')
plt.legend(loc='upper left')
plt.subplots_adjust(top=0.95, bottom=0.15, left=0.23, right=0.77)
plt.savefig('../figure/temp_func_all.png', format='png', dpi=300)
