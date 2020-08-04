"""
Created on July 15, 2020

@author: Yi Wang
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mylib.trend_analysis import trend_analysis

filename = '../data/FAOSTAT_US_N_fertilizer.csv'

df = pd.read_csv(filename)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(df['Year'], df['Value']/1e6, '-o')
ax.set_xlabel('Year')
ax.set_ylabel('US nitrogen fertilizer [Tg N a$^{-1}$]')
ax.set_ylim([0, 15.0])

# trend
percentage_scale = 100.0
yy = np.array(df['Value'])
yy = yy / yy[0] * percentage_scale
ta = trend_analysis()
ta.analysis_yearly(yy)
trend_label = '{:.2f}'.format(ta.popt[1]) + \
        u'\u00B1' + '{:.2f}'.format(ta.trend_std) + \
        r' % a$^{-1}$'
print(trend_label)
ax.text(0.05, 0.1, 'Trend relative to year {:}: '.format(df['Year'][0]) + \
        trend_label, transform=ax.transAxes)

plt.subplots_adjust(left=0.20, right=0.80, bottom=0.25, top=0.75)

plt.savefig('../figure/US_N_fertilizer.png', format='png', dpi=300)
