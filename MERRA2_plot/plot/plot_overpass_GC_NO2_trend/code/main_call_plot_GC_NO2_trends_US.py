import os

year_list = ['2005-2019', '2005-2011', '2011-2019']
month_list = ['06-08']

for year in year_list:
    for month in month_list:

        cmd = 'python main_plot_GC_NO2_trends_US.py ' + \
                '--time_range=' + year + '_' + month

        os.system(cmd)
