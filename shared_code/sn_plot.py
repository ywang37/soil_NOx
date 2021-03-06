"""
Created on January 2, 2020

@author: Yi Wang
"""

import cartopy.feature as cfeature
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from sn_io import read_nc_emissions_multifiles

from mylib.cartopy_plot import add_geoaxes
from mylib.cartopy_plot import cartopy_plot
from mylib.colormap.colormap_utility import truncate_colormap
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.layout import multiFigure, h_1_ax, h_2_ax, panel_tick_label
from mylib.layout import right_center_label
from mylib.scatter_plot import scatter

def layout_1(left=0.08, right=0.97,
        top=0.95, bottom=0.15,
        hspace=-0.1, wspace=0.02):
    """
    """
    plt.subplots_adjust(left=left, right=right,
            top=top, bottom=bottom,
            hspace=hspace, wspace=wspace)

def layout_2(hspace=0.2):
    """
    """
    layout_1(hspace=hspace)

def layout_3(top=0.95, bottom=0.10, 
        left=0.06, right=0.92,
        wspace=0.12):
    """
    """
    plt.subplots_adjust(top=top, bottom=bottom, 
            left=left, right=right,
            wspace=wspace)



#
#------------------------------------------------------------------------------
#
def plot_panel_variables(root_dir, scene_tup, month,
        varname='EmisNO_Soil', read_func=read_nc_emissions_multifiles,
        path=1,
        gc_run='geosfp_2x25_tropchem', res='2x25',
        read_func_varname='emi_dict',
        scene_prefix='GEOS-Chem_',
        cl_res='110m',
        cl_color='k',
        emi_flag=False,
        flag_nan=False,
        flag_diff_nan=False,
        ocean_color=None,
        subdir='',
        scale=1.0,
        vmin=None, vmax=None, units='',
        valid_min=None,
        cmap=None, cb_ticks=None, cb_ticklabels=None,
        seperate_cbar=False,
        cb1_extend=None,
        cb2_extend='max',
        vmin_diff=None, vmax_diff=None, units_diff='',
        cmap_diff=None, cb_ticks_diff=None, cb_ticklabels_diff=None,
        lw=None,
        title_list=None,
        region_limit=None,
        xtick=np.arange(-180, 180.0, 60),
        ytick=None,
        diff=False, verbose=True, layout=layout_1,
        hspace=None):
    """
    """

    # read all data
    if isinstance(root_dir, str):
        data_dict = read_func(root_dir, scene_tup, month,
                varname, gc_run=gc_run, res=res,
                scene_prefix=scene_prefix,
                path=path,
                subdir=subdir,
                verbose=verbose)
    else:
        data_dict = root_dir

    # plot
    fig = plt.figure()
    ax_list = []
    xtick_list = [[], [], xtick, xtick]
    ytick_list = [ytick, [], ytick, []]
    for i in range(len(scene_tup)):
        scene = scene_tup[i]
        if title_list is None:
            title=scene + '_' + varname
        else:
            if title_list[i] is None:
                title=scene + '_' + varname
            else:
                title=title_list[i]
        ax = add_geoaxes(fig, int('22'+str(i+1)),
                title=title,
                lw=lw,
                xtick=xtick_list[i], ytick=ytick_list[i])
        ax_list.append(ax)

    lat_e = data_dict['lat_e']
    lon_e = data_dict['lon_e']
    cp_out_list = []
    if cmap is None:
        cmap = truncate_colormap(WhGrYlRd_map, 
                minval=0.05, maxval=1.0, n=200)
    if cmap_diff is None:
        cmap_diff = gbcwpry_map
        #cmap_diff = plt.get_cmap('seismic')
        
    for i in range(len(scene_tup)):

        scene = scene_tup[i]
        ax = ax_list[i]

        emi = data_dict[read_func_varname][scene] * scale
        if ( (i >= 1) and diff ):
            emi0 = data_dict[read_func_varname][scene_tup[0]] * scale
            emi_diff = emi - emi0
            if (valid_min is not None) and emi_flag:
                flag = np.logical_and(emi0<valid_min, emi<valid_min)
                emi_diff = np.ma.masked_array(emi_diff, flag)
            if flag_diff_nan:
                tmp_d = data_dict[read_func_varname+'_diff_nan'][scene_tup[i]]
                emi_diff[tmp_d] = np.nan
            cp_out = cartopy_plot(lon_e, lat_e, emi_diff, ax=ax,
                    cbar=seperate_cbar, cmap=cmap_diff,
                    vmin=vmin_diff, vmax=vmax_diff,
                    )
            cp_out_list.append(cp_out)
        else:
            if flag_nan:
                tmp = data_dict[read_func_varname+'_nan'][scene_tup[i]]
                emi[tmp] = np.nan
            cp_out = cartopy_plot(lon_e, lat_e, emi, ax=ax, 
                    cbar=seperate_cbar, cmap=cmap,
                    vmin=vmin, vmax=vmax, valid_min=valid_min)
            cp_out_list.append(cp_out)

    # layout
    layout()
    if hspace is not None:
        plt.subplots_adjust(hspace=hspace)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

    for ax in ax_list:
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
        if ocean_color is not None:
            ax.add_feature(cfeature.OCEAN, color=ocean_color, zorder=100)
            ax.coastlines(resolution=cl_res, color=cl_color, lw=lw, zorder=300)
        if region_limit is not None:
            ax.set_xlim([region_limit[1], region_limit[3]])
            ax.set_ylim([region_limit[0], region_limit[2]])

    # colorbar axes parameter
    h = 0.03

    # colorbar cax1
    pos21 = ax_list[2].get_position()
    pos22 = ax_list[3].get_position()
    cax1 = fig.add_axes([pos21.x0+pos21.width/2.0,
        pos21.y0-0.09, 
        pos22.x0+pos22.width/2.0-(pos21.x0+pos21.width/2.0),
        h])

    # colorbar cax2
    if diff:
        pos11 = ax_list[0].get_position()
        cax2 = fig.add_axes([pos11.x0, pos11.y0-0.04, pos11.width, h])

    # parameters for cax1
    if diff:
        cb1_mesh = cp_out_list[3]['mesh']
        if cb1_extend is None:
            cb1_extend = 'both'
        cb1_ticks = cb_ticks_diff
        cb1_ticklabels = cb_ticklabels_diff
        cb1_units = units_diff
    else:
        cb1_mesh = cp_out_list[0]['mesh']
        if cb1_extend is None:
            cb1_extend = 'max'
        cb1_ticks = cb_ticks
        cb1_ticklabels = cb_ticklabels
        cb1_units = units

    # colorbar 1
    cb1 = plt.colorbar(cb1_mesh, cax=cax1,
            orientation='horizontal', extend=cb1_extend)
    cax1.yaxis.set_label_position('right')
    cax1.set_ylabel(units, rotation=0, ha='left', va='center')
    cax1.yaxis.set_label_coords(1.01, 0.5)
    if cb1_ticks is not None:
        cb1.set_ticks(cb1_ticks)
        if cb1_ticklabels is not None:
            cb1.set_ticklabels(cb1_ticklabels)

    # colorbar 2
    if diff:
        cb2 = plt.colorbar(cp_out_list[0]['mesh'], cax=cax2,
                orientation='horizontal', extend=cb2_extend)
        cax2.yaxis.set_label_position('right')
        cax2.set_ylabel(units, rotation=0, ha='left', va='center')
        if cb_ticks is not None:
            cb2.set_ticks(cb_ticks)
            if cb_ticklabels is not None:
                cb2.set_ticklabels(cb_ticklabels)

#
#------------------------------------------------------------------------------
#
def plot_compare_4_to_1_new(data_dict,
        portrait=True,
        left=0.08, right=0.97,
        top=0.95, bottom=0.1,
        hspace=0.1, wspace=0.02,
        region_limit=None,
        xticks=np.arange(-180, 180.0, 60),
        yticks=np.arange(-90, 90.1, 30),
        cl_res='110m',
        cl_color='k',
        lw=None,
        ocean_color=None,
        title_dict={},
        vmin1=None, vmax1=None, units1='',
        cmap1=None, cb_ticks1=None, cb_ticklabels1=None,
        cb1_extend='max',
        vmin_diff=None, vmax_diff=None, units_diff='',
        cmap_diff=None, cb_ticks_diff=None, cb_ticklabels_diff=None,
        cb_diff_extend='both',
        cb_ratio = 0.04,
        cb1_y_off = -0.03,
        cb_diff_y_off = -0.06,
        ):
    """
    """

    # layout
    if portrait:
        nrow = 3
        ncol = 2
        figsize = (6.4, 7.2)
        varn_list = ['OMI', 'ori', 'diff_ori',
                'diff_soil_T_ori', 'diff_surf_T_obs', 'diff_soil_T_obs']
    else:
        nrow = 2
        ncol = 3
        figsize = (9.6, 4.8)
        varn_list = ['OMI', 'diff_ori', 'diff_soil_T_ori', 
                'ori', 'diff_surf_T_obs', 'diff_soil_T_obs']
    projPos = list(range(nrow * ncol))
    layout_dict = multiFigure(nrow, ncol,
            left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace,
            figsize=figsize, projPos=projPos)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']
    lat_e = data_dict['Latitude_e']
    lon_e = data_dict['Longitude_e']

    if vmin1 is None:
        vmin1 = np.nanmin(data_dict['OMI'])
    if vmax1 is None:
        vmax1 = np.nanmax(data_dict['OMI'])
    if vmin_diff is None:
        vmin_diff = -vmax1
    if vmax_diff is None:
        vmax_diff = vmax1

    pout_list = []
    for i in range(len(projPos)):

        ax = axes[i]
        varname = varn_list[i]
        data = data_dict[varname]
        title = title_dict.get(varname, '')

        if ('diff' in varname):
            vmin = vmin_diff
            vmax = vmax_diff
            cmap = cmap_diff
        else:
            vmin = vmin1
            vmax = vmax1
            cmap = cmap1

        pout = cartopy_plot(lon_e, lat_e, data, ax=ax, title=title,
                cbar=False, vmin=vmin, vmax=vmax, cmap=cmap)
        pout_list.append(pout)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

    # set limit
    for ax in axes:
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
        ax.coastlines(resolution=cl_res, color=cl_color, lw=lw, zorder=300)
        if ocean_color is not None:
            ax.add_feature(cfeature.OCEAN, color=ocean_color, zorder=100)
            ax.coastlines(resolution=cl_res, color=cl_color, lw=lw, zorder=300)
        if region_limit is not None:
            ax.set_xlim((region_limit[1],region_limit[3]))
            ax.set_ylim((region_limit[0],region_limit[2]))


    # ticks
    panel_tick_label(axes, ncol, xticks=xticks, yticks=yticks)


    # colorbar axes parameter

    if portrait:
        
        pos11 = pout_list[0]['ax']
        pos12 = pout_list[1]['ax']
        cax1 = h_2_ax(fig, pos11, pos12, ratio=cb_ratio, y_off=cb1_y_off)
        cb1 = plt.colorbar(pout_list[0]['mesh'], cax=cax1,
                orientation='horizontal', extend=cb1_extend)
        right_center_label(cax1, units1)

    else:

        pos11 = pout_list[3]['ax']
        cax1 = h_1_ax(fig, pos11, ratio=cb_ratio, y_off=cb1_y_off)
        cb1 = plt.colorbar(pout_list[3]['mesh'], cax=cax1,
                orientation='horizontal', extend=cb1_extend)
        right_center_label(cax1, units1)


    pos21 = pout_list[4]['ax']
    pos22 = pout_list[5]['ax']
    cax2 = h_2_ax(fig, pos21, pos22, ratio=cb_ratio, 
            y_off=cb_diff_y_off)
    cb2 = plt.colorbar(pout_list[4]['mesh'], cax=cax2,
            orientation='horizontal',
            extend=cb_diff_extend)
    cb2.set_ticks(cb_ticks_diff)
    right_center_label(cax2, units_diff)



#
#------------------------------------------------------------------------------
#
def plot_compare_4_to_1(var_dict, sat_varname, 
        mod_varname_list, sat_sim_name,
        mod_sim_name_dict,
        lw=None,
        sat_cv=None,
        soil_ratio=None,
        cmap_VCD=None, min_VCD=None, max_VCD=None,
        cmap_diff=None, min_diff=None, max_diff=None,
        label_VCD='', label_diff='', label_ratio='',
        min_ratio=None, max_ratio=None,
        ):
    """
    """

    ny_fig, nx_fig = 5, 3

    # plot
    fig = plt.figure(figsize=(10, 9))
    ax_list = []
    for i in range(ny_fig):
        for j in range(nx_fig):

            # title 
            if i == 0:
                if j == 0:
                    title = sat_sim_name
                elif j == 1:
                    title = sat_sim_name + ' stddev/ave'
                else:
                    title = 'soil/total emissions'
            else:
                mod_sim_name = mod_sim_name_dict[mod_varname_list[i-1]]
                if j == 0:
                    title = mod_sim_name
                elif j == 1:
                    title = mod_sim_name + ' diff'
                else:
                    title = mod_sim_name + ' ratio'

            # xtick
            if (i == ny_fig-1):
                xtick = None
            else:
                xtick = []

            # ytick
            if (j == 0):
                ytick = None
            else:
                ytick = []

            ax = add_geoaxes(fig, ny_fig, nx_fig, i*nx_fig+j+1,
                    lw=lw,
                    title=title, xtick=xtick, ytick=ytick)
            ax_list.append(ax)

    if cmap_VCD is None:
        cmap_VCD = truncate_colormap(WhGrYlRd_map,
                minval=0.05, maxval=1.0, n=200)

    if cmap_diff is None:
        cmap_diff = deepcopy(gbcwpry_map)

    # lat_e and lon_e
    lat_e = var_dict['Latitude_e']
    lon_e = var_dict['Longitude_e']

    cbar=False

    if sat_cv is not None:
        ax = ax_list[1]
        sc_pout = cartopy_plot(lon_e, lat_e, sat_cv, ax=ax,
                cmap=deepcopy(cmap_VCD), vmin=0, vmax=1, cbar=cbar)

    if soil_ratio is not None:
        ax = ax_list[2]
        srat_pout = cartopy_plot(lon_e, lat_e, soil_ratio, ax=ax,
                cmap=deepcopy(cmap_VCD), vmin=0, vmax=0.5, cbar=cbar)

    # satellite VCD
    sat_var = var_dict[sat_varname]
    ax = ax_list[0]
    cartopy_plot(lon_e, lat_e, sat_var, ax=ax, cmap=deepcopy(cmap_VCD),
            vmin=min_VCD, vmax=max_VCD, cbar=cbar)

    # model VCD, diff ratio
    for i in range(len(mod_varname_list)):

        mod_varname = mod_varname_list[i]

        # VCD
        mod_var = var_dict[mod_varname]
        ax = ax_list[(i+1)*nx_fig]
        pout_VCD = cartopy_plot(lon_e, lat_e, mod_var, ax=ax, 
                cmap=deepcopy(cmap_VCD),
                vmin=min_VCD, vmax=max_VCD, cbar=cbar)

        # VCD diff
        ax = ax_list[(i+1)*nx_fig+1]
        diff = mod_var - sat_var
        pout_diff = cartopy_plot(lon_e, lat_e, diff, ax=ax, 
                cmap=deepcopy(cmap_diff),
                vmin=min_diff, vmax=max_diff, cbar=cbar)

        # VCD ratio
        ax = ax_list[(i+1)*nx_fig+2]
        ratio = diff / mod_var
        pout_ratio = cartopy_plot(lon_e, lat_e, ratio, ax=ax, 
                cmap=deepcopy(cmap_diff),
                vmin=min_ratio, vmax=max_ratio, cbar=cbar)

    h = 0.02

    cv_pos = sc_pout['ax'].get_position()
    cv_cax = fig.add_axes([cv_pos.x0, cv_pos.y0+0.08, cv_pos.width, h])
    cv_cb = plt.colorbar(sc_pout['mesh'], cax=cv_cax, 
            orientation='horizontal', extend='max')

    srat_pout
    srat_pos = srat_pout['ax'].get_position()
    srat_cax = fig.add_axes([srat_pos.x0, srat_pos.y0+0.08, cv_pos.width, h])
    srat_cb = plt.colorbar(srat_pout['mesh'], cax=srat_cax,
            orientation='horizontal', extend='max')

    # colorbar
    pout_list = [pout_VCD, pout_diff, pout_ratio]
    cb_label = [label_VCD, label_diff, label_ratio]
    for j in range(len(pout_list)):
        pout = pout_list[j]
        pout['ax'].reset_position()
        pos = pout['ax'].get_position()
        cax = fig.add_axes([pos.x0, pos.y0-0.06, pos.width, h])
        cb = plt.colorbar(pout['mesh'], cax=cax, orientation='horizontal')
        cb.set_label(cb_label[j])

    layout_3()

def plot_two_month_diff(mon1_dict, mon2_dict, sat_varname, 
        mod_varname, scene_tup,
        vmin=None, vmax=None,
        diff_min=None, diff_max=None,
        left=0.05, right=0.98, top=0.95, bottom=0.1,
        wspace=0.05, hspace=0.05,
        y_off2=-0.03, y_off3=-0.06,
        xticks=np.arange(-180.0, 180.1, 20.0),
        yticks=np.arange(-90.0, 90.1, 10.0),
        region_limit=[-90.0, -180.0, 90.0, 180.0]):
    """ mon2 - mon1
    """


    nrow = 3
    ncol = 5
    figsize=(13, 8)
    projPos = list(range(nrow * ncol))
    layout_dict = multiFigure(nrow, ncol,
            left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace,
            figsize=figsize, projPos=projPos)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']

    lat_e = mon1_dict['Latitude_e']
    lon_e = mon1_dict['Longitude_e']

    pout1_list = []
    pout2_list = []
    pout3_list = []
    for j in range(ncol):
       
        # satellite 
        if j == 0:

            data1 = mon1_dict[sat_varname]
            data2 = mon2_dict[sat_varname]
            title = 'OMI'

        # model
        else:

            data1 = mon1_dict[mod_varname+scene_tup[j-1]]
            data2 = mon2_dict[mod_varname+scene_tup[j-1]]
            title = scene_tup[j-1]

        diff = data2 - data1

        # plot

        # month 1
        ax = axes[j]
        ax.set_title(title)
        pout = cartopy_plot(lon_e, lat_e, data1, ax=ax, vmin=vmin, vmax=vmax,
                cmap=deepcopy(WhGrYlRd_map), cbar=False)
        pout1_list.append(pout)

        # month 2
        ax = axes[ncol+j]
        pout = cartopy_plot(lon_e, lat_e, data2, ax=ax, vmin=vmin, vmax=vmax,
                cmap=deepcopy(WhGrYlRd_map), cbar=False)
        pout2_list.append(pout)

        # month2 - month1
        ax = axes[ncol*2+j]
        pout = cartopy_plot(lon_e, lat_e, diff, ax=ax, 
                cmap=plt.get_cmap('seismic'), cbar=False,
                vmin=diff_min, vmax=diff_max)
        pout3_list.append(pout)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

    # ticks
    panel_tick_label(axes, ncol, xticks=xticks, yticks=yticks)

    # set limit
    for ax in axes:
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, zorder=200)
        #ax.add_feature(cfeature.OCEAN, color='w', zorder=100)
        ax.set_xlim((region_limit[1],region_limit[3]))
        ax.set_ylim((region_limit[0],region_limit[2]))

    # VCD colorbar
    cax2 = h_2_ax(fig, pout2_list[1]['ax'], pout2_list[3]['ax'],
            y_off=y_off2)
    cb2 = plt.colorbar(pout2_list[1]['mesh'], cax=cax2, \
            orientation='horizontal')
    right_center_label(cax2, r'[molec cm$^{-2}$]')

    # diff colorbar
    cax3 = h_2_ax(fig, pout3_list[1]['ax'], pout3_list[3]['ax'],
            y_off=y_off3)
    cb3 = plt.colorbar(pout3_list[1]['mesh'], cax=cax3, \
            orientation='horizontal')
    right_center_label(cax3, r'[molec cm$^{-2}$]')
#
#------------------------------------------------------------------------------
#
def plot_NO2_T_scatter(data_dict, sat_varname, mod_varname, 
        scene_tup, T_name, T_thre = 20.0,
        left=0.1, right=0.95, top=0.9, bottom=0.15,
        wspace=0.4, hspace=0.4,
        ):
    """
    (ywang, 03/31/20)
    """

    nrow = 2
    ncol = 3
    figsize = (9, 6)
    projPos = [3]
    layout_dict = multiFigure(nrow, ncol,
            left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace,
            figsize=figsize, projPos=projPos)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']

    lat_e = data_dict['Latitude_e']
    lon_e = data_dict['Longitude_e']

#    print(data_dict[T_name]-293.15)
#    print(data_dict[T_name].shape)
#    print(lon_e.shape)
#    exit()

    # dict index for scatter plot
    mod_dict = {1:0, 2:1, 4:2, 5:3}

    T = data_dict[T_name] - 273.15
    flag = T >= T_thre
    #flag = np.logical_and(T<30)
    #flag = T < 30
    x = deepcopy(data_dict[sat_varname])
    x = x[flag]
    for i in [0,1,2,4,5]:

        # satellite
        if i == 0:

            y = data_dict[sat_varname]

        # model
        else:

            y = data_dict[mod_varname+scene_tup[mod_dict[i]]]

        y = y[flag]
        ax = axes[i]
        scatter(ax, x, y)

    # plot satellite NO2
    sat_data = deepcopy(data_dict[sat_varname])
    sat_data_ave = np.nanmean(sat_data, axis=0)
    ax = axes[3]
    pout = cartopy_plot(lon_e, lat_e, sat_data_ave, ax=ax, 
            vmin=None, vmax=None,
            cmap=deepcopy(WhGrYlRd_map), cbar=False)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

    # set limit
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, zorder=200)
    ax.set_xlim([lon_e[0,0],lon_e[0,-1]])
    ax.set_ylim([lat_e[0,0],lat_e[-1,0]])

    # colorbar
    cax = h_1_ax(fig, pout['ax'])
    cb = plt.colorbar(pout['mesh'], cax=cax, \
            orientation='horizontal')

#
#------------------------------------------------------------------------------
#
def plot_NO2_CV_emi_ratio(data_dict, sat_varname='sat_ColumnAmountNO2Trop',
        emi_dict=None,
        left=0.1, right=0.95, top=0.9, bottom=0.2,
        wspace=0.2, hspace=0.5,
        NO2_VCD_unit=r'[molec cm$^{-2}$]',
        y_off1=-0.03,
        xticks=np.arange(-180.0, 180.1, 10.0),
        yticks=np.arange(-90.0, 90.1, 5.0),
        ):
    """
    (ywang 03/31/20)
    """

    nrow = 2
    ncol = 2
    figsize = (8, 6)
    projPos = [0, 1, 2, 3]
    layout_dict = multiFigure(nrow, ncol,
            left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace,
            figsize=figsize, projPos=projPos)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']

    lat_e = data_dict['Latitude_e']
    lon_e = data_dict['Longitude_e']

    # ave
    ave_pout = cartopy_plot(lon_e, lat_e, data_dict[sat_varname + '_ave'], 
            ax=axes[0], vmin=0.0, vmax=None,
            cmap=deepcopy(WhGrYlRd_map), cbar=False)
    axes[0].set_title('ave')

    # std
    std_pout = cartopy_plot(lon_e, lat_e, data_dict[sat_varname + '_std'],
            ax=axes[1], vmin=0.0, vmax=None,
            cmap=deepcopy(WhGrYlRd_map), cbar=False)
    axes[1].set_title('std')

    # cv
    cv_pout = cartopy_plot(lon_e, lat_e, data_dict[sat_varname + '_cv'],
            ax=axes[2], vmin=0.0, vmax=1.0,
            cmap=deepcopy(WhGrYlRd_map), cbar=False)
    axes[2].set_title('cv')

    # emi ratio
    if emi_dict is not None:
        emi_r_pout = cartopy_plot(emi_dict['lon_e'], emi_dict['lat_e'],
                emi_dict['emi_ratio'], 
                ax=axes[3], vmin=0, vmax=1.0,
                cmap=deepcopy(WhGrYlRd_map), cbar=False)
        axes[3].set_title(emi_dict.get('title', ''))

    # ticks
    panel_tick_label(axes, ncol, xticks=xticks, yticks=yticks)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

    # set limit
    for ax in axes:
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, zorder=200)
        ax.set_xlim([lon_e[0,0],lon_e[0,-1]])
        ax.set_ylim([lat_e[0,0],lat_e[-1,0]])


    # colorbar
    ave_cax = h_1_ax(fig, ave_pout['ax'], y_off=y_off1)
    ave_cb = plt.colorbar(ave_pout['mesh'], cax=ave_cax, \
            orientation='horizontal')
    ave_cb.set_label(NO2_VCD_unit)
    std_cax = h_1_ax(fig, std_pout['ax'], y_off=y_off1)
    std_cb = plt.colorbar(std_pout['mesh'], cax=std_cax, \
            orientation='horizontal')
    std_cb.set_label(NO2_VCD_unit)
    cv_cax = h_1_ax(fig, cv_pout['ax'])
    cv_cb = plt.colorbar(cv_pout['mesh'], cax=cv_cax, \
            orientation='horizontal')
    if emi_dict is not None:
        emi_r_cax = h_1_ax(fig, emi_r_pout['ax'])
        emi_r_cb = plt.colorbar(emi_r_pout['mesh'], cax=emi_r_cax, \
                orientation='horizontal')

#
#------------------------------------------------------------------------------
#
def plot_NO2_T_series(data_dict, sat_varname, mod_varname, 
        scene_tup, T_name, T_thre = 20.0,
        left=0.1, right=0.95, top=0.9, bottom=0.15,
        wspace=0.4, hspace=0.4,
        ):
    """
    (ywang, 03/31/20)
    """

    nrow = 2
    ncol = 1
    figsize = (8, 8)
    projPos = [0]
    layout_dict = multiFigure(nrow, ncol,
            left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace,
            figsize=figsize, projPos=projPos)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']

    lat_e = data_dict['Latitude_e']
    lon_e = data_dict['Longitude_e']


    T = data_dict[T_name] - 273.15
    flag = T >= T_thre
    T[np.logical_not(flag)] = np.nan
    T_series = np.nanmean(T, axis=(1,2))


    data_list = []

    # satellite data
    data_list.append(deepcopy(data_dict[sat_varname]))

    # model data
    for i in range(len(scene_tup)):
        data_list.append(deepcopy(data_dict[mod_varname+scene_tup[i]]))

    # 
    for i in range(len(data_list)):
        data_list[i][np.logical_not(flag)] = np.nan


    # plot satellite NO2
    sat_data_ave = np.nanmean(data_list[0], axis=0)
    ax = axes[0]
    pout = cartopy_plot(lon_e, lat_e, sat_data_ave, ax=ax, 
            vmin=None, vmax=None,
            cmap=deepcopy(WhGrYlRd_map), cbar=False)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

    # set limit
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, zorder=200)
    ax.set_xlim([lon_e[0,0],lon_e[0,-1]])
    ax.set_ylim([lat_e[0,0],lat_e[-1,0]])

    # colorbar
    cax = h_1_ax(fig, pout['ax'])
    cb = plt.colorbar(pout['mesh'], cax=cax, \
            orientation='horizontal')

    # time series
    ax = axes[1]
    for i in range(len(data_list)):
        ax.plot(np.nanmean(data_list[i], axis=(1,2)))
    ax_t = ax.twinx()
    ax_t.plot(T_series, 'k')

#
#------------------------------------------------------------------------------
#
def plot_NO2_VS_T(data_dict, sat_varname, mod_varname, 
        scene_tup, T_name,
        emi_ratio=None, emi_ratio_thre=None,
        xticks=np.arange(-180.0, 180.1, 10),
        yticks=np.arange(-90.0, 90.1, 5),
        T_edge = [20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0],
        left=0.1, right=0.95, top=0.9, bottom=0.1,
        wspace=0.3, hspace=0.3, map_region_limit=None,
        color_list=None,
        old_labels=None,
        ax1_ylim=None,
        ):
    """
    (ywang, 03/31/20)
    """
    out_dict = {}

    nrow = 2
    ncol = 3
    figsize = (12, 8)
    projPos = [0]
    layout_dict = multiFigure(nrow, ncol,
            left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace,
            figsize=figsize, projPos=projPos)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']

    lat_e = data_dict['Latitude_e']
    lon_e = data_dict['Longitude_e']

    T_edge = np.linspace(20, 40, 21)


    # soil temperature
    T = data_dict[T_name] - 273.15

    data_list = []

    # satellite data
    data_list.append(deepcopy(data_dict[sat_varname]))

    # model data
    for i in range(len(scene_tup)):
        data_list.append(deepcopy(data_dict[mod_varname+scene_tup[i]]))

    # filter
    if (emi_ratio is not None) and (emi_ratio_thre is not None):
        ratio_flag = emi_ratio > emi_ratio_thre
        for i in range(len(data_list)):
            data_list[i][:,np.logical_not(ratio_flag)] = np.nan

    # plot satellite NO2
    sat_data_ave = np.nanmean(data_list[0], axis=0)
    ax = axes[0]
    ax.set_title(r'OMI NO$_2$ VCD')
    pout = cartopy_plot(lon_e, lat_e, sat_data_ave, ax=ax, 
            vmin=0, vmax=None,
            cmap=deepcopy(WhGrYlRd_map), cbar=False)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

    # latitude and longitude
    panel_tick_label([axes[0]], 1, xticks=xticks, yticks=yticks)

    # set limit
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    ax.add_feature(cfeature.COASTLINE, zorder=200)
    if map_region_limit is None:
        ax.set_xlim([lon_e[0,0],lon_e[0,-1]])
        ax.set_ylim([lat_e[0,0],lat_e[-1,0]])
    else:
        ax.set_xlim([map_region_limit[1], map_region_limit[3]])
        ax.set_ylim([map_region_limit[0], map_region_limit[2]])

    # colorbar
    cax = h_1_ax(fig, pout['ax'])
    cb = plt.colorbar(pout['mesh'], cax=cax, \
            orientation='horizontal')
    #right_center_label(cax, r'[molec cm$^{-2}$]')
    cb.set_label(r'[molec cm$^{-2}$]')


    # T_dict
    T_dict = {}
    T_dict['data'] = []
    T_dict['ave'] = []

    # NO2_list_dict
    NO2_list_dict = []
    for i in range(len(data_list)):
        NO2_list_dict.append( {} )

    flag1 = np.logical_and(T>T_edge[0], T<=T_edge[-1])
    for i in range(len(NO2_list_dict)):

        NO2_list_dict[i]['data'] = []
        NO2_list_dict[i]['ave']  = []
        NO2_list_dict[i]['all_ave'] = np.nanmean(data_list[i])

        # data in the temperature range
        NO2_list_dict[i]['t_range_data'] = data_list[i][flag1]

    # R
    corr_t_range = []
    flag2 = np.logical_not( np.isnan(NO2_list_dict[0]['t_range_data'] ) )
    for i in range(len(NO2_list_dict)):
        tmp2 = pearsonr(NO2_list_dict[i]['t_range_data'][flag2], 
                        NO2_list_dict[0]['t_range_data'][flag2])
        corr_t_range.append( tmp2 )



    for i in range(len(T_edge)-1):

        flag = np.logical_and(T>T_edge[i], T<=T_edge[i+1])

        if np.sum(flag) > 20:

            # soil temperature
            T_dict['data'].append(T[flag])
            T_dict['ave'].append(np.nanmean(T[flag]))

            # NO2
            for j in range(len(data_list)):
                NO2_list_dict[j]['data'].append(data_list[j][flag])
                NO2_list_dict[j]['ave'].append(np.nanmean(data_list[j][flag]))


    
    for i in range(len(data_list)):

        NO2_list_dict[i]['ave'] = np.array(NO2_list_dict[i]['ave'])

        # diff normal
        diff = NO2_list_dict[i]['all_ave'] - NO2_list_dict[0]['all_ave']
        NO2_list_dict[i]['ave_diff_nor'] = \
                NO2_list_dict[i]['ave'] - diff

        # ratio normal
        ratio = NO2_list_dict[i]['all_ave'] / NO2_list_dict[0]['all_ave']
        NO2_list_dict[i]['ave_ratio_nor'] = \
                NO2_list_dict[i]['ave'] / ratio

        # bias
        NO2_list_dict[i]['bias'] = NO2_list_dict[i]['ave'] - \
                NO2_list_dict[0]['ave']

        # nmb
        NO2_list_dict[i]['nmb'] = NO2_list_dict[i]['bias'] / \
                NO2_list_dict[0]['ave'] * 100.0
        

    if old_labels is None:
        old_labels = ['OMI', 'Control', r'T$_{soil}$_S$_{old}$',
                r'T$_{air}$_S$_{new}$', r'T$_{soil}$_S$_{new}$']
    labels = deepcopy(old_labels)

    xlabel = u'Soil temperature [\u00B0C]'
    ylabel = r'NO$_2$ VCD [molec cm$^{-2}$]'

    flag = np.logical_not( np.isnan( NO2_list_dict[0]['ave'] ) )
#    print(flag)
#    print(type(flag))
#    print(NO2_list_dict[0]['ave'])
#    print(type(NO2_list_dict[0]['ave']))
#    exit()

    if color_list is None:
        color_list = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']

    # NO2 VS T
    ax = axes[1]
    corr_no_nor = []
    for i in range(len(NO2_list_dict)):
        zorder = 10 - i
        corr_tmp = pearsonr(NO2_list_dict[i]['ave'][flag],
                NO2_list_dict[0]['ave'][flag])
        if i == 0:
            label = labels[i]
        else:
            label = labels[i] + '({:.2f})'.format(corr_tmp[0])
        ax.plot(T_dict['ave'], NO2_list_dict[i]['ave'],
                c=color_list[i],
                marker='o', label=label, zorder=zorder)
        corr_no_nor.append(corr_tmp )
    if ax1_ylim is not None:
        ax.set_ylim(ax1_ylim)
    ax.legend(loc='best')
    ax.set_title('No normalization')
    print('--- No normalization ---')
    for i in range(len(labels)):
        print(labels[i], corr_no_nor[i])

    # create pandas DataFrame
    df = {}
    df['soil_temp'] = T_dict['ave']
    df['OMI'] = NO2_list_dict[0]['ave']
    for i in range(len(scene_tup)):
        df[scene_tup[i]] = NO2_list_dict[i+1]['ave']
    df = pd.DataFrame(df)
    out_dict['df'] = df

    # NO2 VS T (diff normal)
    ax = axes[2]
    corr_diff_nor = []
    for i in range(len(NO2_list_dict)):
        zorder = 10 - i
        ax.plot(T_dict['ave'], NO2_list_dict[i]['ave_diff_nor'],
                c=color_list[i],
                marker='o', label=labels[i], zorder=zorder)
        corr_diff_nor.append( pearsonr(NO2_list_dict[i]['ave_diff_nor'][flag],
                                      NO2_list_dict[0]['ave_diff_nor'][flag]) )      
    ax.set_title('Difference normalization')
    print('--- Diff normalization ---')
    for i in range(len(labels)):
        print(labels[i], corr_diff_nor[i])

    # NO2 VS T (ratio normal)
    ax = axes[3]
    corr_ratio_nor = []
    for i in range(len(NO2_list_dict)):
        zorder = 10 - i
        ax.plot(T_dict['ave'], NO2_list_dict[i]['ave_ratio_nor'],
                c=color_list[i],
                marker='o', label=labels[i], zorder=zorder)
        corr_ratio_nor.append(pearsonr(NO2_list_dict[i]['ave_ratio_nor'][flag],
                                      NO2_list_dict[0]['ave_ratio_nor'][flag]) )
    ax.set_title('Ratio normalization')
    print('--- Ratio normalization ---')

    print('--- corr_t_range ---')
    for i in range(len(labels)):
        print(labels[i], corr_t_range[i])

    # bias VS T
    ax = axes[4]
    for i in range(1,len(NO2_list_dict)):
        zorder = 10 - i
        ax.plot(T_dict['ave'], NO2_list_dict[i]['bias'],
                marker='o', c=color_list[i], zorder=zorder,
                label=old_labels[i])
    ax.plot([T_dict['ave'][0], T_dict['ave'][-1]], [0, 0], 'k--')
    ax.set_title('Bias')
    ax.legend(loc='best')

    # nmb
    ax = axes[5]
    for i in range(1,len(NO2_list_dict)):
        zorder = 10 - i
        ax.plot(T_dict['ave'], NO2_list_dict[i]['nmb'],
                marker='o', c=color_list[i], zorder=zorder,
                label=old_labels[i])
    ax.plot([T_dict['ave'][0], T_dict['ave'][-1]], [0, 0], 'k--')
    ax.set_title('Relative bias')
    ax.legend(loc='best')

    for i in range(1, len(axes)):
        axes[i].set_xlabel(xlabel)
        axes[i].set_ylabel(ylabel)

    axes[4].set_ylabel(r'GC NO$_2$ VCD bias [molec cm$^{-2}$]')
    axes[5].set_ylabel(r'GC NO$_2$ VCD relative bias [%]')

    return out_dict
#
#------------------------------------------------------------------------------
#
def plot_NOx_emi_ratio(data_dict, 
        left=0.1, right=0.95, top=0.95, bottom=0.15,
        wspace=0.2, hspace=0.6,
        NO2_VCD_unit=r'[molec cm$^{-2}$]',
        ratio=0.03,
        y_off1=-0.02,
        y_off2=-0.04,
        xticks=np.arange(-180.0, 180.1, 60.0),
        yticks=np.arange(-90.0, 90.1, 30.0),
        ):
    """
    (ywang 03/31/20)
    """

    nrow = 3
    ncol = 1
    figsize = (6, 9)
    projPos = [0, 1, 2]
    layout_dict = multiFigure(nrow, ncol,
            left=left, right=right, top=top, bottom=bottom,
            wspace=wspace, hspace=hspace,
            figsize=figsize, projPos=projPos)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']

    lat_e = data_dict['Latitude_e']
    lon_e = data_dict['Longitude_e']

    # totol
    total_pout = cartopy_plot(lon_e, lat_e, data_dict['EmisNO_Total'], 
            ax=axes[0], vmin=0.0, vmax=None,
            cmap=deepcopy(WhGrYlRd_map), cbar=False)
    axes[0].set_title('Total')

    # soil
    soil_pout = cartopy_plot(lon_e, lat_e, data_dict['EmisNO_Soil'],
            ax=axes[1], vmin=0.0, vmax=None,
            cmap=deepcopy(WhGrYlRd_map), cbar=False)
    axes[1].set_title('Soil')

    # ratio
    ratio_pout = cartopy_plot(lon_e, lat_e, data_dict['EmisNO_Soil_ratio'],
            ax=axes[2], vmin=0.0, vmax=1.0,
            cmap=deepcopy(WhGrYlRd_map), cbar=False)
    axes[2].set_title('Soil / Total')

    # ticks
    panel_tick_label(axes, ncol, xticks=xticks, yticks=yticks)

    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')

    # set limit
    for ax in axes:
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, zorder=200)
        ax.set_xlim([lon_e[0,0],lon_e[0,-1]])
        ax.set_ylim([lat_e[0,0],lat_e[-1,0]])


    # colorbar
    total_cax = h_1_ax(fig, total_pout['ax'], y_off=y_off1, ratio=ratio)
    total_cb = plt.colorbar(total_pout['mesh'], cax=total_cax, \
            orientation='horizontal')
    total_cb.set_label(NO2_VCD_unit)

    soil_cax = h_1_ax(fig, soil_pout['ax'], y_off=y_off1, ratio=ratio)
    soil_cb = plt.colorbar(soil_pout['mesh'], cax=soil_cax, \
            orientation='horizontal')
    soil_cb.set_label(NO2_VCD_unit)

    ratio_cax = h_1_ax(fig, ratio_pout['ax'], y_off=y_off2, ratio=ratio)
    ratio_cb = plt.colorbar(ratio_pout['mesh'], cax=ratio_cax, \
            orientation='horizontal')
#
#------------------------------------------------------------------------------
#
def plot_ave_series(data, area, lat_e, lon_e, units='',
        m_vmin=None,
        xticks=np.arange(-180.0, 180.1, 10.0),
        yticks=np.arange(-90.0, 90.1, 5.0),
        time_ticks=None,
        rotation=30.0,
        ):
    """
    """

    # ave
    ave = np.nanmean(data, axis=0)

    # series
    series = np.zeros((data.shape[0],))
    for i in range(data.shape[0]):
        series[i] = np.nansum(data[i,:,:] * area) / np.nansum(area)

    # begin plot
    nrow = 2
    ncol = 1
    figsize = (8, 6)
    projPos = [0]
    layout_dict = multiFigure(nrow, ncol,
            figsize=figsize, projPos=projPos)
    fig  = layout_dict['fig']
    axes = layout_dict['axes']

    # plot ave
    pout = cartopy_plot(lon_e, lat_e, ave,
            ax=axes[0], vmin=m_vmin, cbar=True)
    cb = pout['cb']
    panel_tick_label([axes[0]], ncol, xticks=xticks, yticks=yticks)
    cb.set_label(units)
    states_provinces = cfeature.NaturalEarthFeature(
            category='cultural',
            name='admin_1_states_provinces_lines',
            scale='50m',
            facecolor='none')
    axes[0].add_feature(cfeature.BORDERS)
    axes[0].add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    axes[0].add_feature(cfeature.COASTLINE, zorder=200)
    axes[0].set_xlim((lon_e[0,0],lon_e[0,-1]))
    axes[0].set_ylim((lat_e[0,0],lat_e[-1,0]))

    # plot series
    x = range(len(series))
    axes[1].plot(x, series, '-o', markersize=3)
    axes[1].set_xlabel('Month')
    axes[1].set_ylabel(units)

    if time_ticks is not None:
        axes[1].set_xticks(time_ticks[0])
        axes[1].set_xticklabels(time_ticks[1], rotation=rotation)
#
#------------------------------------------------------------------------------
#
def plot_time_series(data_dict, varnames, xx, label_dict={},
        ratio=False, xlabel='', ylabel='', xlim=None, ylim=None,
        scale=1.0):
    """
    (ywang, 07/08/2020)
    """

    out_dict = {}

    fig = plt.figure()
    ax = fig.add_subplot(111)

    out_dict['fig'] = fig
    out_dict['ax']  = ax

    for varn in varnames:

        yy = deepcopy(data_dict[varn])
        if ratio:
            yy = yy / yy[0]

        label = label_dict.get(varn, '')

        ax.plot(xx, yy * scale, marker="o", label=label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)


    plt.subplots_adjust(left=0.20, right=0.80, bottom=0.25, top=0.75)

    plt.legend(loc='best')

    return out_dict
#
#------------------------------------------------------------------------------
#
def plot_stackplot(data_dict, varnames, xx, label_dict={},
         xlabel='', ylabel='', xlim=None, ylim=None,
         scale=1.0, xminorticks=None):
    """
    (ywang, 07/10/2020)
    """

    out_dict = {}

    fig = plt.figure()
    ax = fig.add_subplot(111)

    out_dict['fig'] = fig
    out_dict['ax']  = ax

    yy_list = []
    labels = []
    for varn in varnames:

        y = deepcopy(data_dict[varn]) * scale
        yy_list.append(y)

        labels.append(label_dict[varn])

    yy = np.vstack(yy_list)

    ax.stackplot(xx, yy, labels=labels)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    if xminorticks is not None:
        ax.set_xticks(xminorticks, minor=True)

    plt.subplots_adjust(left=0.20, right=0.80, bottom=0.25, top=0.75)

    plt.legend(loc='best')

    return out_dict
#
#------------------------------------------------------------------------------
#
