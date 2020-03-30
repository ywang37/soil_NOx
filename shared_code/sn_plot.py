"""
Created on January 2, 2020

@author: Yi Wang
"""

import cartopy.feature as cfeature
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

from sn_io import read_nc_emissions_multifiles

from mylib.cartopy_plot import add_geoaxes
from mylib.cartopy_plot import cartopy_plot
from mylib.colormap.colormap_utility import truncate_colormap
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map
from mylib.layout import multiFigure, h_2_ax, panel_tick_label
from mylib.layout import right_center_label

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



def plot_panel_variables(root_dir, scene_tup, month,
        varname='EmisNO_Soil', read_func=read_nc_emissions_multifiles,
        gc_run='geosfp_2x25_tropchem', res='2x25',
        read_func_varname='emi_dict',
        vmin=None, vmax=None, units='',
        cmap=None, cb_ticks=None, cb_ticklabels=None,
        seperate_cbar=False,
        vmin_diff=None, vmax_diff=None, units_diff='',
        cmap_diff=None, cb_ticks_diff=None, cb_ticklabels_diff=None,
        lw=None,
        title_list=None,
        diff=False, verbose=True, layout=layout_1):
    """
    """

    # read all data
    data_dict = read_func(root_dir, scene_tup, month,
            varname, gc_run=gc_run, res=res,
            verbose=verbose)

    # plot
    fig = plt.figure()
    ax_list = []
    xtick = [[], [], np.arange(-180, 180.0, 60), None]
    ytick = [None, [], None, []]
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
                xtick=xtick[i], ytick=ytick[i])
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

        emi = data_dict[read_func_varname][scene]
        if ( (i >= 1) and diff ):
            emi0 = data_dict[read_func_varname][scene_tup[0]]
            emi_diff = emi - emi0
            flag = np.logical_and(emi0<1e-6, emi<1e-6)
            emi_diff = np.ma.masked_array(emi_diff, flag)
            cp_out = cartopy_plot(lon_e, lat_e, emi_diff, ax=ax,
                    cbar=seperate_cbar, cmap=cmap_diff,
                    vmin=vmin_diff, vmax=vmax_diff,
                    )
            cp_out_list.append(cp_out)
        else:
            cp_out = cartopy_plot(lon_e, lat_e, emi, ax=ax, 
                    cbar=seperate_cbar, cmap=cmap,
                    vmin=vmin, vmax=vmax, valid_min=1e-6)
            cp_out_list.append(cp_out)

    # layout
    layout()

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
        cb1_extend = 'both'
        cb1_ticks = cb_ticks_diff
        cb1_ticklabels = cb_ticklabels_diff
        cb1_units = units_diff
    else:
        cb1_mesh = cp_out_list[0]['mesh']
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
                orientation='horizontal', extend='max')
        cax2.yaxis.set_label_position('right')
        cax2.set_ylabel(units, rotation=0, ha='left', va='center')
        if cb_ticks is not None:
            cb2.set_ticks(cb_ticks)
            if cb_ticklabels is not None:
                cb2.set_ticklabels(cb_ticklabels)

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
def plot_NO2_T(data_dict, sat_varname, mod_varname, 
        scene_tup, T_name,
        left=0.1, right=0.95, top=0.9, bottom=0.15,
        wspace=0.4, hspace=0.4,
        ):
    """
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










