"""
Created on January 2, 2020

@author: Yi Wang
"""

from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

from sn_io import read_nc_emissions_multifiles

from mylib.cartopy_plot import add_geoaxes
from mylib.cartopy_plot import cartopy_plot
from mylib.colormap.colormap_utility import truncate_colormap
from mylib.colormap.gbcwpry_map import gbcwpry_map
from mylib.colormap.WhGrYlRd_map import WhGrYlRd_map

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


def plot_panel_variables(root_dir, scene_tup, month,
        varname='EmisNO_Soil', read_func=read_nc_emissions_multifiles,
        gc_run='geosfp_2x25_tropchem', res='2x25',
        vmin=None, vmax=None, units='',
        cmap=None, cb_ticks=None, cb_ticklabels=None,
        seperate_cbar=False,
        vmin_diff=None, vmax_diff=None, units_diff='',
        cmap_diff=None, cb_ticks_diff=None, cb_ticklabels_diff=None,
        lw=None,
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
        ax = add_geoaxes(fig, int('22'+str(i+1)),
                title=scene + '_' + varname,
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

        emi = data_dict['emi_dict'][scene]
        if ( (i >= 1) and diff ):
            emi0 = data_dict['emi_dict'][scene_tup[0]]
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











