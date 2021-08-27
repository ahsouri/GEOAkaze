from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from palettable import colorbrewer
import numpy.ma as ma
import matplotlib as mpl
mpl.cm.get_cmap
import matplotlib.cm as cm

def make_kml(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
             figs, colorbar=None, **kw):
    """TODO: LatLon bbox, list of figs, optional colorbar figure,
    and several simplekml kw..."""
    kml = Kml()
    altitude = kw.pop('altitude', 2e7)
    roll = kw.pop('roll', 0)
    tilt = kw.pop('tilt', 0)
    altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)
    camera = Camera(latitude=np.mean([urcrnrlat, llcrnrlat]),
                    longitude=np.mean([urcrnrlon, llcrnrlon]),
                    altitude=altitude, roll=roll, tilt=tilt,
                    altitudemode=altitudemode)
    kml.document.camera = camera
    draworder = 0
    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kw.pop('visibility', 1)
        ground.name = kw.pop('name', 'overlay')
        ground.color = kw.pop('color', '9effffff')
        ground.atomauthor = kw.pop('author', 'ocefpaf')
        ground.latlonbox.rotation = kw.pop('rotation', 0)
        ground.description = kw.pop('description', 'Matplotlib figure')
        ground.gxaltitudemode = kw.pop('gxaltitudemode',
                                       'clampToSeaFloor')
        ground.icon.href = fig
        ground.latlonbox.east = urcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.west = llcrnrlon
    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='ScreenOverlay')
        screen.icon.href = colorbar
        screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1

    kmzfile = kw.pop('kmzfile')

    kml.savekmz(kmzfile)

def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
    """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
    aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
    xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = np.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize
    print(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat)
    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)
    if False:
        plt.ioff()  # Make `True` to prevent the KML components from poping-up.
    fig = plt.figure(figsize=figsize,
                     frameon=False,
                     dpi=pixels//10)
    # KML friendly image.  If using basemap try: `fix_aspect=False`.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    return fig, ax

def safe_convolve(ts,agg_fac):
    # Aggregate
    import scipy as sp
    import scipy.ndimage
    # Define weighting
    weights = np.ones(agg_fac)
    weights = weights / np.sum(weights[:])
    ct = np.ones(ts.shape)
    ct[ts.mask] = 0
    ts[ts.mask] = 0.0
    # Perform convolution
    y = sp.ndimage.filters.convolve(ts, weights, mode='constant')
    y_ct = sp.ndimage.filters.convolve(ct, weights, mode='constant')
    # Return mean
    return y/y_ct

def make_kmz(lon_coadd,lat_coadd,var_coadd,outfile_pre,vmin=None,vmax=None,legend=None,cmap=None):
    import os
    # Pixel dimensions
    pixels = 1024 * 10
    if(cmap is None):
        #cmap = colorbrewer.get_map('Spectral', 'diverging', 11, reverse=True).mpl_colormap
        cmap = cm.bone
    if(vmin is None):
        vmin = var_coadd.min()
    if(vmax is None):
        vmax = var_coadd.max()
    # Plot Figure
    overlay_png = outfile_pre+'_overlay.png'
    fig, ax = gearth_fig(llcrnrlon=lon_coadd.min(),
                         llcrnrlat=lat_coadd.min(),
                         urcrnrlon=lon_coadd.max(),
                         urcrnrlat=lat_coadd.max(),
                         pixels=pixels)

    cs = ax.pcolormesh(lon_coadd, lat_coadd, var_coadd, cmap=cm.bone,vmin=vmin,vmax=vmax)
    ax.set_axis_off()
    fig.savefig(overlay_png, transparent=False, format='png')
    plt.close()
    aspect = np.cos(np.mean([lat_coadd.min(), lat_coadd.max()]) * np.pi/180.0)
    xsize = np.ptp([lon_coadd.max(), lon_coadd.min()]) * aspect
    ysize = np.ptp([lat_coadd.max(), lat_coadd.min()])
    aspect = ysize / xsize
    # Plot Colorbar
    outfile = outfile_pre 
    legend_png = outfile_pre+'_legend.png'

    if(legend):
        # Make the KML File
        make_kml(llcrnrlon=lon_coadd.min(), llcrnrlat=lat_coadd.min(),
                urcrnrlon=lon_coadd.max(), urcrnrlat=lat_coadd.max(),
                figs=[overlay_png],kmzfile = outfile, colorbar='legend.png')
    else:
        # Make the KML File
        make_kml(llcrnrlon=lon_coadd.min(), llcrnrlat=lat_coadd.min(),
                urcrnrlon=lon_coadd.max(), urcrnrlat=lat_coadd.max(),
                figs=[overlay_png],kmzfile = outfile)
    # fig = plt.figure(figsize=(0.5, 2.0), facecolor=None, frameon=False)
    # ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
    # cb = fig.colorbar(cs, cax=ax)
    # fig.savefig(legend_png, transparent=False, format='png')  # Change transparent to True if your colorbar is not on space :)
    # plt.close()
    fig = plt.figure(figsize=(8, 3), facecolor=None, frameon=False)
    ax  = fig.add_axes([0.05, 0.80, 0.9, 0.15])
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    #cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,norm=norm,
    #                               orientation='horizontal')
    fig.savefig(legend_png, transparent=True, format='png')  # Change transparent to True if your colorbar is not on space :)
    plt.close()