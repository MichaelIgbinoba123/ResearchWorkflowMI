# imports (unchanged)
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
from wrf import (getvar, interplevel, to_np, latlon_coords, g_geoht, srhel, g_vorticity)
from metpy.units import units
from scipy.ndimage import gaussian_filter
from metpy.plots import colortables
from math import radians, sin, cos, sqrt, atan2

# read the file (unchanged)
path_file = 'D:\wrfout_d02_2022-09-23_09_00_10'
ncfile = Dataset(path_file)

# get the total number of time indices (unchanged)
num_time_indices = len(ncfile.dimensions['Time'])

# Loop through all time indices (unchanged)
for i in range(num_time_indices):
    # get the variables for the current time index (unchanged)
    p = getvar(ncfile, "pressure", timeidx=i)
    ua = getvar(ncfile, "ua", timeidx=i)
    va = getvar(ncfile, "va", timeidx=i)
    height = g_geoht.get_height(ncfile, timeidx=i, units='m')
    terrain = getvar(ncfile, "HGT", timeidx=i)
    wind = getvar(ncfile, "WSPD10MAX", timeidx=i)
    slp = getvar(ncfile, "AFWA_MSLP", timeidx=i)
    mslp = slp / 100
    ref = getvar(ncfile, "REFD_COM", timeidx=i)
    wind_kt = wind * 1.944
    avo = g_vorticity.get_avo(ncfile, timeidx=i)
    U = getvar(ncfile, "U10", timeidx=i)
    V = getvar(ncfile, "V10", timeidx=i)
    U10 = U * 1.944
    V10 = V * 1.944

    # specify your map boundaries (unchanged)
    lat_min = 14
    lat_max = 19
    lon_min = 121
    lon_max = 132

    # get the lat, lon grid (unchanged)
    lats, lons = latlon_coords(wind)

    # specify your colormap and projection (unchanged)
    cmap = colortables.get_colortable('NWSStormClearReflectivity')
    crs = ccrs.PlateCarree()

    # Apply smoothing to th_200 using a 2D Gaussian filter (unchanged)
    sigma = 0.0  # You can adjust this value for more or less smoothing
    ref_smoothed = gaussian_filter(ref, sigma=sigma)

    time_var = ncfile.variables['Times']
    time_str = ''.join([i.decode('utf-8') for i in time_var[i]])

    # plot (unchanged)
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, facecolor='gray', projection=crs)
    ax.coastlines(resolution='10m', alpha=0.75)
    # Use vmin and vmax to set the color range from 260 to 360 K (unchanged)
    plot_ref = ax.pcolormesh(lons, lats, ref_smoothed, cmap=cmap, vmin=-30, vmax=90)
    cbar = fig.colorbar(plot_ref)
    cbar.ax.set_ylabel('Reflectivity (dBz)')

    # some fancy schmancy grid lines (unchanged)
    gl = ax.gridlines(crs=crs, draw_labels=True, alpha=0.5)
    gl.top_labels = None
    gl.right_labels = None
    xgrid = np.arange(lon_min, lon_max + 1, 1.)
    ygrid = np.arange(lat_min, lat_max + 1, 1.)
    gl.xlocator = mticker.FixedLocator(xgrid.tolist())
    gl.ylocator = mticker.FixedLocator(ygrid.tolist())
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}


# Add the 500 hPa wind barbs, only plotting every 125th data point.
    plt.barbs(to_np(lons[::16, ::16]), to_np(lats[::16, ::16]),
              to_np(U10[::16, ::16]), to_np(V10[::16, ::16]),
              length=5, color='white')

    # Add the time information in the top right corner (unchanged)
    time_var = ncfile.variables['Times']
    time_str = ''.join([i.decode('utf-8') for i in time_var[i]])
    ax.text(1, 1.06, time_str, transform=ax.transAxes,
            ha='right', va='top', fontsize=6, color='black', backgroundcolor='white', weight='bold')

    # Display the minimum sea level pressure and maximum wind speed (new)
    min_mslp = np.nanmin(to_np(mslp))
    max_wspd = np.nanmax(to_np(wind_kt))
    ax.text(0.05, 0.05, f"Min MSLP: {min_mslp:.2f} hPa", transform=ax.transAxes,
            ha='left', va='bottom', fontsize=6, color='black', backgroundcolor='white', weight='bold')
    ax.text(0.95, 0.05, f"Max Wind Speed: {max_wspd:.2f} kt", transform=ax.transAxes,
            ha='right', va='bottom', fontsize=6, color='black', backgroundcolor='white', weight='bold')


    # set other plot parameters
    plt.xlim((lon_min, lon_max))
    plt.ylim((lat_min, lat_max))
    plt.title('Composite Reflectivity + Wind (kt) + MSLP (hPa)')

    # Save the plot as a .png image in the "Downloads" folder
    plt.savefig(f'C:/Users/33017/Downloads/plot_{i:03d}.png', bbox_inches='tight')

    # Close the plot window
    plt.close()

# Close the netCDF file
ncfile.close()
