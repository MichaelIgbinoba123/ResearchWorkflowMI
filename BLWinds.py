# imports
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
from wrf import (getvar, interplevel, to_np, latlon_coords, g_temp)
from metpy.units import units
from scipy.ndimage import gaussian_filter

# read the file
path_file = 'D:\wrfout_d02_2022-09-23_18_00_00'
ncfile = Dataset(path_file)

# get the total number of time indices
num_time_indices = len(ncfile.dimensions['Time'])

# Loop through all time indices
for i in range(num_time_indices):
    # get the variables for the current time index
    h = getvar(ncfile, "z", timeidx=i)
    PBLH = getvar(ncfile, "PBLH", timeidx=i)
    ua = getvar(ncfile, "ua", timeidx=i)
    va = getvar(ncfile, "va", timeidx=i)

    # select your pressure level
    h_level = PBLH

    # interpolate ua and va to pressure level
    u_bl = interplevel(ua, h, h_level)
    v_bl = interplevel(va, h, h_level)

    #Calculate scalar
    wind_bl = (np.sqrt(u_bl**2 + v_bl**2)) * 1.944

    # specify your map boundaries
    lat_min = 14
    lat_max = 19
    lon_min = 121
    lon_max = 133

    # get the lat, lon grid
    lats, lons = latlon_coords(ua)

    # specify your colormap and projection
    cmap = plt.get_cmap('cubehelix')
    crs = ccrs.PlateCarree()

    # Apply smoothing to th_200 using a 2D Gaussian filter
    sigma = 1.0  # You can adjust this value for more or less smoothing
    wind_bl_smoothed = gaussian_filter(wind_bl, sigma=sigma)

    time_var = ncfile.variables['Times']
    time_str = ''.join([i.decode('utf-8') for i in time_var[i]])

    # plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, facecolor='gray', projection=crs)
    ax.coastlines(resolution='10m', alpha=0.5)
    # Use vmin and vmax to set the color range from 260 to 360 K
    plot_wind_bl_smoothed = ax.pcolormesh(lons, lats, wind_bl_smoothed, cmap=cmap,     norm=plt.Normalize(0, 200))
    cbar = fig.colorbar(plot_wind_bl_smoothed)
    cbar.ax.set_ylabel('Wind (kt)')

    # some fancy schmancy grid lines
    gl = ax.gridlines(crs=crs, draw_labels=True, alpha=0.9)
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

    # Add the time information in the top right corner
    time_var = ncfile.variables['Times']
    time_str = ''.join([i.decode('utf-8') for i in time_var[i]])
    ax.text(1, 1.07, time_str, transform=ax.transAxes,
            ha='right', va='top', fontsize=8, color='black', backgroundcolor='white', weight='bold')

    #Find max wind
    max_wspd = np.nanmax(to_np(wind_bl_smoothed))
    ax.text(0.95, 0.05, f"Max Wind Speed: {max_wspd:.2f} kt", transform=ax.transAxes,
            ha='right', va='bottom', fontsize=6, color='black', backgroundcolor='white', weight='bold')

    # set other plot parameters
    plt.xlim((lon_min, lon_max))
    plt.ylim((lat_min, lat_max))
    plt.title('Wind Speed at Top of Boundary Layer (kt)')

    # Save the plot as a .png image in the "Downloads" folder
    plt.savefig(f'C:/Users/33017/Downloads/plot_{i:03d}.png', bbox_inches='tight')

    # Close the plot window
    plt.close()

# Close the netCDF file
ncfile.close()