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
path_file = 'D:\wrfout_d02_2023-04-13_12_05_00'
ncfile = Dataset(path_file)

# select your pressure level
p_level = 850

# get the total number of time indices
num_time_indices = len(ncfile.dimensions['Time'])

# Loop through all time indices
for i in range(num_time_indices):
    # get the variables for the current time index
    p = getvar(ncfile, "pressure", timeidx=i)
    PV = getvar(ncfile, "pvo", timeidx=i)
    slp = getvar(ncfile, "AFWA_MSLP", timeidx=i)
    mslp = slp / 100
    

    theta = g_temp.get_theta(ncfile, timeidx=1, units='K')

    # interpolate ua and va to pressure level
    PV_850 = interplevel(PV, p, p_level)

    # specify your map boundaries
    lat_min = -20
    lat_max = -18
    lon_min = 119
    lon_max = 121

    # get the lat, lon grid
    lats, lons = latlon_coords(PV)

    # specify your colormap and projection
    cmap = plt.get_cmap('coolwarm')
    crs = ccrs.PlateCarree()

    # Apply smoothing to th_200 using a 2D Gaussian filter
    sigma = 1.0  # You can adjust this value for more or less smoothing
    PV_850_smoothed = gaussian_filter(PV_850, sigma=sigma)

    time_var = ncfile.variables['Times']
    time_str = ''.join([i.decode('utf-8') for i in time_var[i]])

    # plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, facecolor='None', projection=crs)
    ax.coastlines(resolution='10m', alpha=0.5)
    # Use vmin and vmax to set the color range from 260 to 360 K
    plot_PV_850_smoothed = ax.pcolormesh(lons, lats, PV_850, cmap=cmap,     norm=plt.Normalize(-25, 25))
    cbar = fig.colorbar(plot_PV_850_smoothed)
    cbar.ax.set_ylabel('Potential Vorticity (PVU)')

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
    ax.text(0.95, 0.95, time_str, transform=ax.transAxes,
            ha='right', va='top', fontsize=8, color='black', backgroundcolor='white', weight='bold')

    # Convert MSLP to numpy array
    mslp_np = to_np(mslp)

    # Specify the contour levels for MSLP
    mslp_contour_levels = np.arange(900, 1000, 2)

    plt.contour(lons, lats, mslp_np, levels=mslp_contour_levels, colors='black', linestyles='dashed', linewidths=1, transform=crs)

    # Find the minimum MSLP value and its corresponding latitude and longitude
    mslp_min = np.min(mslp_np)
    mslp_min_lat, mslp_min_lon = np.unravel_index(np.argmin(mslp_np), mslp_np.shape)

    # Add the text annotation for the minimum MSLP value in the top-left corner
    ax.text(0.05, 0.95, f'Minimum SLP = {mslp_min:.1f} hPa', transform=ax.transAxes,
            ha='left', va='top', fontsize=8, color='black', backgroundcolor='white', weight='bold')

    # set other plot parameters
    plt.xlim((lon_min, lon_max))
    plt.ylim((lat_min, lat_max))
    plt.title('850mb Potential Vorticity (PVU) + MSLP (hPa)')

    # Save the plot as a .png image in the "Downloads" folder
    plt.savefig(f'C:/Users/33017/Downloads/plot_{i:03d}.png', bbox_inches='tight')

    # Close the plot window
    plt.close()

# Close the netCDF file
ncfile.close()

