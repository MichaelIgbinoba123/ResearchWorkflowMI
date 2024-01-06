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
path_file = 'D:\wrfout_d02_2022-09-24_12_05_00'
ncfile = Dataset(path_file)

# select your pressure level
p_level = 200

# get the total number of time indices
num_time_indices = len(ncfile.dimensions['Time'])

# Loop through all time indices
for i in range(num_time_indices):
    # get the variables for the current time index
    p = getvar(ncfile, "pressure", timeidx=i)
    ua = getvar(ncfile, "ua", timeidx=i)
    va = getvar(ncfile, "va", timeidx=i)
    te = getvar(ncfile, "T", timeidx=i)

    theta = g_temp.get_theta(ncfile, timeidx=i, units='K')

    # interpolate ua and va to pressure level
    u_500 = interplevel(ua, p, p_level)
    v_500 = interplevel(va, p, p_level)
    th_200 = interplevel(theta, p, p_level)

    # specify your map boundaries
    lat_min = 14
    lat_max = 19
    lon_min = 121
    lon_max = 133

    # get the lat, lon grid
    lats, lons = latlon_coords(th_200)

    # specify your colormap and projection
    cmap = plt.get_cmap('RdBu_r')
    crs = ccrs.PlateCarree()

    # Compute the mean potential temperature at 200mb
    th_200_mean = np.mean(th_200)

    # Subtract the mean from the actual potential temperature field
    th_200_anomaly = th_200 - th_200_mean

    # Apply smoothing to th_200_anomaly using a 2D Gaussian filter
    sigma = 0.0  # You can adjust this value for more or less smoothing
    th_200_anomaly_smoothed = gaussian_filter(th_200_anomaly, sigma=sigma)

    # Calculate the maximum anomaly
    current_max_anomaly = np.max(th_200_anomaly_smoothed)

    time_var = ncfile.variables['Times']
    time_str = ''.join([i.decode('utf-8') for i in time_var[i]])

    # plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, facecolor='None', projection=crs)
    ax.coastlines(resolution='10m', alpha=0.5)

    # Use vmin and vmax to set the color range for the anomaly plot
    plot_th200_anomaly = ax.pcolormesh(lons, lats, th_200_anomaly_smoothed, cmap=cmap, norm=plt.Normalize(-10, 10))
    cbar = fig.colorbar(plot_th200_anomaly)
    cbar.ax.set_ylabel('Potential Temp Anomaly (K)')

    # some fancy schmancy grid lines
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

    # Create streamlines
    # Spacing for streamlines, you can adjust this value to control the density of streamlines
    x_stride = 10
    y_stride = 10
    ax.streamplot(lons[::y_stride, ::x_stride], lats[::y_stride, ::x_stride],
                  to_np(u_500[::y_stride, ::x_stride]), to_np(v_500[::y_stride, ::x_stride]),
                  linewidth=1.0, density=2, arrowsize=1.5, color='black', transform=ccrs.PlateCarree())

    # Add the time information in the top right corner
    time_var = ncfile.variables['Times']
    time_str = ''.join([i.decode('utf-8') for i in time_var[i]])
    ax.text(0.95, 0.95, time_str, transform=ax.transAxes,
            ha='right', va='top', fontsize=8, color='black', backgroundcolor='white', weight='bold')

    # Add the max anomaly information in the top left corner
    ax.text(0.05, 0.95, f'Max anomaly: {current_max_anomaly:.2f} K', transform=ax.transAxes,
            ha='left', va='top', fontsize=8, color='black', backgroundcolor='white', weight='bold')


    # set other plot parameters
    plt.xlim((lon_min, lon_max))
    plt.ylim((lat_min, lat_max))
    plt.title('200mb Wind Vector + 200mb Potential Temperature Anomaly')

    # Save the plot as a .png image in the "Downloads" folder
    plt.savefig(f'C:/Users/33017/Downloads/plot_{i:03d}.png', bbox_inches='tight')

    # Close the plot window
    plt.close()

# Close the netCDF file
ncfile.close()
