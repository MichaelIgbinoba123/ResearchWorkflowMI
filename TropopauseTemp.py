# imports 
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
from wrf import (getvar, interplevel, to_np, latlon_coords, g_geoht, srhel, g_vorticity, g_dewpoint, tk, g_temp)
from metpy.units import units
from metpy.plots import colortables
from scipy.ndimage import gaussian_filter
from math import radians, sin, cos, sqrt, atan2
import metpy.calc as mpcalc

# read the file 
path_file = 'D:\wrfout_d02_2022-09-23_18_00_00'
ncfile = Dataset(path_file)

# get the total number of time indices 
num_time_indices = len(ncfile.dimensions['Time'])

# Loop through all time indices
for i in range(num_time_indices):

	# get the variables for the current time index (unchanged)
	# get the variables for the current time index (unchanged)
	p = getvar(ncfile, "pressure", timeidx=i)
	pres = p * 100
	ua = getvar(ncfile, "ua", timeidx=i)
	va = getvar(ncfile, "va", timeidx=i)
	wind = getvar(ncfile, "WSPD10MAX", timeidx=i)
	ref = getvar(ncfile, "REFL_10CM", timeidx=i)
	theta = g_temp.get_theta(ncfile, timeidx=i, units='K')
	te = tk(pres, theta, meta=True, units='degC')
	
	# get the lat, lon grid 
	lats, lons = latlon_coords(wind)

# Find the lowest te value at each pixel	
	temp = np.min(te, axis=0)

    # specify your map boundaries
	lat_min = 14
	lat_max = 19
	lon_min = 121
	lon_max = 133

    # get the lat, lon grid
	lats, lons = latlon_coords(ua)

    # specify your colormap and projection
	cmap = colortables.get_colortable('ir_rgbv')
	crs = ccrs.PlateCarree()

    # Apply smoothing to th_200 using a 2D Gaussian filter
	sigma = 1.0  # You can adjust this value for more or less smoothing
	wind_bl_smoothed = gaussian_filter(temp, sigma=sigma)

	time_var = ncfile.variables['Times']
	time_str = ''.join([i.decode('utf-8') for i in time_var[i]])

    # plot
	fig = plt.figure(figsize=(10, 6))
	ax = fig.add_subplot(111, facecolor='gray', projection=crs)
	ax.coastlines(resolution='10m', alpha=0.5)
    # Use vmin and vmax to set the color range from 260 to 360 K
	plot_wind_bl_smoothed = ax.pcolormesh(lons, lats, wind_bl_smoothed, cmap=cmap, norm=plt.Normalize(0, -100))
	cbar = fig.colorbar(plot_wind_bl_smoothed)
	cbar.ax.set_ylabel('Temperature (C)')

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

    # set other plot parameters
	plt.xlim((lon_min, lon_max))
	plt.ylim((lat_min, lat_max))
	plt.title('Tropopause Temperature (C)')

    # Save the plot as a .png image in the "Downloads" folder
	plt.savefig(f'C:/Users/33017/Downloads/plot_{i:03d}.png', bbox_inches='tight')

    # Close the plot window
	plt.close()

# Close the netCDF file
ncfile.close()