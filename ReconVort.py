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
path_file = 'D:\wrfout_d02_2022-09-23_21_00_00'
ncfile = Dataset(path_file)

# get the total number of time indices 
num_time_indices = len(ncfile.dimensions['Time'])

initial_radar_lat = 16.929
initial_radar_lon = 125.5

# Loop through all time indices
for i in range(num_time_indices):

	# get the variables for the current time index (unchanged)
	# get the variables for the current time index (unchanged)
	p = getvar(ncfile, "pressure", timeidx=i)
	ua = getvar(ncfile, "ua", timeidx=i)
	va = getvar(ncfile, "va", timeidx=i)
	wind = getvar(ncfile, "WSPD10MAX", timeidx=i)
	ref = getvar(ncfile, "REFL_10CM", timeidx=i)
	
	# get the lat, lon grid 
	lats, lons = latlon_coords(wind)

    # Calculate the new radar position at each time step (100 km eastward movement)
	radar_lat = initial_radar_lat + 0 * (i / num_time_indices)
	radar_lon = initial_radar_lon + 4 * (i / num_time_indices)  # Gradually move east by 1 degree

	# Define pressure levels for interpolation (adjusted by 2.5 mb for every 1/222 degree away from the radar)
	delta_pressure_mb_per_degree = 2.5
	radar_distance_deg = np.sqrt((lons - radar_lon) ** 2 + (lats - radar_lat) ** 2)
	p_levels = [850 - delta_pressure_mb_per_degree * dist_deg * 85 for dist_deg in radar_distance_deg]

	u_level = np.array(interplevel(ua, p, p_levels), dtype=float)
	v_level = np.array(interplevel(va, p, p_levels), dtype=float)
	ref = interplevel(ref, p, p_levels)

	u_level1 = u_level * units("m/s")
	v_level1 = v_level * units("m/s")
	dx = 2.250 * units("km")
	dy = 2.250 * units("km")

	vort = mpcalc.vorticity(u_level1, v_level1, dx=dx, dy=dy)

	# specify your map boundaries
	lat_min = 14
	lat_max = 19
	lon_min = 121
	lon_max = 133

	# specify your colormap and projection 
	cmap = plt.get_cmap('twilight_shifted')
	crs = ccrs.PlateCarree()

	sigma = 0.0
	ref_smoothed = gaussian_filter(vort, sigma=sigma)
	
	# Create a grid of distances from the radar location
	dlon = lons - radar_lon
	dlat = lats - radar_lat
	radius = np.sqrt(dlon ** 2 + dlat ** 2)

	# Mask out the data outside the 100 km radius
	ref_smoothed = np.ma.masked_where(radius > 4.14, ref_smoothed)

	noise_level = 0
	ref_smoothed = ref_smoothed + np.random.randn(*ref_smoothed.shape) * noise_level

	# Apply masking to ref_smoothed where ref_smoothed2 is less than or equal to 0 dBz
	ref_smoothed = np.ma.masked_where(ref <= 0, ref_smoothed)


	# Find the maximum and minimum velocity/vorticity values, excluding NaN pixels
	max_vorticity = np.nanmax(ref_smoothed)
	min_vorticity = np.nanmin(ref_smoothed)


    # Create a single subplot for vorticity
	fig, ax = plt.subplots(1, 1, figsize=(10, 6), subplot_kw={'projection': crs, 'facecolor': 'black'})

    # Plot the vorticity in the single subplot
	plot_vort = ax.pcolormesh(lons, lats, ref_smoothed, cmap=cmap, transform=crs, vmin=-10, vmax=10)
	ax.coastlines(resolution='10m', color='red', alpha=0.9)
	ax.set_title('Radar Vorticity ($10^{-5} s^{-1}$)')
	cbar = fig.colorbar(plot_vort, ax=ax, fraction=0.046, pad=0.04)
	cbar.ax.set_ylabel('Radar Vorticity ($10^{-5} s^{-1}$)')

    # Set gridlines and other properties for the subplot
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

    # Plot a gold-colored star at the radar coordinates
	ax.plot(radar_lon, radar_lat, marker='4', color='gold', markersize=15,
            markeredgewidth=1.5, markeredgecolor='gold')

    # Set plot limits
	ax.set_xlim((lon_min, lon_max))
	ax.set_ylim((lat_min, lat_max))

    # Add the time information in the top right corner
	time_var = ncfile.variables['Times']
	time_str = ''.join([i.decode('utf-8') for i in time_var[i]])
	ax.text(0.95, 0.95, time_str, transform=ax.transAxes,
            fontsize=12, ha='right', va='top', color='white')

    # Add max/min vorticity information to the plot (excluding NaN pixels)
	ax.text(0.05, 0.05, f'Max/Min Vorticity: ({max_vorticity:.2f}/{min_vorticity:.2f}) $10^{-5} s^{-1}$)',
            transform=ax.transAxes, fontsize=9, ha='left', va='bottom', color='black', fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))

    # Save the plot as a .png image in the "Downloads" folder
	plt.savefig(f'C:/Users/33017/Downloads/plot_{i:03d}.png', bbox_inches='tight')

    # Close the plot window
	plt.close()

# Close the netCDF file
ncfile.close()