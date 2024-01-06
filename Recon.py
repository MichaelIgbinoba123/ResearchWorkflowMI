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
path_file = 'D:\wrfout_d02_2022-09-24_06_00_00'
ncfile = Dataset(path_file)

# get the total number of time indices 
num_time_indices = len(ncfile.dimensions['Time'])

initial_radar_lat = 16.244
initial_radar_lon = 124.5

# Loop through all time indices
for i in range(num_time_indices):

	# get the variables for the current time index (unchanged)
	# get the variables for the current time index (unchanged)
	p = getvar(ncfile, "pressure", timeidx=i)
	ua = getvar(ncfile, "ua", timeidx=i)
	va = getvar(ncfile, "va", timeidx=i)
	wind = getvar(ncfile, "WSPD10MAX", timeidx=i)
	ref = getvar(ncfile, "REFL_10CM", timeidx=i)
	wind_kt = wind * 1.944
	wa = getvar(ncfile, "wa", timeidx=i)
	avo = g_vorticity.get_avo(ncfile, timeidx=i)
	height = g_geoht.get_height(ncfile, timeidx=i, units='m')
	mslp = getvar(ncfile, "AFWA_MSLP", timeidx=i)
	theta = g_temp.get_theta(ncfile, timeidx=i, units='K')
	temp = tk(p * 100, theta, meta=True, units='degC')
	dew = g_dewpoint.get_dp(ncfile, timeidx=i, units='degC')
	

	# get the lat, lon grid 
	lats, lons = latlon_coords(wind)

    # Calculate the new radar position at each time step (100 km eastward movement)
	radar_lat = initial_radar_lat + 0 * (i / num_time_indices)
	radar_lon = initial_radar_lon + 4 * (i / num_time_indices)  # Gradually move east by 1 degree

	# Define pressure levels for interpolation (adjusted by 2.5 mb for every 1/222 degree away from the radar)
	delta_pressure_mb_per_degree = 2.5
	radar_distance_deg = np.sqrt((lons - radar_lon) ** 2 + (lats - radar_lat) ** 2)
	p_levels = [850 - delta_pressure_mb_per_degree * dist_deg * 85 for dist_deg in radar_distance_deg]

	# Calculate the elevation angle (in degrees) for each grid point
	elev = 0.799

	u_level = np.array(interplevel(ua, p, p_levels), dtype=float)
	v_level = np.array(interplevel(va, p, p_levels), dtype=float)
	ref = interplevel(ref, p, p_levels)
	w_level = np.array(interplevel(wa, p, p_levels), dtype=float)


	elev_rad = np.radians(elev)

	# Calculate azimuth angle (in radians) relative to the true north direction
	az1 = np.arctan2(u_level, v_level)
	az_rad = np.radians(az1)

	U = u_level * np.sin(az_rad) * np.cos(elev_rad) 
	V = v_level * np.cos(az_rad) * np.cos(elev_rad) 
	W = w_level * np.sin(elev_rad)

	v_r = U + V + W

	# specify your map boundaries
	lat_min = 14
	lat_max = 19
	lon_min = 121
	lon_max = 133

	# specify your colormap and projection 
	cmap = colortables.get_colortable('NWS8bitVel')
	cmap2 = colortables.get_colortable('NWSStormClearReflectivity')
	crs = ccrs.PlateCarree()

	# Apply smoothing to ref using a 2D Gaussian filter 
	sigma = 0.0  
	ref_smoothed = gaussian_filter(v_r, sigma=sigma)

	sigma2 = 0.0
	ref_smoothed2 = gaussian_filter(ref, sigma=sigma2)

	
	# Create a grid of distances from the radar location
	dlon = lons - radar_lon
	dlat = lats - radar_lat
	radius = np.sqrt(dlon ** 2 + dlat ** 2)

	# Mask out the data outside the 100 km radius
	ref_smoothed = np.ma.masked_where(radius > 4.14, ref_smoothed)
	ref_smoothed2 = np.ma.masked_where(radius > 4.14, ref_smoothed2)


	noise_level = 2.5
	ref_smoothed = ref_smoothed + np.random.randn(*ref_smoothed.shape) * noise_level
	ref_smoothed2 = ref_smoothed2 + np.random.randn(*ref_smoothed2.shape) * noise_level


	# Apply masking to ref_smoothed where ref_smoothed2 is less than or equal to 0 dBz
	ref_smoothed = np.ma.masked_where(ref_smoothed2 <= 0, ref_smoothed)


	# Find the maximum and minimum velocity/vorticity values, excluding NaN pixels
	max_velocity = np.nanmax(ref_smoothed)
	min_velocity = np.nanmin(ref_smoothed)

	# Create two subplots
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 6), subplot_kw={'projection': crs, 'facecolor': 'black'})

	# Plot Reflectivity in the first subplot
	plot_ref = ax1.pcolormesh(lons, lats, ref_smoothed2, cmap=cmap2, transform=crs, vmin=-30, vmax=90)
	ax1.coastlines(resolution='10m', color='red', alpha=0.9)
	ax1.set_title('Reflectivity (dBz)')

	# Add a colorbar for the first subplot
	cbar1 = fig.colorbar(plot_ref, ax=ax1, fraction=0.046, pad=0.04)
	cbar1.ax.set_ylabel('Reflectivity (dBz)')

	# Plot the smoothed reflectivity (ref) in the second subplot
	plot_vr = ax2.pcolormesh(lons, lats, ref_smoothed, cmap=cmap, transform=crs, vmin=-80, vmax=80)
	ax2.coastlines(resolution='10m', color='red', alpha=0.9)
	ax2.set_title('Radar Velocity (m/s\u207B\u00B9)')

	# Add a colorbar for the second subplot
	cbar2 = fig.colorbar(plot_vr, ax=ax2, fraction=0.046, pad=0.04)
	cbar2.ax.set_ylabel('Radar Velocity (m/s\u207B\u00B9)')


	# Set gridlines for both subplots
	for ax in [ax1, ax2]:
		gl = ax.gridlines(crs=crs, draw_labels=True, alpha=0.5)
		gl.top_labels = None
		gl.right_labels = None
		xgrid = np.arange(lon_min, lon_max + 1, 1.)  # Updated to show 1-degree intervals
		ygrid = np.arange(lat_min, lat_max + 1, 1.)  # Updated to show 1-degree intervals
		gl.xlocator = mticker.FixedLocator(xgrid.tolist())
		gl.ylocator = mticker.FixedLocator(ygrid.tolist())
		gl.xformatter = LONGITUDE_FORMATTER
		gl.yformatter = LATITUDE_FORMATTER
		gl.xlabel_style = {'size': 10, 'color': 'black'}
		gl.ylabel_style = {'size': 10, 'color': 'black'}

	# Plot a gold-colored star at the radar coordinates for both subplots
	for ax in [ax1, ax2]:
	    ax.plot(radar_lon, radar_lat, marker='4', color='gold', markersize=15,
            markeredgewidth=1.5, markeredgecolor='gold')

	# Set plot limits for both subplots
	for ax in [ax1, ax2]:
 	   ax.set_xlim((lon_min, lon_max))
 	   ax.set_ylim((lat_min, lat_max))

    # Add the time information in the top right corner
	time_var = ncfile.variables['Times']
	time_str = ''.join([i.decode('utf-8') for i in time_var[i]])
	for ax in [ax1, ax2]:
		ax.text(0.95, 0.95, time_str, transform=ax.transAxes,
            fontsize=12, ha='right', va='top', color='white')

	# Add max/min velocity information to the radar velocity plot (excluding NaN pixels)
	ax2.text(0.05, 0.05, f'Max/Min Velocity: ({max_velocity:.2f}/{min_velocity:.2f}) m/s\u207B\u00B9)',
         transform=ax2.transAxes, fontsize=9, ha='left', va='bottom', color='black', fontweight='bold',
         bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))


    # Save the plot as a .png image in the "Downloads" folder
	plt.savefig(f'C:/Users/33017/Downloads/plot_{i:03d}.png', bbox_inches='tight')

    # Close the plot window
	plt.close()

# Close the netCDF file
ncfile.close()