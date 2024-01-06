# imports
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
from wrf import (getvar, interplevel, to_np, latlon_coords, g_geoht, srhel, g_vorticity, rh, g_temp, tk, slp)
from metpy.units import units
from scipy.ndimage import gaussian_filter
from math import radians, sin, cos, sqrt, atan2
from metpy.plots import colortables


# read the file
path_file = 'D:\wrfout_d02_2022-09-24_06_00_00'
ncfile = Dataset(path_file)

# select your pressure level
p_level = 500
p_level2 = 850
p_level3 = 500
p_level4 = 850


# Create a list of pressure levels from 925 to 700 at intervals of 1
p_levels = [925 - i for i in range(76)]
p_levels2 = [925 - i for i in range(76)]

# get the total number of time indices
num_time_indices = len(ncfile.dimensions['Time'])

# Loop through all time indices
for i in range(num_time_indices):
    # get the variables for the current time index
    p = getvar(ncfile, "pressure", timeidx=i)
    P = getvar(ncfile, "P", timeidx=i)
    PB = getvar(ncfile, "PB", timeidx=i)
    theta = g_temp.get_theta(ncfile, timeidx=i, units='K')
    pres = P + PB
    ua = getvar(ncfile, "ua", timeidx=i)
    va = getvar(ncfile, "va", timeidx=i)
    height = g_geoht.get_height(ncfile, timeidx=i, units='m')
    tkel1 = tk(pres, theta, units='K')
    qv = getvar(ncfile, "QVAPOR", timeidx=i)
    terrain = getvar(ncfile, "HGT", timeidx=i)
    wind = getvar(ncfile, "WSPD10MAX", timeidx=i)
    mslp = getvar(ncfile, "AFWA_MSLP", timeidx=i)
    tkel = getvar(ncfile, "T", timeidx=i)
    RH = rh(qv, pres, tkel1)
    wind_kt = wind * 1.944
    avo = g_vorticity.get_avo(ncfile, timeidx=i)

        # Interpolate avo and height to pressure levels
    z_levels = [interplevel(height, p, p_level) for p_level in p_levels]
    a_levels = [interplevel(avo, p, p_level) for p_level in p_levels2]

    # Find the latitude and longitude of the lowest geopotential height at each pressure level
    z_min_lat = []
    z_min_lon = []
    for z_level in z_levels:
        z_min_lat_idx, z_min_lon_idx = np.unravel_index(np.argmin(to_np(z_level)), z_level.shape)
        z_min_lat_arr, z_min_lon_arr = latlon_coords(z_level)
        z_min_lat.append(float(to_np(z_min_lat_arr[z_min_lat_idx, z_min_lon_idx])))
        z_min_lon.append(float(to_np(z_min_lon_arr[z_min_lat_idx, z_min_lon_idx])))

    # Find the latitude and longitude of the minimum MSLP
    mslp_min_lat_idx, mslp_min_lon_idx = np.unravel_index(np.argmin(to_np(mslp)), mslp.shape)
    mslp_min_lat_arr, mslp_min_lon_arr = latlon_coords(mslp)
    mslp_min_lat = float(to_np(mslp_min_lat_arr[mslp_min_lat_idx, mslp_min_lon_idx]))
    mslp_min_lon = float(to_np(mslp_min_lon_arr[mslp_min_lat_idx, mslp_min_lon_idx]))

    # Find the latitude and longitude of the lowest absolute vorticity at each pressure level
    a_min_lat = []
    a_min_lon = []
    for a_level in a_levels:
        a_min_lat_idx, a_min_lon_idx = np.unravel_index(np.argmax(to_np(a_level)), a_level.shape)
        a_min_lat_arr, a_min_lon_arr = latlon_coords(a_level)
        a_min_lat.append(float(to_np(a_min_lat_arr[a_min_lat_idx, a_min_lon_idx])))
        a_min_lon.append(float(to_np(a_min_lon_arr[a_min_lat_idx, a_min_lon_idx])))

    # Take the average of the above locations to get the point for the cross section
    avg_lat = (sum(z_min_lat) + (mslp_min_lat * 50) + sum(a_min_lat)) / (len(z_min_lat) + len(a_min_lat) + 50)
    avg_lon = (sum(z_min_lon) + (mslp_min_lon * 50) + sum(a_min_lon)) / (len(z_min_lon) + len(a_min_lon) + 50)


    # interpolate ua and va to pressure level
    u_200 = interplevel(ua, p, p_level)
    v_200 = interplevel(va, p, p_level)
    u_850 = interplevel(ua, p, p_level2)
    v_850 = interplevel(ua, p, p_level2)
    RH_450 = interplevel(RH, p, p_level3)
    RH_750 = interplevel(RH, p, p_level4)

    #Create shear vector
    u200_850 = (u_200 - u_850)
    v200_850 = (v_200 - v_850)

    shear_mag = np.sqrt(u200_850**2 + v200_850**2)

    RH450_750 = (RH_450 + RH_750) / 2

    # specify your map boundaries
    lat_min = 14
    lat_max = 19
    lon_min = 121
    lon_max = 133

    # get the lat, lon grid
    lats, lons = latlon_coords(wind)

    # specify your colormap and projection
    cmap = colortables.get_colortable('WVCIMSS')
    crs = ccrs.PlateCarree()

    # Apply smoothing to th_200 using a 2D Gaussian filter
    sigma = 0.0  # You can adjust this value for more or less smoothing
    RH450_750_smoothed = gaussian_filter(RH450_750, sigma=sigma)

    time_var = ncfile.variables['Times']
    time_str = ''.join([i.decode('utf-8') for i in time_var[i]])

    # plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, facecolor='gray', projection=crs)
    ax.coastlines(resolution='10m', alpha=0.9)
    # Use vmin and vmax to set the color range from 0 to 100%
    plot_RH450_750 = ax.pcolormesh(lons, lats, RH450_750_smoothed, cmap=cmap, vmin=0, vmax=100)
    cbar = fig.colorbar(plot_RH450_750)
    cbar.ax.set_ylabel('Relative Humidity (%)')

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
                  to_np(u200_850[::y_stride, ::x_stride]), to_np(v200_850[::y_stride, ::x_stride]),
                  linewidth=1.0, density=1.5, arrowsize=1.0, color='black', transform=ccrs.PlateCarree())

    # Add the time information in the top right corner
    time_var = ncfile.variables['Times']
    time_str = ''.join([i.decode('utf-8') for i in time_var[i]])
    ax.text(0.95, 0.95, time_str, transform=ax.transAxes,
            ha='right', va='top', fontsize=6, color='black', backgroundcolor='white', weight='bold')

    # Calculate the average 200-850mb shear in a box centered on the MSLP minimum
    box_size = 2.5  # Adjust the box size as needed
    lat_box_min = int(max(0, avg_lat - box_size))
    lat_box_max = int(min(lats.shape[0], avg_lat + box_size + 1))
    lon_box_min = int(max(0, avg_lon - box_size))
    lon_box_max = int(min(lons.shape[1], avg_lon + box_size + 1))

# Convert the lons and lats DataArrays to NumPy arrays
    lons_np = lons.values
    lats_np = lats.values

# Calculate the average 200-850mb wind shear in the box, excluding values within the specified radius
    avg_shear = np.copy(shear_mag)

# Convert average latitude and longitude to radians
    avg_lat_rad = radians(avg_lat)
    avg_lon_rad = radians(avg_lon)

# Convert lats and lons to radians
    lats_rad = np.radians(lats_np)
    lons_rad = np.radians(lons_np)

# Calculate the Haversine distance for each grid point to the average position
    delta_lat = lats_rad - avg_lat_rad
    delta_lon = lons_rad - avg_lon_rad

    a = np.sin(delta_lat/2) ** 2 + np.cos(avg_lat_rad) * np.cos(lats_rad) * np.sin(delta_lon/2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    distance_deg = np.degrees(c)

# Set all wind shear values within the radius of 1.5 degrees to NaN (exclude from the average calculation)
    avg_shear[distance_deg <= 5] = np.nan

# Calculate the average 200-850mb shear for the remaining non-NaN values
    avg_shear_value = np.nanmean(avg_shear)

# Multiply the average shear by 1.944 to convert it to knots
    avg_shear_knots = avg_shear_value * 1.944

    # Add the average 200-850mb shear information in the top left corner
    ax.text(0.05, 0.95, f'Average 500-850mb shear = {avg_shear_knots:.2f} kts', transform=ax.transAxes,
            ha='left', va='top', fontsize=6, color='black', backgroundcolor='white', weight='bold')


    # Calculate the direction of the wind shear arrow
    avg_shear_dir = np.arctan2(np.nanmean(v200_850[lat_box_min:lat_box_max, lon_box_min:lon_box_max]),
                               np.nanmean(u200_850[lat_box_min:lat_box_max, lon_box_min:lon_box_max]))

    # Add the wind shear arrow in the center of the plot
    arrow_length = 4
    arrow_dx = arrow_length * np.cos(avg_shear_dir)
    arrow_dy = arrow_length * np.sin(avg_shear_dir)
    ax.quiver(lon_min + 11, lat_min + 1, arrow_dx, arrow_dy, angles='xy', scale_units='xy', scale=4,
              color='red', zorder=5, width=0.01, transform=ccrs.PlateCarree())


    # set other plot parameters
    plt.xlim((lon_min, lon_max))
    plt.ylim((lat_min, lat_max))
    plt.title('500-850mb Relative Humidity + 500-850mb Wind Shear')

    # Save the plot as a .png image in the "Downloads" folder
    plt.savefig(f'C:/Users/33017/Downloads/plot_{i:03d}.png', bbox_inches='tight')

    # Close the plot window
    plt.close()

# Close the netCDF file
ncfile.close()