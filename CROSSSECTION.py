import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset
from wrf import to_np, getvar, CoordPair, latlon_coords, vertcross, g_geoht, slp, g_vorticity, interplevel, smooth2d, g_temp, tk
from matplotlib.ticker import (NullFormatter, ScalarFormatter)


# Open the NetCDF file
path_file = 'D:\wrfout_d02_2023-04-12_12_05_00'
ncfile = Dataset(path_file)

# Create a list of pressure levels from 925 to 700 at intervals of 1
p_levels = [925 - i for i in range(76)]
p_levels2 = [925 - i for i in range(76)]

# get the total number of time indices
num_time_indices = len(ncfile.dimensions['Time'])

for i in range(num_time_indices):
    # Extract the model height and wind speed
    z = getvar(ncfile, "z", timeidx=i)
    height = g_geoht.get_height(ncfile, timeidx=i, units='m')
    tkel = getvar(ncfile, "T", timeidx=i)
    P = getvar(ncfile, "p", timeidx=i)
    p = getvar(ncfile, "pressure", timeidx=i)
    PB = getvar(ncfile, "PB", timeidx=i)
    pres = P + PB
    qv = getvar(ncfile, "QVAPOR", timeidx=i)
    mslp = getvar(ncfile, "AFWA_MSLP", timeidx=i)
    wspd = getvar(ncfile, "uvmet_wspd_wdir", units="kt")[0, :]
    avo = g_vorticity.get_avo(ncfile, timeidx=i)
    ua = getvar(ncfile, "ua", timeidx=i)
    va = getvar(ncfile, "va", timeidx=i)
    ref = getvar(ncfile, "REFL_10CM", timeidx=i)

    theta = g_temp.get_theta(ncfile, timeidx=i, units='K')

    temp = tk(pres, theta, units='K')

    wspd_smooth = smooth2d(wspd, 0)

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
        a_min_lat_idx, a_min_lon_idx = np.unravel_index(np.argmin(to_np(a_level)), a_level.shape)
        a_min_lat_arr, a_min_lon_arr = latlon_coords(a_level)
        a_min_lat.append(float(to_np(a_min_lat_arr[a_min_lat_idx, a_min_lon_idx])))
        a_min_lon.append(float(to_np(a_min_lon_arr[a_min_lat_idx, a_min_lon_idx])))

    # Take the average of the above locations to get the point for the cross section
    avg_lat = (sum(z_min_lat) + (mslp_min_lat * 50) + sum(a_min_lat)) / (len(z_min_lat) + len(a_min_lat) + 50)
    avg_lon = (sum(z_min_lon) + (mslp_min_lon * 50) + sum(a_min_lon)) / (len(z_min_lon) + len(a_min_lon) + 50)

    # Adjust the start and end points for the diagonal cross-section (move 1x1 degrees diagonally)
    start_point = CoordPair(lat=avg_lat + 0.5, lon=avg_lon + 0.5)
    end_point = CoordPair(lat=avg_lat - 0.5, lon=avg_lon - 0.5)

    # Compute the vertical cross-section interpolation. Also, include the
    # lat/lon points along the cross-section.
    avo_cross = vertcross(wspd, p, wrfin=ncfile, start_point=start_point,
                           end_point=end_point, latlon=True, meta=True)

    cmap = plt.get_cmap('nipy_spectral')
    levels = np.linspace(0, 140, 500)  # Define the color scale levels

    # Create the figure
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes()

    # Set the background color to gray
    ax.set_facecolor('lightgray')
   

    # Make the contour plot with fixed color table value range
    avo_mesh = ax.contourf(avo_cross, cmap=cmap, levels=levels, vmin=0, vmax=140)

# Add the color bar
    cbar = plt.colorbar(avo_mesh, ax=ax, orientation='vertical', pad=0.05, aspect=40)
    cbar.set_label('Wind Speed (kts)', fontsize=12)  # Add a label to the color bar

    # Set the x-ticks to use latitude and longitude labels.
    coord_pairs = to_np(avo_cross.coords["xy_loc"])
    x_ticks = np.arange(coord_pairs.shape[0])
    x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}") for pair in to_np(coord_pairs)]
    ax.set_xticks(x_ticks[::20])
    ax.set_xticklabels(x_labels[::20], rotation=45, fontsize=8)

    # Set the y-ticks to be height.
    vert_vals = to_np(avo_cross.coords["vertical"])
    v_ticks = np.arange(vert_vals.shape[0])
    ax.set_yticks(v_ticks[::20])
    ax.set_yticklabels(vert_vals[::20], fontsize=8)

    # Set the x-axis and y-axis labels
    ax.set_xlabel("Latitude, Longitude", fontsize=12)
    ax.set_ylabel("Pressure (hPa)", fontsize=12)

    plt.title("Vertical Cross Section of Wind Speed (kts)")

    # Save the plot as a .png image in the "Downloads" folder
    plt.savefig(f'C:/Users/33017/Downloads/plot_{i:03d}.png', bbox_inches='tight')

    # Close the plot window
    plt.close()

# Close the netCDF file
ncfile.close()