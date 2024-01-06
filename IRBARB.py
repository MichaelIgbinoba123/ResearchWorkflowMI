#imports
import numpy as np
import netCDF4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
from wrf import getvar, interplevel, to_np, latlon_coords, g_ctt, ALL_TIMES, g_geoht, g_terrain, g_temp, tk, ctt
from metpy.units import units
from scipy.ndimage import gaussian_filter
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
from cpt_convert import loadCPT
import math


# read the file
path_file = 'D:\wrfout_d02_2022-09-23_21_00_00'
ncfile = Dataset(path_file)

# get the total number of time indices
num_time_indices = len(ncfile.dimensions['Time'])

# Loop through all time indices
for i in range(num_time_indices):

    # get the variables for the current time index
    p = getvar(ncfile, "pressure", timeidx=i)
    theta = g_temp.get_theta(ncfile, timeidx=i, units='K')
    tkel = tk(p * 100, theta, meta=True, units='K')
    height = g_geoht.get_height(ncfile, timeidx=i, units='m')
    terrain = g_terrain.get_terrain(ncfile, timeidx=i, units='m')
    qice1 = getvar(ncfile, "QICE", timeidx=i)
    qv = getvar(ncfile, "QVAPOR", timeidx=i)
    qcld = getvar(ncfile, "QCLOUD", timeidx=i)
    olr = getvar(ncfile, "OLR", timeidx=i)
    th = getvar(ncfile, "TH2", timeidx=i)
    pres_hpa = p / 100
    EIR = np.power(olr, 1/4) * 83.98 - 328.83
    IR = getvar(ncfile, "ctt", timeidx=i)
    SSTSK = getvar(ncfile, "SSTSK", timeidx=i)
    SST = SSTSK - 273.15 
    RMPI = (28.2 + (55.8 * np.exp(0.1813 * (SST - 30)))) * 1.944
    cdck = 0.9
    OL = ctt(p, tkel, qv, qcld, height, terrain, units='degC')
    U = getvar(ncfile, "U10", timeidx=i)
    V = getvar(ncfile, "V10", timeidx=i)
    U10 = U * 1.944
    V10 = V * 1.944

    # specify your map boundaries
    lat_min = 15
    lat_max = 18
    lon_min = 124
    lon_max = 130

    # get the lat, lon grid
    lats, lons = latlon_coords(qv)

    colortable = loadCPT('C:/Users/33017/Desktop/WRF Scripts/IR4AVHRR6.cpt')
    ct = LinearSegmentedColormap('cpt', colortable)
    crs = ccrs.PlateCarree()

    # Apply smoothing to th_200 using a 2D Gaussian filter
    sigma = 1.0  # You can adjust this value for more or less smoothing
    IR_smoothed = gaussian_filter(IR, sigma=sigma)

    time_var = ncfile.variables['Times']
    time_str = ''.join([i.decode('utf-8') for i in time_var[i]])

    # plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, facecolor='gray', projection=crs)
    ax.coastlines(resolution='10m', alpha=0.5)

    # Use vmin and vmax to set the color range from -100 to 50 C
    IR_final = ax.pcolormesh(lons, lats, IR_smoothed, cmap=ct,     norm=plt.Normalize(-100, 40))


    cbar = fig.colorbar(IR_final)
    cbar.ax.set_ylabel('Temperature (C)')

# Add the 500 hPa wind barbs, only plotting every 125th data point.
    plt.barbs(to_np(lons[::14, ::14]), to_np(lats[::14, ::14]),
              to_np(U10[::14, ::14]), to_np(V10[::14, ::14]),
              length=4, color='gold')

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
    ax.text(1.0, 1.07, time_str, transform=ax.transAxes,
            ha='right', va='top', fontsize=8, color='black', backgroundcolor='white', weight='bold')

    # set other plot parameters
    plt.xlim((lon_min, lon_max))
    plt.ylim((lat_min, lat_max))
    plt.title('Infrared Imagery (C) + RMPI Contours (kt)')

    # Save the plot as a .png image in the "Downloads" folder
    plt.savefig(f'C:/Users/33017/Downloads/plot_{i:03d}.png', bbox_inches='tight')

    # Close the plot window
    plt.close()

# Close the netCDF file
ncfile.close()