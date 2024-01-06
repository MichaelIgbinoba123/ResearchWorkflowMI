# imports 
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
from wrf import (getvar, interplevel, interpline, to_np, latlon_coords, g_geoht, srhel, g_vorticity, g_dewpoint, tk, CoordPair, g_temp)
from metpy.units import units
from metpy.plots import colortables
from scipy.ndimage import gaussian_filter
from math import radians, sin, cos, sqrt, atan2


# read the file 
path_file = 'D:\wrfout_d02_2022-09-24_06_00_00'
ncfile = Dataset(path_file)

initial_radar_lat = 16.244
initial_radar_lon = 124.5

# get the variables for the current time index (unchanged)
p = getvar(ncfile, "pressure", timeidx=35)
ua = getvar(ncfile, "ua", timeidx=35)
va = getvar(ncfile, "va", timeidx=35)
wind = getvar(ncfile, "WSPD10MAX", timeidx=35)
ref1 = getvar(ncfile, "REFL_10CM", timeidx=35)
wind_kt = wind * 1.944
wa = getvar(ncfile, "wa", timeidx=35)
avo = g_vorticity.get_avo(ncfile, timeidx=35)
height = g_geoht.get_height(ncfile, timeidx=35, units='m')
mslp = getvar(ncfile, "AFWA_MSLP", timeidx=35)
turbulence = getvar(ncfile, "rh", timeidx=35)
theta1 = g_temp.get_theta(ncfile, timeidx=35, units='K')
temp = tk(p * 100, theta1, meta=True, units='degC')
dew = g_dewpoint.get_dp(ncfile, timeidx=35, units='degC')
XLAT = getvar(ncfile, "XLAT")
XLONG = getvar(ncfile, "XLONG")

lats = XLAT[:, :]
lons = XLONG[:, :]

radar_lat = initial_radar_lat 
radar_lon = initial_radar_lon + 4 

p_level1 = 850

w700 = (interplevel(wa, p, p_level1)) * 1.944
u700 = interplevel(ua, p, p_level1)
v700 = interplevel(va, p, p_level1)
wind700 = (np.sqrt(u700**2 + v700**2)) * 1.944
temp700 = interplevel(temp, p, p_level1)
dew700 = interplevel(dew, p, p_level1)
ref700 = interplevel(ref1, p, p_level1)
rr700 = (0.017) * (10**(0.0714 * ref700))
MSLP = mslp / 100
turb700 = interplevel(turbulence, p, p_level1)

init = CoordPair(lat=initial_radar_lat, lon=initial_radar_lon)
final = CoordPair(lat=radar_lat, lon=radar_lon)

w700interp = interpline(w700, ncfile, timeidx=35, start_point=init, end_point=final)

wind700interp = interpline(wind700, ncfile, timeidx=35, start_point=init, end_point=final)

temp700interp = interpline(temp700, ncfile, timeidx=35, start_point=init, end_point=final)

dew700interp = interpline(dew700, ncfile, timeidx=35, start_point=init, end_point=final)

rr700interp = interpline(rr700, ncfile, timeidx=35, start_point=init, end_point=final)

MSLPinterp = interpline(MSLP, ncfile, timeidx=35, start_point=init, end_point=final)

turb700interp = interpline(turb700, ncfile, timeidx=35, start_point=init, end_point=final)

SFCWindinterp = interpline(wind_kt, ncfile, timeidx=35, start_point=init, end_point=final)

# Add noise to wind700, temp700, dew700, rr700, MSLP, and theta700
noise_level = 0.0  # Adjust this value to control the amount of noise

np.random.seed(42)  # Set a seed for reproducibility

wind700interp_noisy = wind700interp + np.random.randn(*wind700interp.shape) * noise_level
temp700interp_noisy = temp700interp + np.random.randn(*temp700interp.shape) * noise_level
dew700interp_noisy = dew700interp + np.random.randn(*dew700interp.shape) * noise_level
rr700interp_noisy = rr700interp + np.random.randn(*rr700interp.shape) * noise_level
MSLPinterp_noisy = MSLPinterp + np.random.randn(*MSLPinterp.shape) * noise_level
turb700interp_noisy = turb700interp + np.random.randn(*turb700interp.shape) * noise_level
SFCWindinterp_noisy = SFCWindinterp + np.random.randn(*SFCWindinterp.shape) * noise_level
w700interp_noisy = w700interp + np.random.randn(*w700interp.shape) * noise_level


# Interpolate lons to match the shape of MSLPinterp
lons_interp = np.interp(np.linspace(0, 1, len(MSLPinterp_noisy)), np.linspace(0, 1, len(lons[0])), lons[0])

# Create the first plot with MSLPinterp and wind700interp
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

# Custom function to extract and filter the desired longitude range
def extract_and_filter(lon, var):
    idx = np.where((lon >= initial_radar_lon) & (lon <= radar_lon))
    return lon[idx], var[idx]

# Extract and filter the interpolated values within the desired longitude range
lons_interp_filtered, MSLPinterp_filtered = extract_and_filter(lons_interp, MSLPinterp_noisy)
lons_interp_filtered, wind700interp_filtered = extract_and_filter(lons_interp, wind700interp_noisy)
lons_interp_filtered, rr700interp_filtered = extract_and_filter(lons_interp, rr700interp_noisy)
lons_interp_filtered, temp700interp_filtered = extract_and_filter(lons_interp, temp700interp_noisy)
lons_interp_filtered, w700interp_filtered = extract_and_filter(lons_interp, w700interp_noisy)
lons_interp_filtered, turb700interp_filtered = extract_and_filter(lons_interp, turb700interp_noisy)

# Create the first plot with MSLPinterp and wind700interp
ax1 = axs[0, 1]
ax2 = ax1.twinx()

# Plot the filtered data
ax1.plot(lons_interp_filtered, MSLPinterp_filtered, 'r', label='Sea Level Pressure (hPa)')
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Sea Level Pressure (hPa)')
ax1.tick_params(axis='y')
ax1.grid(True)  # Add gridlines to ax1
ax2.plot(lons_interp_filtered, wind700interp_filtered, 'b', label='850mb Wind (kt)')
ax2.set_ylabel('850mb Wind (kt)')
ax2.tick_params(axis='y')
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

# Extract and filter the interpolated values within the desired longitude range
lons_interp_filtered, SFCWindinterp_filtered = extract_and_filter(lons_interp, SFCWindinterp_noisy)
lons_interp_filtered, dew700interp_filtered = extract_and_filter(lons_interp, dew700interp_noisy)

# Second plot - [1, 1]
ax1 = axs[1, 1]
ax2 = ax1.twinx()

# Plot the filtered data
ax1.plot(lons_interp_filtered, rr700interp_filtered, 'r', label='Rain Rate (mm/hr)')
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Rain Rate (mm/hr)')
ax1.tick_params(axis='y')
ax1.grid(True)  # Add gridlines to ax1
ax2.plot(lons_interp_filtered, SFCWindinterp_filtered, 'b', label='Surface Wind (kt)')
ax2.set_ylabel('Surface Wind (kt)')
ax2.tick_params(axis='y')
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

# Third plot - [0, 0]
ax1 = axs[0, 0]
ax2 = ax1.twinx()

# Plot the filtered data
ax1.plot(lons_interp_filtered, temp700interp_filtered, 'orange', label='850mb Temp (C)')
ax1.set_xlabel('Longitude')
ax1.set_ylabel('850mb Temp (C)')
ax1.set_ylim(-40, 40)
ax1.tick_params(axis='y')
ax1.grid(True)  # Add gridlines to ax1
ax2.plot(lons_interp_filtered, dew700interp_filtered, 'g', label='850mb Dewpoint (C)')
ax2.set_ylabel('850mb Dewpoint (C)')
ax2.set_ylim(-40, 40)
ax2.tick_params(axis='y')
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

# Fourth plot - [1, 0]
ax1 = axs[1, 0]
ax2 = ax1.twinx()

# Plot the filtered data
ax1.plot(lons_interp_filtered, w700interp_filtered, 'blue', label='850mb Vert. Motion (m/s)')
ax1.set_xlabel('Longitude')
ax1.set_ylabel('850mb Vertical Motion (kt)')
ax1.tick_params(axis='y')
ax1.grid(True)  # Add gridlines to ax1
ax2.plot(lons_interp_filtered, turb700interp_filtered, 'r', label='850mb Rel. Humidity (%)')
ax2.set_ylabel('850mb Rel. Humidity (%)')
ax2.tick_params(axis='y')
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

# Adjust the layout
fig.tight_layout()

plt.show()