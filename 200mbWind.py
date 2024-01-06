import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
from metpy.plots import colortables

# Specify the path to the netCDF file
file_path = r'D:\Katrina2005.nc'

# Open the netCDF file using xarray
ds = xr.open_dataset(file_path)

# Extract relevant variables for the first time step
latitude = ds['latitude'].values
longitude = ds['longitude'].values
u_component = ds['u'].sel(level=200).isel(time=75).values
v_component = ds['v'].sel(level=200).isel(time=75).values
geopotential_850 = ds['z'].sel(level=850).isel(time=66).values

# Calculate the vector magnitude
magnitude = np.sqrt(u_component**2 + v_component**2)

# Specify latitude and longitude bounds
lat_min, lat_max = 15, 45
lon_min, lon_max = -110, -60

# Create a 2D grid of latitude and longitude
lon_grid, lat_grid = np.meshgrid(longitude, latitude)

# Create a mask for the specified bounds
mask = (lat_grid >= lat_min) & (lat_grid <= lat_max) & (lon_grid >= lon_min) & (lon_grid <= lon_max)

# Extract data within the specified bounds
latitude = lat_grid[mask]
longitude = lon_grid[mask]
u_component = u_component[mask]
v_component = v_component[mask]
magnitude = magnitude[mask]

# Find the location of the minimum geopotential height at 850mb
min_index = np.unravel_index(np.argmin(geopotential_850), geopotential_850.shape)
min_lat, min_lon = latitude[min_index[0]], longitude[min_index[1]]

# Create a wind plot with streamlines and vector magnitude pcolormesh
fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

# Plot streamlines
stream = ax.streamplot(longitude, latitude, u_component, v_component,
                       color='black', linewidth=1, density=2)

cmap = colortables.get_colortable('NWSStormClearReflectivity')

# Add vector magnitude pcolormesh
pcolormesh = ax.pcolormesh(lon_grid, lat_grid, magnitude.reshape(lon_grid.shape),
                           cmap=cmap, shading='auto', alpha=0.8, vmin=0, vmax=90)

# Add red dot for the location of the lowest geopotential height at 850mb
ax.plot(min_lon, min_lat, 'ro', markersize=8, transform=ccrs.PlateCarree())

# Add colorbar for vector magnitude
cbar_pcolormesh = plt.colorbar(pcolormesh, ax=ax, orientation='vertical', pad=0.08, shrink=0.75)
cbar_pcolormesh.set_label('Wind Speed (m/s)')

# Add map features
ax.coastlines()
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

# Set plot title
plt.title('200mb Wind Streamlines with Vector Magnitude (Red Dot = Min 850mb GPH)')

# Show the plot
plt.show()