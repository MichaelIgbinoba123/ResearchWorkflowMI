import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import numpy as np
from metpy.plots import colortables

# Specify the path to the netCDF file
file_path = r'D:\Noru2022.nc'

# Open the netCDF file using xarray
ds = xr.open_dataset(file_path)

# Extract relevant variables for the first time step
latitude = ds['latitude'].values
longitude = ds['longitude'].values

# Extract u and v components for 200mb and 850mb
u_200 = ds['u'].sel(level=200).isel(time=66).values
v_200 = ds['v'].sel(level=200).isel(time=66).values
u_850 = ds['u'].sel(level=850).isel(time=66).values
v_850 = ds['v'].sel(level=850).isel(time=66).values

# Calculate the vector magnitude for 200mb and 850mb
magnitude_200 = np.sqrt(u_200**2 + v_200**2)
magnitude_850 = np.sqrt(u_850**2 + v_850**2)

# Calculate the difference in vector magnitude
magnitude_diff = magnitude_200 - magnitude_850

# Specify latitude and longitude bounds
lat_min, lat_max = 0, 30
lon_min, lon_max = 110, 180

# Create a 2D grid of latitude and longitude
lon_grid, lat_grid = np.meshgrid(longitude, latitude)

# Create a mask for the specified bounds
mask = (lat_grid >= lat_min) & (lat_grid <= lat_max) & (lon_grid >= lon_min) & (lon_grid <= lon_max)

# Extract data within the specified bounds
latitude = lat_grid[mask]
longitude = lon_grid[mask]
u_200 = u_200[mask]
v_200 = v_200[mask]
u_850 = u_850[mask]
v_850 = v_850[mask]
magnitude_diff = magnitude_diff[mask]

# Create a wind plot with streamlines and vector magnitude difference pcolormesh
fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

# Plot streamlines for 200mb
stream_200 = ax.streamplot(longitude, latitude, u_200-u_850, v_200-v_850,
                           color='black', linewidth=1, density=2)

cmap = colortables.get_colortable('NWSStormClearReflectivity')

# Add vector magnitude difference pcolormesh
pcolormesh = ax.pcolormesh(lon_grid, lat_grid, magnitude_diff.reshape(lon_grid.shape),
                           cmap=cmap, shading='auto', alpha=0.8, vmin=0, vmax=25)

# Add colorbar for vector magnitude difference
cbar_pcolormesh = plt.colorbar(pcolormesh, ax=ax, orientation='vertical', pad=0.08, shrink=0.75)
cbar_pcolormesh.set_label('Shear Magniude (m/s)')

# Add map features
ax.coastlines()
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

# Set plot title
plt.title('200mb - 850mb Wind Shear')

# Show the plot
plt.show()