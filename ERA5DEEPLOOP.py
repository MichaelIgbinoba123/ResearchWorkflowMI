import os
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

# Extract relevant variables
latitude = ds['latitude'].values
longitude = ds['longitude'].values

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

# Loop through all time steps
for time_step in range(len(ds['time'])):
    # Extract u and v components for 200mb and 850mb
    u_200 = ds['u'].sel(level=200).isel(time=time_step).values[mask]
    v_200 = ds['v'].sel(level=200).isel(time=time_step).values[mask]
    u_850 = ds['u'].sel(level=850).isel(time=time_step).values[mask]
    v_850 = ds['v'].sel(level=850).isel(time=time_step).values[mask]

    # Calculate the vector magnitude for 200mb and 850mb
    magnitude_200 = np.sqrt(u_200**2 + v_200**2)
    magnitude_850 = np.sqrt(u_850**2 + v_850**2)

    # Calculate the difference in vector magnitude
    magnitude_diff = magnitude_200 - magnitude_850

    # Create a wind plot with streamlines and vector magnitude difference pcolormesh
    fig, ax = plt.subplots(figsize=(12, 10), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    # Plot streamlines for 200mb
    stream_200 = ax.streamplot(longitude, latitude, u_200-u_850, v_200-v_850,
                               color='black', linewidth=1, density=2)

    cmap = colortables.get_colortable('NWSStormClearReflectivity')

    # Add vector magnitude difference pcolormesh
    pcolormesh = ax.pcolormesh(lon_grid, lat_grid, magnitude_diff.reshape(lon_grid.shape),
                               cmap=cmap, shading='auto', alpha=0.8, vmin=0, vmax=25)

    # Add colorbar for vector magnitude difference
    cbar_pcolormesh = plt.colorbar(pcolormesh, ax=ax, orientation='vertical', pad=0.08, shrink=0.5)
    cbar_pcolormesh.set_label('Shear Magnitude (m/s)')

    # Add map features
    ax.coastlines()
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)

    # Set plot title
    plt.title(f'200mb - 850mb Wind Shear | Time Step: {time_step}')

    # Save the figure in the Downloads folder
    download_folder = os.path.join(os.path.expanduser("~"), "Downloads")
    figure_path = os.path.join(download_folder, f'wind_shear_{time_step}.png')
    plt.savefig(figure_path)
    plt.close()

print("Figures saved in Downloads folder.")