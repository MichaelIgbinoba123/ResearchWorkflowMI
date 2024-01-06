import os
import pyart
import fsspec 
from metpy.plots import USCOUNTIES
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import numpy as np
import warnings

fs = fsspec.filesystem("s3", anon=True)

# Get a list of all files in the specified directory
files = sorted(fs.glob("s3://noaa-nexrad-level2/2020/09/14/KEVX/KEVX20200914*"))

# Exclude files with "_MDM" in their names
files = [file for file in files if "_MDM" not in os.path.basename(file)]

batch_size = 2  # Adjust as needed
for i in range(0, len(files), batch_size):
    batch_files = files[i:i+batch_size]
    for file in batch_files:
        # Read radar data
        radar = pyart.io.read_nexrad_archive(f's3://{file}')

        # Create a gate filter which specifies gates to exclude from dealiasing
        gatefilter = pyart.filters.GateFilter(radar)
        gatefilter.exclude_transition()
        gatefilter.exclude_invalid("velocity")
        gatefilter.exclude_invalid("reflectivity")
        gatefilter.exclude_outside("reflectivity", 0, 80)

        # Perform dealiasing
        dealias_data = pyart.correct.dealias_region_based(radar, gatefilter=gatefilter)
        radar.add_field("corrected_velocity", dealias_data)

        # Get velocity and spectrum width data
        velocity_data = radar.fields['corrected_velocity']['data'].copy()
        spectrum_width_data = radar.fields['spectrum_width']['data'].copy()

        # Mask invalid values
        valid_mask = np.isfinite(velocity_data) & np.isfinite(spectrum_width_data)

        # Multiply pixel values at corresponding locations
        product_data = np.abs(velocity_data) ** ((1/10) * spectrum_width_data)

        # Calculate the mean value of velspec_mask
        mean_velspec = np.mean(product_data[valid_mask])

        # Create anomalies based on the mean value
        velspec_anomalies = product_data - mean_velspec

        anomalies_dict = {
            "data": velspec_anomalies,
            "units": "m/s",
            "long_name": "VELSPEC Anomalies",
            "_FillValue": velspec_anomalies,
            "standard_name": "velspec_anomalies",
        }

        radar.add_field("velspec_anomalies", anomalies_dict, replace_existing=True)

        # Create a figure
        fig = plt.figure(figsize=[18, 6])

        # Display for reflectivity
        display1 = pyart.graph.RadarMapDisplay(radar)
        ax1 = plt.subplot(131, projection=ccrs.PlateCarree())
        ax1.set_facecolor('black')  # Set subplot background color to black
        display1.plot_ppi_map('reflectivity',
                              sweep=0,
                              vmin=-10,
                              vmax=70,
                              cmap="pyart_ChaseSpectral",
                              projection=ccrs.LambertConformal(
                                  central_latitude=radar.latitude["data"][0],
                                  central_longitude=radar.longitude["data"][0]),
                              ax=ax1,)
        
        ax1.add_feature(cfeature.COASTLINE, linestyle=':', edgecolor='white')  # Add coastline
        ax1.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='white')  # Add borders
        ax1.add_feature(cfeature.STATES, linestyle=':', edgecolor='white')  # Add states
        ax1.add_feature(cfeature.LAND, edgecolor='white')  # Add land
        
        plt.xlim(-88, -84)
        plt.ylim(28, 30)

        # Display for velocity
        display2 = pyart.graph.RadarMapDisplay(radar)
        ax2 = plt.subplot(133, projection=ccrs.PlateCarree())
        ax2.set_facecolor('black')  # Set subplot background color to black
        display2.plot_ppi_map('corrected_velocity',
                              sweep=1,
                              vmin=-50,
                              vmax=50,
                              cmap="pyart_NWSVel",
                              projection=ccrs.PlateCarree(),
                              ax=ax2,)
        
        ax2.add_feature(cfeature.COASTLINE, linestyle=':', edgecolor='white')  # Add coastline
        ax2.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='white')  # Add borders
        ax2.add_feature(cfeature.STATES, linestyle=':', edgecolor='white')  # Add states
        ax2.add_feature(cfeature.LAND, edgecolor='white')  # Add land
        
        plt.xlim(-88, -84)
        plt.ylim(28, 30)

        # Display for velspec_mask anomalies
        ax3 = plt.subplot(132, projection=ccrs.PlateCarree())
        ax3.set_facecolor('black')  # Set subplot background color to black
        display3 = pyart.graph.RadarMapDisplay(radar)
        display3.plot_ppi_map('velspec_anomalies',
                              sweep=1,
                              vmin=0,
                              vmax=5,
                              cmap="pyart_StepSeq25_r",
                              projection=ccrs.PlateCarree(),
                              ax=ax3,)
        
        ax3.add_feature(cfeature.COASTLINE, linestyle=':', edgecolor='white')  # Add coastline
        ax3.add_feature(cfeature.BORDERS, linestyle=':', edgecolor='white')  # Add borders
        ax3.add_feature(cfeature.STATES, linestyle=':', edgecolor='white')  # Add states
        ax3.add_feature(cfeature.LAND, edgecolor='white')  # Add land
        
        plt.xlim(-88, -84)
        plt.ylim(28, 30)

        # Save the figure as a PNG file in the "Downloads" folder
        png_filename = os.path.join(os.path.expanduser('~'), 'Downloads', f"{os.path.basename(file)[:-3]}.png")
        fig.savefig(png_filename, facecolor=fig.get_facecolor())  # Maintain black background
        plt.close()  # Close the figure to release resources