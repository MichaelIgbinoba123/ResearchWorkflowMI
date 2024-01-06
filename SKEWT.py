import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from netCDF4 import Dataset, num2date, date2num
from datetime import datetime
from cpt_convert import loadCPT
from scipy.ndimage import gaussian_filter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm


#For Map
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker 
from matplotlib.patches import Rectangle
import matplotlib.colors as colors

#Stuff for Meteorological calculations 
import wrf
import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units
from wrf import (to_np, interplevel, geo_bounds, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, g_geoht, g_terrain, srhel, g_vorticity, g_dewpoint, tk, CoordPair, g_temp)

#Next thing: we need to know open up every file in the Dean Directory and loop through them!
import glob
import pandas as pd

#Color Table Stuff

def loadCPT(path):

    try:
        f = open(path)
    except:
        print ("File ", path, "not found")
        return None

    lines = f.readlines()

    f.close()

    x = np.array([])
    r = np.array([])
    g = np.array([])
    b = np.array([])

    colorModel = 'RGB'

    for l in lines:
        ls = l.split()
        if l[0] == '#':
            if ls[-1] == 'HSV':
                colorModel = 'HSV'
                continue
            else:
                continue
        if ls[0] == 'B' or ls[0] == 'F' or ls[0] == 'N':
            pass
        else:
            x=np.append(x,float(ls[0]))
            r=np.append(r,float(ls[1]))
            g=np.append(g,float(ls[2]))
            b=np.append(b,float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

        x=np.append(x,xtemp)
        r=np.append(r,rtemp)
        g=np.append(g,gtemp)
        b=np.append(b,btemp)

    if colorModel == 'HSV':
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
        r[i] = rr ; g[i] = gg ; b[i] = bb

    if colorModel == 'RGB':
        r = r/255.0
        g = g/255.0
        b = b/255.0

    xNorm = (x - x[0])/(x[-1] - x[0])

    red   = []
    blue  = []
    green = []

    for i in range(len(x)):
        red.append([xNorm[i],r[i],r[i]])
        green.append([xNorm[i],g[i],g[i]])
        blue.append([xNorm[i],b[i],b[i]])

    colorDict = {'red': red, 'green': green, 'blue': blue}

    return colorDict

#Function Below is for a storm following map projection. Credit: Dr. Sharanya Majumdar @The University of Miami
dataproj = ccrs.PlateCarree()

def create_map_background(tc_lon,tc_lat):
    
    crs = ccrs.PlateCarree()

    lat_min = 14
    lat_max = 19
    lon_min = 121
    lon_max = 133
     
    #Specify plot size
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, facecolor='gray', projection=crs)
    ax.coastlines(resolution='10m', alpha=0.9)
    
    #specify boundary extent from TC center. Current specifications give us a 10° x 10° projection
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
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAKES)
    return fig, ax

#the unadjusted domain has 9 more files than the adjusted run due to our lovely HPC
Parent_wrfout_unadjusted_d02 = Dataset(r'D:\wrfout_d02_2022-09-24_09_00_00')

#first, we are going to calculate the pressure center for our center following plots

tc_lon = []                                             #Empty List for TC center longitude
tc_lat = []                                             #Empty List for TC center latitude


#Start looking for pressure minimum in WRF out files

wrf_out_data = Parent_wrfout_unadjusted_d02

    #Next, isolate the variable we care to look at, in this case we will only care about the pressure pertubations
P_pertubation = getvar(wrf_out_data, "AFWA_MSLP", timeidx=35)
surface_P_pertubation = P_pertubation
    
    #Now, to find the location of the lowest surface pressure pertubation
minpressure = np.min(surface_P_pertubation[:,:])
minp = minpressure.values

    #Now we know what the lowest pressure is, we have to find the index value for where it occurs
da = surface_P_pertubation
p_index = np.argwhere(da.where(da == minp,0).values)
s_n = p_index[0][0]
w_e = p_index[0][1]

    #Append lat and longitude centers into a list for plotting
lon = wrf_out_data['XLONG'][0,s_n,w_e]
tc_lon.append(lon)
lat = wrf_out_data['XLAT'][0,s_n,w_e]
tc_lat.append(lat)

X_Y_Sounding_Center = []

#now we need to get the x_y points of each lat longitude coordinate

file = Parent_wrfout_unadjusted_d02
wrfin = file
x_ycenter = wrf.ll_to_xy(wrfin, tc_lat, tc_lon)
X_Y_Sounding_Center.append(x_ycenter)

#constant Needed for Calcualtions
sigma = 5.67 * 10**-8                           #Stefon-Boltzman Constant

#compare sounding location to IR brightness to make sure it is in the correct location
wrf_out_data = wrfin         #Open WRF_Out Data (I sorted glob glob my whole directory)
wrf_out_ctt = getvar(wrf_out_data, "ctt", timeidx=35)

sigma = 0.5  # You can adjust this value for more or less smoothing
IR_smoothed = gaussian_filter(wrf_out_ctt, sigma=sigma)
    
ncfile = file
Time = wrf.extract_times(ncfile, timeidx=35, method='cat', squeeze=True, cache=None, meta=False, do_xtime=False)
timestr = (str(Time))
    
titletime=(timestr[0:10]+' '+timestr[11:16])
filetime=(timestr[0:10]+'_'+timestr[11:13])

colortable = loadCPT('C:/Users/33017/Desktop/WRF Scripts/IR4AVHRR6.cpt')
ct = LinearSegmentedColormap('cpt', colortable)

    #Plot OLR Data
fig, ax = create_map_background(tc_lon,tc_lat) #We are only looking at 11 indices, and our list goes from 0-11, not 24 and up

plot = ax.pcolormesh(wrf_out_data['XLONG'][0,:,:], wrf_out_data['XLAT'][0,:,:], IR_smoothed, cmap=ct, norm=plt.Normalize(-100, 40))

cbar = plt.colorbar(plot, orientation = 'horizontal', pad = .05, shrink = .75, aspect = 30, extend = 'both')
cbar.ax.set_xlabel('WRF Cloud Top Temperature (°C)', fontsize = 12)
 
    #For Sounding Location
lon_i = tc_lon[0]
lat_i = tc_lat[0]

plt.plot(lon_i, lat_i, "x", color = 'black', label = "Center of Area Averaged Sounding")

# Calculate the center of the area-averaged sounding
center_lon = tc_lon[0]
center_lat = tc_lat[0]

# Calculate the new bottom left point based on the center of the sounding
new_bl_lon = center_lon - 2.5
new_bl_lat = center_lat - 2.5

# Calculate the top right point based on the new bottom left point
new_tr_lon = center_lon + 2.5
new_tr_lat = center_lat + 2.5

# Calculate the width and height of the square
width = new_tr_lon - new_bl_lon
height = new_tr_lat - new_bl_lat

# Add the square to the plot using the new bottom left point
ax.add_patch(Rectangle((new_bl_lon, new_bl_lat), width, height, edgecolor='black', facecolor='limegreen', fill=True, alpha=0.75, lw=2))

    #Add Coastlines
ax.coastlines('50m', linewidth=2)
ax.add_feature(cfeature.STATES, linewidth=2)
ax.add_feature(cfeature.BORDERS, linewidth=2)

plt.title(f'DATE: '+' '+titletime+' UTC', loc = 'right', fontsize = 12)
plt.title('TYPHOON Noru', loc = 'left', fontsize = 12)
plt.legend()

#Now that we have specified that the location of interest is indeed where we want it, 
#we can take an area average sounding!
from metpy.plots import Hodograph

sigma = 5.67 * 10**-8                                                    #Stefon-Boltzman Constant
    
wrfin = file
    
wrf_out_data = wrfin          #Open WRF_Out Data (I sorted glob glob my whole directory)
wrf_out_ctt = getvar(ncfile, "ctt", timeidx=35)
    
    #Get pertinent values necessary for skew-T
p1 = wrf.getvar(wrfin,"pressure",timeidx=35)               #Pressure
T1 = wrf.getvar(wrfin,"tc",timeidx=35)                     #Temperature
Td1 = wrf.getvar(wrfin,"td",timeidx=35)                    #Dewpoint
u = wrf.getvar(wrfin,"ua",timeidx=35)                     #Zonal Winds
v = wrf.getvar(wrfin,"va",timeidx=35)                     #Meridional Winds
u1 = u * 1.944
v1 = v * 1.944
pot1 = wrf.getvar(wrfin,"th",timeidx=35)                   #potential Temperature
height1 = wrf.getvar(wrfin,"height_agl",timeidx=35)        #Height above ground level for mass point
height = height = g_geoht.get_height(wrfin, timeidx=35, units='m')
terrain = g_terrain.get_terrain(ncfile, timeidx=35, units='m')

THREESRH = srhel(u, v, height, terrain, top=3000.0)

print(THREESRH.shape)
print(THREESRH.coords)
print(THREESRH.attrs)
    
    #Now that we are done with our textfile, we want to make a nice skew_T:
p_area =  (p1[:,int(new_bl_lon):int(new_tr_lon),
                 int(new_bl_lat):int(new_tr_lat)]   * units.hPa)        #Average Pressure in Hpa 
    
T_area =  (T1[:,int(new_bl_lon):int(new_tr_lon),
                 int(new_bl_lat):int(new_tr_lat)]   * units.degC)       #Average Temperature in Degrees Celsius
  
Td_area = (Td1[:,int(new_bl_lon):int(new_tr_lon),
                 int(new_bl_lat):int(new_tr_lat)]  * units.degC)       #Average Dewpoint Temperature in Degrees Celsius
   
u_area =  (u1[:,int(new_bl_lon):int(new_tr_lon),
                 int(new_bl_lat):int(new_tr_lat)]   * units('kt'))     #Average Zonal Winds in m/s                  

v_area =  (v1[:,int(new_bl_lon):int(new_tr_lon),
                 int(new_bl_lat):int(new_tr_lat)]   * units('kt'))     #Average Meridional winds in m/s

height_area =  (height[:,int(new_bl_lon):int(new_tr_lon),
                 int(new_bl_lat):int(new_tr_lat)]   * units('m'))     #Average Geopotential Height in meters

height1_area =  (height1[:,int(new_bl_lon):int(new_tr_lon),
                 int(new_bl_lat):int(new_tr_lat)]   * units('m'))     #Average Height AGL in meters

SRH_area = THREESRH[int(new_bl_lon):int(new_tr_lon), int(new_bl_lat):int(new_tr_lat)]
SRH_avg = SRH_area.mean()

print(SRH_avg)
    
    #to take the area average sounding, we need to average (mean) about the S_N and W_E arrays
p_avg_we = p_area.mean(dim = 'south_north')      #Average out area sounding across the S_n array
p_avg    = p_avg_we.mean(dim = 'west_east')      #Average our area + S_N average by the W_E so we have a 1d array with pressure levels and values

    #Same for Tempearture
T_avg_we = T_area.mean(dim = 'south_north')
T_avg    = T_avg_we.mean(dim = 'west_east')

    #Same For Dewpopint
Td_avg_we = Td_area.mean(dim = 'south_north')
Td_avg    = Td_avg_we.mean(dim = 'west_east')

    #Same for zonal wind
u_avg_we = u_area.mean(dim = 'south_north')
u_avg    = u_avg_we.mean(dim = 'west_east')

    #Same for meridional wind
v_avg_we = v_area.mean(dim = 'south_north')
v_avg    = v_avg_we.mean(dim = 'west_east')

    #Same for Geopotential Height
height_avg_we = height_area.mean(dim = 'south_north')
height_avg    = height_avg_we.mean(dim = 'west_east')

    #Same for Height AGL
height1_avg_we = height1_area.mean(dim = 'south_north')
height1_avg    = height1_avg_we.mean(dim = 'west_east')
    
    
ncfile = wrfin
Time = wrf.extract_times(ncfile, timeidx=35, method='cat', squeeze=True, cache=None, meta=False, do_xtime=False)
timestr = (str(Time))
    
titletime=(timestr[0:10]+' '+timestr[11:16])
filetime=(timestr[0:10]+'_'+timestr[11:13])

fig = plt.figure(figsize=(9, 9))
skew = SkewT(fig, rotation=45)

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot.
skew.plot(p_avg, T_avg, 'r')
skew.plot(p_avg, Td_avg, 'g')
skew.plot_barbs(p_avg[:-5].values, u_avg[:-5], v_avg[:-5])
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-60, 60)

# Set some better labels than the default
skew.ax.set_xlabel(f'Temperature (C)')
skew.ax.set_ylabel(f'Pressure (hPa)')

# Calculate LCL height and plot as black dot. Because `p`'s first value is
# ~1000 mb and its last value is ~250 mb, the `0` index is selected for
# `p`, `T`, and `Td` to lift the parcel from the surface. If `p` was inverted,
# i.e. start from low value, 250 mb, to a high value, 1000 mb, the `-1` index
# should be selected.
lcl_pressure, lcl_temperature = mpcalc.lcl(p_avg[0], T_avg[0], Td_avg[0])
skew.ax.axhline(lcl_pressure, color='#B8860B', linestyle='--', linewidth=2)

text_offset = 40  # You can adjust this value to control the vertical position of the text
skew.ax.text(-60, lcl_pressure.magnitude - text_offset, 'LCL', color='#B8860B')

#EL
el_pressure, el_temperature = mpcalc.el(p_avg, T_avg, Td_avg)
skew.ax.axhline(el_pressure, color='darkred', linestyle='--', linewidth=2)

text_offset = 10  # You can adjust this value to control the vertical position of the text
skew.ax.text(-120, el_pressure.magnitude - text_offset, 'EL', color='darkred')


# Calculate full parcel profile and add to plot as black line
prof = mpcalc.parcel_profile(p_avg, T_avg[0], Td_avg[0])
prof = prof.values - 273.15
skew.plot(p_avg, prof, 'k', linewidth=2)

# Shade areas of CAPE and CIN
skew.shade_cin(p_avg.values, T_avg.values, prof) 
skew.shade_cape(p_avg.values, T_avg.values, prof)

# An example of a slanted line at constant T -- in this case the 0
# isotherm
skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

# Add the relevant special lines
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()

# Create a hodograph
ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
h = Hodograph(ax_hod, component_range=50.)
h.add_grid(increment=5)
hodo = h.plot_colormapped(u_avg, v_avg, p_avg, cmap=plt.get_cmap('nipy_spectral'), norm=colors.Normalize(vmin=0, vmax=1000))

# Create a color bar on the left side of the plot
cbar_ax = fig.add_axes([0.765, 0.56, 0.005, 0.31])  # Adjust the position and size as needed
cbar = plt.colorbar(hodo, cax=cbar_ax)
cbar.ax.tick_params(axis='y', labelsize=8)

MUCAPE, MUCIN = mpcalc.most_unstable_cape_cin(p_avg, T_avg, Td_avg)
mixed_parcel_profile = mpcalc.mixed_parcel(p_avg, T_avg, Td_avg)
MLCAPE, MLCIN = mpcalc.mixed_layer_cape_cin(p_avg, T_avg, Td_avg)
SBCAPE, SBCIN = mpcalc.surface_based_cape_cin(p_avg, T_avg, Td_avg)
KIndex = mpcalc.k_index(p_avg, T_avg, Td_avg)
print("MUCAPE:", MUCAPE.item())
print("MLCAPE:", MLCAPE.item())
print("SBCAPE:", SBCAPE.item())
print("KIndex:", KIndex.item())



# Show the plot
plt.show()