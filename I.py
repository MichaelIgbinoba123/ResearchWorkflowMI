import tropycal.tracks as tracks
import datetime as dt
import matplotlib.pyplot as plt

basin = tracks.TrackDataset(basin='north_atlantic',source='hurdat',include_btk=False)


basin.gridded_stats(request="maximum wind")

# Let's look at the average change in sustained wind speed over a 24-hour period. By default, the value plotted is for the midpoint of the 24-hour period (so 12 hours preceding and following). We'll use the "prop" keyword argument to set the colormap to "bwr" and set the contour level range:

basin.gridded_stats(request="average wind change in 24 hours",prop={'cmap':'bwr','clevs':[-80,80]})

plt.show()