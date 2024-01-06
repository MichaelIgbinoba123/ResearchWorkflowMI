import tropycal.tracks as tracks
import datetime as dt
import matplotlib.pyplot as plt

basin = tracks.TrackDataset(basin='east_pacific',source='hurdat',include_btk=True)

storm = basin.get_storm(('dora',2023))

storm.plot(map_prop={'figsize':(18,9),'linewidth':1.0})

plt.show()