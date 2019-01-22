"""
Analysis of the cloud fraction with height, collapsing the mcl dataset into a 
mask of cloud vs. no cloud f(t,z,y,x) then to a horizontal mean (a.k.a. cloud 
fraction) f(t,z)

Compare the time series of liquid water path to cloud fraction total, and at 
different heights.
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from analysis_tools import get_CTZ

hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]
# Open the netCDF for lwp (only one for the simulation)
lwp_nc = Dataset('lwp_00.nc', 'r')
lwp_key = u'STASH_m01s30i405'
# Read the data from the lwp netCDF
lwp_data = lwp_nc.variables[lwp_key][:]*1.
# Read the times for the lwp netCDF
times_lwp = lwp_nc.variables['min5_0'][:]*1.
# Read the latitude and longitude coordinates
lat = lwp_nc.variables[u'latitude_t'][:]*1.
lon = lwp_nc.variables[u'longitude_t'][:]*1.
lwp_nc.close()

mcl_key = u'STASH_m01s00i392'
for hour in hours:
    print hour
    # Open the netCDF for mr
    mr_nc = Dataset('mr_'+hour+'.nc', 'r')
    # Read the data for mcl from the mr netCDF
    mcl_data = mr_nc.variables[mcl_key][:]*1.
    # Read the height coordinate
    z = mr_nc.variables[u'thlev_zsea_theta'][:]*1.
    
    if hour == '00':
        # Mask out the mcl_data
        mcl_mask = np.where((mcl_data >= 1e-08), 1., 0.)
        times_mcl = mr_nc.variables[u'min10_0'][:]*1.
    else:
        mcl_mask = np.concatenate((mcl_mask, np.where((mcl_data >= 1e-08), 1., 0.)), axis = 0)
        times_mcl = np.concatenate((times_mcl, mr_nc.variables[u'min10_0'][:]*1.), axis = 0)
    mr_nc.close()

# Collapse mcl_mask with horizontal mean
mcl_fraction = np.nanmean(mcl_mask, axis = (2, 3)) # f(t,z,y,x) -> f(t,z)
lwp_fraction = np.nanmean(np.where((lwp_data >= 1e-08), 1., 0.), axis = (1, 2)) # f(t,y,x) -> f(t)

# heights
iz1 = np.where(np.abs(z - 750.0) == np.min(np.abs(z - 750.0)))[0][0]
iz2 = np.where(np.abs(z - 1500.0) == np.min(np.abs(z - 1500.0)))[0][0]
iz3 = np.where(np.abs(z - 2500.0) == np.min(np.abs(z - 2500.0)))[0][0]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(times_lwp, lwp_fraction, 'k', lw = 2, label = 'Total Cloud Cover')
ax.plot(times_mcl, mcl_fraction[:,iz1], 'r', label = 'Cloud Cover ('+str(round(z[iz1], 2))+'m)')
ax.plot(times_mcl, mcl_fraction[:,iz2], 'purple', label = 'Cloud Cover ('+str(round(z[iz2], 2))+'m)')
ax.plot(times_mcl, mcl_fraction[:,iz3], 'b', label = 'Cloud Cover ('+str(round(z[iz3], 2))+'m)')
ax.legend(loc = 2)

fig.tight_layout()
plt.savefig('cloud_cover_time_series.png', dpi = 100)
plt.show()
plt.close('all')

