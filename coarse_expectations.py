"""
Coarse control plots...
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from coarse_graining import coarse_grain as cg
from STASH_keys import *

# Read some data...
my_path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'
wind_nc = Dataset(my_path + 'wind_09.nc', 'r')

w_data = wind_nc.variables[n_key][:]
z = wind_nc.variables['thlev_zsea_theta'][:]

# Choose a time and height to do the coarse graining so that it doesn't take forever...
it = -1
target_height = 600.
iz = np.where(np.abs(z - target_height) == np.min(np.abs(z - target_height)))[0][0]

# Coarse grain the w_data
X, Y = np.meshgrid(np.arange(0., 116000., 100.)/1000., np.arange(0., 31900., 100.)/1000.)

w_data_800 = cg(input_data = w_data[it, iz, :, :], reduce_y = 8, reduce_x = 8, operation = np.nanmean, l_periodic = True)
X800, Y800 = np.meshgrid(np.linspace(0., 116000.1, w_data_800.shape[1])/1000., np.linspace(0., 31900.1, w_data_800.shape[0])/1000.)

w_data_1600 = cg(input_data = w_data[it, iz, :, :], reduce_y = 16, reduce_x = 16, operation = np.nanmean, l_periodic = True)
X16, Y16 = np.meshgrid(np.linspace(0., 116000.1, w_data_1600.shape[1])/1000., np.linspace(0., 31900.1, w_data_1600.shape[0])/1000.)

# Plot a side-by-side comparison
my_levels = [np.percentile(w_data[it,iz,:,:], pp) for pp in np.arange(10, 100, 10)]
my_cmap = 'bwr'
fig = plt.figure(tight_layout = True)
ax1 = fig.add_subplot(3, 1, 1, adjustable = 'box', aspect = 1)
ax1.contourf(X, Y, w_data[it,iz,:,:], cmap = my_cmap, levels = my_levels, extend = 'both')

ax2 = fig.add_subplot(3, 1, 2, adjustable = 'box', aspect = 1)
ax2.contourf(X800, Y800, w_data_800, cmap = my_cmap, levels = my_levels, extend = 'both')

ax3 = fig.add_subplot(3, 1, 3, adjustable = 'box', aspect = 1)
ax3.contourf(X16, Y16, w_data_1600, cmap = my_cmap, levels = my_levels, extend = 'both')

plt.show()


