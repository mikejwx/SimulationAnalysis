"""
Coarse control plots...
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from coarse_graining import coarse_grain as cg
from coarse_graining import fix_x
from STASH_keys import *

# Read some data...
my_path = '/nerc/n02/n02/xb899100/CloudTrail/'

# My experiments
exp0100 = 'Control/'
exp0800 = 'Control_0800m/'

# side-by-side comparison
lwp_0100_nc = Dataset(my_path + exp0100 + 'lwp_00.nc', 'r')
lwp_0800_nc = Dataset(my_path + exp0800 + 'lwp_00.nc', 'r')

lwp_0100_data = lwp_0100_nc.variables[lwp_key][:]*1000.
lwp_0800_data = lwp_0800_nc.variables[lwp_key][:]*1000.

# Choose a time to do the coarse graining so that it doesn't take forever...
it1 = lwp_0100_data.shape[0]/2
it8 = lwp_0800_data.shape[0]/2

# Coarse grain the w_data
X1, Y1 = np.meshgrid(np.arange(1160)*0.1, np.arange(319)*0.1)
X8, Y8 = np.meshgrid(np.arange(146)*0.8, np.arange(40)*0.8)

lwp_cg_data = fix_x(cg(lwp_0100_data[it1,:,:], reduce_y = 8, reduce_x = 8), axis = 1)

# Plot the data
my_cmap = 'viridis'
my_levels = np.arange(11)*100
my_levels[0] = 1
fig = plt.figure(tight_layout = True)
axa = fig.add_subplot(3, 1, 1, adjustable = 'box', aspect = 1)
axa.set_ylabel('y (km)')
axa.set_title('Control, dx = 100 m')

axb = fig.add_subplot(3, 1, 2, adjustable = 'box', aspect = 1)
axb.set_ylabel('y (km)')
axb.set_title('Control, dx = 100 m coarse-grained to 800 m')

axc = fig.add_subplot(3, 1, 3, adjustable = 'box', aspect = 1)
axc.set_ylabel('y (km)')
axc.set_xlabel('x (km)')
axc.set_title('Control, dx = 800 m (explicit)')

lwp100 = axa.contourf(X1, Y1, lwp_0100_data[it1,:,:], cmap = my_cmap, levels = my_levels, extend = 'max', vmax = 2000)
fig.colorbar(lwp100, ax = axa, label = r'LWP g m$^{-2}$')

# read and plot a 100 m lsm
lsm100m_data = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
axa.contour(X1, Y1, lsm100m_data.variables['lsm'][0,0,:,:], colors = ['darkred'], levels = [0], linewidths = [2])
lsm100m_data.close()

lwp100_cg = axb.contourf(X8, Y8, lwp_cg_data, cmap = my_cmap, levels = my_levels, extend = 'max', vmax = 2000)
fig.colorbar(lwp100_cg, ax = axb, label = r'LWP g m$^{-2}$')

# read and plot a 800 m lsm
lsm800m_data = Dataset('/work/n02/n02/xb899100/island_masks/lsm50_0800m.nc', 'r')
axb.contour(X8, Y8, lsm800m_data.variables['lsm'][0,0,:,:], colors = ['darkred'], levels = [0], linewidths = [2])

lwp800 = axc.contourf(X8, Y8, lwp_0800_data[it8,:,:], cmap = my_cmap, levels = my_levels, extend = 'max', vmax = 2000)
fig.colorbar(lwp800, ax = axc, label = r'LWP g m$^{-2}$')
axc.contour(X8, Y8, lsm800m_data.variables['lsm'][0,0,:,:], colors = ['darkred'], levels = [0], linewidths = [2])

lsm800m_data.close()

plt.savefig('../Resolution_Comparison.png', dpi = 150)
plt.show()

lwp_0100_nc.close()
lwp_0800_nc.close()

