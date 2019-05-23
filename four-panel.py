"""
A four-panel plot originally designed to identify the major regime shift that
occurs for the U05 case.

Panels include: across-flow winds, equivalent potential temperature,
some measure of precipitation, some measure of cloud.
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import n_key, theta_key, temp_key, pthe_key, lwp_key, ls_rain_amt_key
from SkewT_archer import getThetaE

# Read in the wind, bouy, and lwp nc files
exp = 'Control_short'
my_path   = '/nerc/n02/n02/xb899100/CloudTrail/' + exp + '/'
wind_nc   = Dataset(my_path + 'wind_08.nc', 'r')
bouy_nc   = Dataset(my_path + 'bouy_08.nc', 'r')
fluxes_nc = Dataset(my_path + 'fluxes_08.nc', 'r')
lwp_nc    = Dataset(my_path + 'lwp_00.nc', 'r')

# Let's look at one time, e.g. peak heating
target_time = 720. - 240. # subtract 240. because this is a short simulation
wind_time_key = [tkey for tkey in wind_nc.variables.keys() if 'min' in tkey][0]
bouy_time_key = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
rain_time_key = [tkey for tkey in lwp_nc.variables.keys() if 'accum' in tkey][0]
lwp_time_key  = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]

wind_times = wind_nc.variables[wind_time_key][:]
bouy_times = bouy_nc.variables[bouy_time_key][:]
rain_times = lwp_nc.variables[rain_time_key][:]
lwp_times  = lwp_nc.variables[lwp_time_key][:]

windex_t = np.where((np.abs(wind_times - target_time) == np.min(np.abs(wind_times - target_time))))[0][0]
bindex_t = np.where((np.abs(bouy_times - target_time) == np.min(np.abs(bouy_times - target_time))))[0][0]
rindex_t = np.where((np.abs(rain_times - target_time) == np.min(np.abs(rain_times - target_time))))[0][0]
lindex_t = np.where((np.abs(lwp_times - target_time) == np.min(np.abs(lwp_times - target_time))))[0][0]

# Choose some vertical levels to do the plots:
z = wind_nc.variables['thlev_zsea_theta'][:]
iz10 = np.where(np.abs(z - 10.0) == np.min(np.abs(z - 10.0)))[0][0]

# Create coordinate system
X, Y = np.meshgrid(np.arange(0., 116000., 100.)/1000., np.arange(0., 31900., 100.)/1000.)

# Compute equivalent potential temperature
theta_e = getThetaE(bouy_nc.variables[theta_key][bindex_t,0,:,:], bouy_nc.variables[temp_key][bindex_t,0,:,:], fluxes_nc.variables[pthe_key][bindex_t,0,:,:], t_units = 'K', p_units = 'Pa')

# Read land-sea mask
lsm  = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
my_lsm = lsm.variables['lsm'][0,0,:,:]
lsm.close()

# Create the plot
fig = plt.figure(tight_layout = True, figsize = (14, 9))
axa = fig.add_subplot(2, 2, 1, adjustable = 'box', aspect = 1)
axa.contour(X, Y, my_lsm, colors = ['k'], linewidths = [2])
axa.set_xticklabels([])
axa.set_ylabel('y (km)')

axb = fig.add_subplot(2, 2, 2, adjustable = 'box', aspect = 1)
axb.contour(X, Y, my_lsm, colors = ['k'], linewidths = [2])
axb.set_yticklabels([])
axb.set_xticklabels([])

axc = fig.add_subplot(2, 2, 3, adjustable = 'box', aspect = 1)
axc.contour(X, Y, my_lsm, colors = ['r'], linewidths = [2])
axc.set_ylabel('y (km)')
axc.set_xlabel('x (km)')

axd = fig.add_subplot(2, 2, 4, adjustable = 'box', aspect = 1)
axd.contour(X, Y, my_lsm, colors = ['r'], linewidths = [2])
axd.set_yticklabels([])
axd.set_xlabel('x (km)')

wnd = axa.contourf(X, Y, wind_nc.variables[n_key][windex_t, iz10, :, :], cmap = 'bwr', levels = [-4., -3, -2., -1, 1, 2., 3, 4.], extend = 'both')
fig.colorbar(wnd, ax = axa, orientation = 'horizontal', label = 'v$^{\prime}$ (m/s)')

the = axb.contourf(X, Y, theta_e, cmap = 'bwr', levels = [int(np.nanmean(theta_e)) + bloop for bloop in [-4., -3, -2., -1, 1, 2., 3, 4.]], extend = 'both')
fig.colorbar(the, ax = axb,orientation = 'horizontal', label = '$\\theta_{e}$ (K)')

pre = axc.contourf(X, Y, lwp_nc.variables[ls_rain_amt_key][rindex_t,:,:], cmap = 'viridis', levels = [1, 2, 4, 8, 16, 32, 64], extend = 'max')
fig.colorbar(pre, ax = axc, orientation = 'horizontal', label = 'Total Accumulated Precipitation (mm)')

lwp = axd.contourf(X, Y, lwp_nc.variables[lwp_key][lindex_t,:,:]*1000., cmap = 'viridis', levels = [1, 50, 100, 200, 400, 800, 1600, 3200], extend = 'max')
fig.colorbar(lwp, ax = axd, orientation = 'horizontal', label = ' Liquid Water Path (g/m$^{2}$)')

plt.savefig('../four-panel_' + exp + '.png', dpi = 150)
plt.show()


