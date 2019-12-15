import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Read the data
paths = {'BLm25' : '/nerc/n02/n02/xb899100/CloudTrail/RH_BLm25/',
         'FAm25' : '/nerc/n02/n02/xb899100/CloudTrail/RH_FAm25/',
         'Control_short' : '/nerc/n02/n02/xb899100/CloudTrail/Control_short/'}

my_data = {}
for path_key in paths.keys():
    my_data[path_key] = {}
    # Read the n-wind
    with Dataset(paths[path_key] + 'wind_' + ['09' if path_key == 'Control' else '04'][0] + '.nc', 'r') as wind_nc:
        my_data[path_key]['z'] = wind_nc.variables['thlev_zsea_theta'][:]*1.
        iz10 = np.where(np.abs(my_data[path_key]['z'] - 10.0) == np.min(np.abs(my_data[path_key]['z'] - 10.0)))[0][0]
        my_data[path_key][n_key] = wind_nc.variables[n_key][-1,iz10,:,:]*1.
    
    # Read and compute the potential temperature anomaly
    with Dataset(paths[path_key] + 'bouy_' + ['09' if path_key == 'Control' else '04'][0] + '.nc', 'r') as bouy_nc:
        iz2 = np.where(np.abs(my_data[path_key]['z'] - 2.0) == np.min(np.abs(my_data[path_key]['z'] - 2.0)))[0][0]
        my_data[path_key]['theta_p'] = bouy_nc.variables[theta_key][-1,iz2,:,:] - bouy_nc.variables[theta_key][-1,iz2,:,:].mean()
    
    # Read and compute the daytime cloud frequency
    with Dataset(paths[path_key] + 'lwp_00.nc', 'r') as lwp_nc:
        time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
        times = lwp_nc.variables[time_key][:]*1.
        if times.max() < 1400.:
            times += 120.0
        t_idx = [it for it in range(times.size) if (360.0 <= times[it])*(times[it] <= 1080.0)]
        my_data[path_key]['cloud_frequency'] = np.where(lwp_nc.variables[lwp_key][t_idx,:,:] > 0, 1.0, 0.0).mean(axis = 0)

# Make the x and y coordinates
X, Y = np.meshgrid(np.arange(my_data[path_key]['cloud_frequency'].shape[1])*0.1, np.arange(my_data[path_key]['cloud_frequency'].shape[0])*0.1)
island_radius = np.sqrt(50.0/np.pi)
x_c, y_c = [100.0+island_radius, 4*island_radius]
R = np.sqrt((X - x_c)**2 + (Y - y_c)**2)

# Make the figure
my_levels = np.array([level for level in np.arange(-3, 3.1, 0.5) if level != 0.0])
fig = plt.figure()

# Panel showing the control simulation
axa = fig.add_subplot(3, 1, 1, adjustable = 'box', aspect = 1)
im = axa.contourf(X, Y, my_data['Control_short'][n_key], levels = my_levels, cmap = 'bwr', extend = 'both')
im.cmap.set_under('navy')
im.cmap.set_over('firebrick')
axa.contourf(X, Y, my_data['Control_short']['cloud_frequency'], levels = [0.25, 1.001], colors = ['k'], alpha = 1./2.)
axa.contour(X, Y, R, levels = [island_radius], colors = ['darkred'], linewidths = [2])
axa.contour(X, Y, my_data['Control_short']['theta_p'], levels = [0.1], colors = ['k'])
axa.set_title('a) Control (short)')
axa.set_xticklabels([''])
axa.set_ylabel('y (km)')

axb = fig.add_subplot(3, 1, 2, adjustable = 'box', aspect = 1)
im = axb.contourf(X, Y, my_data['BLm25'][n_key], levels = my_levels, cmap = 'bwr', extend = 'both')
im.cmap.set_under('navy')
im.cmap.set_over('firebrick')
axb.contourf(X, Y, my_data['BLm25']['cloud_frequency'], levels = [0.25, 1.001], colors = ['k'], alpha = 1./2.)
axb.contour(X, Y, R, levels = [island_radius], colors = ['darkred'], linewidths = [2])
axb.contour(X, Y, my_data['BLm25']['theta_p'], levels = [0.1], colors = ['k'])
axb.set_title('b) BLm25')
axb.set_xticklabels([''])
axb.set_ylabel('y (km)')

axc = fig.add_subplot(3, 1, 3, adjustable = 'box', aspect = 1)
im = axc.contourf(X, Y, my_data['FAm25'][n_key], levels = my_levels, cmap = 'bwr', extend = 'both')
im.cmap.set_under('navy')
im.cmap.set_over('firebrick')
axc.contourf(X, Y, my_data['FAm25']['cloud_frequency'], levels = [0.25, 1.001], colors = ['k'], alpha = 1./2.)
axc.contour(X, Y, R, levels = [island_radius], colors = ['darkred'], linewidths = [2])
axc.contour(X, Y, my_data['FAm25']['theta_p'], levels = [0.1], colors = ['k'])
axc.set_title('c) FAm25')
axc.set_xlabel('x (km)')
axc.set_ylabel('y (km)')

cbar_ax = inset_axes(axc, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (0, -0.525, 1.0, 2.0), bbox_transform = axc.transAxes, borderpad = 0.0)
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = '$n^{\prime}$ (m s$^{-1}$)')
cbar.set_ticks([-3, -2, -1, 0.0, 1, 2, 3])
plt.subplots_adjust(left = 0.10, bottom = 0.2, right = 1.0, wspace = 0.0, hspace = 0.3)
plt.savefig('../Ch5_Figure17.png', dpi = 250, bbox_inches = 'tight')
plt.show()


