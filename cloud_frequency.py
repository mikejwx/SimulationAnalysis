import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

lsm  = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
my_lsm = lsm.variables['lsm'][0,0,:,:]*1.
lsm.close()

lwp_key = u'STASH_m01s30i405'
control_lwp_nc = Dataset('../lwp_00.nc', 'r')
exp01_lwp_nc   = Dataset('/work/n02/n02/xb899100/cylc-run/u-bg023/share/data/history/lwp_00.nc', 'r')
exp02_lwp_nc   = Dataset('/work/n02/n02/xb899100/cylc-run/u-bg113/share/data/history/lwp_00.nc', 'r')

# Make some dimensions
x = np.arange(0., 116000., 100.)/1000.
y = np.arange(0., 31900., 100.)/1000.
X, Y = np.meshgrid(x, y)

################################################################################
#
# Section 1: Cloud Fraction, i.e. Mean cloud mask
#
################################################################################

control_cf = np.nanmean(np.where((control_lwp_nc.variables[lwp_key][:,:,:] > 5e-03), 1., 0.), axis = 0)
exp01_cf = np.nanmean(np.where((exp01_lwp_nc.variables[lwp_key][:,:,:] > 5e-03), 1., 0.), axis = 0)
exp02_cf = np.nanmean(np.where((exp02_lwp_nc.variables[lwp_key][:,:,:] > 5e-03), 1., 0.), axis = 0)

cf_levels = np.arange(0., 0.26, 0.025)
cf_cmap = 'viridis'

# Start figure
fig = plt.figure(figsize = (22, 7), tight_layout = True)
ax0 = fig.add_subplot(3, 3, 1, adjustable = 'box', aspect = 1)
cf = ax0.contourf(X, Y, control_cf, cmap = cf_cmap, levels = cf_levels)
axins = inset_axes(ax0, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax0.transAxes, borderpad = 0.)
ax0.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
plt.colorbar(cf, cax = axins, label = 'Cloud Fraction')
ax0.set_title('Cloud Frequency (Control, H = 250 W m$^{-2}$, E = 250 W m$^{-2}$)')
ax0.set_xlabel('x (km)')
ax0.set_ylabel('y (km)')

ax1 = fig.add_subplot(3, 3, 2, adjustable = 'box', aspect = 1)
ax1.contourf(X, Y, exp01_cf, cmap = cf_cmap, levels = cf_levels)
ax1.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
ax1.set_title('Cloud Frequency (exp01, H = 125 W m$^{-2}$, E = 375 W m$^{-2}$)')
ax1.set_xlabel('x (km)')
ax1.set_ylabel('y (km)')

ax2 = fig.add_subplot(3, 3, 3, adjustable = 'box', aspect = 1)
ax2.contourf(X, Y, exp02_cf, cmap = cf_cmap, levels = cf_levels)
ax2.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
ax2.set_title('Cloud Frequency (exp02, H = 375 W m$^{-2}$, E = 125 W m$^{-2}$)')
ax2.set_xlabel('x (km)')
ax2.set_ylabel('y (km)')

################################################################################
#
# Section 2: Maximum liquid water path
#
################################################################################

control_lwp_max = np.nanmax(control_lwp_nc.variables[lwp_key][:,:,:], axis = 0)*1000.
exp01_lwp_max = np.nanmax(exp01_lwp_nc.variables[lwp_key][:,:,:], axis = 0)*1000.
exp02_lwp_max = np.nanmax(exp02_lwp_nc.variables[lwp_key][:,:,:], axis = 0)*1000.

lwp_max_levels = [50., 100., 200., 500., 1000., 2000., 5000.]
n_lwp_max_levels = float(len(lwp_max_levels))
lwp_max_cmap = mpl.cm.get_cmap('binary')
lwp_max_cmap = lwp_max_cmap(np.arange(n_lwp_max_levels)/n_lwp_max_levels)

ax3 = fig.add_subplot(3, 3, 4, adjustable = 'box', aspect = 1)
lwp_max = ax3.contourf(X, Y, control_lwp_max, colors = lwp_max_cmap, levels = lwp_max_levels, extend = 'max')
ax3.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
axins = inset_axes(ax3, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax3.transAxes, borderpad = 0.)
plt.colorbar(lwp_max, cax = axins, label = r'LWP (g kg$^{-1}$)')
ax3.set_title('Maximum LWP (Control, H = 250 W m$^{-2}$), E = 250 W m$^{-2}$')
ax3.set_xlabel('x (km)')
ax3.set_ylabel('y (km)')

ax4 = fig.add_subplot(3, 3, 5, adjustable = 'box', aspect = 1)
ax4.contourf(X, Y, exp01_lwp_max, colors = lwp_max_cmap, levels = lwp_max_levels, extend = 'max')
ax4.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
ax4.set_title('Maximum LWP (exp01, H = 125 W m$^{-2}$), E = 375 W m$^{-2}$')
ax4.set_xlabel('x (km)')
ax4.set_ylabel('y (km)')

ax5 = fig.add_subplot(3, 3, 6, adjustable = 'box', aspect = 1)
ax5.contourf(X, Y, exp02_lwp_max, colors = lwp_max_cmap, levels = lwp_max_levels, extend = 'max')
ax5.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
ax5.set_title('Maximum LWP (exp02, H = 375 W m$^{-2}$, E = 125 W m$^{-2}$)')
ax5.set_xlabel('x (km)')
ax5.set_ylabel('y (km)')

################################################################################
#
# Section 3: Mean liquid water path
#
################################################################################

control_lwp_mean = np.nanmean(control_lwp_nc.variables[lwp_key][1:,:,:], axis = 0)*1000.
exp01_lwp_mean = np.nanmean(exp01_lwp_nc.variables[lwp_key][1:,:,:], axis = 0)*1000.
exp02_lwp_mean = np.nanmean(exp02_lwp_nc.variables[lwp_key][1:,:,:], axis = 0)*1000.

lwp_mean_levels = [1., 2., 5., 10., 20., 50., 100.]
n_lwp_mean_levels = float(len(lwp_mean_levels))
lwp_mean_cmap = mpl.cm.get_cmap('binary')
lwp_mean_cmap = lwp_mean_cmap(np.arange(n_lwp_mean_levels)/n_lwp_mean_levels)

ax6 = fig.add_subplot(3, 3, 7, adjustable = 'box', aspect = 1)
lwp_mean = ax6.contourf(X, Y, control_lwp_mean, colors = lwp_mean_cmap, levels = lwp_mean_levels, extend = 'max')
ax6.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
axins = inset_axes(ax6, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax6.transAxes, borderpad = 0.)
plt.colorbar(lwp_mean, cax = axins, label = r'LWP (g kg$^{-1}$)')
ax6.set_title('Mean LWP (Control, H = 250 W m$^{-2}$, E = 250 W m$^{-2}$')
ax6.set_xlabel('x (km)')
ax6.set_ylabel('y (km)')

ax7 = fig.add_subplot(3, 3, 8, adjustable = 'box', aspect = 1)
ax7.contourf(X, Y, exp01_lwp_mean, colors = lwp_mean_cmap, levels = lwp_mean_levels, extend = 'max')
ax7.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
ax7.set_title('Mean LWP (exp01, H = 125 W m$^{-2}$, E = 375 W m$^{-2}$)')
ax7.set_xlabel('x (km)')
ax7.set_ylabel('y (km)')

ax8 = fig.add_subplot(3, 3, 9, adjustable = 'box', aspect = 1)
ax8.contourf(X, Y, exp02_lwp_mean, colors = lwp_mean_cmap, levels = lwp_mean_levels, extend = 'max')
ax8.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
ax8.set_title('Mean LWP (exp02, H = 375 W m$^{-2}$, E = 125 W m$^{-2}$)')
ax8.set_xlabel('x (km)')
ax8.set_ylabel('y (km)')

plt.savefig('../Cloud_Properties.png', dpi = 100)
plt.show()
plt.close('all')

# Finally, close all the netCDF
control_lwp_nc.close()
exp01_lwp_nc.close()
exp02_lwp_nc.close()
