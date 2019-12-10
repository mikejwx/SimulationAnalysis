import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

lsm  = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
my_lsm = lsm.variables['lsm'][0,0,:,:]*1.
lsm.close()

lwp_key = u'STASH_m01s30i405'
control_lwp_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control/lwp_00.nc', 'r')
time_key = 'min5_0'
control_times = control_lwp_nc.variables[time_key][:]*1.
exp_rsu_lwp_nc = Dataset('/work/n02/n02/xb899100/cylc-run/u-bg387/share/data/history/lwp.nc', 'r')
exp_rsu_times = exp_rsu_lwp_nc.variables[time_key][:]*1.

# Make some dimensions
x = np.arange(0., 116000., 100.)/1000.
y = np.arange(0., 31900., 100.)/1000.
X, Y = np.meshgrid(x, y)

################################################################################
#
# Section 1: Cloud Fraction, i.e. Mean cloud mask over the heating period
#
################################################################################

control_indexes = [it for it in xrange(len(control_times)) if 360. <= control_times[it] <= 1080.]
control_cf = np.nanmean(np.where((control_lwp_nc.variables[lwp_key][control_indexes,:,:] > 1e-16), 1., 0.), axis = 0)
exp_rsu_indexes = [it for it in xrange(len(exp_rsu_times)) if 360. <= exp_rsu_times[it] <= 1080.]
exp_rsu_cf = np.nanmean(np.where((exp_rsu_lwp_nc.variables[lwp_key][exp_rsu_indexes,:,:] > 1e-16), 1., 0.), axis = 0)

cf_levels = np.arange(0., 0.51, 0.05)
cf_cmap = 'viridis'

# Start figure
fig = plt.figure(figsize = (22, 7))
ax0 = fig.add_subplot(3, 3, 1, adjustable = 'box', aspect = 1)
cf = ax0.contourf(X, Y, control_cf, cmap = cf_cmap, levels = cf_levels)
axins = inset_axes(ax0, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax0.transAxes, borderpad = 0.)
ax0.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
plt.colorbar(cf, cax = axins, label = 'Cloud Fraction')
ax0.set_title('Cloud Frequency (Control, 6 hours spin up)')
ax0.set_xlabel('x (km)')
ax0.set_ylabel('y (km)')

ax1 = fig.add_subplot(3, 3, 2, adjustable = 'box', aspect = 1)
ax1.contourf(X, Y, exp_rsu_cf, cmap = cf_cmap, levels = cf_levels)
ax1.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
ax1.set_title('Cloud Frequency (exp_rsu, 2 hours spin up)')
ax1.set_xlabel('x (km)')
ax1.set_ylabel('y (km)')

ax2 = fig.add_subplot(3, 3, 3, adjustable = 'box', aspect = 1)
cfd = ax2.contourf(X, Y, exp_rsu_cf - control_cf, cmap = 'bwr', levels = [i for i in np.arange(-0.30, 0.31, 0.06) if i != 0.], extend = 'both')
axins = inset_axes(ax2, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax2.transAxes, borderpad = 0.)
plt.colorbar(cfd, cax = axins, label = 'Cloud Fraction Difference')
ax2.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
ax2.set_title('Difference')
ax2.set_xlabel('x (km)')
ax2.set_ylabel('y (km)')

################################################################################
#
# Section 2: Maximum liquid water path
#
################################################################################

control_lwp_max = np.nanmax(control_lwp_nc.variables[lwp_key][control_indexes,:,:], axis = 0)*1000.
exp_rsu_lwp_max = np.nanmax(exp_rsu_lwp_nc.variables[lwp_key][exp_rsu_indexes,:,:], axis = 0)*1000.

lwp_max_levels = [50., 100., 200., 500., 1000., 2000., 5000.]
n_lwp_max_levels = float(len(lwp_max_levels))
lwp_max_cmap = mpl.cm.get_cmap('binary')
lwp_max_cmap = lwp_max_cmap(np.arange(n_lwp_max_levels)/n_lwp_max_levels)

ax3 = fig.add_subplot(3, 3, 4, adjustable = 'box', aspect = 1)
lwp_max = ax3.contourf(X, Y, control_lwp_max, colors = lwp_max_cmap, levels = lwp_max_levels, extend = 'max')
ax3.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
axins = inset_axes(ax3, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax3.transAxes, borderpad = 0.)
plt.colorbar(lwp_max, cax = axins, label = r'LWP (g kg$^{-1}$)')
ax3.set_title('Maximum LWP (Control, 6 hours spin up)')
ax3.set_xlabel('x (km)')
ax3.set_ylabel('y (km)')

ax4 = fig.add_subplot(3, 3, 5, adjustable = 'box', aspect = 1)
ax4.contourf(X, Y, exp_rsu_lwp_max, colors = lwp_max_cmap, levels = lwp_max_levels, extend = 'max')
ax4.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
ax4.set_title('Maximum LWP (exp_rsu, 2 hours spin up)')
ax4.set_xlabel('x (km)')
ax4.set_ylabel('y (km)')

ax5 = fig.add_subplot(3, 3, 6, adjustable = 'box', aspect = 1)
lwp_max_d = ax5.contourf(X, Y, exp_rsu_lwp_max - control_lwp_max, cmap = 'bwr', levels = [i for i in np.arange(-2000., 2000.1, 200) if i != 0], extend = 'both')
ax5.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
axins = inset_axes(ax5, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax5.transAxes, borderpad = 0.)
plt.colorbar(lwp_max_d, cax = axins, label = 'LWP difference (g kg $^{-1}$)')
ax5.set_title('difference in maximum LWP')
ax5.set_xlabel('x (km)')
ax5.set_ylabel('y (km)')

################################################################################
#
# Section 3: Mean liquid water path
#
################################################################################

control_lwp_mean = np.nanmean(control_lwp_nc.variables[lwp_key][control_indexes,:,:], axis = 0)*1000.
exp_rsu_lwp_mean = np.nanmean(exp_rsu_lwp_nc.variables[lwp_key][exp_rsu_indexes,:,:], axis = 0)*1000.

lwp_mean_levels = [1., 2., 5., 10., 20., 50., 100.]
n_lwp_mean_levels = float(len(lwp_mean_levels))
lwp_mean_cmap = mpl.cm.get_cmap('binary')
lwp_mean_cmap = lwp_mean_cmap(np.arange(n_lwp_mean_levels)/n_lwp_mean_levels)

ax6 = fig.add_subplot(3, 3, 7, adjustable = 'box', aspect = 1)
lwp_mean = ax6.contourf(X, Y, control_lwp_mean, colors = lwp_mean_cmap, levels = lwp_mean_levels, extend = 'max')
ax6.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
axins = inset_axes(ax6, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax6.transAxes, borderpad = 0.)
plt.colorbar(lwp_mean, cax = axins, label = r'LWP (g kg$^{-1}$)')
ax6.set_title('Mean LWP (Control, 6 hours spin up)')
ax6.set_xlabel('x (km)')
ax6.set_ylabel('y (km)')

ax7 = fig.add_subplot(3, 3, 8, adjustable = 'box', aspect = 1)
ax7.contourf(X, Y, exp_rsu_lwp_mean, colors = lwp_mean_cmap, levels = lwp_mean_levels, extend = 'max')
ax7.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
ax7.set_title('Mean LWP (exp_rsu, 2 hours spin up)')
ax7.set_xlabel('x (km)')
ax7.set_ylabel('y (km)')

ax8 = fig.add_subplot(3, 3, 9, adjustable = 'box', aspect = 1)
lwp_mean_d = ax8.contourf(X, Y, exp_rsu_lwp_mean - control_lwp_mean, cmap = 'bwr', levels = [i for i in np.arange(-100., 100.1, 10.) if i != 0], extend = 'both')
ax8.contour(X, Y, my_lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
axins = inset_axes(ax8, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax8.transAxes, borderpad = 0.)
plt.colorbar(lwp_mean_d, cax = axins, label = r'difference LWP (g kg$^{-1}$)')
ax8.set_title('Mean LWP difference')
ax8.set_xlabel('x (km)')
ax8.set_ylabel('y (km)')

plt.savefig('../Comparison_of_Cloud_Properties.png', dpi = 150)
plt.show()
plt.close('all')

# Finally, close all the netCDF
control_lwp_nc.close()
exp_rsu_lwp_nc.close()
