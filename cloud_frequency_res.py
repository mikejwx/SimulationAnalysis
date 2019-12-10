import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from STASH_keys import lwp_key

lsm0100  = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
my_lsm0100 = lsm0100.variables['lsm'][0,0,:,:]*1.
lsm0100.close()

lsm0800  = Dataset('/work/n02/n02/xb899100/island_masks/lsm50_0800m.nc', 'r')
my_lsm0800 = lsm0800.variables['lsm'][0,0,:,:]*1.
lsm0800.close()

lsm1600  = Dataset('/work/n02/n02/xb899100/island_masks/lsm50_1600m.nc', 'r')
my_lsm1600 = lsm1600.variables['lsm'][0,0,:,:]*1.
lsm1600.close()

base_path = '/nerc/n02/n02/xb899100/CloudTrail/'
lwp_0100m_nc = Dataset(base_path + 'Control/lwp_00.nc', 'r')
lwp_0800m_nc   = Dataset(base_path + 'Control_0800m/lwp_00.nc', 'r')
lwp_1600m_nc   = Dataset(base_path + 'Control_1600m/lwp_00.nc', 'r')

# Make some dimensions
X01, Y01 = np.meshgrid(np.arange(1160)*0.1, np.arange(319)*0.1)
X08, Y08 = np.meshgrid(np.arange(146)*0.8, np.arange(40)*0.8)
X16, Y16 = np.meshgrid(np.arange(74)*1.6, np.arange(20)*1.6)

################################################################################
#
# Section 1: Cloud Fraction, i.e. Mean cloud mask
#
################################################################################

# get the indexes for 0100m sim that are between 6AM and 12PM
time_key0100 = [tkey for tkey in lwp_0100m_nc.variables.keys() if 'min' in tkey][0]
its0100 = [it for it in range(lwp_0100m_nc.variables[time_key0100].shape[0]) if 360.0 <= lwp_0100m_nc.variables[time_key0100][it] <= 720.0]

# get the indexes for the 0800m sim that are between 6AM and 12PM
time_key0800 = [tkey for tkey in lwp_0800m_nc.variables.keys() if 'min' in tkey][0]
its0800 = [it for it in range(lwp_0800m_nc.variables[time_key0800].shape[0]) if 120.0 <= lwp_0800m_nc.variables[time_key0800][it] <= 540.0]

# get the indexes for the 1600m sim that are between 6AM and 12PM
time_key1600 = [tkey for tkey in lwp_1600m_nc.variables.keys() if 'min' in tkey][0]
its1600 = [it for it in range(lwp_1600m_nc.variables[time_key1600].shape[0]) if 120.0 <= lwp_1600m_nc.variables[time_key1600][it] <= 540.0]

cf_0100m = np.nanmean(np.where((lwp_0100m_nc.variables[lwp_key][its0100,:,:] > 5e-03), 1., 0.), axis = 0)
cf_0800m = np.nanmean(np.where((lwp_0800m_nc.variables[lwp_key][its0800,:,:] > 5e-03), 1., 0.), axis = 0)
cf_1600m = np.nanmean(np.where((lwp_1600m_nc.variables[lwp_key][its1600,:,:] > 5e-03), 1., 0.), axis = 0)

cf_levels = np.arange(0., 0.26, 0.025)
cf_cmap = 'viridis'

# Start figure
fig = plt.figure(figsize = (12, 9), tight_layout = True)
ax0 = fig.add_subplot(3, 1, 1, adjustable = 'box', aspect = 1)
cf = ax0.contourf(X01, Y01, cf_0100m, cmap = cf_cmap, levels = cf_levels)
axins = inset_axes(ax0, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax0.transAxes, borderpad = 0.)
ax0.contour(X01, Y01, my_lsm0100, levels = [0.0], colors = ['k'], linewidths = [2])
plt.colorbar(cf, cax = axins, label = 'Cloud Fraction')
ax0.set_title('Cloud Frequency (dx = 0100m)')
#ax0.set_xlabel('x (km)')
ax0.set_ylabel('y (km)')

ax1 = fig.add_subplot(3, 1, 2, adjustable = 'box', aspect = 1)
cf = ax1.contourf(X08, Y08, cf_0800m, cmap = cf_cmap, levels = cf_levels)
axins = inset_axes(ax1, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax1.transAxes, borderpad = 0.)
plt.colorbar(cf, cax = axins, label = 'Cloud Fraction')
ax1.contour(X08, Y08, my_lsm0800, levels = [0.0], colors = ['k'], linewidths = [2])
ax1.set_title('Cloud Frequency (dx = 0800m)')
#ax1.set_xlabel('x (km)')
ax1.set_ylabel('y (km)')

ax2 = fig.add_subplot(3, 1, 3, adjustable = 'box', aspect = 1)
cf = ax2.contourf(X16, Y16, cf_1600m, cmap = cf_cmap, levels = cf_levels)
axins = inset_axes(ax2, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax2.transAxes, borderpad = 0.)
plt.colorbar(cf, cax = axins, label = 'Cloud Fraction')
ax2.contour(X16, Y16, my_lsm1600, levels = [0.0], colors = ['k'], linewidths = [2])
ax2.set_title('Cloud Frequency (dx = 1600m)')
ax2.set_xlabel('x (km)')
ax2.set_ylabel('y (km)')

plt.savefig('../output_from_param_analysis/cld_frac.png', dpi = 150)
plt.show()

"""
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
"""
