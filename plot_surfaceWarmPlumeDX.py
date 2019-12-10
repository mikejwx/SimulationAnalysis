import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import theta_key, temp_key, pthe_key
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from SkewT_archer import getThetaE
<<<<<<< HEAD
import matplotlib as mpl
=======
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

path_start = '/nerc/n02/n02/xb899100/CloudTrail/Control_'
path_end   = '_HRIC_INV/'
experiments = ['0200m', '0400m', '0800m', '1600m']
# Island geometry
island_area = 50.0
island_radius = np.sqrt(island_area/np.pi)
x_c, y_c = 108., 16.
R = {}
T_data = {}
for exp in experiments:
    T_data[exp] = {}
    dx = float(exp[:-1])/1000.
    with Dataset(path_start + exp + path_end + 'bouy_09.nc', 'r') as bouy_nc:
        T_data[exp]['theta'] = bouy_nc.variables[theta_key][:]*1.
        T_data[exp]['temp']  = bouy_nc.variables[temp_key][:]*1.
        T_time_key  = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
        T_data[exp]['time'] = bouy_nc.variables[T_time_key][:]*1.
        T_data[exp]['x'], T_data[exp]['y'] = np.meshgrid(np.arange(T_data[exp]['theta'].shape[3])*dx, np.arange(T_data[exp]['theta'].shape[2])*dx)
        z_theta = bouy_nc.variables['thlev_zsea_theta'][:]*1.
    R[exp] = np.sqrt((T_data[exp]['x'] - x_c)**2 + (T_data[exp]['y'] - y_c)**2)
    with Dataset(path_start + exp + path_end + 'fluxes_09.nc', 'r') as fluxes_nc:
        T_data[exp]['pthe'] = fluxes_nc.variables[pthe_key][:]*1.
<<<<<<< HEAD

# define plotting parameters
my_levels = np.array([-2.0, -1.0, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 1.0, 2.0])
n_levels = float(len(my_levels))
my_cmap = mpl.cm.get_cmap('bwr')
my_colors = my_cmap((np.arange(n_levels)+1.0)/(n_levels))
=======
    # compute the equivalent potential temperature
    T_data[exp]['thetae'] = getThetaE(T_data[exp]['theta'], T_data[exp]['temp'], T_data[exp]['pthe'], t_units = 'K', p_units = 'Pa')

# define plotting parameters
my_levels = np.array([level for level in np.arange(-5., 5.1, 0.5) if level != 0.0])
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

# Panel labels
labels = ['a)', 'b)', 'c)', 'd)']
iz = np.where(np.abs(z_theta - 2.0) == np.min(np.abs(z_theta - 2.0)))[0][0]
### Make the plot ###
<<<<<<< HEAD
fig = plt.figure(figsize = (10, 5))
=======
fig = plt.figure(figsize = (8.5, 4))
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
for exp in experiments:
    idx = experiments.index(exp)+1
    it = np.where(T_data[exp]['time'] == 720.)[0][0]
    ax = fig.add_subplot(2, 2, idx, adjustable = 'box', aspect = 1)
<<<<<<< HEAD
    im = ax.contourf(T_data[exp]['x'], T_data[exp]['y'], T_data[exp]['theta'][it,iz,:,:] - T_data[exp]['theta'][it,iz,:,:].mean(), levels = my_levels, colors = my_colors, extend = 'both')
    im.cmap.set_under('navy')
    im.cmap.set_over('firebrick')
    ax.contour(T_data[exp]['x'], T_data[exp]['y'], R[exp], levels = [island_radius], colors = ['k'])
    ax.set_title(labels[idx-1] + ' DX' + exp[:-1])
=======
    im = ax.contourf(T_data[exp]['x'], T_data[exp]['y'], T_data[exp]['thetae'][it,iz,:,:] - T_data[exp]['thetae'][it,iz,:,:].mean(), levels = my_levels, cmap = 'bwr', extend = 'both')
    ax.contour(T_data[exp]['x'], T_data[exp]['y'], R[exp], levels = [island_radius], colors = ['k'])
    ax.set_title('DX' + exp[:-1])
    ax.text(5, 22.5, labels[idx-1], bbox = {'facecolor':'w','edgecolor':'none'})
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
    # do some beautification
    if idx in [1, 3]:
        ax.set_ylabel('y (km)')
    if idx in [1, 2]:
        ax.set_xticklabels([''])
    if idx in [2, 4]:
        ax.set_yticklabels([''])
    if idx in [3, 4]:
        ax.set_xlabel('x (km)')

# add colorbar
<<<<<<< HEAD
axins = inset_axes(ax, width = "100%", height = "10%", loc = 3, bbox_to_anchor = (-1.05, -0.75, 2.05, 2.0), bbox_transform = ax.transAxes, borderpad = 0.0)
cbax = plt.colorbar(im, cax = axins, label = u'$\\theta^{\prime}$ (K)', orientation = 'horizontal')
cbar.set_ticks([-2, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 2])
cbar.ax.set_xticklabels([-2, -1, -0.5, -0.25, -0.1, '', 0.1, 0.25, 0.5, 1, 2])
plt.subplots_adjust(bottom = 0.3, wspace = 0.05, hspace = 0.00)
plt.savefig('../surfaceWarmPlume_dx.png', dpi = 250, bbox_inches = 'tight')
=======
axins = inset_axes(ax, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (-1.05, -0.7, 2.05, 2.0), bbox_transform = ax.transAxes, borderpad = 0.0)
cbax = plt.colorbar(im, cax = axins, label = u'$\\theta_{e,sfc}^{\prime}$ (K)', orientation = 'horizontal')
cbax.set_ticks([-5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
plt.subplots_adjust(bottom = 0.3, wspace = 0.05, hspace = 0.05)
plt.savefig('../surfaceWarmPlume_dx.png', dpi = 150, bbox_inches = 'tight')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
plt.show()

