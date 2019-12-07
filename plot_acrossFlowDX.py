import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import n_key
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

path_start = '/nerc/n02/n02/xb899100/CloudTrail/Control_'
path_end   = '_HRIC_INV/'
experiments = ['0200m', '0400m', '0800m', '1600m']
# Island geometry
island_area = 50.0
island_radius = np.sqrt(island_area/np.pi)
x_c, y_c = 108., 16.
R = {}
n_data = {}
for exp in experiments:
    n_data[exp] = {}
    dx = float(exp[:-1])/1000.
    with Dataset(path_start + exp + path_end + 'wind_09.nc', 'r') as wind_nc:
        n_data[exp]['n'] = wind_nc.variables[n_key][:]*1.
        n_time_key  = [tkey for tkey in wind_nc.variables.keys() if 'min' in tkey][0]
        n_data[exp]['time'] = wind_nc.variables[n_time_key][:]*1.
        n_data[exp]['x'], n_data[exp]['y'] = np.meshgrid(np.arange(n_data[exp]['n'].shape[3])*dx, np.arange(n_data[exp]['n'].shape[2])*dx)
        z_theta = wind_nc.variables['thlev_zsea_theta'][:]*1.
    R[exp] = np.sqrt((n_data[exp]['x'] - x_c)**2 + (n_data[exp]['y'] - y_c)**2)

# define plotting parameters
my_levels = np.array([level for level in np.arange(-3.5, 3.6, 0.5) if level != 0.0])

# Panel labels
labels = ['a)', 'b)', 'c)', 'd)']
iz = np.where(np.abs(z_theta - 10.0) == np.min(np.abs(z_theta - 10.0)))[0][0]
### Make the plot ###
fig = plt.figure(figsize = (8.5, 4))
for exp in experiments:
    idx = experiments.index(exp)+1
    it = np.where(n_data[exp]['time'] == 720.)[0][0]
    ax = fig.add_subplot(2, 2, idx, adjustable = 'box', aspect = 1)
    im = ax.contourf(n_data[exp]['x'], n_data[exp]['y'], n_data[exp]['n'][it,iz,:,:], levels = my_levels, cmap = 'bwr', extend = 'both')
    ax.contour(n_data[exp]['x'], n_data[exp]['y'], R[exp], levels = [island_radius], colors = ['k'])
    ax.set_title('DX' + exp[:-1])
    ax.text(5, 22.5, labels[idx-1], bbox = {'facecolor':'w','edgecolor':'none'})
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
axins = inset_axes(ax, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (-1.05, -0.7, 2.05, 2.0), bbox_transform = ax.transAxes, borderpad = 0.0)
cbax = plt.colorbar(im, cax = axins, label = u'Across-Flow Wind (m s$^{-1}$)', orientation = 'horizontal')
cbax.set_ticks([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0])
plt.subplots_adjust(bottom = 0.3, wspace = 0.05, hspace = 0.05)
plt.savefig('../acrossFlow_dx.png', dpi = 150, bbox_inches = 'tight')
plt.show()

