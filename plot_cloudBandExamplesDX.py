import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import lwp_key
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

path_start = '/nerc/n02/n02/xb899100/CloudTrail/Control_'
path_end   = '_HRIC_INV/'
<<<<<<< HEAD
experiments = ['DX0200', 'DX0400', 'DX0800', 'DX1600']
=======
experiments = ['0200m', '0400m', '0800m', '1600m']
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
# Island geometry
island_area = 50.0
island_radius = np.sqrt(island_area/np.pi)
x_c, y_c = 108., 16.
R = {}
lwp_data = {}
for exp in experiments:
    lwp_data[exp] = {}
<<<<<<< HEAD
    dx = float(exp[2:])/1000.
    with Dataset(path_start + exp[2:] + 'm' + path_end + 'lwp_00.nc', 'r') as lwp_nc:
=======
    dx = float(exp[:-1])/1000.
    with Dataset(path_start + exp + path_end + 'lwp_00.nc', 'r') as lwp_nc:
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
        lwp_data[exp]['lwp'] = lwp_nc.variables[lwp_key][:]*1.
        lwp_time_key  = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
        lwp_data[exp]['time'] = lwp_nc.variables[lwp_time_key][:]*1.
        lwp_data[exp]['x'], lwp_data[exp]['y'] = np.meshgrid(np.arange(lwp_data[exp]['lwp'].shape[2])*dx, np.arange(lwp_data[exp]['lwp'].shape[1])*dx)
    
    R[exp] = np.sqrt((lwp_data[exp]['x'] - x_c)**2 + (lwp_data[exp]['y'] - y_c)**2)

# define plotting parameters
my_levels = np.array([1., 5., 10., 50., 100., 500., 1000.])
n_levels = float(len(my_levels))
my_cmap = mpl.cm.get_cmap('Greys')
my_colors = my_cmap((np.arange(n_levels)+1.0)/(n_levels))

# Panel labels
labels = ['a)', 'b)', 'c)', 'd)']
### Make the plot ###
<<<<<<< HEAD
fig = plt.figure(figsize = (10, 5))
=======
fig = plt.figure(figsize = (8.5, 4))
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
for exp in experiments:
    idx = experiments.index(exp)+1
    it = np.where(lwp_data[exp]['time'] == 720.)[0][0]
    ax = fig.add_subplot(2, 2, idx, adjustable = 'box', aspect = 1)
    im = ax.contourf(lwp_data[exp]['x'], lwp_data[exp]['y'], lwp_data[exp]['lwp'][it,:,:]*1000.0, levels = my_levels, colors = my_colors, extend = 'max')
    ax.contour(lwp_data[exp]['x'], lwp_data[exp]['y'], np.where(lwp_data[exp]['lwp'][it,:,:] > 0.0, 1.0, 0.0), levels = [0.5], colors = ['b'])
    ax.contour(lwp_data[exp]['x'], lwp_data[exp]['y'], R[exp], levels = [island_radius], colors = ['r'])
<<<<<<< HEAD
    ax.set_title(labels[idx-1] + ' ' + exp)
=======
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
plt.colorbar(im, cax = axins, label = u'Liquid Water Path (g kg$^{-1}$)', orientation = 'horizontal')
plt.subplots_adjust(bottom = 0.3, wspace = 0.05, hspace = 0.00)
plt.savefig('../cloudBandExamples_dx.png', dpi = 250, bbox_inches = 'tight')
=======
axins = inset_axes(ax, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (-1.05, -0.7, 2.05, 2.0), bbox_transform = ax.transAxes, borderpad = 0.0)
plt.colorbar(im, cax = axins, label = u'Liquid Water Path (g kg$^{-1}$)', orientation = 'horizontal')
plt.subplots_adjust(bottom = 0.3, wspace = 0.05, hspace = 0.05)
plt.savefig('../cloudBandExamples_dx.png', dpi = 150, bbox_inches = 'tight')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
plt.show()

