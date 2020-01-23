import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import ls_rain_amt_key
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

paths = {'U10' : '/nerc/n02/n02/xb899100/CloudTrail/Control2/',
         'U09' : '/nerc/n02/n02/xb899100/CloudTrail/U09/',
         'U08' : '/nerc/n02/n02/xb899100/CloudTrail/U08/',
         'U07' : '/nerc/n02/n02/xb899100/CloudTrail/U07/',
         'U06' : '/nerc/n02/n02/xb899100/CloudTrail/U06/',
         'U05' : '/nerc/n02/n02/xb899100/CloudTrail/U05/'}

keys = ['U05', 'U06', 'U07', 'U08', 'U09', 'U10']
# Island geometry
island_area = 50.0
island_radius = np.sqrt(island_area/np.pi)

R = {}
lwp_data = {}
for exp in keys:
    lwp_data[exp] = {}
    dx = 0.1
    with Dataset(paths[exp] + 'lwp_00.nc', 'r') as lwp_nc:
        lwp_data[exp]['rain'] = lwp_nc.variables[ls_rain_amt_key][:]*1.
        lwp_time_key  = [tkey for tkey in lwp_nc.variables.keys() if 'accum' in tkey][0]
        lwp_data[exp]['time'] = lwp_nc.variables[lwp_time_key][:]*1. + [240.0 if exp in ['U05'] else 0.0][0]
        lwp_data[exp]['x'], lwp_data[exp]['y'] = np.meshgrid(np.arange(lwp_data[exp]['rain'].shape[2])*dx, np.arange(lwp_data[exp]['rain'].shape[1])*dx)
    x_c, y_c = [100. + island_radius if exp in ['U05', 'U10'] else 108.][0], [4*island_radius if exp in ['U05', 'U10'] else 16.][0]
    R[exp] = np.sqrt((lwp_data[exp]['x'] - x_c)**2 + (lwp_data[exp]['y'] - y_c)**2)

# define plotting parameters
my_levels = np.array([0.01, 0.5, 1., 2., 4., 8., 16., 32.])
n_levels = float(len(my_levels))
my_cmap = mpl.cm.get_cmap('YlGnBu')
my_colors = my_cmap((np.arange(n_levels)+1.0)/(n_levels))

# Panel labels
labels = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']
### Make the plot ###
fig = plt.figure()
for exp in keys:
    idx = keys.index(exp)+1
    it0 = np.where(lwp_data[exp]['time'] == 360.)[0][0]
    it1 = np.where(lwp_data[exp]['time'] == 1080.)[0][0]
    ax = fig.add_subplot(3, 2, idx, adjustable = 'box', aspect = 1)
    precip = lwp_data[exp]['rain'][it1,:,:] - lwp_data[exp]['rain'][it0,:,:]
    im = ax.contourf(lwp_data[exp]['x'], lwp_data[exp]['y'], precip, levels = my_levels, colors = my_colors, extend = 'max')
    ax.contour(lwp_data[exp]['x'], lwp_data[exp]['y'], R[exp], levels = [island_radius], colors = ['k'])
    ax.set_title(labels[idx-1] + ' ' + exp + ', P$_{max}$ = ' + str(round(precip.max(), 2)) + ' mm')
    iy, ix = np.where(precip == precip.max())
    ax.plot(lwp_data[exp]['x'][iy[0],ix[0]], lwp_data[exp]['y'][iy[0],ix[0]], color = 'purple', marker = 'x', mew = 2.5, ms = 7.5)
    # do some beautification
    if idx in [1, 3, 5]:
        ax.set_ylabel('y (km)')
    if idx in [1, 2, 3, 4]:
        ax.set_xticklabels([''])
    if idx in [2, 4, 6]:
        ax.set_yticklabels([''])
    if idx in [5, 6]:
        ax.set_xlabel('x (km)')

# add colorbar
axins = inset_axes(ax, width = "100%", height = "10%", loc = 3, bbox_to_anchor = (-1.05, -0.8, 2.05, 2.0), bbox_transform = ax.transAxes, borderpad = 0.0)
cbax = plt.colorbar(im, cax = axins, label = u'Daytime Precipitation Total (mm)', orientation = 'horizontal')
cbax.set_ticklabels([tick if tick < 1 else int(tick) for tick in my_levels])
plt.subplots_adjust(bottom = 0.2, wspace = 0.05, hspace = -0.25)
plt.savefig('../Ch5_Figure22.png', dpi = 250, bbox_inches = 'tight')
plt.show()

