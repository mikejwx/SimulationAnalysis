import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import lwp_key
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

path = '/nerc/n02/n02/xb899100/CloudTrail/Control2/'

data_nc = Dataset(path + 'lwp_00.nc', 'r')
lwp_data = data_nc.variables[lwp_key][:]*1.
times_key = [tkey for tkey in data_nc.variables.keys() if 'min' in tkey][0]
times = data_nc.variables[times_key][:]*1.
data_nc.close()

x, y = np.meshgrid(np.arange(lwp_data.shape[2])*0.1, np.arange(lwp_data.shape[1])*0.1)

times_of_interest = [3, 5, 7, 9, 11]
times_of_interest = [time*60.0 for time in times_of_interest]

# define plotting parameters
my_levels = np.array([1., 5., 10., 50., 100., 500., 1000.])
n_levels = float(len(my_levels))
my_cmap = mpl.cm.get_cmap('Greys')
my_colors = my_cmap((np.arange(n_levels)+1.0)/(n_levels))

island_area = 50.0
island_radius = np.sqrt(island_area/np.pi)
island_x = 100.0 + island_radius
island_y = 4*island_radius
R = np.sqrt((x-island_x)**2 + (y-island_y)**2)
island_mask = np.where(R > island_radius, 1.0, 0.0)

# island surface fluxes
H0 = 250.
E0 = 250.
t0 = 720.
dt = 720.
H = np.array([np.nanmax([0, h]) for h in H0*np.cos((np.pi/2)*(t0-times)/(dt/2.))**1.5])
E = np.array([np.nanmax([0, e]) for e in E0*np.cos((np.pi/2)*(t0-times)/(dt/2.))**1.3])

labels = ['a)', 'b)', 'c)', 'd)', 'e)']
fig = plt.figure(figsize = (6, 10))
for time in times_of_interest:
    it = np.where(np.abs(times - time) == np.min(np.abs(times - time)))[0][0]
    ax = fig.add_subplot(5, 1, times_of_interest.index(time)+1, adjustable = 'box', aspect = 1)
    im = ax.contourf(x, y, lwp_data[it,:,:]*1000., levels = my_levels, colors = my_colors, extend = 'max')
    ax.contour(x, y, np.where(lwp_data[it,:,:] > 0, 1.0, 0.0), levels = [0.5], colors = ['b'], linewidths = [0.25])
    ax.contour(x, y, island_mask, levels = [0.5], colors = ['r'])
    ax.set_title(labels[times_of_interest.index(time)] + ' T+' + str(int(time/60.)) + 'hrs, H = ' + str(int(H[np.where(times == time)[0][0]])) + u'W m$^{-2}$, E = ' + str(int(E[np.where(times == time)[0][0]])) + u'W m$^{-2}$', fontsize = 12)
    # axes labels
    ax.set_ylabel('y (km)')
    if (times_of_interest.index(time)+1) != 5:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')

cbar_ax = inset_axes(ax, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (0, -0.7, 1.0, 2.0), bbox_transform = ax.transAxes, borderpad = 0.0)
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'LWP (g m$^{-2}$)')
cbar.ax.set_xticklabels([str(int(level)) for level in my_levels])
plt.subplots_adjust(left = 0.10, bottom = 0.175, right = 1.0, wspace = 0.0, hspace = 0.225)
plt.savefig('../Ch5_Figure02.png', dpi = 250, bbox_inches = 'tight')
plt.close('all')

