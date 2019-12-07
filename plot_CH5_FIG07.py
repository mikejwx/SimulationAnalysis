import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import theta_key, zi_new_key
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

path = '/nerc/n02/n02/xb899100/CloudTrail/Control2/'

hours = ["{0:02d}".format(hour) for hour in range(0, 13, 3)]
for hour in hours:
    bouy_nc = Dataset(path + 'bouy_' + hour + '.nc', 'r')
    zi_nc   = Dataset(path + 'zi_' + hour + '.nc', 'r')
    if hour == hours[0]:
        zi_data      = zi_nc.variables[zi_new_key][:]*1.
        theta_data   = bouy_nc.variables[theta_key][1:,:,:,:]*1.
        times_key    = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
        times        = bouy_nc.variables[times_key][1:]*1.
        z            = bouy_nc.variables['thlev_zsea_theta'][:]*1.
    else:
        zi_data    = np.concatenate((zi_data, zi_nc.variables[zi_new_key][:]*1.), axis = 0)
        theta_data = np.concatenate((theta_data, bouy_nc.variables[theta_key][:]*1.), axis = 0)
        times      = np.concatenate((times, bouy_nc.variables[times_key][:]*1.), axis = 0)
        
    bouy_nc.close()

x, y = np.meshgrid(np.arange(theta_data.shape[3])*0.1, np.arange(theta_data.shape[2])*0.1)

# compute the theta anomalies at the horizontal mean height of the boundary layer at a given time
# choose times of interest
times_of_interest = [6, 9, 12]
# convert from hours to minutes since simulation start
times_of_interest = [time*60.0 for time in times_of_interest]
# find the horizontal mean boundary layer height at these times
theta_anom = []
for time in times_of_interest:
    it = np.where(time == times)[0][0]
    zi_mean = zi_data[it,:,:].mean()
    print zi_mean
    iz = np.where(np.abs(z - zi_mean) == np.min(np.abs(z - zi_mean)))[0][0]+1
    theta_anom.append(theta_data[it,iz,:,:] - np.nanmean(theta_data[it,iz,:,:]))

theta_anom = np.array(theta_anom)

# define plotting parameters
my_levels = np.array([-2.0, -1.0, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 1.0, 2.0])
n_levels = float(len(my_levels))
my_cmap = mpl.cm.get_cmap('bwr')
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

labels = ['a)', 'b)', 'c)']
fig = plt.figure(figsize = (6, 10))
for time in times_of_interest:
    it = times_of_interest.index(time)
    ax = fig.add_subplot(5, 1, times_of_interest.index(time)+1, adjustable = 'box', aspect = 1)
    im = ax.contourf(x, y, theta_anom[it,:,:], levels = my_levels, colors = my_colors, extend = 'both')
    im.cmap.set_over('firebrick')
    im.cmap.set_under('navy')
    ax.contour(x, y, island_mask, levels = [0.5], colors = ['k'])
    ax.set_title(labels[times_of_interest.index(time)] + ' T+' + str(int(time/60.)) + 'hrs, H = ' + str(int(H[np.where(times == time)[0][0]])) + u'W m$^{-2}$, E = ' + str(int(E[np.where(times == time)[0][0]])) + u'W m$^{-2}$', fontsize = 12)
    # axes labels
    ax.set_ylabel('y (km)')
    if (times_of_interest.index(time)+1) != 3:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')

ax.arrow(island_x - island_radius, island_y, dx = -35.0, dy = -15.0, head_width = 3.0, color = 'k', ls = '--', lw = 2)
ax.arrow(island_x - island_radius, island_y, dx = -40.0, dy = -1.0, head_width = 3.0, color = 'k', ls = '--', lw = 2)
ax.arrow(70.0, island_y - 1.5*island_radius, dx = -35.0, dy = -8.0, head_width = 3.0, color = 'k', ls = ':', lw = 2)
cbar_ax = inset_axes(ax, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (0, -0.5, 1.0, 2.0), bbox_transform = ax.transAxes, borderpad = 0.0)
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = '$\\theta^{\prime}$ (K)')
cbar.set_ticks([-2, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 2])
cbar.ax.set_xticklabels([-2, -1, -0.5, -0.25, -0.1, '', 0.1, 0.25, 0.5, 1, 2])
plt.subplots_adjust(left = 0.10, bottom = 0.175, right = 1.0, wspace = 0.0, hspace = 0.225)
plt.savefig('../Ch5_Figure07.png', dpi = 250, bbox_inches = 'tight')
plt.show()

