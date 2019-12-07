import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import u_key, v_key
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from analysis_tools import fromComponents
"""
Figure 5.8: Snapshots of the 10 m wind speed at 6, 9, and 12.
"""

path = '/nerc/n02/n02/xb899100/CloudTrail/Control2/'

hours = ["{0:02d}".format(hour) for hour in range(0, 13, 3)]
for hour in hours:
    with Dataset(path + 'wind_' + hour + '.nc', 'r') as wind_nc:
        if hour == hours[0]:
            u_data       = wind_nc.variables[u_key][1:,:,:,:]
            v_data       = wind_nc.variables[v_key][1:,:,:,:]
            times_key    = [tkey for tkey in wind_nc.variables.keys() if 'min' in tkey][0]
            times        = wind_nc.variables[times_key][1:]*1.
            z            = wind_nc.variables['thlev_zsea_theta'][:]*1.
        else:
            u_data     = np.concatenate((u_data, wind_nc.variables[u_key][:]), axis = 0)
            v_data     = np.concatenate((v_data, wind_nc.variables[v_key][:]), axis = 0)
            times      = np.concatenate((times, wind_nc.variables[times_key][:]*1.), axis = 0)

x, y = np.meshgrid(np.arange(u_data.shape[3])*0.1, np.arange(u_data.shape[2])*0.1)

# compute the theta anomalies at the horizontal mean height of the boundary layer at a given time
# choose times of interest
times_of_interest = [6, 9, 12]
# convert from hours to minutes since simulation start
times_of_interest = [time*60.0 for time in times_of_interest]
# find the horizontal mean boundary layer height at these times
wind_select = []
iz = np.where(np.abs(z - 10.0) == np.min(np.abs(z - 10.0)))[0][0]+1
for time in times_of_interest:
    it = np.where(time == times)[0][0]
    wind_spd, wind_dir = fromComponents(u_data[it,iz,:,:], v_data[it,iz,:,:])
    wind_select.append(wind_spd)

wind_select = np.array(wind_select)

# define plotting parameters
my_levels = np.array([level for level in np.arange(-3.0, 3.1, 0.5) if level != 0.0]) + 8.0#wind_spd

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
    im = ax.contourf(x, y, wind_select[it,:,:], levels = my_levels, cmap = 'bwr', extend = 'both')
    im.cmap.set_under('navy')
    im.cmap.set_over('firebrick')
    ax.contour(x, y, island_mask, levels = [0.5], colors = ['k'])
    ax.set_title(labels[times_of_interest.index(time)] + ' T+' + str(int(time/60.)) + 'hrs, H = ' + str(int(H[np.where(times == time)[0][0]])) + u'W m$^{-2}$, E = ' + str(int(E[np.where(times == time)[0][0]])) + u'W m$^{-2}$', fontsize = 12)
    # axes labels
    ax.set_ylabel('y (km)')
    if (times_of_interest.index(time)+1) != 3:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')

ax.arrow(island_x, island_y, dx = u_data[-1,iz,:,:].mean()*3.6, dy = v_data[-1,iz,:,:].mean()*3.6, color = 'k', head_width = 3, length_includes_head = True)
cbar_ax = inset_axes(ax, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (0, -0.5, 1.0, 2.0), bbox_transform = ax.transAxes, borderpad = 0.0)
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = '$U$ (m s$^{-1}$)')
cbar.set_ticks(np.array([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0])+8.0)
plt.subplots_adjust(left = 0.10, bottom = 0.175, right = 1.0, wspace = 0.0, hspace = 0.225)
plt.savefig('../Ch5_Figure08.png', dpi = 250, bbox_inches = 'tight')
plt.show()

