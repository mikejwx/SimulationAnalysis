import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import theta_key
import matplotlib as mpl

path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'

hours = ["{0:02d}".format(hour) for hour in range(0, 13, 3)]
for hour in hours:
    bouy_nc = Dataset(path + 'bouy_' + hour + '.nc', 'r')
    if hour == hours[0]:
        theta_data = bouy_nc.variables[theta_key][:]*1.
        times_key = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
        times = bouy_nc.variables[times_key][:]*1.
        z = bouy_nc.variables['thlev_zsea_theta'][:]*1.
    else:
        theta_data = np.concatenate((theta_data, bouy_nc.variables[theta_key][:]*1.), axis = 0)
        times = np.concatenate((times, bouy_nc.variables[times_key][:]*1.), axis = 0)
    
    bouy_nc.close()

x, y = np.meshgrid(np.arange(theta_data.shape[3])*0.1, np.arange(theta_data.shape[2])*0.1)

times_of_interest = [6, 7, 8, 9, 10, 11]
times_of_interest = [time*60.0 for time in times_of_interest]

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

labels = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']
iz_max = np.where(np.abs(z - 1000) == np.min(np.abs(z - 1000)))[0][0]+1
for iz in range(iz_max):
    fig, axes = plt.subplots(nrows = 3, ncols = 2)
    for time in times_of_interest:
        it = np.where(np.abs(times - time) == np.min(np.abs(times - time)))[0][0]
        ax = axes.flat[times_of_interest.index(time)]
        ax.set_adjustable('box')
        ax.set_aspect(1)
        im = ax.contourf(x, y, theta_data[it,iz,:,:] - np.nanmean(theta_data[it,iz,:,:]), levels = my_levels, colors = my_colors, extend = 'both')
        ax.contour(x, y, island_mask, levels = [0.5], colors = ['k'])
        ax.set_title('T+' + str(int(time/60.)) + 'hrs, H = ' + str(int(H[np.where(times == time)[0][0]])) + u'W m$^{-2}$, E = ' + str(int(E[np.where(times == time)[0][0]])) + u'W m$^{-2}$', fontsize = 12)
        ax.text(5, 25, labels[times_of_interest.index(time)], bbox=dict(facecolor='w', edgecolor='None'))
        
        if (times_of_interest.index(time)+1)%2 == 0:
            ax.set_yticklabels([''])
        else:
            ax.set_ylabel('y (km)')
        if (times_of_interest.index(time)+1) in [1, 2, 3, 4]:
            ax.set_xticklabels([''])
        else:
            ax.set_xlabel('x (km)')
    
    cbar_ax = fig.add_axes([0.1, 0.1, 0.85, 0.02])
    cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = '$\\theta^{\prime}$ (K)')
    cbar.ax.set_xticklabels(['-2', '-1', '-0.5', '-0.25', '-0.1', '0.1', '0.25', '0.5', '1', '2'])
    plt.subplots_adjust(left = 0.1, bottom = 0.2, right = 0.95, top = 0.9, wspace = 0.1, hspace = 0.1)
    plt.savefig('../warm_plume_figure' + ['_' + str(int(z[iz])) + 'm'][0] + '.png', dpi = 150, bbox_inches = 'tight')
    plt.show()

