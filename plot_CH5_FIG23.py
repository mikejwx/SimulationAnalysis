import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from STASH_keys import theta_key, lwp_key

paths = {'U05' : '/nerc/n02/n02/xb899100/CloudTrail/U05/',
         'U05_NoRain' : '/nerc/n02/n02/xb899100/CloudTrail/U05_NoRain/',
         'U05_H125' : '/nerc/n02/n02/xb899100/CloudTrail/U05_H125/'}

keys = ['U05', 'U05_NoRain', 'U05_H125']
my_T_levels = np.array([-2, -1, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 1, 2])
n_T_levels = float(len(my_T_levels))
my_T_cmap = mpl.cm.get_cmap('bwr')
my_T_colors = my_T_cmap((np.arange(n_T_levels)+1.0)/(n_T_levels))
my_C_levels = np.arange(0, 0.61, 0.075)

# Read the data
my_data = {}
for key in keys:
    print key
    my_data[key] = {}
    with Dataset(paths[key] + 'bouy_04.nc', 'r') as bouy_nc:
        time_key = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
        times = bouy_nc.variables[time_key][:]*1. + 240.0
        idx = np.where(times == 720.0)[0][0]
        my_data[key]['thetap'] = bouy_nc.variables[theta_key][idx,1,:,:] - bouy_nc.variables[theta_key][idx,1,:,:].mean()
    
    with Dataset(paths[key] + 'lwp_00.nc', 'r') as lwp_nc:
        time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
        times = lwp_nc.variables[time_key][:]*1. + 240.0
        idx = [it for it in range(len(times)) if (360.0 <= times[it]) and (times[it] <= 1080.0)]
        my_data[key]['cloudfrequency'] = np.nanmean(np.where(lwp_nc.variables[lwp_key][idx,:,:] > 0, 1.0, 0.0), axis = 0)

letters = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']
# define the domain and island
X, Y = np.meshgrid(np.arange(my_data[key]['thetap'].shape[1])*100.0, np.arange(my_data[key]['thetap'].shape[0])*100.0)
R_i = 1000*np.sqrt(50.0/np.pi)
x_c = 100000.0 + R_i
y_c = 4*R_i
R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)

# Make the plots
panel_no = 1
fig = plt.figure()
for key in keys:
    print key
    # Do the cloud frequency panel
    ax1 = fig.add_subplot(3, 2, panel_no, adjustable = 'box', aspect = 1)
    im1 = ax1.contourf(X/1000., Y/1000., my_data[key]['cloudfrequency'], levels = my_C_levels, cmap = 'Greys_r', extend = 'max')
    ax1.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax1.set_title(letters[panel_no-1] + ' ' + key + ' cloud frequency')
    if panel_no in [1, 2, 3, 4]:
        ax1.set_xticklabels([''])
    else:
        ax1.set_xlabel('x (km)')
    if panel_no in [2, 4, 6]:
        ax1.set_yticklabels([''])
    else:
        ax1.set_ylabel('y (km)')
    
    # Do the theta prime panel
    panel_no += 1
    ax2 = fig.add_subplot(3, 2, panel_no, adjustable = 'box', aspect = 1)
    im2 = ax2.contourf(X/1000., Y/1000., my_data[key]['thetap'], colors = my_T_colors, levels = my_T_levels, extend = 'both')
    im2.cmap.set_under('navy')
    im2.cmap.set_over('firebrick')
    ax2.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['k'])
    ax2.set_title(letters[panel_no-1] + ' ' + key + ' $\\theta^{\prime}$')
    if panel_no in [1, 2, 3, 4]:
        ax2.set_xticklabels([''])
    else:
        ax2.set_xlabel('x (km)')
    if panel_no in [2, 4, 6]:
        ax2.set_yticklabels([''])
    else:
        ax2.set_ylabel('y (km)')
    panel_no += 1

cbar_ax = inset_axes(ax1, width = "100%", height = "10%", loc = 3, bbox_to_anchor = (0.0, -0.8, 1., 1.), bbox_transform = ax1.transAxes, borderpad = 0.0)
cbar = fig.colorbar(im1, cax = cbar_ax, orientation = 'horizontal', label = 'Cloud Frequency')
cbar.set_ticks([level for level in my_C_levels[::2]])

cbar_ax = inset_axes(ax2, width = "100%", height = "10%", loc = 3, bbox_to_anchor = (0.0, -0.8, 1., 1.), bbox_transform = ax2.transAxes, borderpad = 0.0)
cbar = fig.colorbar(im2, cax = cbar_ax, orientation = 'horizontal', label = '$\\theta^{\prime}$ (K)')
cbar.set_ticks([-2, -1, -0.25, 0, 0.1, 0.5, 2])
cbar.ax.set_xticklabels([-2, -1, -0.25, '', 0.1, 0.5, 2])
plt.subplots_adjust(bottom = 0.2, wspace = 0.05, hspace = -0.25)
plt.savefig('../Ch5_Figure23.png', dpi = 250, bbox_inches = 'tight')
plt.show()

