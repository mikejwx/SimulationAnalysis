import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from STASH_keys import lwp_key

paths = {'H250E250' : '/nerc/n02/n02/xb899100/CloudTrail/Control/',
         'H400E250' : '/nerc/n02/n02/xb899100/CloudTrail/H400E250/',
         'H350E250' : '/nerc/n02/n02/xb899100/CloudTrail/H350E250/',
         'H300E250' : '/nerc/n02/n02/xb899100/CloudTrail/H300E250/',
         'H200E250' : '/nerc/n02/n02/xb899100/CloudTrail/H200E250/',
         'H150E250' : '/nerc/n02/n02/xb899100/CloudTrail/H150E250/',
         'H100E250' : '/nerc/n02/n02/xb899100/CloudTrail/H100E250/',
         'H050E250' : '/nerc/n02/n02/xb899100/CloudTrail/H050E250/'}

keys = ['H050E250', 'H100E250', 'H150E250', 'H200E250', 'H250E250', 'H300E250', 'H350E250', 'H400E250']
my_levels = np.arange(0, 0.51, 0.05)

fig, axes = plt.subplots(nrows = 4, ncols = 2)
fig.set_size_inches(10, 12)
for key in keys:
    print key
    # open the lwp_00.nc file
    lwp_nc = Dataset(paths[key] + 'lwp_00.nc', 'r')
    time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
    times = lwp_nc.variables[time_key][:]*1.
    idx = [it for it in range(len(times)) if (360.0 <= times[it]) and (times[it] <= 1080.0)]
    cloud_frequency = np.nanmean(np.where(lwp_nc.variables[lwp_key][idx,:,:] > 0, 1.0, 0.0), axis = 0)
    lwp_nc.close()
    
    X, Y = np.meshgrid(np.arange(cloud_frequency.shape[1])*100.0, np.arange(cloud_frequency.shape[0])*100.0)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['H250E250']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)
    ax = axes.flat[keys.index(key)]
    ax.set_adjustable('box')
    ax.set_aspect(1)
    im = ax.contourf(X/1000., Y/1000., cloud_frequency, cmap = 'Greys_r', levels = my_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax.set_title(key + ' at T+' + str(times[idx[0]]) + ' to ' + str(times[idx[-1]]) + 'mins')
    if (keys.index(key) + 1) in [1, 2, 3, 4]:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    if (keys.index(key) + 1) in [2, 4, 6]:
        ax.set_yticklabels([''])
    else:
        ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.1, 0.1, 0.85, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Cloud Frequency')
cbar.ax.set_xticklabels([str(level) for level in my_levels])
plt.subplots_adjust(left = 0.1, bottom = 0.2, right = 0.95, top = 0.9, wspace = 0.1, hspace = 0.1)
plt.savefig('../heat_comparison_cloud_frequency.png', dpi = 150, bbox_inches = 'tight')
plt.show()

