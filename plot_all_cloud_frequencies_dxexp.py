import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from STASH_keys import lwp_key
from coarse_graining import coarse_grain as cg

paths = {'DX0100' : '/nerc/n02/n02/xb899100/CloudTrail/Control2/',
         'DX0200' : '/nerc/n02/n02/xb899100/CloudTrail/Control_0200m_HRIC_INV/',
         'DX0400' : '/nerc/n02/n02/xb899100/CloudTrail/Control_0400m_HRIC_INV/',
         'DX0800' : '/nerc/n02/n02/xb899100/CloudTrail/Control_0800m_HRIC_INV/',
         'DX1600' : '/nerc/n02/n02/xb899100/CloudTrail/Control_1600m_HRIC_INV/'}

keys = ['DX0100', 'DX0200', 'DX0400', 'DX0800', 'DX1600']
my_levels = np.arange(0, 0.51, 0.05)
scales = [1, 2, 4, 8, 16]

fig, axes = plt.subplots(nrows = 5, ncols = 2)
fig.set_size_inches(8, 12)
for key in keys:
    print key
    # open the lwp_00.nc file
    lwp_nc = Dataset(paths[key] + 'lwp_00.nc', 'r')
    time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
    times = lwp_nc.variables[time_key][:]*1.
    idx = [it for it in range(len(times)) if (360.0 <= times[it]) and (times[it] <= 1080.0)]
    cloud_frequency = np.nanmean(np.where(lwp_nc.variables[lwp_key][idx,:,:] > 0, 1.0, 0.0), axis = 0)
    lwp_nc.close()
    
    X, Y = np.meshgrid(np.arange(cloud_frequency.shape[1])*float(key[2:]), np.arange(cloud_frequency.shape[0])*float(key[2:]))
    R_i = 1000*np.sqrt(50.0/np.pi)
    x_c = 108000.0
    y_c = 16000.0
    
    R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)
    ax = axes.flat[keys.index(key)*2]
    ax.set_adjustable('box')
    ax.set_aspect(1)
    im = ax.contourf(X/1000., Y/1000., cloud_frequency, cmap = 'Greys_r', levels = my_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax.set_title(key)
    if (keys.index(key) + 1) in [1, 2, 3, 4]:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    if (keys.index(key) + 1) in [2, 4, 6]:
        ax.set_yticklabels([''])
    else:
        ax.set_ylabel('y (km)')
    
    # if DX100 also coarse grain
    if key == 'DX0100':
        for scale in scales:
            print 'scale = ' + str(scale)
            ax = axes.flat[scales.index(scale)*2 + 1]
            ax.set_adjustable('box')
            ax.set_aspect(1)
            # coarse grain
            cfcg = cg(cloud_frequency, scale, scale)
            xcg, ycg = np.meshgrid(np.arange(cfcg.shape[1])*float(key[2:])*scale, np.arange(cfcg.shape[0])*float(key[2:])*scale)
            # plot
            im = ax.contourf(xcg/1000., ycg/1000., cfcg, cmap = 'Greys_r', levels = my_levels, extend = 'max')
            ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
            ax.set_title(key + ' cg to ' + str(int(key[2:])*scale) + 'm')

cbar_ax = fig.add_axes([0.1, 0.1, 0.85, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Cloud Frequency')
cbar.ax.set_xticklabels([str(level) for level in my_levels])
plt.subplots_adjust(left = 0.1, bottom = 0.2, right = 0.95, top = 0.9, wspace = 0.1, hspace = 0.1)
plt.savefig('../dx_comparison_cloud_frequency.png', dpi = 150, bbox_inches = 'tight')
plt.show()

