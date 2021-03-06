import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from STASH_keys import lwp_key

paths = {'H250E250' : '/nerc/n02/n02/xb899100/CloudTrail/Control2/',
         'H400E250' : '/nerc/n02/n02/xb899100/CloudTrail/H400E250/',
         'H350E250' : '/nerc/n02/n02/xb899100/CloudTrail/H350E250/',
         'H300E250' : '/nerc/n02/n02/xb899100/CloudTrail/H300E250/',
         'H200E250' : '/nerc/n02/n02/xb899100/CloudTrail/H200E250/',
         'H150E250' : '/nerc/n02/n02/xb899100/CloudTrail/H150E250/',
         'H100E250' : '/nerc/n02/n02/xb899100/CloudTrail/H100E250/',
         'H050E250' : '/nerc/n02/n02/xb899100/CloudTrail/H050E250/'}

keys = ['H050E250', 'H100E250', 'H150E250', 'H200E250', 'H250E250', 'H300E250', 'H350E250', 'H400E250']
my_levels = np.arange(0, 0.61, 0.075)
# Read the data
my_data = {}
for key in keys:
    print key
    with Dataset(paths[key] + 'lwp_00.nc', 'r') as lwp_nc:
        time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
        times = lwp_nc.variables[time_key][:]*1.
        idx = [it for it in range(len(times)) if (360.0 <= times[it]) and (times[it] <= 1080.0)]
        my_data[key] = np.nanmean(np.where(lwp_nc.variables[lwp_key][idx,:,:] > 0, 1.0, 0.0), axis = 0)

letters = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)']
# Make the plots
fig = plt.figure(figsize = (9,7.5))
for key in keys:
    print key
    X, Y = np.meshgrid(np.arange(my_data[key].shape[1])*100.0, np.arange(my_data[key].shape[0])*100.0)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['H250E250']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)
    ax = fig.add_subplot(4, 2, keys.index(key)+1, adjustable = 'box', aspect = 1)
    im = ax.contourf(X/1000., Y/1000., my_data[key], cmap = 'Greys_r', levels = my_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax.set_title(letters[keys.index(key)] + ' ' + key)
    if (keys.index(key) + 1) in [1, 2, 3, 4, 5, 6]:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    if (keys.index(key) + 1) in [2, 4, 6, 8]:
        ax.set_yticklabels([''])
    else:
        ax.set_ylabel('y (km)')

cbar_ax = inset_axes(ax, width = "100%", height = "10%", loc = 3, bbox_to_anchor = (-1.05, -0.8, 2.05, 2.0), bbox_transform = ax.transAxes, borderpad = 0.0)
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Cloud Frequency')
cbar.set_ticks([level for level in my_levels[::2]])
plt.subplots_adjust(bottom = 0.2, wspace = 0.05, hspace = -0.20)
plt.savefig('../Ch5_Figure12.png', dpi = 250, bbox_inches = 'tight')
plt.show()

