import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import ls_rain_amt_key as rain_key
base_path = '/nerc/n02/n02/xb899100/CloudTrail/'
paths = ['Control_short/', 'Control_0800m/', 'Control_1600m/']

lsm_base = '/work/n02/n02/xb899100/island_masks/'
lsm = [Dataset(lsm_base + 'lsm50.nc', 'r').variables['lsm'][:], Dataset(lsm_base + 'lsm50_0800m.nc', 'r').variables['lsm'][:], Dataset(lsm_base + 'lsm50_1600m.nc', 'r').variables['lsm'][:]]

res = [0.1, 0.8, 1.6]
my_levels = [0.01, 0.5, 1, 2, 4, 8, 16, 32, 64]
my_colors = ['blue', 'cornflowerblue', 'olive', 'gold', 'orange', 'red', 'magenta', 'aliceblue']
target_time0 = 120.0
target_time1 = 480.0

fig = plt.figure(tight_layout = True, figsize = (9,9))
axa = fig.add_subplot(3, 1, 1, adjustable = 'box', aspect = 1)
axb = fig.add_subplot(3, 1, 2, adjustable = 'box', aspect = 1)
axc = fig.add_subplot(3, 1, 3, adjustable = 'box', aspect = 1)
ax = [axa, axb, axc]
cb = []

for path in paths:
    i = paths.index(path)
    ax[i].set_ylabel('y (km)')
    path = base_path + path
    lwp_nc = Dataset(path + '/lwp_00.nc', 'r')
    # get the index for our target time
    time_key = [tkey for tkey in lwp_nc.variables.keys() if 'accum' in tkey][0]
    time_test0 = np.abs(lwp_nc.variables[time_key][:] - target_time0)
    time_test1 = np.abs(lwp_nc.variables[time_key][:] - target_time1)
    it0 = np.where(time_test0 == np.min(time_test0))[0][0]
    it1 = np.where(time_test1 == np.min(time_test1))[0][0]
    
    # get the coordinates
    X, Y = np.meshgrid(np.arange(lwp_nc.variables[rain_key].shape[2])*res[i], np.arange(lwp_nc.variables[rain_key].shape[1])*res[i])
    rainfall = lwp_nc.variables[rain_key][it1,:,:] - lwp_nc.variables[rain_key][it0,:,:]
    cb = ax[i].contourf(X, Y, rainfall, levels = my_levels, colors = my_colors)
    CB = fig.colorbar(cb, ax = ax[i], label = 'Rainfall Accumulation (mm)')
    CB.ax.set_yticklabels([str(level) for level in my_levels])
    ax[i].contour(X, Y, lsm[i][0,0,:,:], colors = ['k'])
    ax[i].set_title('Total accumulated rainfall by peak heating for dx = ' + str(res[i]) + ' km')
    lwp_nc.close()

ax[i].set_xlabel('x (km)')
plt.savefig('../output_from_param_analysis/rainfall_comparison.png', dpi = 150)
plt.show()


