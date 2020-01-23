import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import fromComponents, bilinear_interpolation, get_cs_coords, lcl as get_lcl, downwind_rectangle
from scipy import integrate, interpolate
from datetime import datetime as dt
from STASH_keys import rho_key, u_key, v_key, w_key, pthe_key, q_key, theta_key, ls_rain_amt_key, lwp_key
from SkewT_archer import PTtoTemp

# you need to re-run all of the wind experiments so that they all output precip
# in the meantime, you can plot the cloud mask, liquid water path, or cloud top heights

paths = {'U05' : '/nerc/n02/n02/xb899100/CloudTrail/U05/',
         'U06' : '/nerc/n02/n02/xb899100/CloudTrail/U06/',
         'U07' : '/nerc/n02/n02/xb899100/CloudTrail/U07/',
         'U08' : '/nerc/n02/n02/xb899100/CloudTrail/U08/',
         'U09' : '/nerc/n02/n02/xb899100/CloudTrail/U09/',
         'U10' : '/nerc/n02/n02/xb899100/CloudTrail/Control_short/'}

keys = ['U05', 'U06', 'U07', 'U08', 'U09', 'U10']
my_levels = [0.01, 0.5, 1., 2., 4., 8., 16., 32.]

fig, axes = plt.subplots(nrows = 3, ncols = 2)
for key in keys:
    print key
    # open the lwp_00.nc file
    lwp_nc = Dataset(paths[key] + 'lwp_00.nc', 'r')
    time_key = [tkey for tkey in lwp_nc.variables.keys() if 'accum' in tkey][0]
    times = lwp_nc.variables[time_key][:]*1.
    if key in ['U05']:
        times += 240.0
    idx = [it for it in range(len(times)) if times[it] == 1080.0][0]
    precipitation = lwp_nc.variables[ls_rain_amt_key][idx,:,:]
    lwp_nc.close()
    
    X, Y = np.meshgrid(np.arange(precipitation.shape[1])*100.0, np.arange(precipitation.shape[0])*100.0)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['U05', 'U10']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)
    ax = axes.flat[keys.index(key)]
    ax.set_adjustable('box')
    ax.set_aspect(1)
    im = ax.contourf(X/1000., Y/1000., precipitation, colors = ['blue', 'cornflowerblue', 'olive', 'gold', 'orange', 'red', 'magenta', 'ghostwhite'], levels = my_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax.set_title(key + ' at P$_{max}$=' + str(round(np.max(precipitation), 2)) + 'mm')
    if (keys.index(key) + 1) in [1, 2, 3, 4]:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    if (keys.index(key) + 1) in [2, 4, 6]:
        ax.set_yticklabels([''])
    else:
        ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.1, 0.1, 0.85, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Total Precipitation (mm)')
cbar.ax.set_xticklabels([str(level) for level in my_levels])
plt.subplots_adjust(left = 0.1, bottom = 0.2, right = 0.95, top = 0.9, wspace = 0.1, hspace = 0.1)
plt.savefig('../wind_comparison_precipitation.png', dpi = 150, bbox_inches = 'tight')
plt.show()

crash

fig, axes = plt.subplots(nrows = 3, ncols = 2)
for key in keys:
    print key
    # open the lwp_00.nc file
    lwp_nc = Dataset(paths[key] + 'lwp_00.nc', 'r')
    time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
    times = lwp_nc.variables[time_key][:]*1.
    if key == 'U05':
        times += 240.0
    idx = [it for it in range(len(times)) if (360.0 <= times[it]) and (times[it] <= 1080.0)]
    cloud_frequency = np.nanmean(np.where(lwp_nc.variables[lwp_key][idx,:,:] > 0, 1.0, 0.0), axis = 0)
    lwp_nc.close()
    
    X, Y = np.meshgrid(np.arange(cloud_frequency.shape[1])*100.0, np.arange(cloud_frequency.shape[0])*100.0)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['U05', 'Control']:
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
plt.savefig('../wind_comparison_cloud_frequency.png', dpi = 150, bbox_inches = 'tight')
plt.show()


