import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import fromComponents, bilinear_interpolation, get_cs_coords, lcl as get_lcl, downwind_rectangle
from scipy import integrate, interpolate
from datetime import datetime as dt
from STASH_keys import rho_key, u_key, v_key, w_key, pthe_key, q_key, theta_key, ls_rain_amt_key, lwp_key, ctz_key
from SkewT_archer import PTtoTemp
import os

# you need to re-run all of the wind experiments so that they all output precip
# in the meantime, you can plot the cloud mask, liquid water path, or cloud top heights

paths = {'DX0100C5' : '/nerc/n02/n02/xb899100/CloudTrail/Control2/',
         'DX0200' : '/nerc/n02/n02/xb899100/CloudTrail/Control_0200m_HRIC_INV/',
         'DX0400' : '/nerc/n02/n02/xb899100/CloudTrail/Control_0400m_HRIC_INV/',
         'DX0800' : '/nerc/n02/n02/xb899100/CloudTrail/Control_0800m_HRIC_INV/',
         'DX1600' : '/nerc/n02/n02/xb899100/CloudTrail/Control_1600m_HRIC_INV/'}

data_dict = {}
keys = ['DX0100C5', 'DX0200', 'DX0400', 'DX0800', 'DX1600']
my_precip_levels = [0.01, 0.5, 1., 2., 4., 8., 16., 32.]
my_precip_colors = ['blue', 'cornflowerblue', 'olive', 'gold', 'orange', 'red', 'magenta', 'ghostwhite']

### get the precipitation data and lwp data fields ###
for key in keys:
    print key
    data_dict[key] = {}
    # precipitation
    lwp_nc = Dataset(paths[key] + 'lwp_00.nc', 'r')
    time_key = [tkey for tkey in lwp_nc.variables.keys() if 'accum' in tkey][0]
    data_dict[key]['precip_times'] = lwp_nc.variables[time_key][:]*1.
    data_dict[key]['precip'] = lwp_nc.variables[ls_rain_amt_key][:]
    # cloud masks for cloud frequency
    time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
    data_dict[key]['lwp_times'] = lwp_nc.variables[time_key][:]*1.
    data_dict[key]['cloud_mask'] = np.where(lwp_nc.variables[lwp_key][:] > 0, 1.0, 0.0)
    lwp_nc.close()
    # cloud top heights
    # get a list of and open the zi_%H.nc files
    zi_files = [ncfile for ncfile in os.listdir(paths[key]) if 'zi_' in ncfile]
    zi_files.sort()
    for ncfile in zi_files:
        with Dataset(paths[key] + ncfile, 'r') as zi_nc:
            if ncfile == zi_files[0]:
                time_key = [tkey for tkey in zi_nc.variables.keys() if 'time' in tkey][0]
                data_dict[key]['zi_times'] = zi_nc.variables[time_key][:]*1.
                data_dict[key]['ctz']      = zi_nc.variables[ctz_key][:]
            else:
                data_dict[key]['zi_times'] = np.concatenate((data_dict[key]['zi_times'], zi_nc.variables[time_key][:]*1.), axis = 0)
                data_dict[key]['ctz']      = np.concatenate((data_dict[key]['ctz'], zi_nc.variables[ctz_key][:]), axis = 0)

################################################################################
# Hours 0 - 6 #
################################################################################
### plot the precipitation data fields ###
"""
fig = plt.figure(figsize = (6, 10))
for key in keys:
    idx0 = [it for it in range(len(data_dict[key]['precip_times'])) if data_dict[key]['precip_times'][it] == 360.0][0]
    precip = data_dict[key]['precip'][idx0,:,:]
    grid_spacing = float(key[2:])
    X, Y = np.meshgrid(np.arange(precip.shape[1])*grid_spacing, np.arange(precip.shape[0])*grid_spacing)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['DX0100']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X - x_c)**2 + (Y - y_c)**2)
    ax = fig.add_subplot(5, 1, keys.index(key)+1, adjustable = 'box', aspect = 1)
    im = ax.contourf(X/1000., Y/1000., precip, colors = my_precip_colors, levels = my_precip_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax.set_title(key + ' at P$_{max}$ =' + str(round(np.max(precip), 2)) + 'mm')
    if (keys.index(key) + 1) in [1, 2, 3, 4]:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.12, 0.1, 0.78, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Total Morning Precipitation (mm)')
cbar.ax.set_xticklabels([str(level) for level in my_precip_levels])
plt.subplots_adjust(bottom = 0.17, wspace = 0.17, hspace = 0.30)
plt.savefig('../DX_comparison_morning_precipitation.png', dpi = 150, bbox_inches = 'tight')
plt.show()
"""
### plot the cloud frequency data ###
letter = {'DX0100C5' : 'a)', 'DX0200' : 'b)', 'DX0400' : 'c)', 'DX0800' : 'd)', 'DX1600' : 'e)'}
my_cldfreq_levels = np.arange(0, 0.61, 0.075)
fig = plt.figure(figsize = (6, 10))
for key in keys:
    idx = [it for it in range(len(data_dict[key]['lwp_times'])) if (360.0 <= data_dict[key]['lwp_times'][it]) and (data_dict[key]['lwp_times'][it] <= 1080.0)]
    cloud_frequency = np.nanmean(data_dict[key]['cloud_mask'][idx,:,:], axis = 0)
    grid_spacing = float(key[2:6])
    X, Y = np.meshgrid(np.arange(cloud_frequency.shape[1])*grid_spacing, np.arange(cloud_frequency.shape[0])*grid_spacing)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['DX0100C5']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)
    # add the panel
    ax = fig.add_subplot(5, 1, keys.index(key) + 1, adjustable = 'box', aspect = 1)
    im = ax.contourf(X/1000., Y/1000., cloud_frequency, cmap = 'Greys_r', levels = my_cldfreq_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax.set_title(letter[key] + ' ' + key)
    if (keys.index(key) + 1) in range(5):
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.12, 0.085, 0.8, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Daytime Cloud Frequency')
cbar.set_ticks(my_cldfreq_levels[::2])
cbar.set_ticklabels([str(level) for level in my_cldfreq_levels[::2]])
plt.subplots_adjust(bottom = 0.17, wspace = 0.17, hspace = 0.30)
plt.savefig('../Ch6_Figure04.png', dpi = 250, bbox_inches = 'tight')
plt.show()

"""
### Mean cloud top heights? ###
my_ctz_levels = np.arange(0., 3000.1, 250.)
fig = plt.figure(figsize = (6, 10))
for key in keys:
    print key                                                               
    idx = [it for it in range(len(data_dict[key]['zi_times'])) if (360.0 >= data_dict[key]['zi_times'][it])]
    ctz_data = np.nanmean(data_dict[key]['ctz'][idx,:,:], axis = 0)
    ctz_data[np.isnan(ctz_data)] = 0.0
    grid_spacing = float(key[2:])
    
    X, Y = np.meshgrid(np.arange(ctz_data.shape[1])*grid_spacing, np.arange(ctz_data.shape[0])*grid_spacing)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['DX0100']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)
    # add the panel
    ax = fig.add_subplot(5, 1, keys.index(key) + 1, adjustable = 'box', aspect = 1)
    im = ax.contourf(X/1000., Y/1000., ctz_data, cmap = 'hot', levels = my_ctz_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['k'], linewidths = [2])
    ax.set_title(key)
    if (keys.index(key) + 1) in range(5):
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.12, 0.1, 0.78, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Morning Mean Cloud Top Height (m)')
cbar.ax.set_xticklabels([str(int(level)) for level in my_ctz_levels[::2]])
plt.subplots_adjust(bottom = 0.17, wspace = 0.17, hspace = 0.30)
plt.savefig('../DX_comparison_morning_ctz.png', dpi = 150, bbox_inches = 'tight')
plt.show()

################################################################################
# Hours 6 - 18 #
################################################################################
### plot the precipitation data fields ###
fig = plt.figure(figsize = (6, 10))
for key in keys:
    idx0 = [it for it in range(len(data_dict[key]['precip_times'])) if data_dict[key]['precip_times'][it] == 360.0][0]
    idx1 = [it for it in range(len(data_dict[key]['precip_times'])) if data_dict[key]['precip_times'][it] == 1080.0][0]
    precip = data_dict[key]['precip'][idx1,:,:] - data_dict[key]['precip'][idx0,:,:]
    grid_spacing = float(key[2:])
    X, Y = np.meshgrid(np.arange(precip.shape[1])*grid_spacing, np.arange(precip.shape[0])*grid_spacing)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['DX0100']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X - x_c)**2 + (Y - y_c)**2)
    ax = fig.add_subplot(5, 1, keys.index(key)+1, adjustable = 'box', aspect = 1)
    im = ax.contourf(X/1000., Y/1000., precip, colors = my_precip_colors, levels = my_precip_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax.set_title(key + ' at P$_{max}$ =' + str(round(np.max(precip), 2)) + 'mm')
    if (keys.index(key) + 1) in [1, 2, 3, 4]:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.12, 0.1, 0.78, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Total Daytime Precipitation (mm)')
cbar.ax.set_xticklabels([str(level) for level in my_precip_levels])
plt.subplots_adjust(bottom = 0.17, wspace = 0.17, hspace = 0.30)
plt.savefig('../DX_comparison_daytime_precipitation.png', dpi = 150, bbox_inches = 'tight')
plt.show()

### plot the cloud frequency data ###
fig = plt.figure(figsize = (6, 10))
for key in keys:
    idx = [it for it in range(len(data_dict[key]['lwp_times'])) if (360.0 <= data_dict[key]['lwp_times'][it]) and (data_dict[key]['lwp_times'][it] <= 1080.0)]
    cloud_frequency = np.nanmean(data_dict[key]['cloud_mask'][idx,:,:], axis = 0)
    grid_spacing = float(key[2:])
    X, Y = np.meshgrid(np.arange(cloud_frequency.shape[1])*grid_spacing, np.arange(cloud_frequency.shape[0])*grid_spacing)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['DX0100']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)
    # add the panel
    ax = fig.add_subplot(5, 1, keys.index(key) + 1, adjustable = 'box', aspect = 1)
    im = ax.contourf(X/1000., Y/1000., cloud_frequency, cmap = 'Greys_r', levels = my_cldfreq_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax.set_title(key)
    if (keys.index(key) + 1) in range(5):
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.12, 0.1, 0.78, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Daytime Cloud Frequency')
cbar.ax.set_xticklabels([str(level) for level in my_cldfreq_levels])
plt.subplots_adjust(bottom = 0.17, wspace = 0.17, hspace = 0.30)
plt.savefig('../DX_comparison_daytime_cloud_frequency.png', dpi = 150, bbox_inches = 'tight')
plt.show()

### Mean cloud top heights? ###
fig = plt.figure(figsize = (6, 10))
for key in keys:
    print key
    idx = [it for it in range(len(data_dict[key]['zi_times'])) if (360.0 <= data_dict[key]['zi_times'][it]) and (data_dict[key]['zi_times'][it] <= 1080.0)]
    ctz_data = np.nanmean(data_dict[key]['ctz'][idx,:,:], axis = 0)
    ctz_data[np.isnan(ctz_data)] = 0.0
    grid_spacing = float(key[2:])
    
    X, Y = np.meshgrid(np.arange(ctz_data.shape[1])*grid_spacing, np.arange(ctz_data.shape[0])*grid_spacing)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['DX0100']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)
    # add the panel
    ax = fig.add_subplot(5, 1, keys.index(key) + 1, adjustable = 'box', aspect = 1)
    im = ax.contourf(X/1000., Y/1000., ctz_data, cmap = 'hot', levels = my_ctz_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['k'], linewidths = [2])
    ax.set_title(key)
    if (keys.index(key) + 1) in range(5):
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.12, 0.1, 0.78, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Daytime Mean Cloud Top Height (m)')
cbar.ax.set_xticklabels([str(int(level)) for level in my_ctz_levels[::2]])
plt.subplots_adjust(bottom = 0.17, wspace = 0.17, hspace = 0.30)
plt.savefig('../DX_comparison_daytime_ctz.png', dpi = 150, bbox_inches = 'tight')
plt.show()

################################################################################
# Hours 18 - 24 #
################################################################################
### plot the precipitation data fields ###
fig = plt.figure(figsize = (6, 10))
for key in keys:
    idx1 = [it for it in range(len(data_dict[key]['precip_times'])) if data_dict[key]['precip_times'][it] == 1080.0][0]
    precip = data_dict[key]['precip'][-1,:,:] - data_dict[key]['precip'][idx1,:,:]
    grid_spacing = float(key[2:])
    X, Y = np.meshgrid(np.arange(precip.shape[1])*grid_spacing, np.arange(precip.shape[0])*grid_spacing)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['DX0100']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X - x_c)**2 + (Y - y_c)**2)
    ax = fig.add_subplot(5, 1, keys.index(key)+1, adjustable = 'box', aspect = 1)
    im = ax.contourf(X/1000., Y/1000., precip, colors = my_precip_colors, levels = my_precip_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax.set_title(key + ' at P$_{max}$ =' + str(round(np.max(precip), 2)) + 'mm')
    if (keys.index(key) + 1) in [1, 2, 3, 4]:
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.12, 0.1, 0.78, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Total Evening Precipitation (mm)')
cbar.ax.set_xticklabels([str(level) for level in my_precip_levels])
plt.subplots_adjust(bottom = 0.17, wspace = 0.17, hspace = 0.30)
plt.savefig('../DX_comparison_evening_precipitation.png', dpi = 150, bbox_inches = 'tight')
plt.show()

### plot the cloud frequency data ###
fig = plt.figure(figsize = (6, 10))
for key in keys:
    idx = [it for it in range(len(data_dict[key]['lwp_times'])) if (data_dict[key]['lwp_times'][it] >= 1080.0)]
    cloud_frequency = np.nanmean(data_dict[key]['cloud_mask'][idx,:,:], axis = 0)
    grid_spacing = float(key[2:])
    X, Y = np.meshgrid(np.arange(cloud_frequency.shape[1])*grid_spacing, np.arange(cloud_frequency.shape[0])*grid_spacing)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['DX0100']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)
    # add the panel
    ax = fig.add_subplot(5, 1, keys.index(key) + 1, adjustable = 'box', aspect = 1)
    im = ax.contourf(X/1000., Y/1000., cloud_frequency, cmap = 'Greys_r', levels = my_cldfreq_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['r'])
    ax.set_title(key)
    if (keys.index(key) + 1) in range(5):
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.12, 0.1, 0.78, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Evening Cloud Frequency')
cbar.ax.set_xticklabels([str(level) for level in my_cldfreq_levels])
plt.subplots_adjust(bottom = 0.17, wspace = 0.17, hspace = 0.30)
plt.savefig('../DX_comparison_evening_cloud_frequency.png', dpi = 150, bbox_inches = 'tight')
plt.show()

### Mean cloud top heights? ###
fig = plt.figure(figsize = (6, 10))
for key in keys:
    print key
    idx = [it for it in range(len(data_dict[key]['zi_times'])) if (data_dict[key]['zi_times'][it] >= 1080.0)]
    ctz_data = np.nanmean(data_dict[key]['ctz'][idx,:,:], axis = 0)
    ctz_data[np.isnan(ctz_data)] = 0.0
    grid_spacing = float(key[2:])
    
    X, Y = np.meshgrid(np.arange(ctz_data.shape[1])*grid_spacing, np.arange(ctz_data.shape[0])*grid_spacing)
    R_i = 1000*np.sqrt(50.0/np.pi)
    if key in ['DX0100']:
        x_c = 100000.0 + R_i
        y_c = 4*R_i
    else:
        x_c = 108000.0
        y_c = 16000.0
    
    R = np.sqrt((X-x_c)**2 + (Y-y_c)**2)
    # add the panel
    ax = fig.add_subplot(5, 1, keys.index(key) + 1, adjustable = 'box', aspect = 1)
    im = ax.contourf(X/1000., Y/1000., ctz_data, cmap = 'hot', levels = my_ctz_levels, extend = 'max')
    ax.contour(X/1000., Y/1000., R, levels = [R_i], colors = ['k'], linewidths = [2])
    ax.set_title(key)
    if (keys.index(key) + 1) in range(5):
        ax.set_xticklabels([''])
    else:
        ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')

cbar_ax = fig.add_axes([0.12, 0.1, 0.78, 0.02])
cbar = fig.colorbar(im, cax = cbar_ax, orientation = 'horizontal', label = 'Evening Mean Cloud Top Height (m)')
cbar.ax.set_xticklabels([str(int(level)) for level in my_ctz_levels[::2]])
plt.subplots_adjust(bottom = 0.17, wspace = 0.17, hspace = 0.30)
plt.savefig('../DX_comparison_evening_ctz.png', dpi = 150, bbox_inches = 'tight')
plt.show()
"""
