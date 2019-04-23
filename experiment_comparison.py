import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import ndimage, signal
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Read in a land-sea mask
lsm = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
lsm_data = lsm.variables['lsm'][0,0,:,:]*1.
lsm.close()

path_keys = ['H125', 'H375', 'RH_BLm25', 'RH_FAm25', 'Control']
paths = {'H125'     : '/nerc/n02/n02/xb899100/CloudTrail/H125/',
         'H375'     : '/nerc/n02/n02/xb899100/CloudTrail/H375/',
         'RH_BLm25' : '/work/n02/n02/xb899100/cylc-run/u-bg665/share/data/history/',
         'RH_FAm25' : '/work/n02/n02/xb899100/cylc-run/u-bg933/share/data/history/',
         'Control'  : '/nerc/n02/n02/xb899100/CloudTrail/Control/'}

theta_key = u'STASH_m01s00i004'
lwp_key = u'STASH_m01s30i405'
n_key = u'Flow-Perpendicular'

data_dict = {'theta' : {},
             'wind'  : {},
             'cloud' : {}}

for path_key in path_keys:
    print path_key
    path = paths[path_key]
    # Read in the netCDF
    # Need the theta contour at noon (this is in hour 09 for H125, and H375, but hour 08 for Control and RH experiments)
    if path_key in ['H125', 'H375', 'Control']:
        hour = '09'
        heat_start = 360.
        heat_end = 1080.
        heat_peak = 720.
    else:
        hour = '08'
        heat_start = 120.
        heat_end = 840.
        heat_peak = 480.
    
    bouy_nc = Dataset(path + 'bouy_' + hour + '.nc', 'r')
    wind_nc = Dataset(path + 'wind_' + hour + '.nc', 'r')
    lwp_nc = Dataset(path + 'lwp_00.nc', 'r')
    
    # Grab the height and time dimensions
    height_key = [key for key in bouy_nc.variables.keys() if 'zsea' in key][0]
    
    z = bouy_nc.variables[height_key][:]*1.
    
    time_key_bouy = [key for key in bouy_nc.variables.keys() if 'min' in key][0]
    time_key_wind = [key for key in wind_nc.variables.keys() if 'min' in key][0]
    time_key_lwp = [key for key in lwp_nc.variables.keys() if 'min' in key][0]
    
    time_bouy = bouy_nc.variables[time_key_bouy][:]*1.
    time_wind = wind_nc.variables[time_key_wind][:]*1.
    time_lwp = lwp_nc.variables[time_key_lwp][:]*1.
    
    # Want to plot the theta contour at 350 m
    iz350 = np.where(np.abs(z - 350.) == np.min(np.abs(z - 350.)))[0][0]
    
    # Want to plot the wind at 10 m
    iz010 = np.where(np.abs(z - 10.) == np.min(np.abs(z - 10.)))[0][0]
    
    # Want to plot at 720 minutes
    itPb = np.where(np.abs(time_bouy - heat_peak) == np.min(np.abs(time_bouy - heat_peak)))[0][0]
    itPw = np.where(np.abs(time_wind - heat_peak) == np.min(np.abs(time_wind - heat_peak)))[0][0]
    
    # Want to get cloud frequency for the whole heating period
    it_heating = [it for it in xrange(len(time_lwp)) if heat_start <= time_lwp[it] <= heat_end]
    
    # Get the theta contour
    theta_data = bouy_nc.variables[theta_key][itPb, iz350, :, :]
    theta_anom = theta_data - np.nanmean(theta_data)
    
    # Define a smoothing structure, circular filter
    smoothing_diam = 6000.
    n_r, n_c = np.array([smoothing_diam, smoothing_diam])/100
    w_c = np.ones((n_r, n_c))
    r, c = np.meshgrid(np.arange(n_r), np.arange(n_c))
    w_c = np.where((np.sqrt((r - n_r/2)**2 + (c - n_c/2)**2) <= n_r/2), w_c, 0.)
    w_c /= np.sum(w_c)
    
    # Smooth theta_anom
    theta_anom_sm = signal.convolve2d(theta_anom, w_c, mode = 'same', boundary = 'wrap')
    
    # Get the winds
    n_data = wind_nc.variables[n_key][itPw, iz010, :, :]
    
    # Get the cloud frequency
    cloud_freq = np.nanmean(np.where((lwp_nc.variables[lwp_key][it_heating,:,:] > 1e-16), 1., 0.), axis = 0)
    
    # Store to data_dict
    data_dict['theta'][path_key] = theta_anom_sm
    data_dict['wind'][path_key]  = n_data
    data_dict['cloud'][path_key] = cloud_freq
    bouy_nc.close()
    wind_nc.close()
    lwp_nc.close()

# Create coordinate system
X, Y = np.meshgrid(np.arange(0, 116., 0.1), np.arange(0., 31.9, 0.1))

# Create plotting args
wind_levels = np.array([-2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8])
wind_lgd = [-2.4, -1.2, -0.6, 0.0, 0.6, 1.2, 2.4]
cloud_levels = [0.25, 1.0]
fs = 12 # Fontsize
fmt = '%1.1f'

# Create title dictionary
title_dict = {'Control' : 'Control, $\\beta = 1.0$',
              'H125' : '$\\beta = 1/3$',
              'H375' : '$\\beta = 3.0$',
              'RH_BLm25' : 'BL RH minus 25%',
              'RH_FAm25' : 'FA RH minus 25%'}

fig = plt.figure(tight_layout = True, figsize = (20, 10))
for path_key in path_keys:
    ax = fig.add_subplot(3, 2, path_keys.index(path_key) + 1, adjustable = 'box', aspect = 1)
    WND = ax.contourf(X, Y, data_dict['wind'][path_key], cmap = 'bwr', levels = wind_levels, extend = 'both')
    ax.contourf(X, Y, data_dict['cloud'][path_key], colors = ['k'], levels = cloud_levels, alpha = 2./3.)
    TA = ax.contour(X, Y, data_dict['theta'][path_key], colors = ['k'], linewidths = [2], levels = [0.1])
    ax.contour(X, Y, lsm_data, levels = [0, 1e-16], colors = 'darkred', linewidths = [2])
    #ax.clabel(TA, fontsize = fs, fmt = fmt)
    if path_key == 'Control':
        axins = inset_axes(ax, width = "2.5%", height = "100%", loc = 7, bbox_to_anchor = (0.05, 0.0, 1, 1), bbox_transform = ax.transAxes, borderpad = 0.)
        plt.colorbar(WND, cax = axins, label = 'Flow-perpendicular wind (m/s)', ticks = wind_lgd)
    ax.set_ylabel('y (km)')
    ax.set_xlabel('x (km)')
    ax.set_title(title_dict[path_key], fontsize = fs+4)

plt.savefig('../Experiment_Comparison.png', dpi = 500)
plt.show()


