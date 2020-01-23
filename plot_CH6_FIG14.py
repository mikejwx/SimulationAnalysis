import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import ctz_key, lwp_key
from scipy import ndimage
import os

"""
Four-panel plot showing the cloud top height box-plot time series
"""
resolutions = ['0200m', '0400m', '0800m', '1600m']

paths = {}
for res in resolutions:
    paths[res] = '/nerc/n02/n02/xb899100/CloudTrail/Control_' + res + '_HRIC_INV/'

# download all of the data
my_data = {}
for res in resolutions:
    path = paths[res]
    path_files = os.listdir(path)
    my_data[res] = {}
    
    # Read the liquid water path data
    lwp_files = [my_file for my_file in path_files if 'lwp_' in my_file]
    lwp_files.sort()
    for lwp_file in lwp_files:
        with Dataset(path + lwp_file, 'r') as lwp_nc:
            if lwp_file == lwp_files[0]:
                my_data[res][lwp_key]      = lwp_nc.variables[lwp_key][:]*1.
                lwp_time_key               = [key for key in lwp_nc.variables.keys() if 'min' in key][0]
                my_data[res][lwp_time_key] = lwp_nc.variables[lwp_time_key][:]*1.
            else:
                my_data[res][lwp_key]  = np.concatenate((my_data[res][lwp_key][:], lwp_nc.variables[lwp_key][:]), axis = 0)
                my_data[res][lwp_time_key] = np.concatenate((my_data[res][lwp_time_key][:], lwp_nc.variables[lwp_time_key][:]), axis = 0)
    
    # Read the mcl data
    ctz_files = [my_file for my_file in path_files if 'zi_' in my_file]
    ctz_files.sort()
    for ctz_file in ctz_files:
        with Dataset(path + ctz_file, 'r') as ctz_nc:
            if ctz_file == ctz_files[0]:
                my_data[res]['z'] = Dataset(path + 'bouy_03.nc', 'r').variables['thlev_zsea_theta'][:]*1.
                ctz_time_key = 'time'
                my_data[res][ctz_time_key] = ctz_nc.variables[ctz_time_key][:]*1.
                my_data[res][ctz_key] = ctz_nc.variables[ctz_key][:]*1.
            else:
                my_data[res][ctz_time_key] = np.concatenate((my_data[res][ctz_time_key][:], ctz_nc.variables[ctz_time_key][:]*1.), axis = 0)
                my_data[res][ctz_key] = np.concatenate((my_data[res][ctz_key][:], ctz_nc.variables[ctz_key][:]*1.), axis = 0)
    
    # Filter out the lwp data to only include times concurrent with mcl data
    lwp_time_idx = [idx for idx in range(len(my_data[res][lwp_time_key])) if my_data[res][lwp_time_key][idx] in my_data[res][ctz_time_key]]
    my_data[res][lwp_key] = my_data[res][lwp_key][lwp_time_idx,:,:]
    
    # Create a cloud mask
    cloud_mask = np.where(my_data[res][lwp_key] > 0, 1.0, 0.0)
    
    # For each time, identify all of the clouds, and find the maximum height in each cloud
    max_cld_top_height = []
    n_cld = []
    for it in range(len(my_data[res][ctz_time_key])):
        # Use scipy.ndimage to count and label the clouds
        clouds, n_clds = ndimage.label(cloud_mask[it,:,:])
        print 'Time = ' + str(int(my_data[res][ctz_time_key][it])) + ', # Clouds = ' + str(n_clds)
        cloud_ctz = np.array([np.nanmax(np.where(clouds == cloud, my_data[res][ctz_key][it,:,:], np.nan)) for cloud in range(n_clds)])
        max_cld_top_height.append(cloud_ctz)
        n_cld.append(n_clds)
    
    my_data[res]['max_ctz'] = np.array(max_cld_top_height)
    my_data[res][ctz_time_key] /= 60.

letter = ['a)', 'b)', 'c)', 'd)']
# plot boxplot timeseries
fig = plt.figure(figsize=(10,6), tight_layout = True)
for idx in range(4):
    res = resolutions[idx]
    axa = fig.add_subplot(2, 2, idx+1)
    # Plot the median
    axa.plot(my_data[resolutions[idx]][ctz_time_key], [np.nanpercentile(ob, 50) for ob in my_data[resolutions[idx]]['max_ctz']], 'k', lw = 2)
    # Plot the interquartile range
    axa.fill_between(my_data[resolutions[idx]][ctz_time_key], [np.nanpercentile(ob, 25) for ob in my_data[resolutions[idx]]['max_ctz']], [np.nanpercentile(ob, 75) for ob in my_data[resolutions[idx]]['max_ctz']], facecolor = 'k', edgecolor = 'None', alpha = 0.5)
    # Plot the 90th percentile
    axa.plot(my_data[resolutions[idx]][ctz_time_key], [np.nanpercentile(ob, 90) for ob in my_data[resolutions[idx]]['max_ctz']], 'k', ls = '--')
    # Plot the 10th percentile
    axa.plot(my_data[resolutions[idx]][ctz_time_key], [np.nanpercentile(ob, 10) for ob in my_data[resolutions[idx]]['max_ctz']], 'k', ls = '--')
    # Plot the maximum
    axa.plot(my_data[resolutions[idx]][ctz_time_key], [np.nanmax(ob) if len(ob) > 0 else np.nan for ob in my_data[resolutions[idx]]['max_ctz']], 'rx')
    # Panel formatting
    if idx in [0, 2]:
        axa.set_ylabel('z$_{top}$ (m)')
    else:
        axa.set_yticklabels([''])
    if idx in [2, 3]:
        axa.set_xlabel('Time (hrs)')
    else:
        axa.set_xticklabels([''])
    axa.set_xticks(range(0, int(my_data[res][ctz_time_key].max() + 1), [inc for inc in [6, 24] if my_data[res][ctz_time_key].max()/inc > 1][0]))
    axa.set_ylim([0, 4000])
    axa.set_xlim([0, my_data[res][ctz_time_key].max()])
    axa.set_title(letter[idx] + ' DX' + res[:-1])

plt.savefig('../Ch6_Figure14.png', dpi = 250, bbox_inches = 'tight')
plt.show()



