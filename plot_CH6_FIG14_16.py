import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import ctz_key, lwp_key
from scipy import ndimage
from analysis_tools import round5
import os

paths = ['/nerc/n02/n02/xb899100/CloudTrail/Control2/']
for path in paths:
    path_files = os.listdir(path)
    label = path.split('/')[-2]
    
    # Read the liquid water path data
    lwp_files = [my_file for my_file in path_files if 'lwp_' in my_file]
    lwp_files.sort()
    for lwp_file in lwp_files:
        with Dataset(path + lwp_file, 'r') as lwp_nc:
            if lwp_file == lwp_files[0]:
                lwp_data     = lwp_nc.variables[lwp_key][:]*1.
                lwp_time_key = [key for key in lwp_nc.variables.keys() if 'min' in key][0]
                lwp_times    = lwp_nc.variables[lwp_time_key][:]*1.
            else:
                lwp_data  = np.concatenate((lwp_data, lwp_nc.variables[lwp_key][:]), axis = 0)
                lwp_times = np.concatenate((lwp_times, lwp_nc.variables[lwp_time_key][:]), axis = 0)
    
    # Read the mcl data
    ctz_files = [my_file for my_file in path_files if 'zi_' in my_file]
    ctz_files.sort()
    for ctz_file in ctz_files:
        with Dataset(path + ctz_file, 'r') as ctz_nc:
            if ctz_file == ctz_files[0]:
                z = Dataset(path + 'bouy_03.nc', 'r').variables['thlev_zsea_theta'][:]*1.
                ctz_time_key = 'time'
                ctz_times = ctz_nc.variables[ctz_time_key][:]*1.
                ctz_data = ctz_nc.variables[ctz_key][:]*1.
            else:
                ctz_times = np.concatenate((ctz_times, ctz_nc.variables[ctz_time_key][:]*1.), axis = 0)
                ctz_data = np.concatenate((ctz_data, ctz_nc.variables[ctz_key][:]*1.), axis = 0)
    
    # Filter out the lwp data to only include times concurrent with mcl data
    lwp_time_idx = [idx for idx in range(len(lwp_times)) if lwp_times[idx] in ctz_times]
    lwp_data = lwp_data[lwp_time_idx,:,:]
    
    # Create a cloud mask
    cloud_mask = np.where(lwp_data > 0, 1.0, 0.0)
    
    # For each time, identify all of the clouds, and find the maximum height in each cloud
    max_cld_top_height = []
    n_cld = []
    ctz_03to06 = np.array([np.nan])
    ctz_09to12 = np.array([np.nan])
    for it in range(len(ctz_times)):
        # Use scipy.ndimage to count and label the clouds
        clouds, n_clds = ndimage.label(cloud_mask[it,:,:])
        print 'Time = ' + str(int(ctz_times[it])) + ', # Clouds = ' + str(n_clds)
        cloud_ctz = np.array([np.nanmax(np.where(clouds == cloud, ctz_data[it,:,:], np.nan)) for cloud in range(n_clds)])
        max_cld_top_height.append(cloud_ctz)
        n_cld.append(n_clds)
        
        if (3*60 <= ctz_times[it]) and (ctz_times[it] <= 6*60):
            ctz_03to06 = np.concatenate((ctz_03to06, max_cld_top_height[it]), axis = 0)
        elif (9*60 <= ctz_times[it]) and (ctz_times[it] <= 12*60):
            ctz_09to12 = np.concatenate((ctz_09to12, max_cld_top_height[it]), axis = 0)
    
    # Remove the nans from arrays
    ctz_03to06 = np.array([ob for ob in ctz_03to06 if ob == ob])
    ctz_09to12 = np.array([ob for ob in ctz_09to12 if ob == ob])
    
    before = plt.hist(ctz_03to06, bins = z, normed = True, cumulative = True, alpha = 0.5, color = 'blue')
    after = plt.hist(ctz_09to12, bins = z, normed = True, cumulative = True, alpha = 0.5, color = 'red')
    # plot histogram
    plt.xlim([400, 4500])
    plt.ylim([0, 1])
    plt.plot([400, 4500], [0.5, 0.5])
    plt.savefig('../CTZ_histograms_'+label+'.png', dpi = 150, bbox_inches = 'tight')
    plt.show()
    
    # plot difference
    plt.plot(0.5*(z[1:] + z[:-1]), (after[0] - before[0]))
    plt.ylim([-0.6, 0.6])
    plt.xlim([400, 4500])
    plt.savefig('../CTZ_histograms_shift_'+label+'.png', dpi = 150, bbox_inches = 'tight')
    plt.show()
    
    # plot boxplot timeseries
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(ctz_times, [np.nanpercentile(ob, 50) for ob in max_cld_top_height], 'k', lw = 2)
    ax.fill_between(ctz_times, [np.nanpercentile(ob, 25) for ob in max_cld_top_height], [np.nanpercentile(ob, 75) for ob in max_cld_top_height], facecolor = 'k', edgecolor = 'None', alpha = 0.5)
    ax.plot(ctz_times, [np.nanpercentile(ob, 90) for ob in max_cld_top_height], 'k', ls = '--')
    ax.plot(ctz_times, [np.nanpercentile(ob, 10) for ob in max_cld_top_height], 'k', ls = '--')
    ax.plot(ctz_times, [np.nanmax(ob) if len(ob) > 0 else np.nan for ob in max_cld_top_height], 'rx')
    ax.set_ylabel('Cloud Top Heights (m)')
    ax.set_xlabel('Time (mins)')
    ax.set_xticks(range(0, int(ctz_times.max() + 1), [inc for inc in [360, 1440] if ctz_times.max()/inc > 1][0]))
    ax.set_ylim([0, 5000])
    ax.set_xlim([0, ctz_times.max()])
    plt.savefig('../CTZ_timeseries_'+label+'.png', dpi = 150, bbox_inches = 'tight')
    plt.show()

    # plot number of clouds timeseries
    plt.plot(ctz_times, n_cld)
    plt.ylabel('Number of clouds')
    plt.xlabel('Time (mins)')
    plt.xticks(range(0, int(ctz_times.max() + 1), [inc for inc in [360, 1440] if ctz_times.max()/inc > 1][0]))
    plt.xlim([0, ctz_times.max()])
    plt.ylim([0, round5(max(n_cld)+1)])
    plt.savefig('../n_clds_'+label+'.png', dpi = 150, bbox_inches = 'tight')
    plt.show()
    
    # plot cloud spacing sqrt(A/N)
    plt.plot(ctz_times, np.sqrt(116*31.9/np.array(n_cld)))
    plt.ylabel(u'Mean Cloud Separation $\\approx \sqrt{\\frac{A}{N}}$ (km)')
    plt.xlabel('Time (mins)')
    plt.xticks(range(0, int(ctz_times.max() + 1), [inc for inc in [360, 1440] if ctz_times.max()/inc > 1][0]))
    plt.xlim([0, ctz_times.max()])
    plt.savefig('../cld_separation_'+label+'.png', dpi = 150, bbox_inches = 'tight')
    plt.show()



