import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import ctz_key, lwp_key
from scipy import ndimage

#paths = ['/nerc/n02/n02/xb899100/CloudTrail/U05/', '/nerc/n02/n02/xb899100/CloudTrail/Control/', 
paths = ['/nerc/n02/n02/xb899100/CloudTrail/Control/']#['/nerc/n02/n02/xb899100/CloudTrail/Control_0800m/', '/work/n02/n02/xb899100/cylc-run/u-bn821/share/data/history/']

for path in paths:
    # Read the liquid water path data
    lwp_nc = Dataset(path + 'lwp_00.nc', 'r')
    lwp_data = lwp_nc.variables[lwp_key][:]*1.
    lwp_time_key = [key for key in lwp_nc.variables.keys() if 'min' in key][0]
    lwp_times = lwp_nc.variables[lwp_time_key][:]*1.
    lwp_nc.close()
    
    # Read the mcl data
    hours = ["{0:02d}".format(hour) for hour in range(0, 13, 3)]#[4 if 'nerc' in path else 3][0])]
    label = 'Control'#['Control_0800m' if 'nerc' in path else 'Control_0800m_INV'][0]

    for hour in hours:
        ctz_nc = Dataset(path + 'zi_' + hour + '.nc', 'r')
        if hour == hours[0]:
            z = Dataset(path + 'bouy_00.nc', 'r').variables['thlev_zsea_theta'][:]*1.
            ctz_time_key = 'time'#[key for key in ctz_nc.variables.keys() if 'min' in key][0]
            ctz_times = ctz_nc.variables[ctz_time_key][:]*1.
            ctz_data = ctz_nc.variables[ctz_key][:]*1.
        else:
            ctz_times = np.concatenate((ctz_times, ctz_nc.variables[ctz_time_key][:]*1.), axis = 0)
            ctz_data = np.concatenate((ctz_data, ctz_nc.variables[ctz_key][:]*1.), axis = 0)
        
        ctz_nc.close()

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

    for it in range(len(ctz_times)):
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
    plt.plot(ctz_times, [np.nanpercentile(ob, 50) for ob in max_cld_top_height], 'k', lw = 2)
    plt.fill_between(ctz_times, [np.nanpercentile(ob, 25) for ob in max_cld_top_height], [np.nanpercentile(ob, 75) for ob in max_cld_top_height], facecolor = 'k', edgecolor = 'None', alpha = 0.5)
    plt.plot(ctz_times, [np.nanpercentile(ob, 90) for ob in max_cld_top_height], 'k', ls = '--')
    plt.plot(ctz_times, [np.nanpercentile(ob, 10) for ob in max_cld_top_height], 'k', ls = '--')
    plt.plot(ctz_times, [np.nanmax(ob) if len(ob) > 0 else np.nan for ob in max_cld_top_height], 'rx')
    plt.ylabel('Cloud Top Heights (m)')
    plt.xlabel('Time (mins)')
    plt.xticks(range(0, 901, 60))
    plt.ylim([0, 10000])
    plt.xlim([0, 900])
    plt.savefig('../CTZ_timeseries_'+label+'.png', dpi = 150, bbox_inches = 'tight')
    plt.show()

    # plot number of clouds timeseries
    plt.plot(ctz_times, n_cld)
    plt.ylabel('Number of clouds')
    plt.xlabel('Time (mins)')
    plt.xticks(range(0, 901, 60))
    plt.xlim([0, 900])
    plt.ylim([0, 2000])
    plt.savefig('../n_clds_'+label+'.png', dpi = 150, bbox_inches = 'tight')
    plt.show()




