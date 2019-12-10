import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import lwp_key, mcl_key, w_key, theta_key, q_key
from scipy import ndimage
from datetime import datetime as dt

### Read the data ###
path_start = '/nerc/n02/n02/xb899100/CloudTrail/Control_'
path_end   = '_HRIC_INV/'
experiments = ['0200m', '0400m', '0800m', '1600m']
hours = ["{0:02d}".format(hour) for hour in range(0, 24, 3)]

my_data = {}
for exp in experiments:
    print '[' + dt.now().strftime('%H:%M:%S') + '] Starting experiment: ' + exp
    my_data[exp] = {}
    dx = float(exp[:-1])/1000.
    with Dataset(path_start + exp + path_end + 'lwp_00.nc', 'r') as lwp_nc:
        print '[' + dt.now().strftime('%H:%M:%S') + '] Reading the liquid water path data'
        my_data[exp][lwp_key] = lwp_nc.variables[lwp_key][:]*1.
        lwp_time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
        my_data[exp]['lwp_times'] = lwp_nc.variables[lwp_time_key][:]*1.
        my_data[exp]['x'], my_data[exp]['y'] = np.meshgrid(np.arange(my_data[exp][lwp_key].shape[2])*dx, np.arange(my_data[exp][lwp_key].shape[1])*dx)
    print '[' + dt.now().strftime('%H:%M:%S') + '] Reading the cloud liquid water content data'
    for hour in hours:
        with Dataset(path_start + exp + path_end + 'mr_' + hour + '.nc', 'r') as mr_nc:
            if hour == hours[0]:
                my_data[exp][mcl_key] = mr_nc.variables[mcl_key][:]*1.
                z = mr_nc.variables['thlev_zsea_theta'][:]*1.
                mcl_time_key = [tkey for tkey in mr_nc.variables.keys() if 'min' in tkey][0]
                my_data[exp]['mcl_times'] = mr_nc.variables[mcl_time_key][:]*1.
            else:
                my_data[exp][mcl_key] = np.concatenate((my_data[exp][mcl_key], mr_nc.variables[mcl_key][:]*1.), axis = 0)
                my_data[exp]['mcl_times'] = np.concatenate((my_data[exp]['mcl_times'], mr_nc[mcl_time_key][:]*1.), axis = 0)
    
    # define the cloud masks
    my_data[exp]['sfc_mask'] = np.where(my_data[exp][lwp_key] > 0, 1.0, 0.0)
    my_data[exp]['mask'] = np.where(my_data[exp][mcl_key] > 0, 1.0, 0.0)
    
    ### Get statistics for the bulk cloud foot print ###
    print '[' + dt.now().strftime('%H:%M:%S') + '] Computing statistics about the bulk cloud footprints'
    N_scale = float(my_data[exp][lwp_key][0,:,:].size)
    my_data[exp]['n_clds'] = np.empty_like(my_data[exp]['lwp_times'])
    my_data[exp]['A_cld']  = []
    my_data[exp]['r_cld']  = []
    for it in range(my_data[exp]['lwp_times'].size):
        # Use scipy.ndimage to count and label the clouds
        clouds, n_clds = ndimage.label(my_data[exp]['sfc_mask'][it,:,:])
        # ^- returns (2D field with numbered cloud objects, number of clouds)
        my_data[exp]['n_clds'][it] = n_clds/N_scale
        my_data[exp]['A_cld'].append(np.array([np.sum(np.where(clouds == cloud, my_data[exp]['sfc_mask'][it,:,:]*dx*dx, 0.0)) for cloud in range(1, n_clds)]))
        my_data[exp]['r_cld'].append(np.sqrt(my_data[exp]['A_cld'][it]/np.pi))
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Computing statistics profiles'
    my_data[exp]['profiles'] = {}
    my_data[exp]['profiles']['n_clds'] = np.empty_like(my_data[exp][mcl_key][:,:,0,0])
    my_data[exp]['profiles']['A_cld']  = []
    my_data[exp]['profiles']['r_cld']  = []
    for it in range(my_data[exp]['mcl_times'].size):
        # append an empty list for each time step
        my_data[exp]['profiles']['A_cld'].append([])
        my_data[exp]['profiles']['r_cld'].append([])
        for iz in range(z.size):
            # Use scipy.ndimage to count and label the clouds
            clouds, n_clds = ndimage.label(my_data[exp]['mask'][it,iz,:,:])
            # ^- returns (2D field with numbered cloud objects, number of clouds)
            my_data[exp]['profiles']['n_clds'][it,iz] = n_clds/N_scale
            # append each height into the empty time list
            my_data[exp]['profiles']['A_cld'][it].append(np.array([np.sum(np.where(clouds == cloud, my_data[exp]['mask'][it,iz,:,:]*dx*dx, 0.0)) for cloud in range(1, n_clds)]))
            my_data[exp]['profiles']['r_cld'][it].append(np.sqrt(my_data[exp]['profiles']['A_cld'][it][iz]/np.pi))

print '[' + dt.now().strftime('%H:%M:%S') + '] Making plots'
### Make some plots ###
exp_labels = {'0200m':'DX0200', '0400m':'DX0400', '0800m':'DX0800','1600m':'DX1600'}
my_colors = {'0200m':'brown', '0400m':'red', '0800m' : 'orange', '1600m':'gold'}
# N clouds
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
for exp in experiments:
    ax.plot(my_data[exp]['lwp_times']/60., my_data[exp]['n_clds'], color = my_colors[exp], label = 'DX' + exp[:-1])

ax.set_ylabel(u'$\hat{N}_{cld}$')
ax.set_xlabel(u'Time (hours)')
ax.set_xlim([0, 24])
ax.set_xticks([int(hour) for hour in hours[1:]])
ax.set_xticklabels(hours[1:])
plt.legend(loc = 0)
plt.show()

# A clouds
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
for exp in experiments:
    lowerQ = np.array([np.nanpercentile(my_data[exp]['A_cld'][it], 25) for it in range(my_data[exp]['lwp_times'].size)])
    median = np.array([np.nanpercentile(my_data[exp]['A_cld'][it], 50) for it in range(my_data[exp]['lwp_times'].size)])
    upperQ = np.array([np.nanpercentile(my_data[exp]['A_cld'][it], 75) for it in range(my_data[exp]['lwp_times'].size)])
    ax.plot(my_data[exp]['lwp_times']/60., median, color = my_colors[exp], label = 'DX' + exp[:-1], lw = 1)
    ax.fill_between(my_data[exp]['lwp_times']/60., lowerQ, upperQ, facecolor = my_colors[exp], edgecolor = 'none', alpha = 0.5)

ax.set_ylabel(u'$\overline{A}_{cld}$ (km$^{2}$)')
ax.set_xlabel(u'Time (hours)')
ax.set_xlim([0, 24])
ax.set_xticks([int(hour) for hour in hours[1:]])
ax.set_xticklabels(hours[1:])
plt.legend(loc = 0)
plt.show()

# r clouds
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
for exp in experiments:
    lowerQ = np.array([np.nanpercentile(my_data[exp]['r_cld'][it], 25) for it in range(my_data[exp]['lwp_times'].size)])
    median = np.array([np.nanpercentile(my_data[exp]['r_cld'][it], 50) for it in range(my_data[exp]['lwp_times'].size)])
    upperQ = np.array([np.nanpercentile(my_data[exp]['r_cld'][it], 75) for it in range(my_data[exp]['lwp_times'].size)])
    ax.plot(my_data[exp]['lwp_times']/60., median, color = my_colors[exp], label = 'DX' + exp[:-1], lw = 1)
    ax.fill_between(my_data[exp]['lwp_times']/60., lowerQ, upperQ, facecolor = my_colors[exp], edgecolor = 'none', alpha = 0.5)

ax.set_ylabel(u'$\overline{r}_{cld}$ (km)')
ax.set_xlabel(u'Time (hours)')
ax.set_xlim([0, 24])
ax.set_xticks([int(hour) for hour in hours[1:]])
ax.set_xticklabels(hours[1:])
plt.legend(loc = 0)
plt.show()

# some profiles of these quantities
# N clouds
fig = plt.figure(figsize = (10,6), tight_layout = True)
ax = fig.add_subplot(1, 2, 1)
for exp in experiments:
    my_it0 = np.where(my_data[exp]['mcl_times'] == 660.)[0][0]
    my_it1 = np.where(my_data[exp]['mcl_times'] == 720.)[0][0] + 1
    ax.plot(np.array([my_data[exp]['profiles']['n_clds'][it,:] for it in range(my_it0, my_it1)]).mean(axis = 0), z/1000., color = my_colors[exp], label = exp_labels[exp], lw = 2)

ax.set_xlabel(u'$\hat{N}_{cld}$')
ax.set_ylabel(u'Height (km)')
ax.set_ylim([0, 3.5])
ax.set_xlim([0, 0.012])
ax.set_title('a) 11AM - 12PM')
ax.legend(loc = 0, frameon = False, fontsize = 12)

ax = fig.add_subplot(1, 2, 2)
for exp in experiments:
    my_it0 = np.where(my_data[exp]['mcl_times'] == 1380.)[0][0]
    my_it1 = np.where(my_data[exp]['mcl_times'] == 1440.)[0][0] + 1
    ax.plot(np.array([my_data[exp]['profiles']['n_clds'][it,:] for it in range(my_it0, my_it1)]).mean(axis = 0), z/1000., color = my_colors[exp], lw = 2)

ax.set_xlabel(u'$\hat{N}_{cld}$')
ax.set_yticklabels([''])
ax.set_ylim([0, 3.5])
ax.set_xlim([0, 0.012])
ax.set_title('b) 11PM - 12AM')
plt.savefig('../Ch6_Figure10.png', dpi = 250, bbox_inches = 'tight')
plt.show()

# A clouds
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
for exp in experiments:
    my_it = np.where(my_data[exp]['mcl_times'] >= 360.)[0][0]
    lowerQ = np.array([np.nanpercentile(my_data[exp]['profiles']['A_cld'][my_it][iz], 25) for iz in range(z.size)])
    median = np.array([np.nanpercentile(my_data[exp]['profiles']['A_cld'][my_it][iz], 50) for iz in range(z.size)])
    upperQ = np.array([np.nanpercentile(my_data[exp]['profiles']['A_cld'][my_it][iz], 75) for iz in range(z.size)])
    ax.plot(median, z, color = my_colors[exp], label = 'DX' + exp[:-1])
    ax.fill_betweenx(z, lowerQ, upperQ, facecolor = my_colors[exp], edgecolor = 'none', alpha = 0.5)

ax.set_xlabel(u'$\overline{A}_{cld}$ (km$^{2}$)')
ax.set_ylabel(u'Height (m)')
ax.set_ylim([0, 3500])
plt.legend(loc = 0)
plt.show()

# r clouds
fig = plt.figure(figsize = (10,6), tight_layout = True)
ax = fig.add_subplot(1, 2, 1)
for exp in experiments:
    my_it0 = np.where(my_data[exp]['mcl_times'] == 660.)[0][0]
    my_it1 = np.where(my_data[exp]['mcl_times'] == 720.)[0][0] + 1
    radii = []
    for it in range(my_it0, my_it1):
        for iz in range(z.size):
            radii.append([])
            if len(my_data[exp]['profiles']['r_cld'][it][iz]) > 0:
                for radius in my_data[exp]['profiles']['r_cld'][it][iz]:
                    radii[iz].append(radius)
    lowerQ = np.array([np.nanpercentile(radii[iz], 25) for iz in range(z.size)])
    median = np.array([np.nanpercentile(radii[iz], 50) for iz in range(z.size)])
    upperQ = np.array([np.nanpercentile(radii[iz], 75) for iz in range(z.size)])
    ax.plot(median, z/1000., color = my_colors[exp], label = exp_labels[exp], lw = 2)
    ax.fill_betweenx(z/1000., lowerQ, upperQ, facecolor = my_colors[exp], edgecolor = 'none', alpha = 0.5)

ax.set_xlabel(u'$\overline{r}_{cld}$ (km)')
ax.set_ylabel(u'Height (m)')
ax.set_ylim([0, 3.5])
ax.set_title('a) 11AM - 12PM')
ax.legend(loc = 0, frameon = False, fontsize = 12)

ax = fig.add_subplot(1, 2, 2)
for exp in experiments:
    my_it0 = np.where(my_data[exp]['mcl_times'] == 1380.)[0][0]
    my_it1 = np.where(my_data[exp]['mcl_times'] == 1440.)[0][0] + 1
    radii = []
    for it in range(my_it0, my_it1):
        for iz in range(z.size):
            radii.append([])
            if len(my_data[exp]['profiles']['r_cld'][it][iz]) > 0:
                for radius in my_data[exp]['profiles']['r_cld'][it][iz]:
                    radii[iz].append(radius)
    lowerQ = np.array([np.nanpercentile(radii[iz], 25) for iz in range(z.size)])
    median = np.array([np.nanpercentile(radii[iz], 50) for iz in range(z.size)])
    upperQ = np.array([np.nanpercentile(radii[iz], 75) for iz in range(z.size)])
    ax.plot(median, z/1000., color = my_colors[exp], lw = 2)
    ax.fill_betweenx(z/1000., lowerQ, upperQ, facecolor = my_colors[exp], edgecolor = 'none', alpha = 0.5)

ax.set_xlabel(u'$\overline{r}_{cld}$ (km)')
ax.set_yticklabels([''])
ax.set_ylim([0, 3.5])
ax.set_title('b) 11PM - 12AM')
plt.savefig('../Ch6_Figure11.png', dpi = 250, bbox_inches = 'tight')
plt.show()

# daytime mean cloud radii
fig = plt.figure(figsize = (6,6), tight_layout = True)
ax = fig.add_subplot(1, 1, 1)
for exp in experiments:
    my_it0 = np.where(my_data[exp]['mcl_times'] == 360.)[0][0]
    my_it1 = np.where(my_data[exp]['mcl_times'] == 1080.)[0][0] + 1
    radii = []
    for it in range(my_it0, my_it1):
        for iz in range(z.size):
            radii.append([])
            if len(my_data[exp]['profiles']['r_cld'][it][iz]) > 0:
                for radius in my_data[exp]['profiles']['r_cld'][it][iz]:
                    radii[iz].append(radius)
    lowerQ = np.array([np.nanpercentile(radii[iz], 25) for iz in range(z.size)])
    median = np.array([np.nanpercentile(radii[iz], 50) for iz in range(z.size)])
    upperQ = np.array([np.nanpercentile(radii[iz], 75) for iz in range(z.size)])
    ax.plot(median, z/1000., color = my_colors[exp], label = exp_labels[exp], lw = 2)
    ax.fill_betweenx(z/1000., lowerQ, upperQ, facecolor = my_colors[exp], edgecolor = 'none', alpha = 0.5)

ax.set_xlabel(u'$\overline{r}_{cld}$ (km)')
ax.set_ylabel(u'Height (m)')
ax.set_ylim([0, 3.5])
ax.set_title('6AM - 6PM')
ax.legend(loc = 0, frameon = False, fontsize = 12)

plt.savefig('../Ch6_Figure12.png', dpi = 250, bbox_inches = 'tight')
plt.show()

