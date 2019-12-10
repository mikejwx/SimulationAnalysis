import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import *
from netCDF4 import Dataset
from STASH_keys import *
from analysis_tools import *
from scipy import integrate, ndimage

"""
Comparing the statistical and spatial distribution of cloud top heights before
and during the cloud trail
"""

needed_vars = [zi_new_key, ctz_key, lwp_key]
path = '/nerc/n02/n02/xb899100/CloudTrail/Control2/'
target_times = [ 360.0, 720.0 ]
my_data = {}
with Dataset(path + 'lwp_00.nc', 'r') as lwp_nc:
    lwp_time_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
    lwp_times = lwp_nc.variables[lwp_time_key][:]*1.
    idx_lwp = [it for it in range(lwp_times.shape[0]) if lwp_times[it] in target_times]
    for needed_var in needed_vars:
        if needed_var in lwp_nc.variables.keys():
            print 'Reading ' + needed_var
            my_data[needed_var] = lwp_nc.variables[needed_var][idx_lwp,:,:]*1.

for hour in ['03', '09']:
    with Dataset(path + 'zi_' + hour + '.nc', 'r') as zi_nc:
        if hour == '03':
            ctz_time_key = 'time'
            ctz_times = zi_nc.variables[ctz_time_key][:]*1.
            idx_ctz = [it for it in range(ctz_times.shape[0]) if ctz_times[it] in target_times]
            for needed_var in needed_vars:
                if needed_var in zi_nc.variables.keys():
                    print 'Reading ' + needed_var
                    my_data[needed_var] = zi_nc.variables[needed_var][idx_ctz,:,:]*1.
        else:
            ctz_times = np.concatenate((ctz_times[:], zi_nc.variables[ctz_time_key][:]*1.), axis = 0)
            for needed_var in needed_vars:
                if needed_var in zi_nc.variables.keys():
                    print 'Reading ' + needed_var
                    my_data[needed_var] = np.concatenate((my_data[needed_var][:], zi_nc.variables[needed_var][idx_ctz,:,:]*1.), axis = 0)

with Dataset(path + 'wind_09.nc', 'r') as wind_nc:
    z = wind_nc.variables['thlev_zsea_theta'][:]*1.
    my_data[u_key] = wind_nc.variables[u_key][-1,:,:,:]*1.
    my_data[v_key] = wind_nc.variables[v_key][-1,:,:,:]*1.

# define domain and island
x, y = np.meshgrid(np.arange(my_data[lwp_key][0,:,:].shape[1])*0.1, np.arange(my_data[lwp_key][0,:,:].shape[0])*0.1)
island_area = 50.0
island_radius = np.sqrt(island_area/np.pi)
island_x, island_y = 100.+island_radius, 4*island_radius
R = np.sqrt((x - island_x)**2 + (y - island_y)**2)

cloud_mask = np.where(my_data[lwp_key] > 0, 1.0, 0.0)

# get cloud top height statistics
max_cld_top_height = []
n_cld = []
for it in range(len(idx_lwp)):
    # Use scipy.ndimage to count and label the clouds
    clouds, n_clds = ndimage.label(cloud_mask[it,:,:])
    print '# Clouds = ' + str(n_clds)
    cloud_ctz = np.array([np.nanmax(np.where(clouds == cloud, my_data[ctz_key][it,:,:], np.nan)) for cloud in range(n_clds)])
    max_cld_top_height.append(cloud_ctz)
    n_cld.append(n_clds)

ctz_before = np.array([ctz for ctz in max_cld_top_height[0] if ctz == ctz])
ctz_during = np.array([ctz for ctz in max_cld_top_height[1] if ctz == ctz])

zi_mean = np.nanmean(my_data[zi_new_key][0,:,:])
iz      = np.where(np.abs(z - zi_mean) == np.min(np.abs(z - zi_mean)))[0][0]
U_mean  = np.nanmean(integrate.trapz(x = z[:iz], y = my_data[u_key][:iz,:,:], axis = 0)/z[iz])
V_mean  = np.nanmean(integrate.trapz(x = z[:iz], y = my_data[v_key][:iz,:,:], axis = 0)/z[iz])

wind_spd, wind_dir = fromComponents(U_mean, V_mean)

# Use the wind direction to define our downwind rectangular region
mask, y_prime, x_prime = downwind_rectangle(wind_dir, island_x, island_y, x, y, island_radius, dist_0 = -12.0, dist_1 = 100.0, half_width = 3.0)

# get the histograms
z = 0.5*(z[1:] + z[:-1])
ctz_before_test = np.array([ob for ob in my_data[ctz_key][0,:,:].flatten() if ob == ob])
ctz_during_test = np.array([ob for ob in my_data[ctz_key][-1,:,:].flatten() if ob == ob])

before = plt.hist(ctz_before, bins = z, normed = True, cumulative = True, histtype = 'step')
plt.close('all')
after  = plt.hist(ctz_during, bins = z, normed = True, cumulative = True, histtype = 'step')
plt.close('all')

# make the plot
fig = plt.figure(figsize = (8, 6))
axa = fig.add_subplot(2, 2, 1)
axa.plot(0.5*(before[1][1:] + before[1][:-1])/1000., before[0], color = 'k', lw = 2, label='T+6hrs')
axa.plot(0.5*(after[1][1:] + after[1][:-1])/1000., after[0], color = 'b', label = 'T+12hrs')
axa.set_xlim([0.5, 3.5])
axa.text(0.7, 0.875, 'a)')
axa.set_ylabel('Fraction of Clouds')
axa.set_xlabel('Cloud Top Height, $z_{top}$ (km)')
axa.legend(loc = 0, frameon = False)

axb = fig.add_subplot(2, 2, 2)
axb.plot(0.5*(after[1][1:] + after[1][:-1])/1000., after[0] - before[0], color = 'k')
axb.set_xlabel('Cloud Top Height, $z_{top}$ (km)')
axb.set_ylabel('Change in Fraction')
axb.set_xlim([0.5, 3.5])
axb.set_ylim([-0.02, 0.1])
axb.text(0.7, 0.0875, 'b)')

axc = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
my_levels = [-0.5, 0.5]
axc.contourf(x, y, cloud_mask[-1,:,:], levels = [0.5, 1.1], colors = ['k'])
axc.contour(x, y, R, levels = [island_radius], colors = ['r'], linewidths = [2])
axc.contour(x, y, np.where(mask == mask, 1.0, 0.0), levels = [0.5], colors = ['b'])
axc.set_title('c) Cloud Mask at T+12hrs')
# annotate the mean cloud top height in the CT region
axc.annotate(u'$\overline{z_{top}}$ = ' + str(round(np.nanmean(mask*my_data[ctz_key][-1,:,:]),1)) + ' m', xy = (60, 5.5), xytext = (90, -7.5), color = 'b', arrowprops=dict(facecolor='b', shrink=0.025,edgecolor='none'))
#annotate the mean cloud top height outside the CT region
axc.annotate(u'$\overline{z_{top}}$ = ' + str(round(np.nanmean(np.where(mask == mask, np.nan, my_data[ctz_key][-1,:,:])),1)) + ' m', xy = (20, 5.5), xytext = (-5, -7.5), color = 'k', arrowprops=dict(facecolor='k', shrink=0.025,edgecolor='none'))
axc.set_xlabel('x (km)')
axc.set_ylabel('y (km)')

plt.subplots_adjust(left = 0.10, bottom = 0.1, right = 0.9, wspace = 0.5, hspace = 0.3)
plt.savefig('../Ch5_Figure11.png', dpi = 250, bbox_inches = 'tight')
plt.show()

