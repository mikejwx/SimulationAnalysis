import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
from netCDF4 import Dataset
from scipy import ndimage, signal
from analysis_tools import send_email
import os

# Read in a land-sea mask
lsm = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
lsm_data = lsm.variables['lsm'][0,0,:,:]*1.
lsm.close()

def plume_analysis(variable, threshold = 0.1, smoothing_radius = 6000.):
    """
    Find the length of the plume.
    
    Take the mean of the input variable.
    Smooth variable to "smoothing_radius" resolution using a circular filter.
    
    Mask out features that are warmer than "threshold" Kelvin.
    Find the biggest feature and assume that it is the plume.
    Find the pixel within the plume furthest from the island centrepoint, define
    this as the length of the plume.
    Repeat for all provided time steps.
    Return an array of the plume lengths "distance"
    ----------------------------------------------------------------------------
    INPUT:
    variable  = a 3D array of data e.g. potential temperature (K)
              = f(t, y, x) e.g. theta
    threshold = the number of standard deviations to consider objects
    
    OUTPUT:
    distance  = a 1D array of plume lengths in time
              = f(t)
    """
    # Get the mean and standard deviation
    var_mean = np.nanmean(variable, axis = (1, 2))
    #var_sd   = np.std(variable, axis = (1, 2))
    
    # Calculate the anomaly
    var_anomaly = np.transpose((np.transpose(variable) - var_mean))
    
    # What time is it?
    time = bouy_nc.variables[time_key][:]*1.
    
    # Define a smoothing structure
    n_r, n_c = np.array([smoothing_radius, smoothing_radius])/100
    w_c = np.ones((n_r, n_c))
    r, c = np.meshgrid(np.arange(n_r), np.arange(n_c))
    w_c = np.where((np.sqrt((r - n_r/2)**2 + (c - n_c/2)**2) <= n_r/2), w_c, 0.)
    w_c /= np.sum(w_c)
    
    distance = np.zeros(variable.shape[0])
    for it in xrange(variable.shape[0]):
        # Smooth the anomalies
        smoothed_anomaly = signal.convolve2d(var_anomaly[it,:,:], w_c, mode = 'same', boundary = 'wrap')
        
        # Make some plots of the smoothed theta anomalies
        fig = plt.figure(tight_layout = True)
        ax = fig.add_subplot(1,1,1, adjustable = 'box', aspect = 1)
        sa = ax.contourf(X/1000., Y/1000., smoothed_anomaly, cmap = 'bwr', levels = [level for level in np.linspace(-1., 1., 21) if level != 0], extend = 'both')
        fig.colorbar(sa, ax = ax, orientation = 'horizontal')
        ax.contour(X/1000., Y/1000., lsm_data, levels = [1e-16], colors = ['k'])
        ax.set_title(experiment + ' T+' + "{0:04d}".format(int(time[it])) + ' mins')
        plt.savefig('../WarmPlume/'+experiment+'_T+' + "{0:04d}".format(int(time[it])) + '.png', dpi = 100)
        plt.close('all')
        #send_email(message = 'Working on hour == ' + hour + ', time == ' + "{0:04d}".format(int(time[it])) + ' for experiment: ' + experiment, subject = 'Plume Lengths Update', attachments = ['../WarmPlume/'+experiment+'_T+' + "{0:04d}".format(int(time[it])) + '.png'], isAttach = True)
        # Get a mask of all the pixels above 3 standard deviations
        mask = np.where((smoothed_anomaly >= threshold), 1., 0.)
        
        # Label all the features returned in the mask
        label_image, count = ndimage.label(mask)
        
        # Keep only the largest feature
        sizes = ndimage.sum(mask, label_image, range(count + 1))
        biggest_one = np.where(sizes == np.max(sizes))[0][0]
        
        # Get the maximum distance from that feature
        if sizes[biggest_one] < 1000:
            # Then indistinguishable from background perturbations, and not a plume
            distance[it] = np.nan
        else:
            distance[it] = np.max(np.where((label_image == biggest_one), R, 0.))/1000.
    
    return distance

################################################################################
#                                                                              #
# Section 1: How big is it?                                                    #
#                                                                              #
################################################################################

paths = {'control' : '..',
         'exp01'   : '/work/n02/n02/xb899100/cylc-run/u-bg023/share/data/history',
         'exp02'   : '/work/n02/n02/xb899100/cylc-run/u-bg113/share/data/history'}
experiments = ['control', 'exp01', 'exp02']
times = {}
my_distance = {}

# Define keys for variables of interest
theta_key = u'STASH_m01s00i004'
z_key     = u'thlev_zsea_theta'

hours = ["{0:02d}".format(h) for h in xrange(3, 24, 3)]

for experiment in experiments:
    # Read the bouy netCDF so that we can see theta data
    hour = '00'
    bouy_nc = Dataset(paths[experiment] + '/bouy_' + hour + '.nc', 'r')
    time_key  = [key for key in bouy_nc.variables.keys() if 'min' in key][0]
    
    z       = bouy_nc.variables[z_key][:]*1.
    times[experiment] = bouy_nc.variables[time_key][:]*1.
    
    ################################################################################
    # Length (this will eventually be expanded to include volume of the plume)     #
    ################################################################################
    iz = np.where(np.abs(z - 350) == np.min(np.abs(z - 350)))[0][0]
    
    ### Define plume as areas that are 2 standard deviations warmer than the 
    # horizontal mean at that height at that time ###
    
    x = np.arange(0., 116000., 100.)
    y = np.arange(0., 31900., 100.)
    X, Y = np.meshgrid(x, y)
    island_area = 50. # square kilometres
    island_radius = np.sqrt(island_area/np.pi)*1000. # metres
    
    x_c = 100000. + 2.*island_radius
    y_c = 4.*island_radius
    R = np.sqrt((X - x_c)**2. + (Y - y_c)**2.)
    
    #print 'Working on hour == 00 for experiment: ' + experiment
    send_email(message = 'Working on hour == ' + hour + ' for experiment: ' + experiment, subject = 'Plume Lengths Update', attachments = [''], isAttach = False)
    my_distance[experiment] = plume_analysis(bouy_nc.variables[theta_key][:,iz,:,:])
    bouy_nc.close()
    
    for hour in hours:
        #print 'Working on hour == ' + hour + ' for experiment: ' + experiment
        send_email(message = 'Working on hour == ' + hour + ' for experiment: ' + experiment, subject = 'Plume Lengths Update', attachments = [''], isAttach = False)
        bouy_nc = Dataset(paths[experiment] + '/bouy_' + hour + '.nc', 'r')
        my_distance[experiment] = np.concatenate((my_distance[experiment], plume_analysis(bouy_nc.variables[theta_key][:,iz,:,:])), axis = 0)
        times[experiment] = np.concatenate((times[experiment], bouy_nc.variables[time_key][:]*1.), axis = 0)
        bouy_nc.close()

plt.plot(times['control']/60., my_distance['control'], 'k', lw = 2, label = 'control: $\\beta = 1$')
plt.plot(times['exp01']/60., my_distance['exp01'], 'b', label = 'exp01: $\\beta = 1/3$')
plt.plot(times['exp02']/60., my_distance['exp02'], 'r', label = 'exp02: $\\beta = 3$')
plt.ylabel('Plume length (km)')
plt.xlabel('Time (hrs)')
plt.title('Warm Plume Length')
plt.legend(loc = 2)
plt.savefig('../Warm_Plume_length.png', dpi = 100)
send_email('Finished calculating the length of the warm plume in the control and the first experiment', subject = 'Plume Lengths Complete', attachments = ['../Warm_Plume_length.png'], isAttach = True)

command1 = 'convert ../WarmPlume/control*.png -delay 30 -loop 0 control_anim.gif'
command2 = 'convert ../WarmPlume/exp01*.png -delay 60 -loop 0 exp01_anim.gif'
command3 = 'convert ../WarmPlume/exp02*.png -delay 60 -loop 0 exp02_anim.gif'
commands = [command1, command2, command3]
[os.system(command) for command in commands]
send_email('Here are some animations of the warm plume', subject = 'Plume Lengths Complete', attachments = ['../WarmPlume/control_anim.gif', '../WarmPlume/exp01_anim.gif', '../WarmPlume/exp02_anim.gif'], isAttach = True)


"""
###Control###
# get the original data
bouy_nc = Dataset('../bouy_12.nc', 'r')
theta = bouy_nc.variables[theta_key][0,iz,:,:]
theta_mean = np.nanmean(theta)
theta_sd = np.std(theta)
bouy_nc.close()

resolutions = [500, 1000, 2000, 4000, 8000]
fig = plt.figure(tight_layout = True)
ax1 = fig.add_subplot(len(resolutions)+1, 2, 1, adjustable = 'box', aspect = 1)
original = ax1.contourf(X, Y, theta - theta_mean, levels = [l for l in np.linspace(-1, 1, 11) if l != 0], cmap = 'bwr', extend = 'both')
fig.colorbar(original, ax = ax1)
ax1.set_title('Original, 100 m')

# for dx = dy = 100 m, so 4 km smoothing = 40 grid points
count = 3
for dx in resolutions:
    n_r, n_c = np.array([dx, dx])/100
    w_s = np.ones((n_r, n_c))
    w_s /= np.sum(w_s)
    
    f = signal.convolve2d(theta, w_s, mode = 'same', boundary = 'wrap')
    
    ax2 = fig.add_subplot(len(resolutions)+1, 2, count, adjustable = 'box', aspect = 1)
    filtered = ax2.contourf(X, Y, (f - theta_mean), levels = [l for l in np.linspace(-1, 1, 11) if l != 0], cmap = 'bwr', extend = 'both')
    ax2.set_xlim([np.min(x), np.max(x)])
    ax2.set_ylim([np.min(y), np.max(y)])
    fig.colorbar(filtered, ax = ax2)
    ax2.set_title('Square, ' + str(dx) + ' m')
    count += 1
    
    n_r, n_c = np.array([dx, dx])/100
    w_c = np.ones((n_r, n_c))
    r, c = np.meshgrid(np.arange(n_r), np.arange(n_c))
    w_c = np.where((np.sqrt((r - n_r/2)**2 + (c - n_c/2)**2) <= n_r/2), w_c, 0.)
    w_c /= np.sum(w_c)
    
    f = signal.convolve2d(theta, w_c, mode = 'same', boundary = 'wrap')
    
    ax3 = fig.add_subplot(len(resolutions)+1, 2, count, adjustable = 'box', aspect = 1)
    filtered = ax3.contourf(X, Y, (f - theta_mean), levels = [l for l in np.linspace(-1, 1, 11) if l != 0], cmap = 'bwr', extend = 'both')
    ax3.set_xlim([np.min(x), np.max(x)])
    ax3.set_ylim([np.min(y), np.max(y)])
    fig.colorbar(filtered, ax = ax3)
    ax3.set_title('Circle, ' + str(dx) + ' m')
    count += 1

plt.show()

### Experiment 1###
# get the original data
bouy_nc = Dataset('~/cylc-run/u-bg023/share/data/history/bouy_12.nc', 'r')
theta = bouy_nc.variables[theta_key][0,iz,:,:]
theta_mean = np.nanmean(theta)
theta_sd = np.std(theta)
bouy_nc.close()

fig = plt.figure(tight_layout = True)
ax1 = fig.add_subplot(len(resolutions)+1, 2, 1, adjustable = 'box', aspect = 1)
original = ax1.contourf(X, Y, (theta - theta_mean), levels = [l for l in np.linspace(-1, 1, 11) if l != 0], cmap = 'bwr', extend = 'both')
fig.colorbar(original, ax = ax1)
ax1.set_title('Original, 100 m')

# for dx = dy = 100 m, so 4 km smoothing = 40 grid points
count = 3
for dx in resolutions:
    n_r, n_c = np.array([dx, dx])/100
    w_s = np.ones((n_r, n_c))
    w_s /= np.sum(w_s)
    
    f = signal.convolve2d(theta, w_s, mode = 'same', boundary = 'wrap')
    
    ax2 = fig.add_subplot(len(resolutions)+1, 2, count, adjustable = 'box', aspect = 1)
    filtered = ax2.contourf(X, Y, (f - theta_mean), levels = [l for l in np.linspace(-1, 1, 11) if l != 0], cmap = 'bwr', extend = 'both')
    ax2.set_xlim([np.min(x), np.max(x)])
    ax2.set_ylim([np.min(y), np.max(y)])
    fig.colorbar(filtered, ax = ax2)
    ax2.set_title('Square, ' + str(dx) + ' m')
    count += 1
    
    n_r, n_c = np.array([dx, dx])/100
    w_c = np.ones((n_r, n_c))
    r, c = np.meshgrid(np.arange(n_r), np.arange(n_c))
    w_c = np.where((np.sqrt((r - n_r/2)**2 + (c - n_c/2)**2) <= n_r/2), w_c, 0.)
    w_c /= np.sum(w_c)
    
    f = signal.convolve2d(theta, w_c, mode = 'same', boundary = 'wrap')
    
    ax3 = fig.add_subplot(len(resolutions)+1, 2, count, adjustable = 'box', aspect = 1)
    filtered = ax3.contourf(X, Y, (f - theta_mean), levels = [l for l in np.linspace(-1, 1, 11) if l != 0], cmap = 'bwr', extend = 'both')
    ax3.set_xlim([np.min(x), np.max(x)])
    ax3.set_ylim([np.min(y), np.max(y)])
    fig.colorbar(filtered, ax = ax3)
    ax3.set_title('Circle, ' + str(dx) + ' m')
    count += 1

plt.show()


"""


