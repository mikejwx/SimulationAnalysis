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

def plume_area(variable, threshold = 0.1, smoothing_radius = 6000.):
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
    returnable  = a 1D array of our plume characteristic of choice, to be edited
                  in below
                = f(t)
    """
    # Get the mean and standard deviation
    var_mean = np.nanmean(variable, axis = (1, 2))
    
    # Calculate the anomaly
    var_anomaly = np.transpose((np.transpose(variable) - var_mean))
    
    # Define a smoothing structure
    n_r, n_c = np.array([smoothing_radius, smoothing_radius])/100
    w_c = np.ones((n_r, n_c))
    r, c = np.meshgrid(np.arange(n_r), np.arange(n_c))
    w_c = np.where((np.sqrt((r - n_r/2)**2 + (c - n_c/2)**2) <= n_r/2), w_c, 0.)
    w_c /= np.sum(w_c)
    
    returnable = np.zeros(variable.shape[0])
    for it in xrange(variable.shape[0]):
        # Smooth the anomalies
        smoothed_anomaly = signal.convolve2d(var_anomaly[it,:,:], w_c, mode = 'same', boundary = 'wrap')
        
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
            returnable[it] = np.nan
        else:
            returnable[it] = sizes[biggest_one]
    
    return returnable

def plume_length(variable, threshold = 0.1, smoothing_radius = 6000.):
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
    returnable  = a 1D array of our plume characteristic of choice, to be edited
                  in below
                = f(t)
    """
    # Get the horizontal mean
    var_mean = np.nanmean(variable, axis = (1, 2))
    
    # Calculate the anomaly
    var_anomaly = np.transpose((np.transpose(variable) - var_mean))
    
    # Define a smoothing structure
    n_r, n_c = np.array([smoothing_radius, smoothing_radius])/100
    w_c = np.ones((n_r, n_c))
    r, c = np.meshgrid(np.arange(n_r), np.arange(n_c))
    w_c = np.where((np.sqrt((r - n_r/2)**2 + (c - n_c/2)**2) <= n_r/2), w_c, 0.)
    w_c /= np.sum(w_c)
    
    returnable = np.zeros(variable.shape[0])
    for it in xrange(variable.shape[0]):
        # Smooth the anomalies
        smoothed_anomaly = signal.convolve2d(var_anomaly[it,:,:], w_c, mode = 'same', boundary = 'wrap')
        
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
            returnable[it] = np.nan
        else:
            returnable[it] = np.nanmax(np.where((label_image == biggest_one), R, np.nan)) # R is defined globally!
    
    return returnable
################################################################################
#                                                                              #
# Section 1: How big is it?                                                    #
#                                                                              #
################################################################################

paths = {'control' : '/nerc/n02/n02/xb899100/CloudTrail/Control/',
         'exp01'   : '/work/n02/n02/xb899100/cylc-run/u-bg023/share/data/history',
         'exp02'   : '/work/n02/n02/xb899100/cylc-run/u-bg113/share/data/history'}
experiments = ['control', 'exp01', 'exp02']
times = {}
my_length = {}
my_area = {}

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
    send_email(message = 'Working on experiment: ' + experiment, subject = 'Plume Lengths Update', attachments = [''], isAttach = False)
    my_length[experiment] = plume_length(bouy_nc.variables[theta_key][:,iz,:,:])
    my_area[experiment] = plume_area(bouy_nc.variables[theta_key][:,iz,:,:])
    bouy_nc.close()
    
    for hour in hours:
        #print 'Working on hour == ' + hour + ' for experiment: ' + experiment
        #send_email(message = 'Working on hour == ' + hour + ' for experiment: ' + experiment, subject = 'Plume Lengths Update', attachments = [''], isAttach = False)
        bouy_nc = Dataset(paths[experiment] + '/bouy_' + hour + '.nc', 'r')
        my_length[experiment] = np.concatenate((my_length[experiment], plume_length(bouy_nc.variables[theta_key][:,iz,:,:])), axis = 0)
        my_area[experiment] = np.concatenate((my_area[experiment], plume_area(bouy_nc.variables[theta_key][:,iz,:,:])), axis = 0)
        times[experiment] = np.concatenate((times[experiment], bouy_nc.variables[time_key][:]*1.), axis = 0)
        bouy_nc.close()

fig = plt.figure(figsize = (15, 6))
ax0 = fig.add_subplot(1, 2, 1)
ax0.plot(times['control']/60., my_length['control']/1000., 'k', lw = 2, label = 'control: $H = 250 W m^{-2}$')
ax0.plot(times['exp01']/60., my_length['exp01']/1000., 'b', label = 'exp01: $H = 125 W m^{-2}$')
ax0.plot(times['exp02']/60., my_length['exp02']/1000., 'r', label = 'exp02: $H = 375 W m^{-2}$')
ax0.set_xlim([0, 24])
ax0.set_xticks(range(0, 25, 3))
ax0.set_ylim([0, 125])
ax0.set_ylabel('Warm Plume Length (km)')
ax0.set_xlabel('Time (hrs)')
ax0.set_title('Warm Plume Length at ' + str(int(z[iz])) + 'm')
plt.legend(loc = 2)

ax1 = fig.add_subplot(1, 2, 2)
ax1.plot(times['control']/60., my_area['control']*100./1000./1000., 'k', lw = 2, label = 'control: $H = 250 W m^{-2}$')
ax1.plot(times['exp01']/60., my_area['exp01']*100./1000./1000., 'b', label = 'exp01: $H = 125 W m^{-2}$')
ax1.plot(times['exp02']/60., my_area['exp02']*100./1000./1000., 'r', label = 'exp02: $H = 375 W m^{-2}$')
ax1.set_xlim([0, 24])
ax1.set_xticks(range(0, 25, 3))
ax1.set_ylim([0, 5])
ax1.set_ylabel('Anomalous Warm Area (km$^{2}$)')
ax1.set_xlabel('Time (hrs)')
ax1.set_title('Anomalous Warm Area at ' + str(int(z[iz])) + 'm')
plt.legend(loc = 2)
plt.savefig('../Warm_Plume_analysis.png', dpi = 100)
send_email('Finished analysing the warm plume', subject = 'Plume Analysis Complete', attachments = ['../Warm_Plume_analysis.png'], isAttach = True)

"""
H0_control = 250.
H0_exp1 = 125.
H0_exp2 = 375.
"""
def H(H0, times, dt, t_max):
    """
    Calculate the diurnal cycle of sensible heat fluxes given a H0 and a set of
    times, length of day, and time of peak heating.
    """
    times /= 60.
    H = H0*np.cos(0.5*np.pi*(t_max - times)/(dt/2.))**1.5
    if type(H) != np.ndarray:
        H = np.array([H])
    H = np.where((H <= 0.) or np.isnan(H), 0., H)
    return H


# Plot the heat flux for when the warm plume forms and dissipates
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# Start
ax.plot(0, H(250., times['control'][np.min(np.where(-1*(np.isnan(my_area['control'])-1))[0])], 12., 12.), 'k^', ms = 10)
ax.text(0.1, H(250., times['control'][np.min(np.where(-1*(np.isnan(my_area['control'])-1))[0])], 12., 12.), str(int(times['control'][np.min(np.where(-1*(np.isnan(my_area['control'])-1))[0])]))+'mins')
ax.plot(1, H(125., times['exp01'][np.min(np.where(-1*(np.isnan(my_area['exp01'])-1))[0])], 12., 12.), 'b^', ms = 10)
ax.text(1.1, H(125., times['exp01'][np.min(np.where(-1*(np.isnan(my_area['exp01'])-1))[0])], 12., 12.), str(int(times['exp01'][np.min(np.where(-1*(np.isnan(my_area['exp01'])-1))[0])]))+'mins')
ax.plot(2, H(375., times['exp02'][np.min(np.where(-1*(np.isnan(my_area['exp02'])-1))[0])], 12., 12.), 'r^', ms = 10)
ax.text(2.1, H(375., times['exp02'][np.min(np.where(-1*(np.isnan(my_area['exp02'])-1))[0])], 12., 12.), str(int(times['exp02'][np.min(np.where(-1*(np.isnan(my_area['exp02'])-1))[0])]))+'mins')
# End
ax.plot(0, H(250., times['control'][np.max(np.where(-1*(np.isnan(my_area['control'])-1))[0])], 12., 12.), 'ko', ms = 10)
ax.text(0.1, H(250., times['control'][np.max(np.where(-1*(np.isnan(my_area['control'])-1))[0])], 12., 12.), str(int(times['control'][np.max(np.where(-1*(np.isnan(my_area['control'])-1))[0])]))+'mins')
ax.plot(1, H(125., times['exp01'][np.max(np.where(-1*(np.isnan(my_area['exp01'])-1))[0])], 12., 12.), 'bo', ms = 10)
ax.text(1.1, H(125., times['exp01'][np.max(np.where(-1*(np.isnan(my_area['exp01'])-1))[0])], 12., 12.), str(int(times['exp01'][np.max(np.where(-1*(np.isnan(my_area['exp01'])-1))[0])]))+'mins')
ax.plot(2, H(375., times['exp02'][np.max(np.where(-1*(np.isnan(my_area['exp02'])-1))[0])], 12., 12.), 'ro', ms = 10)
ax.text(2.1, H(375., times['exp02'][np.max(np.where(-1*(np.isnan(my_area['exp02'])-1))[0])], 12., 12.), str(int(times['exp02'][np.max(np.where(-1*(np.isnan(my_area['exp02'])-1))[0])]))+'mins')
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['Control', 'exp01', 'exp02'])
ax.set_ylabel('Sensible Heat Flux (W m$^{-2}$)')
ax.set_xlim([-0.5, 2.5])
ax.set_ylim([-10, 120])
ax.plot([-0.5, 2.5], [0, 0], 'k--')
ax.set_title('Heat Fluxes at start and end of warm plume')
plt.savefig('../Warm_Plume_analysis2.png', dpi = 100)
plt.show()


#command1 = 'convert ../control*.png -delay 60 -loop 0 ../control_anim.gif'
#command2 = 'convert ../exp01*.png -delay 60 -loop 0 ../exp01_anim.gif'
#command3 = 'convert ../exp02*.png -delay 60 -loop 0 ../exp02_anim.gif'
#commands = [command1, command2, command3]
#[os.system(command) for command in commands]
#send_email('Here are some animations of the warm plume', subject = 'Plume Lengths Complete', attachments = ['../WarmPlume/control_anim.gif', '../WarmPlume/exp01_anim.gif', '../WarmPlume/exp02_anim.gif'], isAttach = True)


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


