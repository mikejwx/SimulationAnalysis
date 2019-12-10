"""
Analysis of the cloud fraction with height, collapsing the mcl dataset into a 
mask of cloud vs. no cloud f(t,z,y,x) then to a horizontal mean (a.k.a. cloud 
fraction) f(t,z)

Compare the time series of liquid water path to cloud fraction total, and at 
different heights.
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from datetime import datetime as dt

#path = '/nerc/n02/n02/xb899100/CloudTrail/Control/' # Control
#path = '/nerc/n02/n02/xb899100/CloudTrail/Control_Spinup/' # Spinup for the control simulation
path = '/nerc/n02/n02/xb899100/CloudTrail/U05_Spinup/' # Spinup for the U05 simulation

ID = 'U05_Spinup'
cloud_thresh = 1e-16
l_spinup = 1

lwp_key = u'STASH_m01s30i405'
mcl_key = u'STASH_m01s00i392'

if l_spinup:
    hours = ["{0:02d}".format(h) for h in xrange(1, 11)]
else:
    hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]
    # Open the netCDF for lwp (only one for the simulation)
    lwp_nc = Dataset(path + 'lwp_00.nc', 'r')
    # Read the data from the lwp netCDF
    lwp_data = lwp_nc.variables[lwp_key][:]*1.
    # Read the times for the lwp netCDF
    times_lwp = lwp_nc.variables['min5_0'][:]*1.
    # Read the latitude and longitude coordinates
    lat = lwp_nc.variables[u'latitude_t'][:]*1.
    lon = lwp_nc.variables[u'longitude_t'][:]*1.
    lwp_nc.close()

for hour in hours:
    print '[' + dt.now().strftime('%H:%M:%S') + '] ' + hour
    # Open the netCDF for mr
    mr_nc = Dataset(path + 'mr_'+hour+'.nc', 'r')
    # Read the data for mcl from the mr netCDF
    mcl_data = mr_nc.variables[mcl_key][:]*1.
    
    if l_spinup:
        lwp_nc = Dataset(path + 'lwp_' + hour + '.nc', 'r')
        if hour == hours[0]:
            lwp_data = lwp_nc.variables[lwp_key][:]*1.
            lwp_times_key = [key for key in lwp_nc.variables.keys() if 'min' in key][0]
            times_lwp = lwp_nc.variables[lwp_times_key][:]*1.
            lat = lwp_nc.variables[u'latitude_t'][:]*1.
            lon = lwp_nc.variables[u'longitude_t'][:]*1.
        else:
            lwp_data = np.concatenate((lwp_data, lwp_nc.variables[lwp_key][:]*1.), axis = 0)
            times_lwp = np.concatenate((times_lwp, lwp_nc.variables[lwp_times_key][:]*1.), axis = 0)
        lwp_nc.close()
    
    if hour == hours[0]:
        # Mask out the mcl_data
        mcl_mask = np.where((mcl_data >= cloud_thresh), 1., 0.)
        mcl_times_key = [key for key in mr_nc.variables.keys() if 'min' in key][0]
        times_mcl = mr_nc.variables[mcl_times_key][:]*1.
        # Read the height coordinate
        z = mr_nc.variables[u'thlev_zsea_theta'][:]*1.
    else:
        mcl_mask = np.concatenate((mcl_mask, np.where((mcl_data >= cloud_thresh), 1., 0.)), axis = 0)
        times_mcl = np.concatenate((times_mcl, mr_nc.variables[mcl_times_key][:]*1.), axis = 0)
    mr_nc.close()

# Collapse mcl_mask with horizontal mean and convert to percentages
mcl_fraction = np.nanmean(mcl_mask, axis = (2, 3))*100. # f(t,z,y,x) -> f(t,z)
lwp_fraction = np.nanmean(np.where((lwp_data >= cloud_thresh), 1., 0.), axis = (1, 2))*100. # f(t,y,x) -> f(t)

# heights
iz1 = np.where(np.abs(z - 750.0) == np.min(np.abs(z - 750.0)))[0][0]
iz2 = np.where(np.abs(z - 1500.0) == np.min(np.abs(z - 1500.0)))[0][0]
iz3 = np.where(np.abs(z - 2500.0) == np.min(np.abs(z - 2500.0)))[0][0]

if l_spinup:
    # Time series at different heights for entire simulation
    fig = plt.figure(tight_layout = True)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(times_lwp/60./24., lwp_fraction, 'k', lw = 2, label = 'Total Cloud Frac')
    ax.plot(times_mcl/60./24., mcl_fraction[:,iz1], 'r', label = 'Cloud Frac ('+str(round(z[iz1], 2))+'m)')
    ax.plot(times_mcl/60./24., mcl_fraction[:,iz2], 'purple', label = 'Cloud Frac ('+str(round(z[iz2], 2))+'m)')
    ax.plot(times_mcl/60./24., mcl_fraction[:,iz3], 'b', label = 'Cloud Frac ('+str(round(z[iz3], 2))+'m)')
    ax.set_xlim([0, 10])
    ax.set_xlabel('Time (days)')
    ax.set_xticks(xrange(0, 11))
    ax.set_ylim([0, 40])
    ax.set_ylabel('Cloud Fraction')
    ax.legend()
    plt.savefig('../cloud_fraction_time_series_' + ID + '.png', dpi = 150)
    plt.show()

    # Hovmoller for entire simulation
    mcl_levels = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]
    fig = plt.figure(tight_layout = True)
    ax = fig.add_subplot(1, 1, 1)
    Hov = ax.contourf(times_mcl/60./24., z/1000., np.transpose(mcl_fraction)*100., cmap = 'bone_r', levels = mcl_levels, extend = 'max', vmax = 10.5)
    plt.colorbar(Hov, ax = ax, label = 'Cloud Fraction (%)')
    ax.contour(times_mcl/60./24., z/1000., np.transpose(mcl_fraction)*100., levels = mcl_levels, colors = ['k'])
    [ax.plot([it,it], [0,5], color = 'grey', ls = '--', lw = 2) for it in np.arange(0.25, 10.76, 0.5)]
    ax.set_xlim([0, 10])
    ax.set_xticks(xrange(0, 11))
    ax.set_xticklabels(hours)
    ax.set_xlabel('Time (days)')
    ax.set_ylim([0, 5])
    ax.set_ylabel('Height (km)')
    ax.set_title('Whole Domain, Horizontally Averaged Cloud Fraction')
    plt.savefig('../cloud_fraction_hovmoller_' + ID + '.png', dpi = 150)
    plt.show()
    
    ### repeat the above plots for the mean of last 4 days which we use for IC+forcing ###
    # Time series at different heights averaged over last 4 days
    dt_lwp = times_lwp[1] - times_lwp[0] # the time output frequency
    dt_mcl = times_mcl[1] - times_mcl[0]
    times_lwp_mean = np.arange(0, 1441, dt_lwp) # synthetic new times for a one day period
    times_mcl_mean = times_mcl[:143]
    lwp_fraction_mean = np.zeros_like(times_lwp_mean)
    mcl_fraction_mean = np.zeros((len(times_mcl_mean), len(z)))
    dt_i_lwp = 1440./dt_lwp # number of time indexes in a day
    dt_i_mcl = 1440./dt_mcl
    
    for i in xrange(len(times_lwp_mean)):
        J = [jt for jt in xrange(len(times_lwp)) if (jt >= 6*dt_i_lwp) and (times_lwp[jt]%1440. == times_lwp_mean[i])] # only the last four days and the current time of day "time[i]"
        lwp_fraction_mean[i] = np.nanmean([lwp_fraction[j] for j in J])
    
    for i in xrange(len(times_mcl_mean)):
        K = [kt for kt in xrange(len(times_mcl)) if (kt >= 6*dt_i_mcl) and (times_mcl[kt]%1440. == times_mcl_mean[i])]
        mcl_fraction_mean[i,:] = np.nanmean(np.array([mcl_fraction[k,:] for k in K]), axis = 0)
    
    fig = plt.figure(tight_layout = True)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(times_lwp_mean/60., lwp_fraction_mean, 'k', lw = 2, label = 'Total Cloud Frac')
    ax.plot(times_mcl_mean/60., mcl_fraction_mean[:,iz1], 'r', label = 'Cloud Frac ('+str(round(z[iz1], 2))+'m)')
    ax.plot(times_mcl_mean/60., mcl_fraction_mean[:,iz2], 'purple', label = 'Cloud Frac ('+str(round(z[iz2], 2))+'m)')
    ax.plot(times_mcl_mean/60., mcl_fraction_mean[:,iz3], 'b', label = 'Cloud Frac ('+str(round(z[iz3], 2))+'m)')
    ax.set_xlim([0, 24])
    ax.set_xlabel('Time (hrs)')
    ax.set_xticks(xrange(0, 25, 3))
    ax.set_ylim([0, 40])
    ax.set_ylabel('Cloud Fraction')
    ax.legend()
    plt.savefig('../cloud_fraction_time_series_mean_' + ID + '.png', dpi = 150)
    plt.show()
    
    # Hovmoller faveraged over last 4 days
    mcl_levels = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]
    fig = plt.figure(tight_layout = True)
    ax = fig.add_subplot(1, 1, 1)
    Hov = ax.contourf(times_mcl_mean/60., z/1000., np.transpose(mcl_fraction_mean), cmap = 'bone_r', levels = mcl_levels, extend = 'max', vmax = 10.5)
    plt.colorbar(Hov, ax = ax, label = 'Cloud Fraction (%)')
    ax.contour(times_mcl_mean/60., z/1000., np.transpose(mcl_fraction_mean), levels = mcl_levels, colors = ['k'])
    [ax.plot([it,it], [0,5], color = 'grey', ls = '--', lw = 2) for it in [6,18]]
    ax.set_xlim([0, 24])
    ax.set_xticks(xrange(0, 25, 3))
    ax.set_xlabel('Time (hrs)')
    ax.set_ylim([0, 5])
    ax.set_ylabel('Height (km)')
    ax.set_title('Whole Domain, Horizontally Averaged Cloud Fraction')
    plt.savefig('../cloud_fraction_hovmoller_mean_' + ID + '.png', dpi = 150)
    plt.show()
else:
    # Time series at different heights
    fig = plt.figure(tight_layout = True)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(times_lwp/60., lwp_fraction, 'k', lw = 2, label = 'Total Cloud Frac')
    ax.plot(times_mcl/60., mcl_fraction[:,iz1], 'r', label = 'Cloud Frac ('+str(round(z[iz1], 2))+'m)')
    ax.plot(times_mcl/60., mcl_fraction[:,iz2], 'purple', label = 'Cloud Frac ('+str(round(z[iz2], 2))+'m)')
    ax.plot(times_mcl/60., mcl_fraction[:,iz3], 'b', label = 'Cloud Frac ('+str(round(z[iz3], 2))+'m)')
    ax.set_xlim([0, 24])
    ax.set_xticks(xrange(0, 25, 3))
    ax.set_xlabel('Time (hrs)')
    ax.set_ylim([0, 40.])
    ax.set_ylabel('Cloud Fraction')
    ax.legend()
    plt.savefig('../cloud_fraction_time_series_' + ID + '.png', dpi = 150)
    plt.show()
    
    # Hovmoller
    mcl_levels = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]
    fig = plt.figure(tight_layout = True)
    ax = fig.add_subplot(1, 1, 1)
    Hov = ax.contourf(times_mcl/60., z/1000., np.transpose(mcl_fraction), cmap = 'bone_r', levels = mcl_levels, extend = 'max', vmax = 10.5)
    plt.colorbar(Hov, ax = ax, label = 'Cloud Fraction (%)')
    ax.contour(times_mcl/60., z/1000., np.transpose(mcl_fraction), levels = mcl_levels, colors = ['k'])
    ax.plot([6,6], [0,5], color = 'grey', ls = '--', lw = 2)
    ax.plot([18,18], [0,5], color = 'grey', ls = '--', lw = 2)
    ax.set_xlabel('Time (hrs)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0, 5])
    ax.set_xlim([0, 24])
    ax.set_xticks(xrange(0, 25, 3))
    ax.set_xticklabels(hours)
    ax.set_title('Whole Domain, Horizontally Averaged Cloud Fraction')
    plt.savefig('../cloud_fraction_hovmoller_' + ID + '.png', dpi = 150)
    plt.show()



