import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from netCDF4 import Dataset
import matplotlib as mpl
from STASH_keys import w_key
import os
from multiprocessing import Pool

# Let's plot the 350 m w
target_z = 350.
hours = ["{0:02}".format(h) for h in xrange(0, 16, 4)]
for hour in hours:
    print 'Start: hour ' + hour
    # Open the netCDFs containing the variables required to compute theta_e
    wind_10_nc   = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control_short/wind_' + hour + '.nc', 'r')
    wind_05_nc   = Dataset('/nerc/n02/n02/xb899100/CloudTrail/U05/wind_' + hour + '.nc', 'r')
    
    # Compute theta_e and store to a variable
    if hour == hours[0]:
        time10_key = [key for key in wind_10_nc.variables.keys() if 'min' in key][0]
        time_10 = wind_10_nc.variables[time10_key][1:]
        time05_key = [key for key in wind_05_nc.variables.keys() if 'min' in key][0]
        time_05 = wind_05_nc.variables[time05_key][1:]
        
        iz = np.where(np.abs(wind_10_nc.variables['thlev_zsea_theta'][:] - target_z) == np.min(np.abs(wind_10_nc.variables['thlev_zsea_theta'][:] - target_z)))[0][0]
        w_10 = wind_10_nc.variables[w_key][1:,iz,:,:]
        w_05 = wind_05_nc.variables[w_key][1:,iz,:,:]
    else:
        time_10 = np.concatenate((time_10, wind_10_nc.variables[time10_key][:]), axis = 0)
        time_05 = np.concatenate((time_05, wind_05_nc.variables[time05_key][:]), axis = 0)
        
        w_10 = np.concatenate((w_10, wind_10_nc.variables[w_key][:,iz,:,:]), axis = 0)
        w_05 = np.concatenate((w_05, wind_05_nc.variables[w_key][:,iz,:,:]), axis = 0)
    
    wind_10_nc.close()
    wind_05_nc.close()

# Define a horizontal coordinate system
X, Y = np.meshgrid(np.arange(0., 116000., 100.)/1000., np.arange(0., 31900., 100.)/1000.)

# Read the land sea mask
lsm_nc   = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
lsm = lsm_nc.variables['lsm'][0,0,:,:]
lsm_nc.close()

# Create synthetic surface flux diurnal cycle plot data
t0 = 720.
t = time_10[1:]*1. + 240.
dt2 = 720./2.
H = np.where(((t > 360.)*(t < 1080.)), 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.5, 0.)
E = np.where(((t > 360.)*(t < 1080.)), 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.3, 0.)

# make a scratch directory
os.mkdir('./scratch_w/')

def create_plot(it0):
    """
    Create a plot of the liquid water path for the two experiments.
    """
    my_levels = [l for l in np.linspace(-2.5, 2.5, 11) if l != 0.]
    # The two experiments aren't at the same time output level, match up the indexes
    it1 = np.where(time_05 == time_10[it0])[0][0]
    
    # Start the plots
    fig = plt.figure()
    
    ax0 = fig.add_subplot(2, 1, 1,adjustable = 'box', aspect = 1)
    pcm = ax0.contourf(X, Y, w_10[it0,:,:], extend = 'both',
                       cmap = 'bwr', levels = my_levels, vmin = -2.5, vmax = 2.5)
    ax0.contour(X, Y, lsm, colors = ['r'])
    ax0.set_title('U10, T+' + "{0:04}".format(int(time_10[it0]+240.)) + ' mins, \nH$_{land}$ = ' + "{0:03d}".format(int(H[it0])) + 'W/m$^{2}$, E$_{land}$ = ' +  "{0:03d}".format(int(E[it0])) + 'W/m$^{2}$')
    ax0.set_ylabel('y (km)')
    fig.colorbar(pcm, ax=ax0, label = 'w (m/s)')
    
    ax1 = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
    pcm = ax1.contourf(X, Y, w_05[it1,:,:], extend = 'both',
                       cmap = 'bwr', levels = my_levels, vmin = -2.5, vmax = 2.5)
    ax1.contour(X, Y, lsm, colors = ['r'])
    ax1.set_title('U05, T+' + "{0:04}".format(int(time_05[it1]+240.)) + ' mins, \nH$_{land}$ = ' + "{0:03d}".format(int(H[it0])) + 'W/m$^{2}$, E$_{land}$ = ' +  "{0:03d}".format(int(E[it0])) + 'W/m$^{2}$')
    ax1.set_ylabel('y (km)')
    ax1.set_xlabel('x (km)')
    fig.colorbar(pcm, ax=ax1, label = 'w (m/s)')
    
    plt.savefig('./scratch_w/' + "{0:04}".format(int(time_10[it0])) + 'comparison.png', dpi = 150)
    plt.close('all')

p = Pool()
p.map(create_plot, xrange(len(time_10[1:])))
p.close()
p.join()

command1 = 'convert -delay 60 -loop 0 ./scratch_w/*.png ../w_animation_U10vU05.gif'
command2 = 'rm -rf ./scratch_w/'
[os.system(command) for command in [command1, command2]]
print 'Complete.'

