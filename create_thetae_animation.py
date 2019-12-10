import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from netCDF4 import Dataset
import matplotlib as mpl
from STASH_keys import theta_key, temp_key, pthe_key
from SkewT_archer import getThetaE
import os
from multiprocessing import Pool

# Let's plot the lowest model level theta_e
target_z = 2.
hours = ["{0:02}".format(h) for h in xrange(0, 16, 4)]
for hour in hours:
    print 'Start: hour ' + hour
    # Open the netCDFs containing the variables required to compute theta_e
    bouy_10_nc   = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control_short/bouy_' + hour + '.nc', 'r')
    fluxes_10_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control_short/fluxes_' + hour + '.nc', 'r')
    bouy_05_nc   = Dataset('/nerc/n02/n02/xb899100/CloudTrail/U05/bouy_' + hour + '.nc', 'r')
    fluxes_05_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/U05/fluxes_' + hour + '.nc', 'r')
    
    # Compute theta_e and store to a variable
    if hour == hours[0]:
        time10_key = [key for key in bouy_10_nc.variables.keys() if 'min' in key][0]
        time_10 = bouy_10_nc.variables[time10_key][1:]
        time05_key = [key for key in bouy_05_nc.variables.keys() if 'min' in key][0]
        time_05 = bouy_05_nc.variables[time05_key][1:]
        
        iz = np.where(np.abs(bouy_10_nc.variables['thlev_zsea_theta'][:] - target_z) == np.min(np.abs(bouy_10_nc.variables['thlev_zsea_theta'][:] - target_z)))[0][0]
        
        theta_e_10 = np.zeros_like(bouy_10_nc.variables[theta_key][1:,iz,:,:])
        theta_e_05 = np.zeros_like(bouy_05_nc.variables[theta_key][1:,iz,:,:])
        for it in xrange(theta_e_10.shape[0]):
            theta_e_10[it,:,:] = getThetaE(bouy_10_nc.variables[theta_key][1+it,iz,:,:], bouy_10_nc.variables[temp_key][1+it,iz,:,:], fluxes_10_nc.variables[pthe_key][1+it,iz,:,:], t_units = 'K', p_units = 'Pa')
        for it in xrange(theta_e_05.shape[0]):
            theta_e_05[it,:,:] = getThetaE(bouy_05_nc.variables[theta_key][1+it,iz,:,:], bouy_05_nc.variables[temp_key][1+it,iz,:,:], fluxes_05_nc.variables[pthe_key][1+it,iz,:,:], t_units = 'K', p_units = 'Pa')
    else:
        time_10 = np.concatenate((time_10, bouy_10_nc.variables[time10_key][:]), axis = 0)
        time_05 = np.concatenate((time_05, bouy_05_nc.variables[time05_key][:]), axis = 0)
        
        theta_e_10_t = np.zeros_like(bouy_10_nc.variables[theta_key][:,iz,:,:])
        theta_e_05_t = np.zeros_like(bouy_05_nc.variables[theta_key][:,iz,:,:])
        for it in xrange(theta_e_10_t.shape[0]):
            theta_e_10_t[it,:,:] = getThetaE(bouy_10_nc.variables[theta_key][it,iz,:,:], bouy_10_nc.variables[temp_key][it,iz,:,:], fluxes_10_nc.variables[pthe_key][it,iz,:,:], t_units = 'K', p_units = 'Pa')
        for it in xrange(theta_e_05_t.shape[0]):
            theta_e_05_t[it,:,:] = getThetaE(bouy_05_nc.variables[theta_key][it,iz,:,:], bouy_05_nc.variables[temp_key][it,iz,:,:], fluxes_05_nc.variables[pthe_key][it,iz,:,:], t_units = 'K', p_units = 'Pa')
        
        theta_e_10 = np.concatenate((theta_e_10, theta_e_10_t), axis = 0)
        theta_e_05 = np.concatenate((theta_e_05, theta_e_05_t), axis = 0)
    
    bouy_10_nc.close()
    fluxes_10_nc.close()
    bouy_05_nc.close()
    fluxes_05_nc.close()


# Define a horizontal coordinate system
X, Y = np.meshgrid(np.arange(0., 116000., 100.)/1000., np.arange(0., 31900., 100.)/1000.)

# Read the land sea mask
lsm_nc   = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
lsm = lsm_nc.variables['lsm'][0,0,:,:]
lsm_nc.close()

# Create synthetic surface flux diurnal cycle plot data
t0 = 720.
t = time_10[:]*1. + 240.
dt2 = 720./2.
H = np.where(((t > 360.)*(t < 1080.)), 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.5, 0.)
E = np.where(((t > 360.)*(t < 1080.)), 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.3, 0.)

# make a scratch directory
os.mkdir('./scratch_thetae/')

mean_te = np.nanmean(theta_e_10)
my_levels = np.array([int(mean_te) + inc for inc in xrange(-5, 6)])
def create_plot(it0):
    """
    Create a plot of the liquid water path for the two experiments.
    """
    my_cmap   = mpl.cm.get_cmap('magma')
    my_colors = np.array([my_cmap(float(i)/len(my_levels)) for i in xrange(len(my_levels))])
    # The two experiments aren't at the same time output level, match up the indexes
    it1 = np.where(time_05 == time_10[it0])[0][0]
    
    # Start the plots
    fig = plt.figure()
    
    ax0 = fig.add_subplot(2, 1, 1,adjustable = 'box', aspect = 1)
    pcm = ax0.contourf(X, Y, theta_e_10[it0,:,:], levels = my_levels,
                       extend = 'both', colors = my_colors)
    ax0.contour(X, Y, lsm, colors = ['r'])
    ax0.set_title('U10, T+' + "{0:04}".format(int(time_10[it0]+240.)) + ' mins, \nH$_{land}$ = ' + "{0:03d}".format(int(H[it0])) + 'W/m$^{2}$, E$_{land}$ = ' +  "{0:03d}".format(int(E[it0])) + 'W/m$^{2}$')
    fig.colorbar(pcm, ax=ax0, label = '$\\theta_{e}$ (K)')
    
    ax1 = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
    pcm = ax1.contourf(X, Y, theta_e_05[it1,:,:],  levels = my_levels,
                       extend = 'both', colors = my_colors)
    ax1.contour(X, Y, lsm, colors = ['r'])
    ax1.set_title('U05, T+' + "{0:04}".format(int(time_05[it1]+240.)) + ' mins, \nH$_{land}$ = ' + "{0:03d}".format(int(H[it0])) + 'W/m$^{2}$, E$_{land}$ = ' +  "{0:03d}".format(int(E[it0])) + 'W/m$^{2}$')
    fig.colorbar(pcm, ax=ax1, label = '$\\theta_{e}$ (K)')
    
    plt.savefig('./scratch_thetae/' + "{0:04}".format(int(time_10[it0])) + 'comparison.png', dpi = 150)

p = Pool()
p.map(create_plot, xrange(len(time_10)))
p.close()
p.join()

command1 = 'convert -delay 60 -loop 0 ./scratch_thetae/*.png ../thetae_animation_U10vU05.gif'
command2 = 'rm -rf ./scratch_thetae/'
[os.system(command) for command in [command1, command2]]
print 'Complete.'

