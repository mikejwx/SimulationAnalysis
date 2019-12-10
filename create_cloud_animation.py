import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from netCDF4 import Dataset
import matplotlib as mpl
from STASH_keys import lwp_key
import os
from multiprocessing import Pool

# Read the lwp data for the U10
lwp_10_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control_short/lwp_00.nc', 'r')
lwp_10 = lwp_10_nc.variables[lwp_key][1:,:,:]*1000. # data every 5 minutes
time10_key = [key for key in lwp_10_nc.variables.keys() if 'min' in key][0]
time_10 = lwp_10_nc.variables[time10_key][:]
lwp_10_nc.close()

# Read the lwp data for the U05
lwp_05_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/U05/lwp_00.nc', 'r')
lwp_05 = lwp_05_nc.variables[lwp_key][1:,:,:]*1000. # data every 1 minute
time05_key = [key for key in lwp_05_nc.variables.keys() if 'min' in key][0]
time_05 = lwp_05_nc.variables[time05_key][:]
lwp_05_nc.close()

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
os.mkdir('./scratch_cloud/')

def create_plot(it0):
    """
    Create a plot of the liquid water path for the two experiments.
    """
    my_cmap   = mpl.cm.get_cmap('Greys_r')
    my_levels = [0., 1., 2.5, 5., 7.5, 10., 25., 50., 75., 100., 250., 500., 750., 1000., 2500., 5000.]
    my_colors = np.array([my_cmap(float(i)/len(my_levels)) for i in xrange(len(my_levels))])
    # The two experiments aren't at the same time output level, match up the indexes
    it1 = np.where(time_05 == time_10[it0])[0][0]
    
    # Start the plots
    fig = plt.figure()
    
    ax0 = fig.add_subplot(2, 1, 1,adjustable = 'box', aspect = 1)
    pcm = ax0.contourf(X, Y, lwp_10[it0,:,:], extend = 'max',
                       colors = my_colors, levels = my_levels)
    ax0.contour(X, Y, lsm, colors = ['r'])
    ax0.set_title('U10, T+' + "{0:04}".format(int(time_10[it0]+240.)) + ' mins, \nH$_{land}$ = ' + "{0:03d}".format(int(H[it0])) + 'W/m$^{2}$, E$_{land}$ = ' +  "{0:03d}".format(int(E[it0])) + 'W/m$^{2}$')
    ax0.set_ylabel('y (km)')
    fig.colorbar(pcm, ax=ax0, label = 'LWP (g/m$^{2}$)')
    
    ax1 = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
    pcm = ax1.contourf(X, Y, lwp_05[it1,:,:], extend = 'max',
                       colors = my_colors, levels = my_levels)
    ax1.contour(X, Y, lsm, colors = ['r'])
    ax1.set_title('U05, T+' + "{0:04}".format(int(time_05[it1]+240.)) + ' mins, \nH$_{land}$ = ' + "{0:03d}".format(int(H[it0])) + 'W/m$^{2}$, E$_{land}$ = ' +  "{0:03d}".format(int(E[it0])) + 'W/m$^{2}$')
    ax1.set_ylabel('y (km)')
    ax1.set_xlabel('x (km)')
    fig.colorbar(pcm, ax=ax1, label = 'LWP (g/m$^{2}$)')
    
<<<<<<< HEAD
    plt.show()
    
    #plt.savefig('./scratch_cloud/' + "{0:04}".format(int(time_10[it0])) + 'comparison.png', dpi = 150)
    #plt.close('all')
"""
=======
    plt.savefig('./scratch_cloud/' + "{0:04}".format(int(time_10[it0])) + 'comparison.png', dpi = 150)
    plt.close('all')

>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
p = Pool()
p.map(create_plot, xrange(len(time_10[1:])))
p.close()
p.join()

command1 = 'convert -delay 10 -loop 0 ./scratch_cloud/*.png ../cloud_animation_U10vU05.gif'
command2 = 'rm -rf ./scratch_cloud/'
[os.system(command) for command in [command1, command2]]
<<<<<<< HEAD
"""
=======

>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
