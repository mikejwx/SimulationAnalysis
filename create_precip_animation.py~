import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from netCDF4 import Dataset
import matplotlib as mpl
from STASH_keys import ls_rain_amt_key
import os
from multiprocessing import Pool

# Read the lwp data for the U10
lwp_10_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control_short/lwp_00.nc', 'r')
lwp_10 = lwp_10_nc.variables[ls_rain_amt_key][1:,:,:] # data every 15 minutes
time_10 = lwp_10_nc.variables['accum15'][:]*1.
lwp_10_nc.close()

# Read the lwp data for the U05
lwp_05_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/U05/lwp_00.nc', 'r')
lwp_05 = lwp_05_nc.variables[ls_rain_amt_key][1:,:,:] # data every 15 minutes
time_05 = lwp_05_nc.variables['accum15'][:]*1.
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
os.mkdir('./scratch_rain/')

def create_plot(it0):
    """
    Create a plot of the rain accumulation for the two experiments.
    """
    my_levels = [0.01, 0.5, 1., 2., 4., 8., 16., 32., 64.]
    my_colors = ['blue', 'cornflowerblue', 'olive', 'gold', 'orange', 'red', 'magenta', 'aliceblue']
    # The two experiments aren't at the same time output level, match up the indexes
    it1 = it0
    
    # Start the plots
    fig = plt.figure()
    
    ax0 = fig.add_subplot(2, 1, 1,adjustable = 'box', aspect = 1)
    pcm = ax0.contourf(X, Y, lwp_10[it0,:,:], extend = 'max',
                       colors = my_colors, levels = my_levels)
    ax0.contour(X, Y, lsm, colors = ['k'])
    ax0.set_title('U10, T+' + "{0:04}".format(int(time_10[it0]+240.)) + ' mins, \nH$_{land}$ = ' + "{0:03d}".format(int(H[it0])) + 'W/m$^{2}$, E$_{land}$ = ' +  "{0:03d}".format(int(E[it0])) + 'W/m$^{2}$')
    ax0.set_ylabel('y (km)')
    fig.colorbar(pcm, ax=ax0, label = 'Total Rain (mm)')
    
    ax1 = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
    pcm = ax1.contourf(X, Y, lwp_05[it1,:,:], extend = 'max',
                       colors = my_colors, levels = my_levels)
    ax1.contour(X, Y, lsm, colors = ['k'])
    ax1.set_title('U05, T+' + "{0:04}".format(int(time_05[it1]+240.)) + ' mins, \nH$_{land}$ = ' + "{0:03d}".format(int(H[it0])) + 'W/m$^{2}$, E$_{land}$ = ' +  "{0:03d}".format(int(E[it0])) + 'W/m$^{2}$')
    ax1.set_ylabel('y (km)')
    ax1.set_xlabel('x (km)')
    fig.colorbar(pcm, ax=ax1, label = 'Total Rain (mm)')
    
    plt.savefig('./scratch_rain/' + "{0:04}".format(int(time_10[it0])) + 'comparison.png', dpi = 150)
    plt.close('all')

p = Pool()
p.map(create_plot, xrange(len(time_10[1:])))
p.close()
p.join()

it0 = len(time_10[1:])-1

command1 = 'convert -delay 30 -loop 0 ./scratch_rain/*.png ../rain_animation_U10vU05.gif'
command2 = 'mv ./scratch_rain/' + "{0:04}".format(int(time_10[it0])) + 'comparison.png ../rain_comparison.png'
command3 = 'rm -rf ./scratch_rain/'
[os.system(command) for command in [command1, command2, command3]]
print 'Complete.'

