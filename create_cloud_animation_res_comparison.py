import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from netCDF4 import Dataset
import matplotlib as mpl
from STASH_keys import lwp_key
import os
from multiprocessing import Pool
from coarse_graining import coarse_grain as cg, fix_x

# Read the lwp data for the Control (100 m)
lwp_0100_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control/lwp_00.nc', 'r')
lwp_0100m = lwp_0100_nc.variables[lwp_key][1:,:,:]*1000.
time0100_key = [tkey for tkey in lwp_0100_nc.variables.keys() if 'min' in tkey][0]
time_0100m = lwp_0100_nc.variables[time0100_key][:]
lwp_0100_nc.close()

# Read the lwp data for the Control (800 m)
lwp_0800_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control_0800m/lwp_00.nc', 'r')
lwp_0800m = lwp_0800_nc.variables[lwp_key][1:,:,:]*1000.
time0800_key = [tkey for tkey in lwp_0800_nc.variables.keys() if 'min' in tkey][0]
time_0800m = lwp_0800_nc.variables[time0800_key][:] + 240.
lwp_0800_nc.close()

# Read the lwp data for the Control (1600 m)
lwp_1600_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control_1600m/lwp_00.nc', 'r')
lwp_1600m = lwp_1600_nc.variables[lwp_key][1:,:,:]*1000.
time1600_key = [tkey for tkey in lwp_1600_nc.variables.keys() if 'min' in tkey][0]
time_1600m = lwp_1600_nc.variables[time1600_key][:] + 240.
lwp_1600_nc.close()

# Define a horizontal coordinate system
X1, Y1 = np.meshgrid(np.arange(1160)*0.1, np.arange(319)*0.1)
X8, Y8 = np.meshgrid(np.arange(146)*0.8, np.arange(40)*0.8)
X16, Y16 = np.meshgrid(np.arange(74)*1.6, np.arange(20)*1.6)

# Read the land sea mask
lsm_0100nc   = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
lsm0100 = lsm_0100nc.variables['lsm'][0,0,:,:]
lsm_0100nc.close()

lsm_0800nc   = Dataset('/work/n02/n02/xb899100/island_masks/lsm50_0800m.nc', 'r')
lsm0800 = lsm_0800nc.variables['lsm'][0,0,:,:]
lsm_0800nc.close()

lsm_1600nc   = Dataset('/work/n02/n02/xb899100/island_masks/lsm50_1600m.nc', 'r')
lsm1600 = lsm_1600nc.variables['lsm'][0,0,:,:]
lsm_1600nc.close()

# Create synthetic surface flux diurnal cycle plot data
t0 = 720.
t = time_0100m[1:]*1.
dt2 = 720./2.
H = np.where(((t > 360.)*(t < 1080.)), 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.5, 0.)
E = np.where(((t > 360.)*(t < 1080.)), 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.3, 0.)

# make a scratch directory
os.mkdir('./scratch_cloud/')

def create_plot(it0):
    """
    Create a plot of the liquid water path for the three experiments.
    """
    print 'Plotting T = ' + "{0:04d}".format(int(time_0100m[it0]))
    my_cmap   = mpl.cm.get_cmap('Greys_r')
    my_levels = [0., 1., 2.5, 5., 7.5, 10., 25., 50., 75., 100., 250., 500., 750., 1000., 2500., 5000.]
    my_colors = np.array([my_cmap(float(i)/len(my_levels)) for i in xrange(len(my_levels))])
    # The two experiments aren't at the same time output level, match up the indexes
    it1 = np.where(time_0800m == time_0100m[it0])[0][0]
    it2 = np.where(time_1600m == time_0100m[it0])[0][0]
    
    # Start the plots
    fig = plt.figure(tight_layout = True, figsize = (20, 12))
    
    # control, dx = 100 m
    axa = fig.add_subplot(3, 2, 1, adjustable = 'box', aspect = 1)
    dx0100 = axa.contourf(X1, Y1, lwp_0100m[it0,:,:], extend = 'max',
                       colors = my_colors, levels = my_levels)
    axa.contour(X1, Y1, lsm0100, colors = ['r'])
    axa.set_title('Control, dx = 100 m, T+' + "{0:04}".format(int(time_0100m[it0])) + ' mins, \nH$_{land}$ = ' + "{0:03d}".format(int(H[it0])) + 'W/m$^{2}$, E$_{land}$ = ' +  "{0:03d}".format(int(E[it0])) + 'W/m$^{2}$')
    axa.set_ylabel('y (km)')
    axa.set_xlabel('x (km)')
    fig.colorbar(dx0100, ax=axa, label = 'LWP (g/m$^{2}$)')
    
    # Control 100 m coarse grained to 800 m
    axc = fig.add_subplot(3, 2, 3, adjustable = 'box', aspect = 1)
    dx0800cg = axc.contourf(X8, Y8, fix_x(cg(lwp_0100m[it0,:,:], 8, 8)), extend = 'max', colors = my_colors, levels = my_levels)
    axc.contour(X8, Y8, lsm0800, colors = ['r'])
    axc.set_title('dx = 800 m (coarse grained from 100 m)')
    axc.set_ylabel('y (km)')
    axc.set_xlabel('x (km)')
    fig.colorbar(dx0800cg, ax=axc, label = 'LWP (g/m$^{2}$)')
    
    # 800 m
    axd = fig.add_subplot(3, 2, 4, adjustable = 'box', aspect = 1)
    dx0800 = axd.contourf(X8, Y8, lwp_0800m[it1,:,:], extend = 'max',
                       colors = my_colors, levels = my_levels)
    axd.contour(X8, Y8, lsm0800, colors = ['r'])
    axd.set_title('dx = 800 m')
    axd.set_ylabel('y (km)')
    axd.set_xlabel('x (km)')
    fig.colorbar(dx0800, ax=axd, label = 'LWP (g/m$^{2}$)')
    
    # Control 100 m coarse grained to 1600 m
    axe = fig.add_subplot(3, 2, 5, adjustable = 'box', aspect = 1)
    dx1600cg = axe.contourf(X16, Y16, fix_x(cg(lwp_0100m[it0,:,:], 16, 16)), extend = 'max', colors = my_colors, levels = my_levels)
    axe.contour(X16, Y16, lsm1600, colors = ['r'])
    axe.set_title('dx = 1600 m (coarse grained from 100 m)')
    axe.set_ylabel('y (km)')
    axe.set_xlabel('x (km)')
    fig.colorbar(dx1600cg, ax=axe, label = 'LWP (g/m$^{2}$)')
    
    # 1600 m
    axf = fig.add_subplot(3, 2, 6, adjustable = 'box', aspect = 1)
    dx1600 = axf.contourf(X16, Y16, lwp_1600m[it2,:,:], extend = 'max',
                       colors = my_colors, levels = my_levels)
    axf.contour(X16, Y16, lsm1600, colors = ['r'])
    axf.set_title('dx = 1600 m')
    axf.set_ylabel('y (km)')
    axf.set_xlabel('x (km)')
    fig.colorbar(dx1600, ax=axf, label = 'LWP (g/m$^{2}$)')
    
    #plt.show()
    
    plt.savefig('./scratch_cloud/' + "{0:04}".format(int(time_0100m[it0])) + 'comparison.png', dpi = 150)
    plt.close('all')


my_its = np.where((300. <= time_0100m)*(time_0100m <= 1080.))[0]
p = Pool()
p.map(create_plot, my_its)
p.close()
p.join()

command1 = 'convert -delay 25 -loop 0 ./scratch_cloud/*.png ../cloud_animation_comparing_resolutions.gif'
command2 = 'rm -rf ./scratch_cloud/'
[os.system(command) for command in [command1, command2]]


