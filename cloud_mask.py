import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from multiprocessing import Pool
import os

"""
Code to plot synthetic albedos for my cloud trail simulation i.e. using liquid 
water path to estimate the albedo.
"""

lwp_key = u'STASH_m01s30i405'

# Read the netCDF
print 'Reading...'
lwp = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control_H250Const/lwp_00.nc', 'r')
lsm = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')

lwp_data = lwp.variables[lwp_key][:]*1000.
lsm_data = lsm.variables['lsm'][:]*1.
lwp_times_key = [key for key in lwp.variables.keys() if 'min' in key][0]
times = lwp.variables[lwp_times_key][:]*1. + 240.
lwp.close()
lsm.close()

X = np.arange(0., 116000., 100.)/1000.
Y = np.arange(0., 31900., 100.)/1000.
X, Y = np.meshgrid(X, Y)
"""
albedo = np.where((lwp_data >= 5.), 17.7*np.log(lwp_data)-20.5, 5.)
albedo = np.where((albedo >= 80.), 80., albedo)/100.
"""

mask = np.where((lwp_data >= 5.), 1., 0.)

### Create synthetic surface flux diurnal cycle plot data ###
t0 = 720.
t = times*1.
dt2 = 720./2.
H = np.where((t >= (t0-dt2))*(t <= (t0+dt2)), 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.5, 0.)
E = np.where((t >= (t0-dt2))*(t <= (t0+dt2)), 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.3, 0.)

x_offset = np.max(X) + 0.1
y_offset = np.max(Y) + 0.1
gs = mpl.gridspec.GridSpec(nrows = 1, ncols = 2, height_ratios = [1], width_ratios = [3,1])
vmin = 0.0
vmax = 1.0
my_levels = [0.0, 1e-16, 1.0]#np.linspace(0.0, 0.25, 11)
alpha = 0.8

def do_plot(it):
    print 'Plotting T+' + "{0:04d}".format(int(times[it])) + ' mins'
    print 'H = ' + str(H[it])
    fig = plt.figure(figsize = (21, 7), tight_layout = True)
    ax = fig.add_subplot(gs[0], adjustable = 'box', aspect = 1)
    #top left
    ax.contourf(X-x_offset, Y+y_offset, mask[it,:,:], colors = ['w', 'k'], levels = my_levels, extend = 'max', vmin = vmin, vmax = vmax, alpha = alpha)
    ax.contour(X-x_offset, Y+y_offset, lsm_data[0,0,:,:], levels = [1e-08], colors = ['darkred'], linewidths = [2])
    #top centre
    ax.contourf(X, Y+y_offset, mask[it,:,:], colors = ['w', 'k'], levels = my_levels, extend = 'max', vmin = vmin, vmax = vmax, alpha = alpha)
    ax.contour(X, Y+y_offset, lsm_data[0,0,:,:], levels = [1e-08], colors = ['darkred'], linewidths = [2])
    #top right
    ax.contourf(X+x_offset, Y+y_offset, mask[it,:,:], colors = ['w', 'k'], levels = my_levels, extend = 'max', vmin = vmin, vmax = vmax, alpha = alpha)
    ax.contour(X+x_offset, Y+y_offset, lsm_data[0,0,:,:], levels = [1e-08], colors = ['darkred'], linewidths = [2])
    #centre left
    ax.contourf(X-x_offset, Y, mask[it,:,:], colors = ['w', 'k'], levels = my_levels, extend = 'max', vmin = vmin, vmax = vmax, alpha = alpha)
    ax.contour(X-x_offset, Y, lsm_data[0,0,:,:], levels = [1e-08], colors = ['darkred'], linewidths = [2])
    #centre centre
    ax.contourf(X, Y, mask[it,:,:], colors = ['w', 'k'], levels = my_levels, extend = 'max', vmin = vmin, vmax = vmax, alpha = alpha)
    ax.contour(X, Y, lsm_data[0,0,:,:], levels = [1e-08], colors = ['darkred'], linewidths = [2])
    #centre right
    ax.contourf(X+x_offset, Y, mask[it,:,:], colors = ['w', 'k'], levels = my_levels, extend = 'max', vmin = vmin, vmax = vmax, alpha = alpha)
    ax.contour(X+x_offset, Y, lsm_data[0,0,:,:], levels = [1e-08], colors = ['darkred'], linewidths = [2])
    #bottom left
    ax.contourf(X-x_offset, Y-y_offset, mask[it,:,:], colors = ['w', 'k'], levels = my_levels, extend = 'max', vmin = vmin, vmax = vmax, alpha = alpha)
    ax.contour(X-x_offset, Y-y_offset, lsm_data[0,0,:,:], levels = [1e-08], colors = ['darkred'], linewidths = [2])
    #bottom centre
    alb = ax.contourf(X, Y-y_offset, mask[it,:,:], colors = ['w', 'k'], levels = my_levels, extend = 'max', vmin = vmin, vmax = vmax, alpha = alpha)
    #axins = inset_axes(ax, width = "50%", height = "5%", loc = 8, bbox_to_anchor = (0.0, -0.2, 1, 1), bbox_transform = ax.transAxes, borderpad = 0.)
    #plt.colorbar(alb, cax = axins, orientation = 'horizontal', label = 'Synthetic Albedo', ticks = my_levels)
    ax.contour(X, Y-y_offset, lsm_data[0,0,:,:], levels = [1e-08], colors = ['darkred'], linewidths = [2])
    #bottom right
    ax.contourf(X+np.max(X), Y-np.max(Y), mask[it,:,:], colors = ['w', 'k'], levels = my_levels, extend = 'max', vmin = vmin, vmax = vmax, alpha = alpha)
    ax.contour(X+np.max(X), Y-np.max(Y), lsm_data[0,0,:,:], levels = [1e-08], colors = ['darkred'], linewidths = [2])
    
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.plot([0, 0, np.max(X), np.max(X), 0], [0, np.max(Y), np.max(Y), 0, 0], color = 'indigo', ls = '--', lw = 2)
    ax.set_title("T+{0:04d}".format(int(times[it]))+'mins\n H = ' + str(int(round(H[it],0))) + ' W m$^{-2}$, E = ' + str(int(round(E[it],0))) + ' W m$^{-2}$')
    
    ax = fig.add_subplot(gs[1], adjustable = 'box', aspect = (2./3.)*24./500.)
    ax.plot(t/60., H, 'r', lw = 2, label = 'H')
    ax.plot(t/60., E, 'b', lw = 1, label = 'E')
    ax.set_xlim([0, 24])
    ax.set_xticks(np.arange(0, 1440.1, 180)/60.)
    ax.plot([times[it]/60., times[it]/60.], [0., 500.], 'k--')
    ax.set_xlabel('Time (hrs)')
    ax.set_ylabel('Land-Surface Heat Flux (W m$^{-2}$)')
    plt.legend(loc='upper center', ncol = 2)
    #fig.subplots_adjust(left = 0.0625, right = 0.87, bottom = 0.1, top = 0.9, wspace = 0.21, hspace = 0.31)
    
    plt.savefig('./scratch_cloud/'+"{0:04d}".format(int(times[it]))+'.png')
    plt.close('all')


# make a scratch directory
if not os.path.exists('./scratch_cloud/'):
    print 'Creating scratch directory'
    os.mkdir('./scratch_cloud/')
else:
    print 'Path already exists\nDeleting scratch directory and creating a clean empty one'
    os.system('rm -rf ./scratch_cloud/')
    os.mkdir('./scratch_cloud/')

p = Pool()
p.map(do_plot, xrange(len(times)))
p.close()
p.join()

command1 = 'convert -delay 15 -loop 0 ./scratch_cloud/*.png ../cloud_animation_Control_H250Const.gif'
command2 = 'rm -rf ./scratch_cloud/'
print 'Creating the animation'
os.system(command1)
print 'Deleting the scratch directory'
os.system(command2)

