import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl

"""
Code to plot cloud masks for my cloud trail simulation
"""

lwp_key = u'STASH_m01s30i405'

# Read the netCDF
print 'Reading...'
lwp = Dataset('lwp_00.nc', 'r')
lsm = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')

lwp_data = lwp.variables[lwp_key][:]*1000.
lsm_data = lsm.variables['lsm'][:]*1.

times = lwp.variables['min5_0'][:]*1.
lwp.close()
lsm.close()

X = np.arange(0., 116000., 100.)/1000.
Y = np.arange(0., 31900., 100.)/1000.
X, Y = np.meshgrid(X, Y)

albedo = np.where((lwp_data >= 5.), 20.5*np.log(lwp_data)-17.7, 5.)
albedo = np.where((albedo >= 80.), 80., albedo)/100.

### Create synthetic surface flux diurnal cycle plot data ###
t0 = 720.
t = np.arange(0., 1440.1, 10.)
dt2 = 720./2.
H = 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.5
E = 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.3


for it in xrange(times.shape[0]):
    print 'Plotting T+' + "{0:04d}".format(int(times[it])) + ' mins'
    fig = plt.figure(figsize = (18, 6))
    gs = mpl.gridspec.GridSpec(nrows = 1, ncols = 2, height_ratios = [1], width_ratios = [2,1])
    ax = fig.add_subplot(gs[0], adjustable = 'box', aspect = 'equal')
    #top left
    ax.contourf(X-np.max(X), Y+np.max(Y), albedo[it,:,:], cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
    #top centre
    ax.contourf(X, Y+np.max(Y), albedo[it,:,:], cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
    #top right
    ax.contourf(X+np.max(X), Y+np.max(Y), albedo[it,:,:], cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
    ax.contour(X+np.max(X), Y+np.max(Y), lsm_data[0,0,:,:], levels = [1e-08], colors = ['r'], lw = [5])
    #centre left
    ax.contourf(X-np.max(X), Y, albedo[it,:,:], cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
    #centre centre
    ax.contourf(X, Y, albedo[it,:,:], cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
    ax.contour(X, Y, lsm_data[0,0,:,:], levels = [1e-08], colors = ['r'], lw = [5])
    #centre right
    ax.contourf(X+np.max(X), Y, albedo[it,:,:], cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
    #bottom left
    ax.contourf(X-np.max(X), Y-np.max(Y), albedo[it,:,:], cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
    ax.contour(X-np.max(X), Y-np.max(Y), lsm_data[0,0,:,:], levels = [1e-08], colors = ['r'])
    #bottom centre
    alb = ax.contourf(X, Y-np.max(Y), albedo[it,:,:], cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
    axins = inset_axes(ax, width = "50%", height = "5%", loc = 8, bbox_to_anchor = (0.0, -0.26, 1, 1), bbox_transform = ax.transAxes, borderpad = 0.)
    plt.colorbar(alb, cax = axins, orientation = 'horizontal', label = 'synthetic albedo')
    #bottom right
    ax.contourf(X+np.max(X), Y-np.max(Y), albedo[it,:,:], cmap = 'Greys_r', levels = np.linspace(0., 1., 21))
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.plot([0, 0, np.max(X), np.max(X), 0], [0, np.max(Y), np.max(Y), 0, 0], 'r--')
    ax.set_title("{0:04d}".format(int(times[it]))+'mins')

    ax = fig.add_subplot(gs[1], adjustable = 'box', aspect = 1)
    ax.plot(t, H, 'r', lw = 2, label = 'Sensible HF')
    ax.plot(t, E, 'b', lw = 1, label = 'Latent HF')
    ax.set_xlim([0, 1440])
    ax.set_xticks(np.arange(0, 1440.1, 180))
    ax.plot([times[it], times[it]], [0., 500.], 'k--')
    ax.set_xlabel('Time (mins)')
    ax.set_ylabel('Land-Surface Heat Flux (W m$^{-2}$)')
    lgd=plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.05), prop={'size': 12})
    fig.subplots_adjust(left = 0.0625, right = 0.87, bottom = 0.1, top = 0.9, wspace = 0.21, hspace = 0.31)
    plt.savefig('./cloudmasks/'+"{0:04d}".format(int(times[it]))+'.png', dpi = 100, bbox_extra_artists=(lgd,))
    if it == 1:
        plt.show()
    plt.close()


