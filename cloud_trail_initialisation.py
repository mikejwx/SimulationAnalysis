import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
from netCDF4 import Dataset
from scipy import interpolate, integrate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from analysis_tools import send_email, summary, transform_winds
from datetime import datetime as dt

## First read some data
# STASH codes
print '[' + dt.now().strftime('%H:%M:%S') + '] Generating STASH codes...'
u_key     = u'STASH_m01s00i002'
v_key     = u'STASH_m01s00i003'
theta_key = u'STASH_m01s00i004'
w_key     = u'STASH_m01s00i150'
q_key     = u'STASH_m01s00i010'
lwp_key   = u'STASH_m01s30i405'
z_key     = u'thlev_zsea_theta'
z_key_rho = u'rholev_zsea_rho'
lsm_key   = u'lsm'

lsm  = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
print '[' + dt.now().strftime('%H:%M:%S') + '] Getting the land-sea mask...'
my_lsm = lsm.variables[lsm_key][0,0,:,:]*1.
lsm.close()

### Create synthetic surface flux diurnal cycle plot data ###
t0 = 720.
t = np.arange(0., 1440.1, 10.)
dt2 = 720./2.
H = 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.5
E = 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.3

for hour in ['06', '09']:

    # Open the netCDFs
    print '[' + dt.now().strftime('%H:%M:%S') +'] Opening our netCDFs...'
    wind_nc = Dataset('../wind_'+hour+'.nc', 'r')
    lwp_nc  = Dataset('../lwp_00.nc', 'r')

    # Get an array of heights
    print '[' + dt.now().strftime('%H:%M:%S') +'] Get array of heights...'
    z = wind_nc.variables[z_key][:]

    print '[' + dt.now().strftime('%H:%M:%S') +'] Get array of times...'
    times     = wind_nc.variables[u'min10_0'][:]
    times_lwp = lwp_nc.variables[u'min5_0'][:]

    # wind components
    print '[' + dt.now().strftime('%H:%M:%S') +'] Reading winds...'
    my_u = wind_nc.variables[u_key][:]
    my_v = wind_nc.variables[v_key][:]
    my_w = wind_nc.variables[w_key][:]

    # lwp
    print '[' + dt.now().strftime('%H:%M:%S') + '] Reading liquid water path...'
    its = np.array([i for i in xrange(len(times_lwp)) if times_lwp[i] in times])
    my_lwp = lwp_nc.variables[lwp_key][its,:,:]

    height_1 = 10.0 # Height used for near surface winds
    height_2 = 375.0 # Height used for vertical velocity (approx = z_i/2)
    iz1 = np.where(np.min(np.abs(z - height_1)) == np.abs(z - height_1))[0][0]
    iz2 = np.where(np.min(np.abs(z - height_2)) == np.abs(z - height_2))[0][0]

    # Get the flow-relative winds and the heights at which we want to plot them
    print '[' + dt.now().strftime('%H:%M:%S') +'] Transforming the winds to the rotated coordinate system'
    my_s, my_n = transform_winds(my_u, my_v)

    # manufacture some horizontal grid coordinates
    X = np.arange(0., 116000., 100.)/1000.
    Y = np.arange(0., 31900., 100.)/1000.
    X, Y = np.meshgrid(X, Y)

    #==============================================================================#
    #                                                                              #
    # Plot flow-relative wind anomalies                                            #
    #                                                                              #
    #==============================================================================#
    my_cmap = mpl.cm.get_cmap('Greys')
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="w")
    for it in xrange(len(times)):
        
        print '[' + dt.now().strftime('%H:%M:%S') +'] Plotting time ' + str(times[it])
        
        fig = plt.figure(figsize = (10, 11))
        
        ax = fig.add_subplot(4,1,1, adjustable = 'box', aspect = 1.)
        CS1 = ax.contourf(X, Y, my_s[it,iz1,:,:] - my_s[it,iz1,:].mean(), cmap = 'bwr', levels = [x for x in np.linspace(-2.5, 2.5, 11) if x != 0.], extend = 'both')
        ax.contourf(X, Y, my_lwp[it,:,:], colors = ['k'], levels = [1e-16, 1e03], extend = 'max', alpha = 0.5)
        fig.colorbar(CS1, ax = ax)
        ax.contour(X, Y, my_lsm, levels = [1e-08], colors = ['k'])
        ax.set_title('Along-flow Anomaly (m/s) at '+str(round(z[iz1], 2))+'m')
        ax.set_ylabel('y (km)')
        ax.set_xticklabels([''])
        ax.annotate('A)', (0, 1), xytext=(5,-5),xycoords='axes fraction',
            textcoords='offset points',ha='left', va='top', weight = 'bold', bbox=bbox_props)
        
        ax = fig.add_subplot(4,1,2, adjustable = 'box', aspect = 1.)
        CN1 = ax.contourf(X, Y, my_n[it,iz1,:,:] - my_n[it,iz1,:,:].mean(), cmap = 'bwr', levels = [x for x in np.linspace(-2.5, 2.5, 11) if x != 0.], extend = 'both')
        ax.contourf(X, Y, my_lwp[it,:,:]*1000., colors = ['k'], levels = [1e-16, 1e03], extend = 'max', alpha = 0.5)
        fig.colorbar(CN1, ax = ax)
        ax.contour(X, Y, my_lsm, levels = [1e-08], colors = ['k'])
        ax.set_title('Cross-flow Anomaly (m/s) at '+str(round(z[iz1], 2))+'m')
        ax.set_ylabel('y (km)')
        ax.set_xticklabels([''])
        ax.annotate('B)', (0, 1), xytext=(5,-5),xycoords='axes fraction',
            textcoords='offset points',ha='left', va='top', weight = 'bold', bbox=bbox_props)
        
        ax = fig.add_subplot(4,1,3, adjustable = 'box', aspect = 1.)
        W1 = ax.contourf(X, Y, my_w[it,iz2,:,:], cmap = 'bwr', levels = [x for x in np.linspace(-2.5, 2.5, 11) if x != 0.], extend = 'both')
        ax.contourf(X, Y, my_lwp[it,:,:]*1000., colors = ['k'], levels = [1e-16, 1e03], extend = 'max', alpha = 0.5)
        fig.colorbar(W1, ax = ax)
        ax.contour(X, Y, my_lsm, levels = [1e-08], colors = ['k'])
        ax.set_title('Vertical Velocity (m/s) at '+str(round(z[iz2], 2))+'m')
        ax.set_ylabel('y (km)')
        ax.set_xlabel('x (km)')
        ax.annotate('C)', (0, 1), xytext=(5,-5),xycoords='axes fraction',
            textcoords='offset points',ha='left', va='top', weight = 'bold', bbox=bbox_props)
        
        ax = fig.add_subplot(4,1,4, adjustable = 'box')
        ax.plot(t, H, 'r', lw = 2, label = 'Sensible HF')
        ax.plot(t, E, 'b', lw = 1, label = 'Latent HF')
        ax.set_xlim([0, 1440])
        ax.set_xticks(np.arange(0, 1440.1, 180))
        ax.plot([times[it], times[it]], [0., 500.], 'k--', lw = 2)
        ax.set_xlabel('Time (mins)')
        ax.set_ylabel('Land-Surface Heat Flux (W m$^{-2}$)')
        plt.legend(loc='upper right',prop={'size': 12})
        ax.annotate('D)', (0, 1), xytext=(5,-5),xycoords='axes fraction',
            textcoords='offset points',ha='left', va='top', weight = 'bold', bbox=bbox_props)
        
        fig.suptitle('Flow Relative Anomalies at T+'+str(int(times[it]))+'mins', fontsize = 20)
        fig.subplots_adjust(top = 0.88)
        plt.savefig('../CT_initialisation/FRWA_lwp'+"{0:04d}".format(int(times[it]))+'mins.png', dpi = 100)
        plt.close('all')

    wind_nc.close()
    lwp_nc.close()


