"""
Want to plot arrows to represent the wind speed and filled contour for fields
to discuss buoyancy such as theta, w, and qv.

Want to do this for several heights:
e.g. in the boundary layer, in the lower cloud layer, and in the upper cloud
layer.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from netCDF4 import Dataset
from scipy import interpolate, integrate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from analysis_tools import send_email, regrid, summary, transform_winds

## First read some data
# STASH codes
print 'Generating STASH codes...'
u_key = u'STASH_m01s00i002'
v_key = u'STASH_m01s00i003'
theta_key = u'STASH_m01s00i004'
w_key = u'STASH_m01s00i150'
mv_key = u'STASH_m01s00i391'
z_key = u'thlev_zsea_theta'
z_key_rho = u'rholev_zsea_rho'
lsm_key = u'lsm'

lsm  = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
print ' Getting the land-sea mask...'
my_lsm = lsm.variables[lsm_key][0,0,:,:]*1.
lsm.close()

hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]

# Choose target heights for the thermodynamic slices
my_heights = [1200.0, 700.0, 200.0] # m

# Arguments for plotting the thermodynamics
color_maps = {'theta' : 'hot',
              'w'     : 'bwr',
              'mv'    : 'BrBG'}
label = {'theta' : '$\\theta$ (K)',
         'w' : 'w (m s$^{-1}$)',
         'mv' : 'm$_v$ (g kg$^{-1}$)'}

### Create synthetic surface flux diurnal cycle plot data ###
t0 = 720.
t = np.arange(0., 1440.1, 10.)
dt2 = 720./2.
H = 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.5
E = 250.*np.cos(0.5*np.pi*(t0 - t)/dt2)**1.3

for hour in hours:
    # Open the netCDFs
    print 'Opening our netCDFs...'
    u    = Dataset('u_'+hour+'.nc', 'r')
    v    = Dataset('v_'+hour+'.nc', 'r')
    bouy = Dataset('bouy_'+hour+'.nc', 'r')
    mr   = Dataset('mr_'+hour+'.nc', 'r')
    
    # Get an array of heights
    print 'Get array of heights...'
    z = bouy.variables[z_key][:]*1.
    
    print 'Get array of times...'
    times = bouy.variables[u'min10_0'][:]*1.
    
    # u-wind, needs to be regridded
    print 'Regridding u...'
    my_u = regrid(bouy, u, u_key)[:]*1.
    # v-wind
    print 'Regridding v...'
    my_v = regrid(bouy, v, v_key)[:]*1.
    
    # theta
    print 'Reading theta...'
    my_theta = bouy.variables[theta_key][:]*1.
    # w-wind
    print 'Reading w...'
    my_w = bouy.variables[w_key][:]*1.
    
    # mixing ratio
    print 'Regridding mv...'
    my_mv = regrid(bouy, mr, mv_key)[:]*1.
    
    if hour == '00':
        # Get the indexes for the heights at which we want to plot the flow-relative winds
        height_1 = 200.0
        height_2 = 700.0
        height_3 = 500.0
        iz1 = np.where(np.min(np.abs(z - height_1)) == np.abs(z - height_1))[0][0]
        iz2 = np.where(np.min(np.abs(z - height_2)) == np.abs(z - height_2))[0][0]
        iz3 = np.where(np.min(np.abs(z - height_3)) == np.abs(z - height_3))[0][0]
    # List my variables
    variables = {'theta' : my_theta,
                 'w'     : my_w,
                 'mv'    : my_mv*1000.}
    
    # Get the flow-relative winds and the heights at which we want to plot them
    print 'Transforming the winds to the rotated coordinate system'
    my_s, my_n = transform_winds(my_u, my_v)
    
    # loop through all the times in this netCDF on the hour
    its = np.where(times%60. == 0)[0]
    its = [it for it in its if it != 0]
    for it in its:
        # manufacture some horizontal grid coordinates
        X = np.arange(0., 116000., 100.)/1000.
        Y = np.arange(0., 31900., 100.)/1000.
        X, Y = np.meshgrid(X, Y)
        
        for key in variables.keys():
            
            # Plot the individual snapshots individually
            fig = plt.figure(tight_layout=True, figsize=(8, 8))
            for height in my_heights:
                
                # Find the index for that height
                iz = np.where(np.min(np.abs(z - height)) == np.abs(z - height))[0][0]
                
                # Scale the colorbar so that the middle is the mean of that field
                # and the scale extends symmetrically so that the extremes are the 
                # maximum distance from the mean
                dv = 0.5*np.subtract(*np.percentile(variables[key][:,iz,:,:], [99.99, 0.01]))
                vmean = np.nanmean(variables[key][:,iz,:,:])
                
                titles = {'theta' : 'Potential Temperature (K) at T+' + str(round(times[it],2)) + 'mins',
                    'w'     : 'Vertical Velocity (m s$^{-1}$) at T+' + str(round(times[it],2)) + 'mins',
                    'mv'    : 'Water Vapour Mixing Ratio (g kg$^{-1}$) at T+' + str(round(times[it],2)) + 'mins'}
                
                pp = my_heights.index(height)
                ax = fig.add_subplot(len(my_heights), 1, pp+1, adjustable = 'box', aspect = 1.)
                T = ax.contourf(X, Y, variables[key][it,iz,:,:], levels = np.linspace(vmean-dv, vmean+dv, 9), vmin = vmean-dv, vmax = vmean+dv, extend = 'both', cmap = color_maps[key])
                fig.colorbar(T, ax = ax, label = label[key])
                ax.contour(X, Y, my_lsm, levels = [1e-08], colors = ['k'], lw = [10])
                ax.set_title(titles[key] + ' at ' + str(round(z[iz],2)) + ' m')
                ax.set_xlabel('x (km)')
                ax.set_ylabel('y (km)')
                plt.tight_layout()
        plt.savefig(key + '_'+"{0:04d}".format(int(times[it]))+'mins.png', dpi = 100)
        plt.close('all')
        
        #==============================================================================#
        #                                                                              #
        # Plot flow-relative wind anomalies                                            #
        #                                                                              #
        #==============================================================================#
        
        print 'Plotting flow-relative wind anomalies and vertical velocities'
                
        fig = plt.figure(figsize = (16, 8))
        gs = mpl.gridspec.GridSpec(nrows = 3, ncols = 2, height_ratios = [1,1,1], width_ratios = [1,1])

        ax = fig.add_subplot(gs[0,0], adjustable = 'box', aspect = 1.)
        CS1 = ax.contourf(X, Y, my_s[it,iz1,:,:] - my_s[it,iz1,:].mean(), cmap = 'bwr', levels = np.linspace(-2., 2., 9), extend = 'both')
        axins = inset_axes(ax, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax.transAxes, borderpad = 0.)
        plt.colorbar(CS1, cax = axins)
        ax.contour(X, Y, my_lsm, levels = [1e-08], colors = ['k'], lw = [10])
        ax.set_title('Along-flow Anomaly (m/s) at '+str(round(z[iz1], 2))+'m')
        ax.set_ylabel('y (km)')
        ax.set_xticklabels([''])
        ax.annotate('A)', (0, 1), xytext=(5,-5),xycoords='axes fraction',
            textcoords='offset points',ha='left', va='top')

        ax = fig.add_subplot(gs[1,0], adjustable = 'box', aspect = 1.)
        CN1 = ax.contourf(X, Y, my_n[it,iz1,:,:] - my_n[it,iz1,:,:].mean(), cmap = 'bwr', levels = np.linspace(-2., 2., 9), extend = 'both')
        axins = inset_axes(ax, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax.transAxes, borderpad = 0.)
        plt.colorbar(CN1, cax = axins)
        ax.contour(X, Y, my_lsm, levels = [1e-08], colors = ['k'], lw = [10])
        ax.set_title('Cross-flow Anomaly (m/s) at '+str(round(z[iz1], 2))+'m')
        ax.set_ylabel('y (km)')
        ax.set_xticklabels([''])

        ax = fig.add_subplot(gs[0,1], adjustable = 'box', aspect = 1.)
        CS2 = ax.contourf(X, Y, my_s[it,iz2,:,:] - my_s[it,iz2,:,:].mean(), cmap = 'bwr', levels = np.linspace(-2., 2., 9), extend = 'both')
        axins = inset_axes(ax, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax.transAxes, borderpad = 0.)
        plt.colorbar(CS2, cax = axins)
        ax.contour(X, Y, my_lsm, levels = [1e-08], colors = ['k'], lw = [10])
        ax.set_title('Along-flow Anomaly (m/s) at '+str(round(z[iz2], 2))+'m')
        ax.set_yticklabels([''])
        ax.set_xticklabels([''])

        ax = fig.add_subplot(gs[1,1], adjustable = 'box', aspect = 1.)
        CN2 = ax.contourf(X, Y, my_n[it,iz2,:,:] - my_n[it,iz2,:,:].mean(), cmap = 'bwr', levels = np.linspace(-2., 2., 9), extend = 'both')
        axins = inset_axes(ax, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax.transAxes, borderpad = 0.)
        plt.colorbar(CN2, cax = axins)
        ax.contour(X, Y, my_lsm, levels = [1e-08], colors = ['k'], lw = [10])
        ax.set_title('Cross-flow Anomaly (m/s) at '+str(round(z[iz2], 2))+'m')
        ax.set_xlabel('x (km)')
        ax.set_yticklabels([''])

        ax = fig.add_subplot(gs[2,0], adjustable = 'box', aspect = 1.)
        W1 = ax.contourf(X, Y, my_w[it,iz3,:,:], cmap = 'bwr', levels = np.linspace(-2., 2., 9), extend = 'both')
        axins = inset_axes(ax, width = "5%", height = "100%", loc = 6, bbox_to_anchor = (1.05, 0., 1, 1), bbox_transform = ax.transAxes, borderpad = 0.)
        plt.colorbar(W1, cax = axins)
        ax.contour(X, Y, my_lsm, levels = [1e-08], colors = ['k'], lw = [10])
        ax.set_title('Vertical Velocity (m/s) at '+str(round(z[iz3], 2))+'m')
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')

        ax = fig.add_subplot(gs[2,1])
        ax.plot(t, H, 'r', lw = 2, label = 'Sensible HF')
        ax.plot(t, E, 'b', lw = 1, label = 'Latent HF')
        ax.set_xlim([0, 1440])
        ax.set_xticks(np.arange(0, 1440.1, 180))
        ax.plot([times[it], times[it]], [0., 500.], 'k--', lw = 2)
        ax.set_xlabel('Time (mins)')
        ax.set_ylabel('Land-Surface Heat Flux (W m$^{-2}$)')
        lgd=plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.05), prop={'size': 12})

        fig.suptitle('Flow Relative Anomalies at T+'+str(int(times[it]))+'mins', fontsize = 20)
        fig.subplots_adjust(left = 0.0625, right = 0.87, bottom = 0.1, top = 0.9, wspace = 0.31, hspace = 0.31)
        plt.savefig('FlowRelativeWindAnoms_'+"{0:04d}".format(int(times[it]))+'mins.png', dpi = 100, bbox_extra_artists=(lgd,))
        plt.close('all')

    u.close()
    v.close()
    bouy.close()
    mr.close()
send_email(message = 'Finished wind_buoyancy_analysis.py', subject = 'xb899100@ARCHER updates', attachments = ['theta_0720mins.png', 'mv_0720mins.png', 'w_0720mins.png','FlowRelativeWindAnoms_0720mins.png'], isAttach = True)











