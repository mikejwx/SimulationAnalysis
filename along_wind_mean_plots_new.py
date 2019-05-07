"""
What does the structure of a chunk of the cloud trail region look like?
Look at the along flow mean cloud liquid water mixing ratio, vertical velocity,
 potential temperature and in plane wind arrows/streamlines?

See Kirshbaum and Fairman (2015) figure 9(b) for inspiration
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import bilinear_interpolation, get_cs_coords, transform_winds, send_email
import os
from scipy import interpolate, integrate
from datetime import datetime as dt
from multiprocessing import Pool

print '[' + dt.now().strftime('%H:%M:%S') + '] Importing variable keys'
from STASH_keys import w_key, theta_key, theta_anom_key, mcl_key, lwp_key, np_key

hours = ['09']
wanted_times = [600., 660., 720.]

IDs = ['Control', 'H125', 'H375']
paths = ['/nerc/n02/n02/xb899100/CloudTrail/'+ID+'/' for ID in IDs]

for exp_idx in xrange(len(paths)):
    path      = paths[exp_idx]
    ID        = IDs[exp_idx]
    lwp_nc    = Dataset(path + 'lwp_00.nc', 'r')
    lwp_times = lwp_nc.variables[[key for key in lwp_nc.variables.keys() if 'min' in key][0]][:]*1.
    my_cmap   = mpl.cm.get_cmap('Greys')
    
    for hour in hours:
        # read in some data 
        w_swath_nc     = Dataset(path + 'w_swath_' + hour + '.nc', 'r')
        theta_swath_nc = Dataset(path + 'theta_swath_' + hour + '.nc', 'r')
        mcl_swath_nc   = Dataset(path + 'mcl_swath_' + hour + '.nc', 'r')
        n_swath_nc     = Dataset(path + 'n_swath_' + hour + '.nc', 'r')
        
        # determine which wanted_times are in this netCDF
        time_key = [key for key in w_swath_nc.variables.keys() if 'min' in key][0]
        times = w_swath_nc.variables[time_key][:]
        its = [np.where(wanted_time == times)[0][0] for wanted_time in wanted_times if wanted_time in times]
        
        # read the variables
        print '[' + dt.now().strftime('%H:%M:%S') + '] Reading data from the netCDFs'
        z     = w_swath_nc.variables['thlev_zsea_theta'][:]
        i_max = np.where(np.abs(z - 3000.) == np.min(np.abs(z - 3000.)))[0][0] # Find the index of the level nearest 3000 m
        z     = z[:i_max] # Clip to just the lowest 3000 m
        
        # We know from our domain definition where the centre of the island should be
        R_i = 1000.0*(50.0/np.pi)**0.5 # island radius
        x_c = 100000.0 + R_i
        y_c = 4*R_i
        
        # Define faux cartesian coordinates for horizontal slices of the whole domain
        X, Y = np.meshgrid(np.arange(0., 116000., 100.), np.arange(0., 31900., 100.))
        
        # Read in the interpolated land sea mask
        lsm = theta_swath_nc.variables['lsm'][:]*1.
        
        # Define the distance away from the island centre
        x_prime = w_swath_nc.variables['x_prime'][:]
        y_prime = w_swath_nc.variables['y_prime'][:]
        x_R = x_prime - x_c
        y_R = y_prime - y_c
        R   = - np.sign(x_R)*np.sqrt(x_R**2 + y_R**2)
        
        # Define the coordinate system in terms of along the wind and across the wind
        i_along, i_across = np.where(R == np.min(np.abs(R)))
        R_along  = R[i_along[0],:]
        R_across = R[:,i_across[0]]
        
        # Select a chunk from the 'along the wind' direction
        r_target0 = R_i # i.e. 40 km downwind is + 40000.
        r_target1 = R_i + 80000.
        
        idx0 = np.where(np.abs(R_along - r_target0) == np.min(np.abs(R_along - r_target0)))[0][0]
        idx1 = np.where(np.abs(R_along - r_target1) == np.min(np.abs(R_along - r_target1)))[0][0] + 1
        if idx1 == x_prime.shape[1]:
            idx1 = -1
        
        for it in its:
            theta   = theta_swath_nc.variables[theta_key][it,:i_max,:,:]
            theta_p = theta_swath_nc.variables[theta_anom_key][it,:i_max,:,:]
            w       = w_swath_nc.variables[w_key][it,:i_max,:,:]
            mcl     = mcl_swath_nc.variables[mcl_key][it,:i_max,:,:]
            n_prime = n_swath_nc.variables[np_key][it,:i_max,:,:]
            
            print '[' + dt.now().strftime('%H:%M:%S') + '] Collapsing variable chunks in the along-wind direction'
            theta_mean   = np.nanmean(theta, axis = 2)
            theta_p_mean = np.nanmean(theta_p, axis = 2)
            w_mean       = np.nanmean(w, axis = 2)
            n_mean       = np.nanmean(n_prime, axis = 2)
            mcl_mean     = np.nanmean(mcl, axis = 2)
            
            # Find the time index in lwp that matched the time index in the other variables
            lwp_it = np.where(np.min(np.abs(lwp_times - times[it])) == np.abs(lwp_times - times[it]))[0][0]
            lwp_data = lwp_nc.variables[lwp_key][lwp_it,:,:]
            ### Plot the liquid water path and the bounding box for our chunk ###
            print '[' + dt.now().strftime('%H:%M:%S') + '] Plot where the chunk is'
            fig = plt.figure(tight_layout = True, figsize=(15, 5))
            ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
            LWP_plt = ax.contourf(X/1000., Y/1000., lwp_data[:]*1000., colors = my_cmap(np.arange(0, 9.)/8.), levels = [0., 10., 20., 50., 100., 200., 300., 400., 500., 600.], extend = 'max')
            fig.colorbar(LWP_plt, ax = ax, label = r'LWP (g m$^{-2}$)')
            # Plot the island coastline
            ax.contour(x_prime/1000., y_prime/1000., lsm, colors = ['r'], linewidths = 2)
            x_corners = np.array([x_prime[0,idx0], x_prime[0,idx1], x_prime[-1, idx1], x_prime[-1,idx0], x_prime[0,idx0]])/1000.
            y_corners = np.array([y_prime[0,idx0], y_prime[0,idx1], y_prime[-1, idx1], y_prime[-1,idx0], y_prime[0,idx0]])/1000.
            # Plot the chunk box
            ax.plot(x_corners, y_corners, 'r')
            ax.set_xlabel('x (km)')
            ax.set_ylabel('y (km)')
            ax.set_title(ID + ': ' + str(int((r_target0 - R_i)/1000.)) + ' to ' + str(int((r_target1 - R_i)/1000.)) + ' km downwind of island, T+' + "{0:04d} mins".format(int(times[it])))
            plt.savefig('../' + ID + "_T+{0:04d}".format(int(times[it])) + '_' + str(int((r_target0 - R_i)/1000.)) + 'to' + str(int((r_target1 - R_i)/1000.)) + 'km.png', dpi = 150)
            plt.close('all')
            
            ### Plot the along-wind chunk mean as in KF15 ###
            print '[' + dt.now().strftime('%H:%M:%S') + '] Plot along-wind chunk mean following KF15'
            fig = plt.figure(tight_layout = True, figsize = (20,6))
            ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
            # potential temperature
            tc = ax.contour(R_across/1000., z/1000., theta_mean, colors = ['b'], levels = range(298, 310), linewidths = [2])
            ax.clabel(tc, inline=1, fmt= '%1d')
            # vertical velocity
            ax.contourf(R_across/1000., z/1000., w_mean, cmap = 'bwr', vmin = -2., vmax = 2., levels = [x for x in np.arange(-1.0, 2.01, 0.2) if round(x,1) != 0], extend = 'max')
            ax.contour(R_across/1000., z/1000., w_mean, colors = ['k'], levels = [x for x in np.arange(-1.0, 2.01, 0.2) if round(x,1) != 0])
            # cloud liquid
            ax.contourf(R_across/1000., z/1000., mcl_mean, levels = [0.01e-3, 1e10], hatches = ['////'], colors = 'none')
            # winds
            q = ax.quiver(R_across/1000., z[::3]/1000., n_mean[::3,:], w_mean[::3,:], scale = 50., width = 0.001)
            ax.quiverkey(q, X = 0.55, Y = -0.05, U = 2, label = '2 m/s', labelpos = 'E')
            ax.set_xlabel("y$^{\prime}$ (km)")
            ax.set_ylabel('Height (km)')
            ax.set_title(ID + ': ' + str(int((r_target0 - R_i)/1000.)) + ' to ' + str(int((r_target1 - R_i)/1000.)) + ' km downwind of island, T+' + "{0:04d} mins".format(int(times[it])))
            plt.savefig('../along_wind_mean_' + ID + "_T+{0:04d}".format(int(times[it])) + '_' + str(int((r_target0 - R_i)/1000.)) + 'to' + str(int((r_target1 - R_i)/1000.)) + 'km.png', dpi = 150)
            plt.close('all')
            
            ### Plot the along wind chunk theta anomaly to inspect the warm plume ###
            print '[' + dt.now().strftime('%H:%M:%S') + '] Plot along-wind chunk mean of potential temperature anomalies'
            fig = plt.figure(tight_layout = True,figsize=(20,6))
            ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
            anom_plt = ax.contourf(R_across/1000., z/1000., theta_p_mean, cmap = 'bwr', levels = [x for x in np.arange(-1., 1.1, 0.1) if round(x,1) != 0], extend = 'both')
            fig.colorbar(anom_plt, ax = ax, ticks = np.arange(-1, 1.1, 0.1), label = '$\\theta^{\prime}$ (K)')
            q = ax.quiver(R_across/1000., z[::3]/1000., n_mean[::3,:], w_mean[::3,:], scale = 50., width = 0.001)
            ax.quiverkey(q, X = 0.55, Y = -0.05, U = 2, label = '2 m/s', labelpos = 'E')
            ax.contourf(R_across/1000., z/1000., mcl_mean, levels = [0.01e-3, 1e10], hatches = ['////'], colors = 'none')
            ax.set_ylabel('Height (km)')
            ax.set_xlabel('y$^{\prime}$ (km)')
            ax.set_title(ID + ': ' + '$\\theta^{\prime}$ ' + str(int((r_target0 - R_i)/1000.)) + ' to ' + str(int((r_target1 - R_i)/1000.)) + ' km downwind of island, T+' + "{0:04d} mins".format(int(times[it])))
            plt.savefig('../along_wind_mean_warm_plume' + ID + "_T+{0:04d}".format(int(times[it])) + '_' + str(int((r_target0 - R_i)/1000.)) + 'to' + str(int((r_target1 - R_i)/1000.)) + 'km.png', dpi = 150)
            plt.close('all')
    lwp_nc.close()
    print '[' + dt.now().strftime('%H:%M:%S') + '] Complete.'

