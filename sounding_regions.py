import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate, integrate
from SkewT_archer import *
from analysis_tools import bilinear_interpolation
from mpl_toolkits.axes_grid.inset_locator import inset_axes

restart = 'y'
while restart == 'y':
    # Let the user tellus what they want to look at:
    experiment_input = input('Which experiment would you like to inspect?\nControl, H=250W/m2 (option = 0)\n  exp01, H=125W/m2 (option = 1)\n  exp02, H=375W/m2 (option = 2)\noption = ')
    experiments = ['Control', 'H125', 'H375']
    path = '/nerc/n02/n02/xb899100/CloudTrail/'+experiments[experiment_input]+'/'
    
    l_skewT = input('Plot (option = 0) or SkewT (option = 1)?\noption = ')
    
    hours = ["{0:02d}".format(hour) for hour in xrange(0, 24, 3)]
    prompt = 'Which time chunk would you like to inspect?\n'
    for hour in hours:
        prompt += hour + ' to ' + "{0:02d}".format(int(hour) + 3) + '  (option = ' + str(hours.index(hour)) + ')\n'
    prompt += 'option = '
    hour_idx = input(prompt)
    
    # Grab all the data required to plot the skewT or the plot
    print 'Opening the netCDF needed...'
    bouy_nc   = Dataset(path + 'bouy_' + hours[hour_idx] + '.nc', 'r')
    wind_nc   = Dataset(path + 'wind_' + hours[hour_idx] + '.nc', 'r')
    fluxes_nc = Dataset(path + 'fluxes_' + hours[hour_idx] + '.nc', 'r')
    
    # define keys for our variables
    temp_key = u'STASH_m01s16i004'
    q_key    = u'STASH_m01s00i010'
    u_key    = u'STASH_m01s00i002'
    v_key    = u'STASH_m01s00i003'
    p_key    = u'STASH_m01s00i408'
    
    time_key = [key for key in bouy_nc.variables.keys() if 'min' in key][0]
    
    # Get the height and time dimensions
    z     = bouy_nc.variables['thlev_zsea_theta'][:]*1.
    times = bouy_nc.variables[time_key][:]*1.
    
    # Read in our landsea mask to help us pick a location
    landseamask = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
    lsm = landseamask.variables['lsm'][0,0,:,:]*1.
    y = landseamask.variables['latitude'][:]/1000.
    x = landseamask.variables['longitude'][:]/1000.

    # create 2D coordinate mesh
    X, Y = np.meshgrid(x, y)
    landseamask.close()

    # Also read in liquid water path to help us pick a location
    lwp_nc       = Dataset(path + 'lwp_00.nc', 'r')
    lwp_key      = u'STASH_m01s30i405'
    lwp_data     = lwp_nc.variables[lwp_key][:]*1.
    lwp_time_key = [key for key in lwp_nc.variables.keys() if 'min' in key][0]
    lwp_times    = lwp_nc.variables[lwp_time_key][:]*1.
    
    # Get the indexes for the lwp times that have corresponding (thermo)dynamic data
    its = [it for it in xrange(len(lwp_times)) if lwp_times[it] in times]
    
    lwp_nc.close()
    
    # Need more user input, do we want a point or an area?
    l_area = input('Plot data at a point (option = 0), or the average of an area (option = 1)?\noption = ')
    if l_area:
        ########################################################################
        # Grab the mean of an area
        ########################################################################
        # Start routine to display a scene and allow user to input the corners of polygon
        end_points = []
        point = [(np.nan, np.nan)]
        n_endpoints = 0
        it = 0
        while n_endpoints < 3:
            n_endpoints = input('How many vertices in your polygon?\n')
            if n_endpoints < 3:
                print 'A polygon must have more than 3 vertices...'
        
        fig = plt.figure(tight_layout = True)
        ax0 = fig.add_subplot(1, 3, 1, adjustable = 'box', aspect = 1)
        ax  = fig.add_subplot(1, 3, 2, adjustable = 'box', aspect = 1)
        ax1 = fig.add_subplot(1, 3, 3, adjustable = 'box', aspect = 1)
        ax0.set_xlim(-11, -9)
        ax0.set_ylim(-11, -9)
        ax0.text(-10, -10, 'Previous')
        ax0.set_yticklabels([])
        ax0.set_xticklabels([])
        ax1.set_xlim(-101, -99)
        ax1.set_ylim(-101, -99)
        ax1.text(-100, -100, 'Next')
        ax1.set_yticklabels([])
        ax1.set_xticklabels([])
        while len(end_points) < n_endpoints:
            ax.cla()
            ax.contourf(X, Y, lwp_data[its[it],:,:]*1000., levels = np.arange(0., 1000.1, 50.), cmap = 'Greys_r')
            ax.contour(X, Y, lsm, levels = [0.0, 1e-5], colors = ['r'])
            [ax.plot(point[0][0], point[0][1], 'bo', markersize = 10) for point in end_points]
            ax.plot([point[0][0] for point in end_points], [point[0][1] for point in end_points], color = 'b', lw = 2)
            ax.set_xlabel('x (km)')
            ax.set_ylabel('y (km)')
            ax.set_title('LWP for T+' + "{0:04d}".format(int(lwp_times[its[it]])) + 'mins')
            plt.pause(1e-4)
            point = plt.ginput(1)
            if point[0][0] < 0:
                if abs(point[0][0] + 10) < abs(point[0][0] + 100):
                    it -= 1
                else:
                    it += 1
                it = it%len(its)
            else:
                end_points.append(point)
                print end_points[-1]
            
        # Connect the polygon
        end_points = end_points + [end_points[0]]
        xn = [point[0][0] for point in end_points]
        yn = [point[0][1] for point in end_points]
        
        # Plot the final polygon
        ax.cla()
        ax.contourf(X, Y, lwp_data[its[it],:,:]*1000., levels = np.arange(0., 1000.1, 50.), cmap = 'Greys_r')
        ax.contour(X, Y, lsm, levels = [0.0, 1e-5], colors = ['r'])
        ax.plot(xn, yn, color = 'b', lw = 2)
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        ax.set_title('LWP for T+' + "{0:04d}".format(int(lwp_times[its[it]])) + 'mins')
        plt.pause(1)
        plt.close('all')
        
        # Get a mask of the area enclosed by that polygon
        # Get the mid-point of the polygon as the mean of the xn and yn coordinates
        x_mid = np.mean(xn)
        y_mid = np.mean(yn)

        mask = np.ones_like(X)
        fig = plt.figure(tight_layout = True)
        ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
        for side in xrange(n_endpoints):
            print 'Working on side ' + str(side+1)
            slope_y = (yn[side+1] - yn[side])/(xn[side+1] - xn[side])
            intercept_y = yn[side]
            
            slope_x = (xn[side+1] - xn[side])/(yn[side+1] - yn[side])
            intercept_x = xn[side]
            
            x_side = 0.5*(xn[side+1] + xn[side])
            y_side = 0.5*(yn[side+1] + yn[side])
            ### Assuming points are going clockwise ###
            if abs(slope_y) >= abs(slope_x):
                # If slope_y is large, then approximate it as a vertical line in the y direction and we want points to the left or right of it
                if x_mid > x_side:
                    # If mid-point of the polygon (x_mid) is greater than the midpoint of the side (x_side) then we want points that are greater than the coordinates of the side
                    mask *= np.where((X > (slope_x*(Y-yn[side]) + intercept_x)), 1., np.nan)
                else:
                    mask *= np.where((X < (slope_x*(Y-yn[side]) + intercept_x)), 1., np.nan)
            else:
                if y_mid >= y_side:
                    mask *= np.where((Y > (slope_y*(X-xn[side]) + intercept_y)), 1., np.nan)
                else:
                    mask *= np.where((Y < (slope_y*(X-xn[side]) + intercept_y)), 1., np.nan)
            
            ax.cla()
            ax.contourf(X, Y, mask*lwp_data[its[it],:,:]*1000., levels = np.arange(0., 1000.1, 50.), cmap = 'Greys_r')
            ax.contour(X, Y, lsm, levels = [0, 1e-16], colors = ['r'])
            ax.plot(xn, yn, color = 'b', ls = '--')
            plt.pause(0.25)
        plt.close('all')
        
        print 'Calculating area mean temperatures...'
        temperature_prof = np.nanmean(mask*bouy_nc.variables[temp_key][it,:,:,:], axis = (1, 2))
        temperature_25th = np.nanpercentile(mask*bouy_nc.variables[temp_key][it,:,:,:], 10, axis = (1, 2))
        temperature_75th = np.nanpercentile(mask*bouy_nc.variables[temp_key][it,:,:,:], 90, axis = (1, 2))
        print 'Calculating area mean dewpoints...'
        dewpoint_data    = getDew(bouy_nc.variables[q_key][it,:,:,:], fluxes_nc.variables[p_key][it,:,:,:], q_units = 'kg/kg', p_units = 'Pa')
        dewpoint_prof    = np.nanmean(mask*dewpoint_data, axis = (1, 2))
        dewpoint_25th    = np.nanpercentile(mask*dewpoint_data, 10, axis = (1, 2))
        dewpoint_75th    = np.nanpercentile(mask*dewpoint_data, 90, axis = (1, 2))
        print 'Calculating area mean winds...'
        u_prof = np.nanmean(mask*wind_nc.variables[u_key][it,:,:,:], axis = (1, 2))
        v_prof = np.nanmean(mask*wind_nc.variables[v_key][it,:,:,:], axis = (1, 2))
        p_prof = np.nanmean(mask*fluxes_nc.variables[p_key][it,:,:,:], axis = (1, 2))
    else:
        ########################################################################
        # Grab data at a point
        ########################################################################
        print 'Choose the point...'
        # At a point rather than an average of an area
        point = [(np.nan, np.nan)]
        it = 0
        fig = plt.figure(tight_layout = True)
        ax0 = fig.add_subplot(1, 3, 1, adjustable = 'box', aspect = 1)
        ax  = fig.add_subplot(1, 3, 2, adjustable = 'box', aspect = 1)
        ax1 = fig.add_subplot(1, 3, 3, adjustable = 'box', aspect = 1)
        ax0.set_xlim(-11, -9)
        ax0.set_ylim(-11, -9)
        ax0.text(-10, -10, 'Previous')
        ax0.set_yticklabels([])
        ax0.set_xticklabels([])
        ax1.set_xlim(-101, -99)
        ax1.set_ylim(-101, -99)
        ax1.text(-100, -100, 'Next')
        ax1.set_yticklabels([])
        ax1.set_xticklabels([])
        while point == [(np.nan, np.nan)]:
            ax.cla()
            ax.contourf(X, Y, lwp_data[its[it],:,:]*1000., levels = np.arange(0., 1000.1, 50.), cmap = 'Greys_r', extend = 'max')
            ax.contour(X, Y, lsm, levels = [0,1e-16], colors = ['r'])
            ax.plot(point[0][0], point[0][1], 'r*', markersize = 20)
            ax.set_title('T+' + "{0:04d}".format(int(lwp_times[its[it]])) + ' mins')
            plt.pause(1e-4)
            point = plt.ginput(1)
            if point[0][0] < 0:
                if abs(point[0][0] + 10) < abs(point[0][0] + 100):
                    it -= 1
                else:
                    it += 1
                point = [(np.nan, np.nan)]
                it = it%len(its)
            print point
        
        ax.cla()
        ax.contourf(X, Y, lwp_data[its[it],:,:]*1000., levels = np.arange(0., 1000.1, 50.), cmap = 'Greys_r', extend = 'max')
        ax.contour(X, Y, lsm, levels = [0,1e-16], colors = ['r'])
        ax.plot(point[0][0], point[0][1], 'r*', markersize = 20)
        ax.set_title('T+' + "{0:04d}".format(int(lwp_times[its[it]])) + ' mins')
        plt.pause(1)
        plt.close('all')
        
        print 'Calculating temperature at this point...'
        temperature_prof = np.array([obs[0] for obs in bilinear_interpolation(X, Y, bouy_nc.variables[temp_key][it,:,:,:], [point[0][0]], [point[0][1]], kind = 2)])
        print 'Calculating dewpoint at this point...'
        dewpoint_prof    = np.array([obs[0] for obs in bilinear_interpolation(X, Y, getDew(bouy_nc.variables[q_key][it,:,:,:], fluxes_nc.variables[p_key][it,:,:,:], q_units = 'kg/kg', p_units = 'Pa'), [point[0][0]], [point[0][1]], kind = 2)])
        print 'Calculating winds at this point...'
        u_prof = np.array([obs[0] for obs in bilinear_interpolation(X, Y, wind_nc.variables[u_key][it,:,:,:], [point[0][0]], [point[0][1]], kind = 2)])
        v_prof = np.array([obs[0] for obs in bilinear_interpolation(X, Y, wind_nc.variables[v_key][it,:,:,:], [point[0][0]], [point[0][1]], kind = 2)])
        p_prof = np.array([obs[0] for obs in bilinear_interpolation(X, Y, fluxes_nc.variables[p_key][it,:,:,:], [point[0][0]], [point[0][1]], kind = 2)])
    
    if l_skewT:
        ########################################################################
        # Plot a skewT
        ########################################################################
        print 'Plotting data on a SkewT...'
        fig = plt.figure(figsize = (5, 7.5))
        ax = fig.add_subplot(1, 2, 1, adjustable = 'box', aspect = 1)
        plotSkewT(temperature_prof[:-1]-273.15, dewpoint_prof[:-1]-273.15, p_prof[:-1]/100., u_prof[:-1], v_prof[:-1])
        axins = inset_axes(ax, width = "100%", height = "100%", loc = 9, bbox_to_anchor = (0.0, 2, 2, 2*319./1160.), bbox_transform = ax.transAxes, borderpad = 0.)
        axins.contourf(X, Y, lwp_data[its[it],:,:]*1000., levels = np.arange(0, 1000.1, 50.), cmap = 'Greys_r', extend = 'max')
        axins.contour(X, Y, lsm, levels = [0,1e-16], colors = ['r'])
        if l_area:
            axins.plot(xn, yn, 'b--', lw = 2)
        else:
            axins.plot(point[0][0], point[0][1], 'r*', markersize = 20)
        axins.set_xlabel('x (km)')
        axins.set_ylabel('y (km)')
        axins.set_title(experiments[experiment_input] + ' simulation at T+' "{0:04d}".format(int(lwp_times[its[it]])) + 'mins')
        fig.subplots_adjust(top = 1.5)
        plt.show()
        
    else:
        print 'Plotting data...'
        # Plot the time and space mean of the temperature increments in this region
        fig = plt.figure(tight_layout = True)
        ax  = fig.add_subplot(2, 1, 1, adjustable = 'box', aspect = 1)
        ax.contourf(X, Y, lwp_data[its[it],:,:]*1000., levels = np.arange(0., 1000.1, 50.), cmap = 'Greys_r', extend = 'max')
        ax.contour(X, Y, lsm, levels = [0,1e-16], colors = ['r'])
        if l_area:
            ax.plot(xn, yn, 'b--', lw = 2)
        else:
            ax.plot(point[0][0], point[0][1], 'r*', markersize = 20)
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        # Plot the temperature and dewpoint
        ax1 = fig.add_subplot(2, 2, 3)
        ax1.plot(temperature_prof, z, 'r', lw = 2, label = 'Temperature')
        ax1.plot(dewpoint_prof, z, 'b', label = 'Dewpoint')
        if l_area:
            # Plot a interquartile range
            ax1.fill_betweenx(z, temperature_25th, temperature_75th, alpha = 0.5, color = 'r', edgecolor = '')
            ax1.fill_betweenx(z, dewpoint_25th, dewpoint_75th, alpha = 0.5, color = 'b', edgecolor = '')
        ax1.set_xlim([250, 310])
        ax1.set_ylim([0, 4000])
        ax1.set_ylabel('Height (m)')
        ax1.set_xlabel('Temperature (K)')
        
        # Plot u and v winds
        ax2 = fig.add_subplot(2, 2, 4, sharey = ax1)
        ax2.plot(u_prof, z, 'k', lw = 2, label = 'u-wind')
        ax2.plot(v_prof, z, 'k--', lw = 2, label = 'v-wind')
        ax2.set_ylim([0, 4000])
        ax2.set_xlabel('Wind (m/s)')
        
        ax.set_title(experiments[experiment_input] + ' simulation at T+' "{0:04d}".format(int(lwp_times[its[it]])) + 'mins')
        plt.show()
        bouy_nc.close()
        wind_nc.close()
        fluxes_nc.close()
        
    restart = raw_input('Would you like to make another selection? (y or n)\n')

