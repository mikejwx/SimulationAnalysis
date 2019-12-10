import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate, integrate

restart = 'y'
while restart == 'y':
    ### Read in the tempinc data ###
    path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'
    filename_temp = 'tempinc_09.nc'
    filename_q = 'qinc_09.nc'
    filenames = [filename_temp, filename_q]

    choice = input('Which do we want, tempinc (0) or qinc (1)?\n')
    filename_chosen = filenames[choice]
    suffix = str(181 + choice)
    inc_nc = Dataset(path + filename_chosen, 'r')

    z = inc_nc.variables[u'thlev_zsea_theta'][:]*1.
    inc_times = inc_nc.variables[u'min30_0'][:]*1.

    ### Define keys for the variables of interest ###
    # These are the STASH items I requested, but they might not all exist
    advection     = u'STASH_m01s12i'+suffix
    bdylayr_lscld = u'STASH_m01s09i'+suffix
    bdylayr       = u'STASH_m01s03i'+suffix
    diffusion     = u'STASH_m01s13i'+suffix
    idealised     = u'STASH_m01s53i'+suffix
    lsrain        = u'STASH_m01s04i'+suffix
    QTbalcld      = u'STASH_m01s15i'+suffix
    total         = u'STASH_m01s30i'+suffix

    increments = {'advection'     : advection,
                  'bdylayr_lscld' : bdylayr_lscld,
                  'bdylayr'       : bdylayr,
                  'diffusion'     : diffusion,
                  'idealised'     : idealised,
                  'lsrain'        : lsrain,
                  'QTbalcld'      : QTbalcld,
                  'total'         : total}

    ### Read in our landsea mask for our reference ###
    landseamask = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
    lsm = landseamask.variables['lsm'][0,0,:,:]*1.
    y = landseamask.variables['latitude'][:]/1000.
    x = landseamask.variables['longitude'][:]/1000.

    # create 2D coordinate mesh
    X, Y = np.meshgrid(x, y)
    landseamask.close()

    ### Also read in our liquid water path for reference, where to choose the box ###
    filename = 'lwp_00.nc'
    lwp_nc = Dataset(path + filename, 'r')
    lwp_key = u'STASH_m01s30i405'
    lwp_data = lwp_nc.variables[lwp_key][:]*1.
    lwp_times = lwp_nc.variables['min5_0'][:]*1.
    its = [it for it in xrange(len(lwp_times)) if lwp_times[it] in inc_times]
    lwp_mean = np.nanmean(np.where((lwp_data[min(its):(max(its)+1),:,:] > 1e-5), 1., 0.), axis = 0)
    lwp_nc.close()

    ### Start routine to display a scene and allow user to input the corners of polygon ###
    end_points = []
    n_endpoints = input('How many vertices in your polygon?\n')
    fig = plt.figure(tight_layout = True)
    ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
    while len(end_points) < n_endpoints:
        ax.cla()
        ax.contourf(X, Y, lwp_mean, levels = np.arange(0., 0.30001, 0.05), cmap = 'Greys_r')
        ax.contour(X, Y, lsm, levels = [0.0, 1e-5], colors = ['r'])
        [ax.plot(point[0][0], point[0][1], 'bo', markersize = 10) for point in end_points]
        ax.plot([point[0][0] for point in end_points], [point[0][1] for point in end_points], color = 'b', lw = 2)
        ax.set_xlabel('x (km)')
        ax.set_ylabel('y (km)')
        ax.set_title('Cloud frequency for T+' + str(int(np.nanmean(lwp_times[min(its):(max(its)+1)]))) + 'mins')
        plt.pause(1e-4)
        end_points.append(plt.ginput(1))
        print end_points[-1]
    
    # Connect the polygon
    end_points = end_points + [end_points[0]]
    xn = [point[0][0] for point in end_points]
    yn = [point[0][1] for point in end_points]

    # Plot the final polygon
    ax.cla()
    ax.contourf(X, Y, lwp_mean, levels = np.arange(0., 0.30001, 0.05), cmap = 'Greys_r')
    ax.contour(X, Y, lsm, levels = [0.0, 1e-5], colors = ['r'])
    ax.plot(xn, yn, color = 'b', lw = 2)
    ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    ax.set_title('Cloud mask for T+' + str(int(np.nanmean(lwp_times[min(its):(max(its)+1)]))) + 'mins')
    plt.pause(1)
    plt.close('all')

    ### Get a mask of the area enclosed by that polygon ###
    ## get the indices enclosed by that rectangle
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
        ax.contourf(X, Y, mask*lwp_mean, levels = np.arange(0., 0.30001, 0.05), cmap = 'Greys_r')
        ax.contour(X, Y, mask, levels = [0, 1e-16], colors = ['r'])
        ax.plot(xn, yn, color = 'b', ls = '--')
        plt.pause(0.1)
    plt.close('all')


    # Plot the time and space mean of the temperature increments in this region
    fig = plt.figure(tight_layout = True)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylim([0, 4000])
    ax.set_ylabel('Height (m)')
    if choice:
        ax.set_xlabel('Increment (kg/kg/day)')
    else:
        ax.set_xlabel('Increment (K/day)')
    for key in increments.keys():
        if increments[key] in inc_nc.variables.keys():
            print 'Working on ' + key
            mean_var_profile = np.nanmean(mask*np.nanmean(inc_nc.variables[increments[key]][:]*1., axis = 0), axis = (1, 2))*86400/3.
            ax.plot(mean_var_profile, z, label = key)
            ax.legend()
            plt.pause(1)
    ax.set_title('All Increments')
    plt.show()
    inc_nc.close()
    restart = raw_input('Would you like to make another selection? (y/n)\n')

