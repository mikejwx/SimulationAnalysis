"""
1. Pick a representative wind direction.
2. Find the cross sectional slice through the domain intersecting the middle of 
   the island.
3. Interpolate to that cross sectional slice path at 100 m resolution?
4. Plot the anomalies with respect to the horizontal mean of quantities of 
   interest along our cross section.
"""
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from netCDF4 import Dataset
from datetime import datetime as dt
from scipy import interpolate, integrate
from SkewT_archer import *
from analysis_tools import bilinear_interpolation, zi, find_h

def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 80, fill = '='):
    """
    Call in a loop to create terminal progress bar
    @params:
      iteration - Required : current iteration (Int)
      total     - Required : total iterations (Int)
      prefix    - Optional : prefix string (Str)
      suffix    - Optional : suffix string (Str)
      decimals  - Optional : positive number of decimals in percent complete (Int)
      length    - Optional : character length of bar (Int)
      fill      - Optional : bar fill character (Str)
    """
    from __future__ import print_function
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r['+dt.now().strftime("%H:%M:%S")+'] %s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

testing = True

f = open("output.txt", "w+")

### TO DO: Put into a loop to view changes in time ###

#================================= Section 1 ==================================#
#                                                                              #
# Pick a representative wind direction                                         #
#                                                                              #
#==============================================================================#

#print('Starting section 1')
f = open("output.txt", "a")
f.write('[' + dt.now().strftime("%H:%M:%S") + '] Starting section 1\n')
f.close()
# Assume that the appropriate representative wind direction is equal to the 
# initial conditions for the wind

# Initial conditions taken directly from namelist
u_0 = np.array([-6.09,-7.02,-7.53,-7.89,-8.15,-8.36,-8.53,-8.68,-8.79,-8.89,
                -8.97,-9.02,-9.07,-9.1,-9.12,-9.13,-9.14,-9.15,-9.16,-9.16,
                -9.17,-9.19,-9.21,-9.24,-9.29,-9.35,-9.43,-9.54,-9.68,-9.84,
                -9.98,-10.09,-10.14,-10.15,-10.13,-10.1,-10.08,-10.06,-10.05,
                -10.04,-10.01,-10.0,-10.0,-10.0,-9.66,-0.09,0.0,0.0])
v_0 = np.array([-1.2,-1.38,-1.48,-1.55,-1.6,-1.64,-1.67,-1.7,-1.72,-1.74,
                -1.76,-1.78,-1.8,-1.81,-1.83,-1.84,-1.85,-1.86,-1.86,-1.85,
                -1.84,-1.82,-1.79,-1.76,-1.73,-1.69,-1.65,-1.61,-1.55,-1.45,
                -1.3,-1.08,-0.83,-0.62,-0.46,-0.34,-0.26,-0.2,-0.15,-0.11,-0.03,
                -0.01,0.0,0.0,0.0,0.0,0.0,0.0])
z_0 = np.array([1.0000004,3.6666676,7.666668,13.000004,19.666672,27.666672,
                37.000008,47.66668,59.66668,73.00004,87.66668,103.66668,
                121.00004,139.66672,159.66672,181.00004,203.66672,227.66672,
                337.00008,367.66676,399.66676,433.0,467.6668,503.6668,541.0,
                579.6668,619.6668,661.0,703.6668,747.6668,793.0004,839.6668,
                887.6668,937.0004,987.6668,1039.6668,1093.0004,1147.6672,
                1203.6668,1261.0004,1503.6672,1567.6668,1633.0004,8955.796,
                9205.932,14947.828,15802.464,15802.464])

# Calculate the mean wind direction in the boundary layer (lowest ~ 850 m)
z_1 = np.arange(0., 850.1, 1.)
u_1 = interpolate.interp1d(y = u_0, x = z_0, fill_value = 'extrapolate')(z_1)
v_1 = interpolate.interp1d(y = v_0, x = z_0, fill_value = 'extrapolate')(z_1)

U_0 = integrate.trapz(y = u_1, x = z_1)/850.
V_0 = integrate.trapz(y = v_1, x = z_1)/850.

wind_speed_0 = np.sqrt(U_0**2 + V_0**2)
# Wind speed should be ~ 9.5 m/s
wind_dir_0 = 360.*np.arctan(U_0/V_0)/(2.*np.pi)
# Wind direction should be ~ 80 deg.
f = open("output.txt", "a")
f.write('[' + dt.now().strftime("%H:%M:%S") + '] Wind direction is: ' + str(round(wind_dir_0)) + " deg\n")
f.close()
#================================= Section 2 ==================================#
#                                                                              #
# Find the cross sectional slice through the domain                            #
#                                                                              #
#==============================================================================#

#print('Starting section 2')
f = open("output.txt", "a")
f.write('[' + dt.now().strftime("%H:%M:%S") + '] Starting section 2\n')
f.close()
# Read in a land-sea mask to find the island and center the cross section on the
# island mid-point
f = open("output.txt", "a")
f.write('[' + dt.now().strftime("%H:%M:%S") + '] Reading land-sea mask...\n')
f.close()
landseamask = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')

lsm = landseamask.variables['lsm'][0,0,:,:]*1.
y = landseamask.variables['latitude'][:]*1.
x = landseamask.variables['longitude'][:]*1.

# create 2D coordinate mesh
X, Y = np.meshgrid(x, y)

landseamask.close()

# We know from our domain definition where the centre of the island should be
R_i = 1000.0*(50.0/np.pi)**0.5 # island radius
x_c = 100000.0 + R_i
y_c = 4*R_i

# From this center point, draw line in either direction parallel to wind_dir_0
# get the coordinates of this line at 100 m resolution.
x_cs, y_cs = get_cs_coords(x_c, y_c, wind_dir_0, h = 100.0)

### Testing ###
if testing:
    f = open("output.txt", "a")
    f.write('[' + dt.now().strftime("%H:%M:%S") + '] Testing cross section...\n')
    f.close()
    lwp_nc   = Dataset('../lwp_00.nc', 'r')
    lwp_data = lwp_nc.variables['STASH_m01s30i405'][144,:,:]*1.
    lwp_nc.close()
    
    fig = plt.figure()
    my_x = 12.
    fig.set_size_inches(my_x, my_x*32./116.)
    plt.contourf(X/1000., Y/1000., lwp_data, levels = [0, 1e-8, 1], colors = ['k', 'w'])
    plt.contourf(X/1000., Y/1000., lsm, levels = [0, 1e-16, 1], colors = ['w', 'b'], alpha = 0.5)
    plt.plot(x_cs/1000., y_cs/1000., 'r', lw = 3)
    plt.ylabel('y (km)')
    plt.xlabel('x (km)')
    plt.text(x_cs[0], y_cs[0], 'A', fontsize = 20)
    plt.text(x_cs[-1], y_cs[-1], 'B', fontsize = 20)
    plt.title('Cross Section Path')
    plt.tight_layout()
    plt.savefig('Cross_Section_location.png', dpi = 100)
    plt.close('all')

#================================= Section 3 ==================================#
#                                                                              #
# Interpolate our data to the cross section coordinates                        #
#                                                                              #
#==============================================================================#

#print('Starting section 3')
f = open("output.txt", "a")
f.write('[' + dt.now().strftime("%H:%M:%S") + '] Starting section 3\n')
f.close()

# time slices

# Read some data
#print('Reading data')
f = open("output.txt", "a")
f.write('[' + dt.now().strftime("%H:%M:%S") + '] Reading a netCDF for theta, q, and w...\n')
f.close()

def get_cs(it):
    # only need data to the height just above 3000 m
    iz3000 = np.where(np.abs(z - 3000) == np.min(np.abs(z - 3000)))[0][0]+1
    # leave the rest of the levels = 0

    f = open("output.txt", "a")
    f.write('[' + dt.now().strftime("%H:%M:%S") + '] Interpolating to the cross section coordinates...\n')
    f.close()

    START = dt.now()
    print('[' + dt.now().strftime("%H:%M:%S") + '] Starting interpolation...\n')
    start = dt.now()
    theta_cs = bilinear_interpolation(X, Y, bouy_nc.variables[theta_key][it,:iz3000,:,:]*1., x_cs, y_cs, kind = 3)
    end = dt.now()
    print('Finished interpolating Theta, elapsed time: ' + str(end - start))

    start = dt.now()
    q_cs     = bilinear_interpolation(X, Y, bouy_nc.variables[q_key][it,:iz3000,:,:]*1., x_cs, y_cs, kind = 3)
    end = dt.now()
    print('Finished interpolating q, elapsed time: ' + str(end - start))

    start = dt.now()
    w_cs     = bilinear_interpolation(X, Y, bouy_nc.variables[w_key][it,:iz3000,:,:]*1., x_cs, y_cs, kind = 3)
    end = dt.now()
    print('Finished interpolating w, elapsed time: ' + str(end - start))
    
    END = dt.now()
    f = open('output.txt', 'a')
    f.write('[' + dt.now().strftime("%H:%M:%S") + '] Completed interpolation in ' + str(END - START) + "...\n")
    f.close()

    #================================= Section 4 ==================================#
    #                                                                              #
    # Plot the cross section and anomalies                                         #
    #                                                                              #
    #==============================================================================#

    #print('Starting section 4')

    # Translate the cross section coordinates into the distance downwind of the island
    R = - np.sign(x_c - x_cs)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)

    # Calculate the anomaly from the horizontal mean
    theta_cs_anom = np.zeros_like(theta_cs)
    q_cs_anom = np.zeros_like(q_cs)

    theta_hm = np.nanmean(bouy_nc.variables[theta_key][it,:iz3000,:,:], axis = (1, 2))
    q_hm     = np.nanmean(bouy_nc.variables[q_key][it,:iz3000,:,:], axis = (1, 2))
    for i in xrange(theta_cs.shape[1]):
        theta_cs_anom[:,i] = theta_cs[:,i] - theta_hm
        q_cs_anom[:,i]     = q_cs[:,i] - q_hm

    # Calculate the boundary layer height
    theta_v_cs = theta_cs*(1.0 + 0.608*q_cs)
    zi_cs = zi(theta_v_cs, z[:iz3000])

    # Plotting params
    text_pos_x = -R_i/1000. + R_i*0.1/1000.
    text_pos_y = 100.
    fig_width = 10 # inches

    fig = plt.figure()
    fig.set_size_inches(fig_width, fig_width/1.5)
    ax = plt.subplot(311)
    plt.contourf(R/1000., z[:iz3000], theta_cs_anom, cmap = 'bwr', levels = np.arange(-1.0, 1.01, 0.2), extend = 'both')
    plt.colorbar(label = r'$\theta$ Anomaly (K)')
    plt.title('Potential Temperature Anomaly (K) Cross Section')
    plt.ylim([0, 3000])
    r_i = np.array([-R_i, R_i])/1000.
    plt.plot(r_i, [0, 0], 'k', lw = 5)
    plt.text(text_pos_x, text_pos_y, 'Island')
    plt.plot(R/1000., zi_cs, 'k', lw = 2)
    plt.ylabel('Height (m)')
    plt.setp(ax.get_xticklabels(), visible=False)

    ax1 = plt.subplot(312, sharex = ax)
    plt.contourf(R/1000., z[:iz3000], q_cs_anom*1000., cmap = 'BrBG', levels = np.arange(-2.0, 2.01, 0.2), extend = 'both')
    plt.colorbar(label = r'q Anomaly (g kg$^{-1}$)')
    plt.title('Specific Humidity Anomaly (g kg$^{-1}$) Cross Section')
    plt.ylim([0, 3000])
    plt.plot(r_i, [0, 0], 'k', lw = 5)
    plt.text(text_pos_x, text_pos_y, 'Island')
    plt.plot(R/1000., zi_cs, 'k', lw = 2)
    plt.ylabel('Height (m)')
    plt.setp(ax1.get_xticklabels(), visible=False)

    plt.subplot(313, sharex = ax)
    plt.contourf(R/1000., z[:iz3000], w_cs, cmap = 'bwr', levels = np.arange(-2.0, 2.01, 0.2), extend = 'both')
    plt.colorbar(label = r'w (m s$^{-1}$)')
    plt.title('Vertical Velocity (m s$^{-1}$) Cross Section')
    plt.ylim([0, 3000])
    plt.plot(r_i, [0, 0], 'k', lw = 5)
    plt.text(text_pos_x, text_pos_y, 'Island')
    plt.plot(R/1000., zi_cs, 'k', lw = 2)
    plt.ylabel('Height (m)')
    plt.xlabel('Distance from Island mid-point (km)')
    plt.suptitle('Along-wind Cross Section, at T+'+str(int(times[it])) + 'mins', fontsize = 15)
    plt.tight_layout()
    fig.subplots_adjust(top=0.88)
    plt.savefig('AlongWind_smoothed_T'+"{0:04d}".format(int(times[it]))+'.png', dpi = 100)
    #plt.show()
    
    # For the skewT plots
    # Get a sounding every 5 km downwind
    h  = 100.0 # copied from above because for some reason, we cannot see this here.
    i_cs = np.where(R%5000. < h)[0]
    temp_cs_soundings = bilinear_interpolation(X, Y, bouy_nc.variables[temp_key][it,:,:,:]*1., x_cs[i_cs], y_cs[i_cs], kind = 0) # nearest neighbor
    pres_cs_soundings = bilinear_interpolation(X, Y, fluxes_nc.variables[pres_key][it,:,:,:]*1., x_cs[i_cs], y_cs[i_cs], kind = 0)
    q_cs_soundings    = bilinear_interpolation(X, Y, bouy_nc.variables[q_key][it,:,:,:]*1., x_cs[i_cs], y_cs[i_cs], kind = 0)
    dew_cs_soundings  = getDew(q_cs_soundings, pres_cs_soundings/100.)
    
    #Plot the soundings as well
    for i in xrange(len(i_cs)):
        plotSkewT(temp_cs_soundings[:,i]-273.15, dew_cs_soundings[:,i]-273.15, pres_cs_soundings[:,i]/100., my_title = 'T+' + str(int(times[it])) + 'mins, ' + str(int(R[i_cs[i]]/1000.)) + 'km downwind')
        plt.savefig('AlongWind_sounding_T'+"{0:04d}".format(int(times[it]))+'_R'+"{0:02d}".format(int(R[i_cs[i]]/1000.))+'.png', dpi = 100)
        if it == 0:
            plt.show()
        plt.close('all')
    #================================= Section 5 ==================================#
    #                                                                              #
    # Calculate cross sections perpendicular to the flow and plot them             #
    #                                                                              #
    #==============================================================================#

    # Choose distance from island mid-point to do the cross section
    # e.g. 20 km
    r_targets = [10000.0, 20000.0, 40000.0, 60000.0]
    for r_target in r_targets:
        # Find where the radial distance from the island is nearest 25 km, this will do for now
        R_pos = np.where((R > 0), R, np.nan)
        iR = np.where(np.abs(R_pos - r_target) == np.nanmin(np.abs(R_pos - r_target)))[0][0]

        # From this coordinate, extend a line perpendicular to the previous cross section
        # get the coordinates of this line at h m resolution.
        h  = 100.0
        dx = h*np.sin(np.pi*(wind_dir_0+90.)/180.0)
        dy = h*np.cos(np.pi*(wind_dir_0+90.)/180.0)

        x_cs2 = [x_cs[iR]]
        y_cs2 = [y_cs[iR]]

        # Populate my list of x and y coordinates along the cross section
        while (x_cs2[-1] < np.max(x))*(y_cs2[-1] < np.max(y)):
            x_cs2.append(x_cs2[-1] + dx)
            y_cs2.append(y_cs2[-1] + dy)

        x_cs2 = x_cs2[::-1]
        y_cs2 = y_cs2[::-1]

        while (x_cs2[-1] > np.min(x))*(y_cs2[-1] > np.min(y)):
            x_cs2.append(x_cs2[-1] - dx)
            y_cs2.append(y_cs2[-1] - dy)

        # check that all the coordinates along the cross section are within the domain
        i = 0
        while i < len(x_cs2):
            if (x_cs2[i] > np.max(x)) or (x_cs2[i] < np.min(x)) or (y_cs2[i] > np.max(y)) or (y_cs2[i] < np.min(y)):
                del x_cs2[i]
                del y_cs2[i]
            else:
                i += 1

        # store to arrays
        x_cs2 = np.array(x_cs2)
        y_cs2 = np.array(y_cs2)

        ### Testing ###
        if testing:
            fig = plt.figure()
            my_x = 12.
            fig.set_size_inches(my_x, my_x*32./116.)
            plt.contourf(X/1000., Y/1000., lwp_data, levels = [0, 1e-8, 1], colors = ['k', 'w'])
            plt.contourf(X/1000., Y/1000., lsm, levels = [0, 1e-16, 1], colors = ['w', 'b'], alpha = 0.5)
            plt.plot(x_cs/1000., y_cs/1000., 'r', lw = 3)
            plt.plot(x_cs2/1000., y_cs2/1000., 'b', lw = 3)
            plt.ylabel('y (km)')
            plt.xlabel('x (km)')
            plt.text(x_cs[0], y_cs[0], 'A', fontsize = 20)
            plt.text(x_cs[-1], y_cs[-1], 'B', fontsize = 20)
            plt.text(x_cs2[0], y_cs2[0], 'C', fontsize = 20)
            plt.text(x_cs2[-1], y_cs2[-1], 'D', fontsize = 20)
            plt.title('Cross Section Path')
            plt.tight_layout()
            plt.savefig('Cross_Section_location2.png', dpi = 100)
            #plt.show()
            plt.close('all')

        # Initialise arrays for our cross section
        theta_cs2 = np.zeros((len(z), len(x_cs2)))
        q_cs2 = np.zeros_like(theta_cs2)
        w_cs2 = np.zeros_like(theta_cs2)

        f = open("output.txt", "a")
        f.write('[' + dt.now().strftime("%H:%M:%S") + '] Starting interpolation...\n')
        f.close()
        START = dt.now()

        # Interpolate to cs2
        start = dt.now()
        theta_cs2 = bilinear_interpolation(X, Y, bouy_nc.variables[theta_key][it,:iz3000,:,:], x_cs2, y_cs2, kind = 3)
        end = dt.now()
        print('Finished interpolating Theta, elapsed time: ' + str(end - start))

        start = dt.now()
        q_cs2     = bilinear_interpolation(X, Y, bouy_nc.variables[q_key][it,:iz3000,:,:], x_cs2, y_cs2, kind = 3)
        end = dt.now()
        print('Finished interpolating q, elapsed time: ' + str(end - start))

        start = dt.now()
        w_cs2     = bilinear_interpolation(X, Y, bouy_nc.variables[w_key][it,:iz3000,:,:], x_cs2, y_cs2, kind = 3)
        end = dt.now()
        print('Finished interpolating w, elapsed time: ' + str(end - start))
        END = dt.now()

        f = open('output.txt', 'a')
        f.write('[' + dt.now().strftime("%H:%M:%S") + '] Completed interpolation in ' + str(END - START) + "...\n")
        f.close()

        # Translate the cross section coordinates into the distance away from the island downwind mid-point
        R2 = - np.sign(x_cs[iR] - x_cs2)*np.sqrt((x_cs2 - x_cs[iR])**2 + (y_cs2 - y_cs[iR])**2)

        # Calculate the anomaly from the horizontal mean
        theta_cs2_anom = np.zeros_like(theta_cs2)
        q_cs2_anom = np.zeros_like(q_cs2)

        for i in xrange(theta_cs2.shape[1]):
            theta_cs2_anom[:,i] = theta_cs2[:,i] - theta_hm
            q_cs2_anom[:,i]     = q_cs2[:,i] - q_hm

        # Calculate the boundary layer height
        theta_v_cs2 = theta_cs2*(1.0 + 0.608*q_cs2)
        zi_cs2 = zi(theta_v_cs2, z[:iz3000])

        # Plotting params
        text_pos_x = -R_i/1000. + R_i*0.1/1000.
        text_pos_y = 100.
        fig = plt.figure()

        ax = plt.subplot(311)
        plt.contourf(R2/1000., z[:iz3000], theta_cs2_anom, cmap = 'bwr', levels = np.arange(-1.0, 1.01, 0.2), extend = 'both')
        plt.colorbar()
        plt.title('Potential Temperature Anomaly (K) Cross Section')
        plt.ylim([0, 3000])
        r_i = np.array([-R_i, R_i])/1000.
        plt.plot(r_i, [0, 0], 'k', lw = 10)
        plt.text(text_pos_x, text_pos_y, 'Island')
        plt.plot(R2/1000., zi_cs2, 'k', lw = 2)
        plt.ylabel('Height (m)')
        plt.setp(ax.get_xticklabels(), visible=False)

        ax1 = plt.subplot(312, sharex = ax)
        plt.contourf(R2/1000., z[:iz3000], q_cs2_anom*1000., cmap = 'BrBG', levels = np.arange(-2.0, 2.01, 0.2), extend = 'both')
        plt.colorbar()
        plt.title('Specific Humidity Anomaly (g kg$^{-1}$) Cross Section')
        plt.ylim([0, 3000])
        plt.plot(r_i, [0, 0], 'k', lw = 10)
        plt.text(text_pos_x, text_pos_y, 'Island')
        plt.plot(R2/1000., zi_cs2, 'k', lw = 2)
        plt.ylabel('Height (m)')
        plt.setp(ax1.get_xticklabels(), visible=False)

        plt.subplot(313, sharex = ax)
        plt.contourf(R2/1000., z[:iz3000], w_cs2, cmap = 'bwr', levels = np.arange(-2.0, 2.01, 0.2), extend = 'both')
        plt.colorbar()
        plt.title('Vertical Velocity (m s$^{-1}$) Cross Section')
        plt.ylim([0, 3000])
        plt.plot(r_i, [0, 0], 'k', lw = 10)
        plt.text(text_pos_x, text_pos_y, 'Island')
        plt.plot(R2/1000., zi_cs2, 'k', lw = 2)
        plt.ylabel('Height (m)')
        plt.xlabel('Distance from Island mid-point (km)')
        plt.suptitle('Cross-wind Cross Section ('+str(int(r_target/1000.))+'km downwind), at T+'+str(int(times[it])) + 'mins', fontsize = 15)
        plt.tight_layout()
        fig.subplots_adjust(top=0.88)
        plt.savefig('CrossWind_smoothed_'+str(int(r_target/1000.))+'km_T'+"{0:04d}".format(int(times[it]))+'.png', dpi = 100)
        #plt.show()

hours = ["{0:02d}".format(h) for h in xrange(0, 24,3)]
for hour in hours:
    # e.g. the v-component winds
    bouy_nc   = Dataset('../bouy_' + hour + '.nc', 'r')
    fluxes_nc = Dataset('../fluxes_' + hour + '.nc', 'r')
    theta_key = u'STASH_m01s00i004' #bouy.nc
    temp_key  = u'STASH_m01s16i004' #bouy.nc
    pres_key  = u'STASH_m01s00i408' #fluxes.nc
    q_key     = u'STASH_m01s00i010' #bouy.nc
    w_key     = u'STASH_m01s00i150' #bouy.nc
    z         = bouy_nc.variables[u'thlev_zsea_theta'][:]*1.
    times     = bouy_nc.variables[u'min10_0'][:]*1.

    for I in xrange(len(times)):
        print('Hour = ' + hour + ', I = ' + str(I))
        get_cs(I)

    bouy_nc.close()
    fluxes_nc.close()
