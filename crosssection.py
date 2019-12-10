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
#plt.switch_backend('agg')
from netCDF4 import Dataset
from datetime import datetime as dt
from scipy import interpolate, integrate
from SkewT_archer import *
from analysis_tools import bilinear_interpolation, get_cs_coords


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
z_1 = np.arange(0., 500.1, 1.)
u_1 = interpolate.interp1d(y = u_0, x = z_0, fill_value = 'extrapolate')(z_1)
v_1 = interpolate.interp1d(y = v_0, x = z_0, fill_value = 'extrapolate')(z_1)

U_0 = integrate.trapz(y = u_1, x = z_1)/500.
V_0 = integrate.trapz(y = v_1, x = z_1)/500.

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
x_cs, y_cs = get_cs_coords(x_c, y_c, wind_dir_0, X, Y, h = 100.0)

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
    plt.savefig('../Cross_Section_location.png', dpi = 100)
    plt.show()

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
    iz3000 = np.where(np.abs(z - 5000) == np.min(np.abs(z - 5000)))[0][0]+1
    # leave the rest of the levels = 0
    my_kind = 2
    d = 2000.
    print('[' + dt.now().strftime("%H:%M:%S") + '] Interpolating to the cross section coordinates...')
    
    START = dt.now()
    print('[' + dt.now().strftime("%H:%M:%S") + '] Starting interpolation...')
    start = dt.now()
    theta_cs = bilinear_interpolation(X, Y, bouy_nc.variables[theta_key][it,:iz3000,:,:]*1., x_cs, y_cs, kind = my_kind, d = d)
    end = dt.now()
    print('Finished interpolating Theta, elapsed time: ' + str(end - start))
    
    start = dt.now()
    q_cs     = bilinear_interpolation(X, Y, bouy_nc.variables[q_key][it,:iz3000,:,:]*1., x_cs, y_cs, kind = my_kind, d = d)
    end = dt.now()
    print('Finished interpolating q, elapsed time: ' + str(end - start))
    
    start = dt.now()
    w_cs     = bilinear_interpolation(X, Y, bouy_nc.variables[w_key][it,:iz3000,:,:]*1., x_cs, y_cs, kind = my_kind, d = d)
    end = dt.now()
    print('Finished interpolating w, elapsed time: ' + str(end - start))
    
    start = dt.now()
    zi_cs    = bilinear_interpolation(X, Y, zi_nc.variables[zi_key][:], x_cs, y_cs, kind = my_kind, d = d)
    end = dt.now()
    print('Finished interpolating zi, elapsed time: ' + str(end - start))
    
    start = dt.now()
    lcl_cs    = bilinear_interpolation(X, Y, zi_nc.variables[lcl_key][:], x_cs, y_cs, kind = my_kind, d = d)
    end = dt.now()
    print('Finished interpolating lcl, elapsed time: ' + str(end - start))
    
    END = dt.now()
    print('[' + dt.now().strftime("%H:%M:%S") + '] Completed interpolation in ' + str(END - START) + "...")
    
    #================================= Section 4 ==================================#
    #                                                                              #
    # Plot the cross section and anomalies                                         #
    #                                                                              #
    #==============================================================================#
    
    #print('Starting section 4')
    
    # Translate the cross section coordinates into the distance downwind of the island
    R = np.sign(x_c - x_cs)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)
    
    # Calculate the anomaly from the horizontal mean
    theta_cs_anom = np.zeros_like(theta_cs)
    q_cs_anom = np.zeros_like(q_cs)
    
    theta_hm = np.nanmean(bouy_nc.variables[theta_key][it,:iz3000,:,:], axis = (1, 2))
    q_hm     = np.nanmean(bouy_nc.variables[q_key][it,:iz3000,:,:], axis = (1, 2))
    for i in xrange(theta_cs.shape[1]):
        theta_cs_anom[:,i] = theta_cs[:,i] - theta_hm
        q_cs_anom[:,i]     = q_cs[:,i] - q_hm
    
    # Plotting params
    text_pos_x = -R_i/1000. + R_i*0.1/1000.
    text_pos_y = 100.
    fig_width = 10 # inches
    
    fig = plt.figure()
    fig.set_size_inches(fig_width, fig_width/1.5)
    ax = plt.subplot(311)
    plt.contourf(R/1000., z[:iz3000]/1000., theta_cs_anom, cmap = 'bwr', levels = [l for l in np.linspace(-1., 1., 11) if l != 0.], extend = 'both')
    plt.colorbar(label = r'$\theta$ Anomaly (K)')
    plt.title('Potential Temperature Anomaly (K) Cross Section')
    plt.ylim([0, 5])
    r_i = np.array([-R_i, R_i])/1000.
    plt.plot(r_i, [0, 0], 'k', lw = 5)
    plt.text(text_pos_x, text_pos_y, 'Island')
    plt.plot(R/1000., zi_cs[it,:]/1000., 'k', lw = 2)
    plt.plot(R/1000., lcl_cs[it,:]/1000., 'grey', lw = 2)
    plt.ylabel('Height (km)')
    plt.setp(ax.get_xticklabels(), visible=False)
    
    ax1 = plt.subplot(312, sharex = ax)
    plt.contourf(R/1000., z[:iz3000]/1000., q_cs_anom*1000., cmap = 'BrBG', levels = [l for l in np.linspace(-2., 2., 11) if l != 0.], extend = 'both')
    plt.colorbar(label = r'q Anomaly (g kg$^{-1}$)')
    plt.title('Specific Humidity Anomaly (g kg$^{-1}$) Cross Section')
    plt.ylim([0, 5])
    plt.plot(r_i, [0, 0], 'k', lw = 5)
    plt.text(text_pos_x, text_pos_y, 'Island')
    plt.plot(R/1000., zi_cs[it,:]/1000., 'k', lw = 2)
    plt.plot(R/1000., lcl_cs[it,:]/1000., 'grey', lw = 2)
    plt.ylabel('Height (km)')
    plt.setp(ax1.get_xticklabels(), visible=False)
    
    plt.subplot(313, sharex = ax)
    plt.contourf(R/1000., z[:iz3000]/1000., w_cs, cmap = 'bwr', levels = [l for l in np.linspace(-2.5, 2.5, 11) if l != 0.], extend = 'both')
    plt.colorbar(label = r'w (m s$^{-1}$)')
    plt.title('Vertical Velocity (m s$^{-1}$) Cross Section')
    plt.ylim([0, 5])
    plt.plot(r_i, [0, 0], 'k', lw = 5)
    plt.text(text_pos_x, text_pos_y, 'Island')
    plt.plot(R/1000., zi_cs[it,:]/1000., 'k', lw = 2)
    plt.plot(R/1000., lcl_cs[it,:]/1000., 'grey', lw = 2)
    plt.ylabel('Height (km)')
    plt.xlabel('Distance from Island mid-point (km)')
    plt.suptitle('Along-wind Cross Section, at T+'+str(int(times[it])) + 'mins', fontsize = 15)
    plt.tight_layout()
    fig.subplots_adjust(top=0.88)
    plt.savefig('../AlongWind_interpolated_T'+"{0:04d}".format(int(times[it]))+'.png', dpi = 100)
    plt.close('all')
    #plt.show()
    
    #================================= Section 5 ==================================#
    #                                                                              #
    # Calculate cross sections perpendicular to the flow and plot them             #
    #                                                                              #
    #==============================================================================#

    # Choose distance from island mid-point to do the cross section
    # e.g. 20 km
    r_targets = [10000.0, 20000.0, 40000.0, 60000.0]
    for r_target in r_targets:
        # Find where the distance to the island is nearest to the target distance
        R_pos = np.where((R > 0), R, np.nan)
        iR = np.where(np.abs(R_pos - r_target) == np.nanmin(np.abs(R_pos - r_target)))[0][0]
        
        x_cs2, y_cs2 = get_cs_coords(x_cs[iR], y_cs[iR], wind_dir_0+90., X, Y, h = 100.0)
        
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
            plt.savefig('../Cross_Section_location2.png', dpi = 100)
            plt.show()

        f = open("output.txt", "a")
        f.write('[' + dt.now().strftime("%H:%M:%S") + '] Starting interpolation...\n')
        f.close()
        START = dt.now()

        # Interpolate to cs2
        start = dt.now()
        theta_cs2 = bilinear_interpolation(X, Y, bouy_nc.variables[theta_key][it,:iz3000,:,:], x_cs2, y_cs2, kind = my_kind, d = d)
        end = dt.now()
        print('Finished interpolating Theta, elapsed time: ' + str(end - start))

        start = dt.now()
        q_cs2     = bilinear_interpolation(X, Y, bouy_nc.variables[q_key][it,:iz3000,:,:], x_cs2, y_cs2, kind = my_kind, d = d)
        end = dt.now()
        print('Finished interpolating q, elapsed time: ' + str(end - start))

        start = dt.now()
        w_cs2     = bilinear_interpolation(X, Y, bouy_nc.variables[w_key][it,:iz3000,:,:], x_cs2, y_cs2, kind = my_kind, d = d)
        end = dt.now()
        print('Finished interpolating w, elapsed time: ' + str(end - start))
        END = dt.now()
        
        start = dt.now()
        zi_cs2     = bilinear_interpolation(X, Y, zi_nc.variables[zi_key][:], x_cs2, y_cs2, kind = my_kind, d = d)
        end = dt.now()
        print('Finished interpolating zi, elapsed time: ' + str(end - start))
        END = dt.now()
        
        start = dt.now()
        lcl_cs2     = bilinear_interpolation(X, Y, zi_nc.variables[lcl_key][:], x_cs2, y_cs2, kind = my_kind, d = d)
        end = dt.now()
        print('Finished interpolating lcl, elapsed time: ' + str(end - start))
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

        # Plotting params
        text_pos_x = -R_i/1000. + R_i*0.1/1000.
        text_pos_y = 100.
        fig = plt.figure()
        
        ax = plt.subplot(311)
        plt.contourf(R2/1000., z[:iz3000]/1000., theta_cs2_anom, cmap = 'bwr', levels = [l for l in np.linspace(-1., 1., 11) if l != 0.], extend = 'both')
        plt.colorbar()
        plt.title('Potential Temperature Anomaly (K) Cross Section')
        plt.ylim([0, 5])
        r_i = np.array([-R_i, R_i])/1000.
        plt.plot(r_i, [0, 0], 'k', lw = 10)
        plt.text(text_pos_x, text_pos_y, 'Island')
        plt.plot(R2/1000., zi_cs2[it,:]/1000., 'k', lw = 2)
        plt.plot(R2/1000., lcl_cs2[it,:]/1000., 'grey', lw = 2)
        plt.ylabel('Height (m)')
        plt.setp(ax.get_xticklabels(), visible=False)
        
        ax1 = plt.subplot(312, sharex = ax)
        plt.contourf(R2/1000., z[:iz3000]/1000., q_cs2_anom*1000., cmap = 'BrBG', levels = [l for l in np.linspace(-2., 2., 11) if l != 0.], extend = 'both')
        plt.colorbar()
        plt.title('Specific Humidity Anomaly (g kg$^{-1}$) Cross Section')
        plt.ylim([0, 5])
        plt.plot(r_i, [0, 0], 'k', lw = 10)
        plt.text(text_pos_x, text_pos_y, 'Island')
        plt.plot(R2/1000., zi_cs2[it,:]/1000., 'k', lw = 2)
        plt.plot(R2/1000., lcl_cs2[it,:]/1000., 'grey', lw = 2)
        plt.ylabel('Height (m)')
        plt.setp(ax1.get_xticklabels(), visible=False)
        
        plt.subplot(313, sharex = ax)
        plt.contourf(R2/1000., z[:iz3000]/1000., w_cs2, cmap = 'bwr', levels = [l for l in np.linspace(-2.5, 2.5, 11) if l != 0.], extend = 'both')
        plt.colorbar()
        plt.title('Vertical Velocity (m s$^{-1}$) Cross Section')
        plt.ylim([0, 5])
        plt.plot(r_i, [0, 0], 'k', lw = 10)
        plt.text(text_pos_x, text_pos_y, 'Island')
        plt.plot(R2/1000., zi_cs2[it,:]/1000., 'k', lw = 2)
        plt.plot(R2/1000., lcl_cs2[it,:]/1000., 'grey', lw = 2)
        plt.ylabel('Height (m)')
        plt.xlabel('Distance from Island mid-point (km)')
        plt.suptitle('Cross-wind Cross Section ('+str(int(r_target/1000.))+'km downwind), at T+'+str(int(times[it])) + 'mins', fontsize = 15)
        plt.tight_layout()
        fig.subplots_adjust(top=0.88)
        plt.savefig('../CrossWind_smoothed_'+str(int(r_target/1000.))+'km_T'+"{0:04d}".format(int(times[it]))+'.png', dpi = 100)
        plt.close('all')
        #plt.show()

#hours = ["{0:02d}".format(h) for h in xrange(0, 24,3)]
hours = ['06', '09']
for hour in hours:
    # e.g. the v-component winds
    bouy_nc   = Dataset('../bouy_' + hour + '.nc', 'r')
    fluxes_nc = Dataset('../fluxes_' + hour + '.nc', 'r')
    zi_nc     = Dataset('../zi_' + hour + '.nc', 'r')
    theta_key = u'STASH_m01s00i004' #bouy.nc
    temp_key  = u'STASH_m01s16i004' #bouy.nc
    pres_key  = u'STASH_m01s00i408' #fluxes.nc
    q_key     = u'STASH_m01s00i010' #bouy.nc
    w_key     = u'STASH_m01s00i150' #bouy.nc
    zi_key    = u'new boundary layer depth' #zi.nc
    lcl_key   = u'lifting condensation level' #zi.nc
    z         = bouy_nc.variables[u'thlev_zsea_theta'][:]*1.
    times     = bouy_nc.variables[u'min10_0'][:]*1.

    for I in xrange(len(times)):
        print('Hour = ' + hour + ', I = ' + str(I))
        get_cs(I)

    bouy_nc.close()
    fluxes_nc.close()
    zi_nc.close()



