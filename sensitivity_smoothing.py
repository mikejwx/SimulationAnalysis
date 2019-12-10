"""
Want to test the sensitivity of cross sections and hovmollers to the radius of averaging.
"""

import numpy as np
import matplotlib.pyplot as plt
from analysis_tools import bilinear_interpolation, get_cs_coords
from netCDF4 import Dataset
from scipy import interpolate, integrate
from datetime import datetime as dt

# Calculate the points along the cross section.
# Assume that the appropriate representative wind direction is equal to the 
# initial conditions for the wind

# Initial conditions taken directly from namelist
print '[' + dt.now().strftime('%H:%M:%S') + '] Defining the mean boundary layer winds'
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

print '[' + dt.now().strftime('%H:%M:%S') + '] Reading the land-sea mask'
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

print '[' + dt.now().strftime('%H:%M:%S') + '] Generating cross section coordinates'
x_cs, y_cs = get_cs_coords(x_c, y_c, wind_dir_0, X, Y, h = 100.)

# define the hours from which to draw data
hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]

# Consider specific humidity
q_key  = u'STASH_m01s00i010'
zi_key = u'new boundary layer depth'
z_key  = u'thlev_zsea_theta'

dict_keys = ['0100', '0250', '0500', '1000', '2000', '4000']
hovmoller_dict    = {}
crosssection_dict = {}
zi_dict           = {}
selected_pts_dict = {}

for key in dict_keys:
    hovmoller_dict[key]    = np.zeros((18*len(hours), len(x_cs)))
    selected_pts_dict[key] = np.zeros_like(lsm)

for hour in hours:
    print '[' + dt.now().strftime('%H:%M:%S') + '] Reading the data for hour ' + hour
    # Read the data
    bouy_nc = Dataset('../bouy_' + hour + '.nc', 'r')
    if hour == '00':
        z = bouy_nc.variables[z_key][:]
        # skip the first time step to avoid issue with moisture resetting
        it_start = 1
        iz = np.where(np.abs(z - 360.) == np.min(np.abs(z - 360.)))[0][0] # these is a height from Matthews et al (2007)
    else:
        it_start = 0
    
    # Calculate the hovmoller and cross section
    for key in dict_keys:
        print '[' + dt.now().strftime('%H:%M:%S') + '] Working on ' + key + 'm smoothing Radius'
        hovmoller_dict[key][18*hours.index(hour):(18*hours.index(hour) + 18),:] = bilinear_interpolation(X, Y, bouy_nc.variables[q_key][it_start:,iz,:,:], x_cs, y_cs, kind = 3, d = int(key))
        
        if hour == '09':
            print '[' + dt.now().strftime('%H:%M:%S') + '] Doing the selected points, boundary layer depth, and cross section smoothing'
            # Offline, get the selected points
            for i in xrange(len(x_cs)):
                r = np.sqrt((X - x_cs[i])**2 + (Y - y_cs[i])**2)
                iy, ix = np.where((r < int(key)))
                selected_pts_dict[key][iy,ix] = 1.0
            
            # Additionally read in the boundary layer height
            zi_nc = Dataset('../zi_' + hour + '.nc', 'r')
            zi_data = zi_nc.variables[zi_key][:]
            zi_nc.close()
            
            # choose to consider hour == '09' for the cross section test
            crosssection_dict[key] = bilinear_interpolation(X, Y, bouy_nc.variables[q_key][-1,:, :, :], x_cs, y_cs, kind = 3, d = int(key))
            zi_dict[key]           = bilinear_interpolation(X, Y, zi_data, x_cs, y_cs, kind = 3, d = int(key))[-1,:]
    bouy_nc.close()

# convert island radius and along wind distance from m into km
R_i /= 1000.
R = -np.sign(x_cs - x_c)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)/1000.

# create a times array in units of hours
times = np.arange(10., 1440.1, 10.)/60.

fig = plt.figure(tight_layout = True, figsize = (16, 9))
for key in dict_keys:
    ax = fig.add_subplot(3, 2, dict_keys.index(key)+1, adjustable = 'box', aspect = 'equal')
    my_plt = ax.contourf(R, times, hovmoller_dict[key]*1000., cmap = 'BrBG', levels = np.linspace(14.5, 18.5, 11), extend = 'both')
    fig.colorbar(my_plt, ax = ax, label = 'Specific Humidity (g kg$^{-1}$)')
    ax.plot([R_i, R_i], [0, 24], 'k--', lw = 2)
    ax.plot([-R_i, -R_i], [0, 24], 'k--', lw = 2)
    #ax.plot([0, 102.6], [6, 9], 'b--', lw = 2)
    ax.set_xlabel('Distance Downwind of Island Center (km)')
    ax.set_ylabel('Time (hrs)')
    ax.set_yticks(xrange(0, 25, 3))
    ax.set_xlim([np.min(R), np.max(R)])
    ax.set_title('Hovmoller w/ smoothing radius=' + str(int(key)) + ' m')

fig.suptitle('Specific Humidity at ' + str(round(z[iz], 2)) + ' m', fontsize = 20)
fig.subplots_adjust(top = 0.75)
plt.savefig('../Sensitivity_to_smoothing_Hovmoller.png', dpi = 100)
plt.close('all')

fig = plt.figure(tight_layout = True, figsize = (16, 9))
for key in dict_keys:
    ax = fig.add_subplot(3, 2, dict_keys.index(key)+1, adjustable = 'box')
    my_plt = ax.contourf(R, z/1000., crosssection_dict[key]*1000., cmap = 'BrBG', levels = np.linspace(0., 20., 11), extend = 'max')
    ax.plot(R, zi_dict[key]/1000., 'k', lw = 2)
    ax.plot([-R_i, R_i], [0,0], 'k', lw = 5)
    fig.colorbar(my_plt, ax = ax, label = 'Specific Humidity (g kg$^{-1}$)')
    ax.set_xlabel('Distance Downwind of Island Center (km)')
    ax.set_ylabel('Height (km)')
    ax.set_ylim([0, 3])
    ax.set_title('Cross Section w/ smoothing radius=' + str(int(key)) + ' m')

fig.suptitle('Specific Humidity at T+' + str(round(times[71], 2)) + ' hr', fontsize = 20)
fig.subplots_adjust(top = 0.5)
plt.savefig('../Sensitivity_to_smoothing_cross_section.png', dpi = 100)
plt.close('all')

lwp_key = u'STASH_m01s30i405'
lwp_nc = Dataset('../lwp_00.nc', 'r')
lwp = lwp_nc.variables[lwp_key][144,:,:]
lwp_nc.close()

fig = plt.figure(figsize = (16, 9))
for key in dict_keys:
    ax = fig.add_subplot(3, 2, dict_keys.index(key)+1, adjustable = 'box', aspect = 'equal')
    my_plt = ax.contourf(X/1000., Y/1000., lwp, colors = ['w', 'k'], levels = [0.0, 1e-08, 1.0])
    ax.contour(X/1000., Y/1000., lsm, levels = [1e-16], colors = ['k'], linewidths = [2])
    ax.contourf(X/1000., Y/1000., selected_pts_dict[key], colors = ['w', 'r'], levels = [0.0, 1e-16, 1.0], alpha = 0.5)
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_title('Selected points w/ smoothing radius=' + str(int(key)) + ' m')

fig.suptitle('Cloud Mask at T+' + str(round(times[71], 2)) + ' hr', fontsize = 20)
#fig.subplots_adjust(top = 0.85)
plt.savefig('../Sensitivity_to_smoothing_smoothed_area.png', dpi = 100)
plt.close('all')

