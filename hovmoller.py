"""
We are going to do a hovmoller diagram of different variables at different levels
along the length of the cloud trail
"""
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from netCDF4 import Dataset
from datetime import datetime as dt
from scipy import interpolate, integrate
from analysis_tools import bilinear_interpolation, find_h, zi, get_cs_coords
import os

# Calculate the points along the cross section.
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

x_cs, y_cs = get_cs_coords(x_c, y_c, wind_dir_0, X, Y, h = 100.)

hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]
# use dictionaries to hold the stash codes for the variables
var_keys = {'theta' : u'STASH_m01s00i004',
            'qv'    : u'STASH_m01s00i010',
            'w'     : u'STASH_m01s00i150'}
var_cmaps = {'theta' : 'hot',
             'qv'    : 'BrBG',
             'w'     : 'bwr'}
var_factor = {'theta' : 1.,
              'qv'    : 1000.,
              'w'     : 1.}
var_units = {'theta' : '(K)',
             'qv'    : r'(g kg$^{-1}$)',
             'w'     : r'(m s^{-1}$)'}
# use dictionary to hold the hovmollers
hovmollers = {}
for key in var_keys.keys():
    hovmollers[key] = np.zeros((18*len(hours), len(x_cs)))

# convert island radius and along wind distance from m into km
R_i /= 1000.
R = -np.sign(x_cs - x_c)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)/1000.

# create a times array in units of hours
times = np.arange(10., 1440.1, 10.)/60.

# get an array of the heights
bouy_nc = Dataset('../bouy_00.nc', 'r')
z = bouy_nc.variables['thlev_zsea_theta'][:]*1.
bouy_nc.close()

# get the indexes for target heights
izs = [np.where(np.abs(z - i) == np.min(np.abs(z - i)))[0][0] for i in [90, 180, 360, 710]] # these heights are from Matthews et al (2007)

for iz in izs:
    print 'Working on height ' + str(int(z[iz])) + ' m'
    for hour in hours:
        # read in the data
        print 'Reading hour ' + hour
        bouy_nc = Dataset('../bouy_'+hour+'.nc', 'r')
        if hour == '00':
            # skip the first time step to avoid issue with moisture resetting
            it_start = 1
        else:
            it_start = 0
        for key in var_keys.keys():
            hovmollers[key][18*hours.index(hour):(18*hours.index(hour) + 18),:] = bilinear_interpolation(X, Y, bouy_nc.variables[var_keys[key]][it_start:,iz,:,:], x_cs, y_cs, kind = 3)
        bouy_nc.close()

    for key in var_keys.keys():
        print 'Making Hovmoller for ' + key + ' at height ' + str(int(z[iz])) + ' m'
        cbar_min = round(var_factor[key]*(np.nanmean(hovmollers[key]) - np.max(2*np.std(hovmollers[key]), 0.5)), 1)
        cbar_max = round(var_factor[key]*(np.nanmean(hovmollers[key]) + np.max(2*np.std(hovmollers[key]), 0.5)), 1)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        hov = ax.contourf(R, times, hovmollers[key]*var_factor[key], cmap = var_cmaps[key], levels = np.linspace(cbar_min, cbar_max, 21), extend = 'both')
        plt.colorbar(hov, ax = ax, label = key + ' ' + var_units[key])
        ax.plot([-R_i, -R_i], [0, 24], 'k--')
        ax.plot([R_i,R_i], [0, 24], 'k--')
        ax.set_yticks(np.arange(0., 24.1, 3.))
        ax.set_ylabel('Time (hrs)')
        ax.set_xlabel('Distance Downwind (km)')
        ax.set_title(r'Hovmoller of ' + key + ' at ' + "{0:04d}".format(int(z[iz])) + 'm')
        if not os.path.isdir('../Hovmoller_' + key):
            os.mkdir('../Hovmoller_' + key)
        plt.savefig('../Hovmoller_' + key + '/Hovmoller_' + key + '_' + "{0:04d}".format(int(z[iz])) + 'm.png', dpi = 100)
        plt.show()


