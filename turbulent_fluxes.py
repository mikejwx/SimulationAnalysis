"""
Scripts to analyse turbulent fluxes in our simulations.

See the Matthews et al. (2007) NCAR/COMET module schematic
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import bilinear_interpolation, get_cs_coords
import os
from scipy import interpolate, integrate
from datetime import datetime as dt

# Calculate the points along the cross section.
# Assume that the appropriate representative wind direction is equal to the 
# initial conditions for the wind

print '[' + dt.now().strftime('%H:%M:%S') + '] Generating the cross section'
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

# convert island radius and along wind distance from m into km
R_i /= 1000.
R = -np.sign(x_cs - x_c)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)/1000.

# Make use of our netCDF naming convention
hours = ["{0:02d}".format(h) for h in xrange(0, 24,3)]

print '[' + dt.now().strftime('%H:%M:%S') + '] Defining physical constants and STASH keys'
# Define some constants
Lv = 2.501e06

# Define some keys
shf_key = u'STASH_m01s03i216'
lhf_key = u'STASH_m01s03i222'
u_key   = u'STASH_m01s00i002'
v_key   = u'STASH_m01s00i003'
w_key   = u'STASH_m01s00i150'
zi_key  = u'new boundary layer depth'

for hour in ['09', '12']:#hours:
    print '[' + dt.now().strftime('%H:%M:%S') + '] Starting files for hour ' + hour
    # Read netcdf
    fluxes_nc = Dataset('../fluxes_' + hour + '.nc', 'r')
    z_rho     = fluxes_nc.variables['rholev_zsea_rho'][:]
    z_rho[0]  = 0. # the fluxes are on a weird grid where the first level is actually the surface
    
    wind_nc   = Dataset('../wind_' + hour + '.nc', 'r')
    z_theta   = wind_nc.variables['thlev_zsea_theta'][:]
    
    zi_nc     = Dataset('../zi_' + hour + '.nc', 'r')
    times = wind_nc.variables['min10_0'][:]
    
    if hour == '09':
        # Initialise the arrays to store our analysis
        shf_cs = np.zeros((len(z_rho), len(x_cs)))
        lhf_cs = np.zeros((len(z_rho), len(x_cs)))
        tke_cs = np.zeros((len(z_theta), len(x_cs)))
        zi_cs  = np.zeros(len(x_cs))
        nt = 0.
    """
    Cross sections and horizontal slices within the boundary layer of turbulent 
    heat and moisture fluxes and resolved turbulent kinetic energy.
    """
    
    for it in xrange(len(times)):
        print '[' + dt.now().strftime('%H:%M:%S') + '] time ' + str(int(times[it]))
        print '[' + dt.now().strftime('%H:%M:%S') + '] Calculating cross sections'
        shf_cs += bilinear_interpolation(X, Y, fluxes_nc.variables[shf_key][it,:,:,:], x_cs, y_cs, kind = 2)
        print '[' + dt.now().strftime('%H:%M:%S') + '] Completed cross section for shf'
        lhf_cs += bilinear_interpolation(X, Y, fluxes_nc.variables[lhf_key][it,:,:,:], x_cs, y_cs, kind = 2)
        print '[' + dt.now().strftime('%H:%M:%S') + '] Completed cross section for lhf'
        # calculate the tke data
        u_p = (wind_nc.variables[u_key][it,:,:,:].transpose() - np.nanmean(wind_nc.variables[u_key][it,:,:,:], axis = (1, 2))).transpose()
        v_p = (wind_nc.variables[v_key][it,:,:,:].transpose() - np.nanmean(wind_nc.variables[v_key][it,:,:,:], axis = (1, 2))).transpose()
        w_p = (wind_nc.variables[w_key][it,:,:,:].transpose() - np.nanmean(wind_nc.variables[w_key][it,:,:,:], axis = (1, 2))).transpose()
        tke_data = 0.5*(u_p**2 + v_p**2 + w_p**2)
        tke_cs += bilinear_interpolation(X, Y, tke_data, x_cs, y_cs, kind = 2)
        print '[' + dt.now().strftime('%H:%M:%S') + '] Completed cross section for tke'
    zi_cs += np.nanmean(bilinear_interpolation(X, Y, zi_nc.variables[zi_key][:], x_cs, y_cs, kind = 2), axis = 0)
    
    nt += len(times)
    fluxes_nc.close()
    wind_nc.close()

shf_cs /= nt
lhf_cs /= nt
tke_cs /= nt

# Plot the Cross Section
fig = plt.figure(tight_layout = True)
# Sensible heat flux cross section
ax = fig.add_subplot(3, 1, 1)
SHF = ax.contourf(R, z_rho/1000., shf_cs, cmap = 'bwr', levels = np.linspace(-100., 100., 11), extend = 'both')
fig.colorbar(SHF, ax = ax)
ax.plot(R, zi_cs/1000., 'k', lw = 2)
ax.set_xlabel('Distance Downwind (km)')
ax.set_ylabel('Height (km)')
ax.set_ylim([0, 5])
ax.set_title('Sensible Heat Flux (W m$^{-2}$)')

# Latent heat flux cross section
ax = fig.add_subplot(3, 1, 2)
LHF = ax.contourf(R, z_rho/1000., Lv*lhf_cs, cmap = 'bwr', levels = np.linspace(-100., 100., 11), extend = 'both')
fig.colorbar(LHF, ax = ax)
ax.plot(R, zi_cs/1000., 'k', lw = 2)
ax.set_xlabel('Distance Downwind (km)')
ax.set_ylabel('Height (km)')
ax.set_ylim([0, 5])
ax.set_title('Latent Heat Flux (W m$^{-2}$)')

# Resolved TKE
ax = fig.add_subplot(3, 1, 3)
TKE = ax.contourf(R, z_theta/1000., tke_cs, cmap = 'viridis', levels = np.linspace(0., 5., 11), extend = 'max')
fig.colorbar(TKE, ax = ax)
ax.plot(R, zi_cs/1000., 'k', lw = 2)
ax.set_xlabel('Distance Downwind (km)')
ax.set_ylabel('Height (km)')
ax.set_ylim([0, 5])
ax.set_title('Turbulence Kinetic Energy (m$^{2}$ s$^{-2}$)')

plt.suptitle('Time-mean values for the period T+550 to T+900 mins', fontsize = 20)
plt.savefig('../turbulence_analysis.png', dpi = 100)
plt.close('all')


