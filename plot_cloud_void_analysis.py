import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate, integrate
from analysis_tools import downwind_rectangle, fromComponents
from STASH_keys import theta_key, q_key, w_key, u_key, v_key, zi_new_key, lcl_key

# We want to show that the warm plume is capping the boundary layer
# We suspect that this is responsible for the void downwind of the termination 
# point of the surface warm plume

### Read some data from the reference experiment ###
path = '/nerc/n02/n02/xb899100/CloudTrail/Control_0800m_HRIC_INV/'
my_data = {}
# Read thermodynamic data
with Dataset(path + 'bouy_09.nc', 'r') as bouy_nc:
    my_data[theta_key] = bouy_nc.variables[theta_key][:]*1.
    my_data[q_key]     = bouy_nc.variables[q_key][:]*1.
    my_data['z']       = bouy_nc.variables['thlev_zsea_theta'][:]*1.
    timeKey            = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
    my_data['times']   = bouy_nc.variables[timeKey][:]*1.
    
# Read wind data
with Dataset(path + 'wind_09.nc', 'r') as wind_nc:
    my_data[w_key] = wind_nc.variables[w_key][:]*1.
    my_data[u_key] = wind_nc.variables[u_key][:]*1.
    my_data[v_key] = wind_nc.variables[v_key][:]*1.

# Read the boundary layer heights
with Dataset(path + 'zi_09.nc', 'r') as zi_nc:
    my_data[zi_new_key] = zi_nc.variables[zi_new_key][:]*1.
    my_data[lcl_key]    = zi_nc.variables[lcl_key][:]*1.

# Create the horizontal coordinate system and define the island
dx = 800.0
x, y = np.meshgrid(np.arange(my_data[theta_key].shape[3])*dx, np.arange(my_data[theta_key].shape[2])*dx)

# Radius of the island, m
island_area = 50.0 # km^2
R_i = 1000.0*np.sqrt(island_area/np.pi)
x_c = 108000.0#100000.0 + R_i
y_c = 16000.0#4*R_i

# Distance from the island, m
R = np.sqrt((x - x_c)**2 + (y - y_c)**2)

# Find the mean boundary layer wind speed and direction
zi_mean = np.nanmean(my_data[zi_new_key][0,:,:])
iz      = np.where(np.abs(my_data['z'] - zi_mean) == np.min(np.abs(my_data['z'] - zi_mean)))[0][0]
U_mean  = np.nanmean(integrate.trapz(x = my_data['z'][:iz], y = my_data[u_key][0,:iz,:,:], axis = 0)/my_data['z'][iz])
V_mean  = np.nanmean(integrate.trapz(x = my_data['z'][:iz], y = my_data[v_key][0,:iz,:,:], axis = 0)/my_data['z'][iz])

wind_spd, wind_dir = fromComponents(U_mean, V_mean)

# Use the wind direction to define our downwind rectangular region
mask, y_prime, x_prime = downwind_rectangle(wind_dir, x_c, y_c, x, y, R_i, dist_0 = -4*R_i, dist_1 = 80000.0, half_width = 5000.0)

# Compute the potential temperature cross section

x_prime_new  = np.arange(-4*R_i, 80000.1, dx)
theta_cs     = np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[theta_key][it,:,:,:], np.nan), axis = (1, 2)) for x_0 in x_prime_new] for it in [-1]]) # array[t, x', z]
w_cs         = np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[w_key][it,:,:,:], np.nan), axis = (1, 2)) for x_0 in x_prime_new] for it in [-1]])
zi_cs       = np.transpose(np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[zi_new_key][it,:,:], np.nan)) for x_0 in x_prime_new] for it in [-1]]))
lcl_cs      = np.transpose(np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[lcl_key][it,:,:], np.nan)) for x_0 in x_prime_new] for it in [-1]]))

# Compute the anomaly from the surface parcel ascent
theta_parcel_cs = np.transpose(theta_cs[0,:,:]) - theta_cs[0,:,0]
theta_parcel_anom_cs = np.transpose(np.transpose(theta_parcel_cs) - np.nanmean(theta_parcel_cs, axis = 1))
xp, z = np.meshgrid(x_prime_new, my_data['z'])

# make the plot
fig = plt.figure(tight_layout = True)
axa = fig.add_subplot(2, 1, 1, adjustable = 'box', aspect = 10)
my_plt = axa.contourf(xp/1000., z/1000., theta_parcel_cs, cmap = 'bwr', vmin = -1, vmax = 1, levels = [lev for lev in np.arange(-1., 1.01, 0.1) if lev != 0.0], extend = 'both')
plt.colorbar(my_plt, ax = axa, orientation = 'horizontal', label = u'$(\\theta - \\theta_{sfc,p})$', pad = 0.25, shrink = 0.75)
axa.plot(x_prime_new/1000., lcl_cs/1000., 'k', lw = 2)
axa.plot(x_prime_new/1000., zi_cs/1000., 'k')
axa.set_ylim([0, 2.5])
axa.set_ylabel('Height (km)')


axb = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 10)
my_plt = axb.contourf(xp/1000., z/1000., theta_parcel_anom_cs, cmap = 'bwr', vmin = -1, vmax = 1, levels = [lev for lev in np.arange(-1, 1.01, 0.1) if lev != 0.0], extend = 'both')
plt.colorbar(my_plt, ax = axb, orientation = 'horizontal', label = u'$(\\theta - \\theta_{sfc,p}) - \overline{(\\theta - \\theta_{sfc,p})}$', pad = 0.25, shrink = 0.75)
axb.plot(x_prime_new/1000., lcl_cs/1000., 'k', lw = 2)
axb.plot(x_prime_new/1000., zi_cs/1000., 'k')
axb.set_ylim([0, 2.5])
axb.set_ylabel('Height (km)')
axb.set_xlabel('Distance Downwind (km)')
plt.savefig('../cloud_void_analysis_0800m_HRIC_INV.png', dpi = 150, bbox_inches = 'tight')
plt.show()

# make a separate plot to show the vertical velocities
fig = plt.figure(tight_layout = True)
ax0 = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 10)
my_plt = ax0.contourf(xp/1000., z/1000., np.transpose(w_cs[0]), cmap = 'bwr', vmin = -1, vmax = 1, levels = [lev for lev in np.arange(-1., 1.01, 0.1) if lev != 0.0], extend = 'both')
plt.colorbar(my_plt, ax = ax0, orientation = 'horizontal', label = u'w (m s$^{-1}$)', pad = 0.25, shrink = 0.75)
ax0.plot(x_prime_new/1000., lcl_cs/1000., 'k', lw = 2)
ax0.plot(x_prime_new/1000., zi_cs/1000., 'k')
ax0.set_ylim([0, 2.5])
ax0.set_ylabel('Height (km)')
plt.savefig('../cloud_void_analysis2.png', dpi = 150, bbox_inches = 'tight')
plt.show()

distance_downwind = wind_spd*zi_mean/(integrate.trapz(x = my_data['z'][:iz], y = mask*my_data[w_key][it,:iz,:,:], axis = 0)/my_data['z'][iz-1])
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
my_plt = ax1.contourf(x_prime, y_prime, distance_downwind, cmap = 'bwr', levels = np.arange(-20000, 20000.1, 1000), vmin = -20000, vmax = 20000, extend = 'both')
ax1.contour(x_prime, y_prime, R, levels = [R_i], colors = ['k'])
ax1.set_ylim([-5000, 5000])
ax1.set_xlim([-16000, 80000])
plt.colorbar(my_plt, ax = ax1, orientation = 'horizontal')
plt.show()

distance_downwind_x  = np.array([np.nanmin(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx))*(distance_downwind > 0), distance_downwind, np.nan)) for x_0 in np.arange(-4*R_i-dx/2., 80000.1, dx)])
plt.plot(np.arange(-4*R_i-dx/2, 80000.1, dx), distance_downwind_x)
plt.show()
distance_downwind_xy = np.nanmean(distance_downwind_x)
print str(round(distance_downwind_xy/1000., 3)) + ' km'
print str(round(distance_downwind_xy/wind_spd, 0)) + ' seconds'
print str(round(zi_mean/(distance_downwind_xy/wind_spd), 1)) + ' m/s'

