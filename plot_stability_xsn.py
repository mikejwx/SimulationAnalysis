import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate, integrate
from analysis_tools import downwind_rectangle, fromComponents, lcl, find_h
from STASH_keys import theta_key, q_key, w_key, u_key, v_key, zi_new_key, lcl_key, pthe_key
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl

# We want to show that the warm plume is capping the boundary layer
# We suspect that this is responsible for the void downwind of the termination 
# point of the surface warm plume

### Read some data from the reference experiment ###
path = '/gws/nopw/j04/paracon_rdg/users/mcjohnston/experiments/Control/'

my_data = {}
# Read thermodynamic data
with Dataset(path + 'bouy/bouy_09.nc', 'r') as bouy_nc:
    my_data[theta_key] = bouy_nc.variables[theta_key][:]*1.
    my_data[q_key]     = bouy_nc.variables[q_key][:]*1.
    my_data['z']       = bouy_nc.variables['thlev_zsea_theta'][:]*1.
    timeKey            = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
    my_data['times']   = bouy_nc.variables[timeKey][:]*1.
    my_data[w_key]     = bouy_nc.variables[w_key][:]*1.
    # Compute lcl
    my_data[lcl_key]   = lcl(bouy_nc.variables[temp_key][:]*1.,
                             bouy_nc.variables[q_key][:]*1.,
                             bouy_nc.variables[pthe_key][:]*1.,
                             my_data['z'])
                             
# Read wind data
with Dataset(path + 'u/u_09.nc', 'r') as u_nc:
    my_data[u_key] = u_nc.variables[u_key][:]*1.

with Dataset(path + 'v/v_09.nc', 'r') as v_nc:
    my_data[v_key] = v_nc.variables[v_key][:]*1.

# Compute zi
my_data['thetav'] = my_data[theta_key]*(1 + 0.608*my_data[q_key])
my_data[zi_new_key] = find_h(my_data['thetav'], my_data[u_key], my_data[v_key], my_data['z'])

# Create the horizontal coordinate system and define the island
dx = 100.0
x, y = np.meshgrid(np.arange(my_data[theta_key].shape[3])*dx, np.arange(my_data[theta_key].shape[2])*dx)

# Radius of the island, m
island_area = 50.0 # km^2
R_i = 1000.0*np.sqrt(island_area/np.pi)
x_c = 100000.0 + R_i
y_c = 4*R_i

# Distance from the island, m
R = np.sqrt((x - x_c)**2 + (y - y_c)**2)

# Find the mean boundary layer wind speed and direction
zi_mean = np.nanmean(my_data[zi_new_key][0,:,:])
iz      = np.where(np.abs(my_data['z'] - zi_mean) == np.min(np.abs(my_data['z'] - zi_mean)))[0][0]
U_mean  = np.nanmean(integrate.trapz(x = my_data['z'][:iz], y = my_data[u_key][0,:iz,:,:], axis = 0)/my_data['z'][iz])
V_mean  = np.nanmean(integrate.trapz(x = my_data['z'][:iz], y = my_data[v_key][0,:iz,:,:], axis = 0)/my_data['z'][iz])

wind_spd, wind_dir = fromComponents(U_mean, V_mean)

# Use the wind direction to define our downwind rectangular region
mask, y_prime, x_prime = downwind_rectangle(wind_dir, x_c, y_c, x, y, R_i, dist_0 = -12000.0, dist_1 = 100000.0, half_width = 3000.0)

# Compute the potential temperature cross section

x_prime_new  = np.arange(-12000.0, 100000.1, dx)
theta_cs     = np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[theta_key][it,:,:,:], np.nan), axis = (1, 2)) for x_0 in x_prime_new] for it in [-1]]) # array[t, x', z]
w_cs         = np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[w_key][it,:,:,:], np.nan), axis = (1, 2)) for x_0 in x_prime_new] for it in [-1]])
zi_cs       = np.transpose(np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[zi_new_key][it,:,:], np.nan)) for x_0 in x_prime_new] for it in [-1]]))
lcl_cs      = np.transpose(np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[lcl_key][it,:,:], np.nan)) for x_0 in x_prime_new] for it in [-1]]))

# Compute the anomaly from the surface parcel ascent
theta_parcel_cs = np.transpose(theta_cs[0,:,:]) - theta_cs[0,:,0]
theta_parcel_anom_cs = np.transpose(np.transpose(theta_parcel_cs) - np.nanmean(theta_parcel_cs, axis = 1))
xp, z = np.meshgrid(x_prime_new, my_data['z'])

### Make the plot ###
# Define some plotting parameters
my_levels = np.array([-2.0, -1.0, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 1.0, 2.0])
n_levels = float(len(my_levels))
my_cmap = mpl.cm.get_cmap('bwr')
my_colors = my_cmap((np.arange(n_levels)+1.0)/(n_levels))

fig = plt.figure(figsize = (6, 6))
# Panel A, the difference between the surface parcel and the environment theta
axa = fig.add_subplot(2, 1, 1, adjustable = 'box', aspect = 10)
cf = axa.contourf(xp/1000., z/1000., theta_parcel_cs, cmap = 'bwr', levels = my_levels, extend = 'both')
cf.cmap.set_over('firebrick')
cf.cmap.set_under('navy')
cbar_ax = inset_axes(axa, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (0, -0.25, 1.0, 2.0), bbox_transform = axa.transAxes, borderpad = 0.0)
cbar = fig.colorbar(cf, cax = cbar_ax, orientation = 'horizontal', label = u'$(\\theta - \\theta_{sfc,p})$ (K)')
cbar.set_ticks([-2, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 2])
cbar.ax.set_xticklabels([-2, -1, -0.5, -0.25, -0.1, '', 0.1, 0.25, 0.5, 1, 2])
axa.plot(x_prime_new/1000., lcl_cs/1000., 'k', lw = 2)
axa.plot(x_prime_new/1000., zi_cs/1000., 'k')
axa.set_yticks(np.arange(0, 1.6, 0.5))
axa.set_ylim([0, 1.5])
axa.set_xlim([-12, 50])
axa.set_xticklabels([''])
axa.set_ylabel('Height (km)')
axa.text(-10, 1.275, 'a)', bbox = {'edgecolor':'none', 'facecolor':'white','pad':0.0})

# Panel B, the anomaly of A from the horizontal mean
axb = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 10)
cf = axb.contourf(xp/1000., z/1000., theta_parcel_anom_cs, cmap = 'bwr', levels = my_levels, extend = 'both')
cf.cmap.set_over('firebrick')
cf.cmap.set_under('navy')
cbar_ax = inset_axes(axb, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (0, -0.5, 1.0, 2.0), bbox_transform = axb.transAxes, borderpad = 0.0)
cbar = fig.colorbar(cf, cax = cbar_ax, orientation = 'horizontal', label = u'$(\\theta - \\theta_{sfc,p}) - \overline{(\\theta - \\theta_{sfc,p})}$ (K)')
cbar.set_ticks([-2, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 2])
cbar.ax.set_xticklabels([-2, -1, -0.5, -0.25, -0.1, '', 0.1, 0.25, 0.5, 1, 2])
axb.plot(x_prime_new/1000., lcl_cs/1000., 'k', lw = 2)
axb.plot(x_prime_new/1000., zi_cs/1000., 'k')
axb.set_yticks(np.arange(0, 1.6, 0.5))
axb.set_ylim([0, 1.5])
axb.set_xlim([-12, 50])
axb.set_ylabel('Height (km)')
axb.text(-10, 1.275, 'b)', bbox = {'edgecolor':'none', 'facecolor':'white','pad':0.0})
axb.set_xlabel('Distance Downwind (km)')
plt.subplots_adjust(left = 0.10, bottom = 0.175, right = 0.90, wspace = 0.0, hspace = -0.10)
#plt.savefig('../Ch5_Figure06.png', dpi = 250, bbox_inches = 'tight')
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
print(str(round(distance_downwind_xy/1000., 3)) + ' km')
print(str(round(distance_downwind_xy/wind_spd, 0)) + ' seconds')
print(str(round(zi_mean/(distance_downwind_xy/wind_spd), 1)) + ' m/s')

