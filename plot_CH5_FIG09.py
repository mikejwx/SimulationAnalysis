import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate, integrate
from analysis_tools import downwind_rectangle, fromComponents
from STASH_keys import w_key, u_key, v_key, zi_new_key, lcl_key, s_key, n_key
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl

# We want to show that the warm plume is capping the boundary layer
# We suspect that this is responsible for the void downwind of the termination 
# point of the surface warm plume

### Read some data from the reference experiment ###
path = '/nerc/n02/n02/xb899100/CloudTrail/Control2/'

my_data = {}
# Read wind data
with Dataset(path + 'wind_09.nc', 'r') as wind_nc:
    my_data['z']   = wind_nc.variables['thlev_zsea_theta'][:]*1.
    iz_max = np.where(np.abs(my_data['z'] - 3000.0) == np.min(np.abs(my_data['z'] - 3000.0)))[0][0] + 1
    my_data['z']   = my_data['z'][:iz_max]
    my_data[u_key] = wind_nc.variables[u_key][0,:iz_max,:,:]*1.
    my_data[v_key] = wind_nc.variables[v_key][0,:iz_max,:,:]*1.
    my_data[s_key] = wind_nc.variables[s_key][:,:iz_max,:,:]*1.
    my_data[n_key] = wind_nc.variables[n_key][:,:iz_max,:,:]*1.

# Read the boundary layer heights
with Dataset(path + 'zi_09.nc', 'r') as zi_nc:
    my_data[zi_new_key] = zi_nc.variables[zi_new_key][:]*1.
    my_data[lcl_key]    = zi_nc.variables[lcl_key][:]*1.

# Create the horizontal coordinate system and define the island
dx = 100.0
x, y = np.meshgrid(np.arange(my_data[s_key].shape[3])*dx, np.arange(my_data[s_key].shape[2])*dx)

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
U_mean  = np.nanmean(integrate.trapz(x = my_data['z'][:iz], y = my_data[u_key][:iz,:,:], axis = 0)/my_data['z'][iz])
V_mean  = np.nanmean(integrate.trapz(x = my_data['z'][:iz], y = my_data[v_key][:iz,:,:], axis = 0)/my_data['z'][iz])

wind_spd, wind_dir = fromComponents(U_mean, V_mean)

# Use the wind direction to define our downwind rectangular region
d_min = -12000.0
d_max = 50000.0
d_mid = 3000.0
mask, y_prime, x_prime = downwind_rectangle(wind_dir, x_c, y_c, x, y, R_i, dist_0 = d_min, dist_1 = d_max, half_width = d_mid)

# Compute the s-wind cross section
x_prime_new  = np.arange(d_min, d_max + 0.1, dx)
s_xcs        = np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[s_key][it,:,:,:], np.nan), axis = (1, 2)) for x_0 in x_prime_new] for it in [-1]]) # array[t, x', z]
zi_xcs       = np.transpose(np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[zi_new_key][it,:,:], np.nan)) for x_0 in x_prime_new] for it in [-1]]))
lcl_xcs      = np.transpose(np.array([[np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*my_data[lcl_key][it,:,:], np.nan)) for x_0 in x_prime_new] for it in [-1]]))

# Compute the s-wind anomaly from the surface parcel ascent
s_anom_xcs    = s_xcs[it,:,:] - my_data[s_key][it,:,:,:].mean(axis = (1, 2))

# Compute the n-wind cross section
y_prime_new   = np.arange(-d_mid, d_mid + 0.1, dx)
n_ycs         = np.array([[np.nanmean(np.where((y_0 <= y_prime)*(y_prime <= (y_0 + dx)), mask*my_data[n_key][it,:,:,:], np.nan), axis = (1, 2)) for y_0 in y_prime_new] for it in [-1]])
zi_ycs        = np.transpose(np.array([[np.nanmean(np.where((y_0 <= y_prime)*(y_prime <= (y_0 + dx)), mask*my_data[zi_new_key][it,:,:], np.nan)) for y_0 in y_prime_new] for it in [-1]]))
lcl_ycs       = np.transpose(np.array([[np.nanmean(np.where((y_0 <= y_prime)*(y_prime <= (y_0 + dx)), mask*my_data[lcl_key][it,:,:], np.nan)) for y_0 in y_prime_new] for it in [-1]]))

xp, xz = np.meshgrid(x_prime_new, my_data['z'])
yp, yz = np.meshgrid(y_prime_new, my_data['z'])

### Make the plot ###
# Define some plotting parameters
my_levels = np.array([level for level in np.arange(-1.0, 1.1, 0.25) if level != 0.0])

fig = plt.figure(figsize = (6, 8))
# Panel A, the difference between the surface parcel and the environment theta
axa = fig.add_subplot(2, 1, 1, adjustable = 'box', aspect = 10)
cf = axa.contourf(xp/1000., xz/1000., np.transpose(s_anom_xcs), cmap = 'bwr', levels = my_levels, extend = 'both')
cf.cmap.set_over('firebrick')
cf.cmap.set_under('navy')
cbar_ax = inset_axes(axa, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (0, -0.4, 1.0, 2.0), bbox_transform = axa.transAxes, borderpad = 0.0)
cbar = fig.colorbar(cf, cax = cbar_ax, orientation = 'horizontal', label = u'$s^{\prime}$ (m s$^{-1}$)')
cbar.set_ticks(np.array([-1, -0.5, 0, 0.5, 1]))
axa.plot(x_prime_new/1000., lcl_xcs/1000., 'k', lw = 2)
axa.plot(x_prime_new/1000., zi_xcs/1000., 'k')
axa.set_yticks(np.arange(0, 3.1, 0.5))
axa.set_ylim([0, 2.5])
axa.set_xlim([-12, 50])
axa.set_xlabel('Distance Downwind (km)')
axa.set_ylabel('Height (km)')
axa.text(-10, 2.125, 'a)', bbox = {'edgecolor':'none', 'facecolor':'white','pad':0.0})

# Panel B, the anomaly of A from the horizontal mean
axb = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
cf = axb.contourf(yp/1000., yz/1000., np.transpose(n_ycs[0,:,:]), cmap = 'bwr', levels = my_levels, extend = 'both')
cf.cmap.set_over('firebrick')
cf.cmap.set_under('navy')
cbar_ax = inset_axes(axb, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (0, -0.4, 1.0, 2.0), bbox_transform = axb.transAxes, borderpad = 0.0)
cbar = fig.colorbar(cf, cax = cbar_ax, orientation = 'horizontal', label = u'$n^{\prime}$ (m s$^{-1}$)')
cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
axb.plot(y_prime_new/1000., lcl_ycs/1000., 'k', lw = 2)
axb.plot(y_prime_new/1000., zi_ycs/1000., 'k')
axb.set_yticks(np.arange(0, 3.1, 0.5))
axb.set_ylim([0, 2.5])
axb.set_xlim([-3, 3])
axb.set_ylabel('Height (km)')
axb.text(-2.8, 2.125, 'b)', bbox = {'edgecolor':'none', 'facecolor':'white','pad':0.0})
axb.set_xlabel('Distance Across-wind (km)')
plt.subplots_adjust(left = 0.10, bottom = 0.175, right = 0.90, wspace = 0.0, hspace = 0.3)
plt.savefig('../Ch5_Figure09.png', dpi = 250, bbox_inches = 'tight')
plt.show()

