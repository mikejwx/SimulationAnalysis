import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import lwp_key, u_key, v_key, zi_new_key
from analysis_tools import downwind_rectangle, fromComponents
from scipy import interpolate, integrate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Read the data
path = '/nerc/n02/n02/xb899100/CloudTrail/Control2/'
label = path.split('/')[-2]
my_data = {}
with Dataset(path + 'lwp_00.nc', 'r') as lwp_nc:
    my_data[lwp_key] = lwp_nc.variables[lwp_key][:]*1.
    times_key = [tkey for tkey in lwp_nc.variables.keys() if 'min' in tkey][0]
    my_data['lwp_times'] = lwp_nc.variables[times_key][:]*1.

with Dataset(path + 'wind_00.nc', 'r') as wind_nc:
    my_data[u_key] = wind_nc.variables[u_key][1,:,:,:]*1.
    my_data[v_key] = wind_nc.variables[v_key][1,:,:,:]*1.
    my_data['z']   = wind_nc.variables['thlev_zsea_theta'][:]*1.

with Dataset(path + 'zi_00.nc', 'r') as zi_nc:
    my_data[zi_new_key] = zi_nc.variables[zi_new_key][0,:,:]*1.

# Grab the mean initial wind speed and direction from the boundary layer
zi_mean = np.nanmean(my_data[zi_new_key])
z_int = np.arange(0., zi_mean, 1.)
U_fun = interpolate.interp1d(x = my_data['z'], y = np.nanmean(my_data[u_key], axis = (1, 2)), fill_value = 'extrapolate')(z_int)
V_fun = interpolate.interp1d(x = my_data['z'], y = np.nanmean(my_data[v_key], axis = (1, 2)), fill_value = 'extrapolate')(z_int)

U_mean = integrate.trapz(x = z_int, y = U_fun)/z_int.max()
V_mean = integrate.trapz(x = z_int, y = V_fun)/z_int.max()

wind_spd, wind_dir = fromComponents(U_mean, V_mean)

# Define the original coordinate system
dx = 100.0
x, y = np.meshgrid(np.arange(my_data[lwp_key].shape[2])*dx, np.arange(my_data[lwp_key].shape[1])*dx)
island_area = 50.0
island_radius = np.sqrt(island_area/np.pi)*1000.0
x_c = 100000.0 + island_radius
y_c = 4*island_radius
R = np.sqrt((x - x_c)**2 + (y - y_c)**2)

# Define the downwind rectangle
mask, y_prime, x_prime = downwind_rectangle(wind_dir, x_c, y_c, x, y, island_radius, dist_0 = -12000.0, dist_1 = 100000.0, half_width = 3000.0)

# Compute the daytime cloud frequency
start = np.where(my_data['lwp_times'] == 360.)[0][0]
end   = np.where(my_data['lwp_times'] == 1080.)[0][0] + 1

cloud_frequency = np.nanmean(np.where(my_data[lwp_key][start:end,:,:] > 0, 1.0, 0.0), axis = 0)
cloud_frequency_before = np.nanmean(np.where(my_data[lwp_key][:(start+1),:,:] > 0, 1.0, 0.0), axis = 0)

# Collapse cloud frequency from array[x_prime, y_prime] -> array[x_prime]
cloud_frequency_x = np.array([np.nanmean(np.where((x_0 <= x_prime)*(x_prime <= (x_0 + dx)), mask*cloud_frequency, np.nan)) for x_0 in np.arange(-12000.0, 100000.1, dx)])
x_prime_new = np.arange(-12000.0, 100000.1, dx)

new_mask = np.where(mask==mask, 1.0, 0.0)

fig = plt.figure(figsize = (6, 6))
axa = fig.add_subplot(2, 1, 1, adjustable = 'box', aspect = 75)
axa.plot(x_prime_new/1000., cloud_frequency_x, 'k', lw = 2)
axa.plot([-20., 100.], np.nanmean(cloud_frequency*(1-new_mask))*np.ones(2), 'grey', lw = 0.5)
axa.plot([-2*island_radius/1000., 0.0], [0.005, 0.005], 'red', lw = 5)
axa.set_ylim([0, 0.5])
axa.set_ylabel('Cloud Frequency')
axa.set_xlim([-12, 100.])
axa.set_xlabel('Distance downwind (km)')
axa.text(25, 0.15, '"void"', verticalalignment='bottom', horizontalalignment='center')
axa.text(-10, 0.425, 'a)')
axa.text((-island_radius/1000.), 0.00875, 'island', color = 'red', verticalalignment='bottom', horizontalalignment='center')

axb = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
cf = axb.contourf(x/1000., y/1000., cloud_frequency, cmap = 'Greys_r', levels = np.arange(0., 0.61, 0.075), extend = 'max')
cbar_ax = inset_axes(axb, width = "100%", height = "5%", loc = 3, bbox_to_anchor = (0, -0.5, 1.0, 2.0), bbox_transform = axb.transAxes, borderpad = 0.0)
cbar = fig.colorbar(cf, cax = cbar_ax, orientation = 'horizontal', label = 'Cloud Frequency')
cbar.set_ticks([0, 0.15, 0.30, 0.45, 0.6])
cbar.set_ticklabels(['0.0', '0.15', '0.3', '0.45', '0.6'])
axb.contour(x/1000., y/1000., new_mask, levels = [0.5], colors = ['b'])
axb.contour(x/1000., y/1000., R, levels = [island_radius], colors = ['red'], linewidths = [2])
axb.text((5./115.)*116, (0.425/0.5)*32, 'b)', bbox = dict(ec = 'None',fc = 'w', pad = 0.2))
axb.set_ylabel('y (km)')
axb.set_xlabel('x (km)')
plt.subplots_adjust(left = 0.10, bottom = 0.175, right = 0.90, wspace = 0.0, hspace = -0.10)
plt.savefig('../Ch5_Figure3.png', dpi = 250, bbox_inches = 'tight')
plt.show()

