import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate,integrate
from netCDF4 import Dataset

Lv = 2.501e6
#hours = ["{0.02d}".format(h) for h in xrange(0, 24, 3)]
hours = ['06', '09']
blf_key = u'STASH_m01s03i222'
cld_key = u'STASH_m01s09i182'
rho_key = u'STASH_m01s00i389'
zth_key = 'thlev_zsea_theta' # bl+cld is on theta levels
zrh_key = 'rholev_zsea_rho' # bl is on rho levels. ugh.

u_key = u'STASH_m01s00i002'
v_key = u'STASH_m01s00i003'

# Define a coordinate system
x = np.arange(0., 116000., 100.)
y = np.arange(0., 31900., 100.)
X, Y = np.meshgrid(x, y)
# define a period coordinate system
X_per = np.concatenate((X-116000., X, X+116000.), axis = 1)
X_per = np.concatenate((X_per, X_per, X_per), axis = 0)
Y_per = np.concatenate((Y, Y, Y), axis = 1)
Y_per = np.concatenate((Y_per - 31900., Y_per, Y_per +  31900.), axis = 0)

# Define the coordinates of the island centrepoint
x_c, y_c = [108000., 15950.]

base_path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'

### Should take about three hours to flow through the whole domain ###
# read the data
qinc_nc = Dataset(base_path + 'qinc_' + hours[1] + '.nc', 'r')
cld_times_key = [key for key in qinc_nc.variables.keys() if 'min' in key][0]
fluxes_nc = Dataset(base_path + 'fluxes_' + hours[1] + '.nc', 'r')
blf_times_key = [key for key in fluxes_nc.variables.keys() if 'min' in key][0]
wind_nc = Dataset(base_path + 'wind_' + hours[1] + '.nc', 'r')

cld_times = qinc_nc.variables[cld_times_key][:]*1.
blf_times = fluxes_nc.variables[blf_times_key][:]*1.

z_rho = fluxes_nc.variables[zrh_key][:]
z_the = qinc_nc.variables[zth_key][:]
"""
We want to compare the total boundary layer moisture fluxes with the ls cloud moisture fluxes
-> in qinc there is a variable that contains the bl and ls cld
    This is the tendency due to the boundary layer scheme and the large scale cloud scheme kg/kg/timestep
-> in fluxes there is a variable that contains bl
    This is the boundary layer moisture fluxes kg/m2/s
To get the ls cloud, simply subtract the bl from bl+ls cld, we want to consider vertically integrated quantities, and follow the flow (Lagrangian perspective)
"""

# Pick a point over the island and extrapolate forward and backward
# e.g. pick the island centrepoint

BLF = np.zeros_like(cld_times)
CLD = np.zeros_like(cld_times)
x_p, y_p = [x_c, y_c]
x_ps = [x_p]
y_ps = [y_p]
R = np.sqrt((X_per-x_p)**2. + (Y_per-y_p)**2.)
mask = np.where((R < 2000.), 1., 0.) # only consider grid cells that are less than 2 km away from the centrepoint
mask = np.max([mask[(0+j*319):(319+j*319),(0+i*1160):(1160+i*1160)] for j in xrange(3) for i in xrange(3)], axis = 0)

#advect with the boundary layer winds lowest ~ 500 m
iz500 = np.where(np.abs(z_the - 500.) == np.min(np.abs(z_the - 500.)))[0][0]+1
for it in xrange(len(cld_times)):
    # Define the current time
    target_time = cld_times[it]

    blf_it = np.where(blf_times == target_time)[0][0]
    cld_it = np.where(cld_times == target_time)[0][0]
    
    # Define our point at this time, starting from the island centre point
    if it != 0:
        x_p = (x_p + np.nanmean(np.where(mask, integrate.trapz(y = wind_nc.variables[u_key][blf_it,1:(iz500+1),:,:], x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))*(cld_times[it] - cld_times[it-1])*60.)%116000. # wind nc starts at z = 0, while other nc on theta levels start at z = 2 m
        y_p = (y_p + np.nanmean(np.where(mask, integrate.trapz(y = wind_nc.variables[v_key][blf_it,1:(iz500+1),:,:], x = z_the[:iz500], axis = 0)/np.max(z_the[:iz500]), np.nan))*(cld_times[it] - cld_times[it-1])*60.)%31900.
        x_ps.append(x_p)
        y_ps.append(y_p)
    print 'x_p = ' + str(int(x_p)) + ', y_p = ' + str(int(y_p))
    # Redefine our mask
    R = np.sqrt((X_per - x_p)**2. + (Y_per - y_p)**2.)
    mask = np.where((R < 2000.), 1., 0.)
    mask = np.max([mask[(0+j*319):(319+j*319),(0+i*1160):(1160+i*1160)] for j in xrange(3) for i in xrange(3)], axis = 0)
    
    # Compute the fluxes in that ring.
    BLF[it] = np.nanmean(np.where(mask, Lv*fluxes_nc.variables[blf_key][blf_it,0,:,:], np.nan))
    rho_data = interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][blf_it,:,:,:]*1., axis = 0, fill_value = 'extrapolate')(z_the)
    CLD[it] = np.nanmean(np.where(mask, integrate.trapz(y = Lv*rho_data*qinc_nc.variables[cld_key][cld_it,:,:,:]/3., x = z_the, axis = 0), np.nan))

fig = plt.figure(tight_layout = True, figsize = (12,6))
ax = fig.add_subplot(1, 2, 1)
ax.plot(cld_times/60., BLF, 'b--', label = '$\\rho L_{v} (\overline{w^{\prime} q_{t}^{\prime}})_{sfc}$')
ax.plot(cld_times/60., CLD, 'b', label = 'BL + LS CLD = $L_{v} \int_{0}^{z_{max}}-\\frac{\partial (w^{\prime} q_{v}^{\prime})}{\partial z} + \\frac{\partial q_{v}}{\partial t}_{ls rain} dz$ ?')
plt.legend(loc = 'lower right')
ax.set_xlabel('Time (hrs)')
ax.set_ylabel('Flux (W/m$^{2}$)')

ax = fig.add_subplot(1, 2, 2, adjustable = 'box', aspect = 1)
R = np.sqrt((X-x_c)**2. + (Y-y_c)**2.)
ax.contourf(X/1000., Y/1000., np.where((R < 1000*np.sqrt(50./np.pi)), 1., 0), levels = [0., 1e-16, 1.], colors = ['w', 'k'])
for i in xrange(len(x_ps)):
    R = np.sqrt((X_per - x_ps[i])**2. + (Y_per - y_ps[i])**2.)
    mask = np.where((R < 2000.), 1., 0.)
    mask = np.max([mask[(0+j*319):(319+j*319),(0+i*1160):(1160+i*1160)] for j in xrange(3) for i in xrange(3)], axis = 0)
    ax.contourf(X/1000., Y/1000., mask, levels = [1e-16, 1.], colors = ['r'], alpha = 0.5)

ax.plot(np.array(x_ps)/1000., np.array(y_ps)/1000., '--r*')
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
plt.savefig('../lagrangian_latent_heat_fluxes.png', dpi = 150)
plt.show()


