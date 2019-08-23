"""
Cheap and nasty way to get lwp from convection scheme?
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import integrate, interpolate
from STASH_keys import rho_key

path = '/nerc/n02/n02/xb899100/CloudTrail/Control_0800m_param/'
qlconv_key = u'STASH_m01s05i213'
# read the convection nc file
days  = ["{0:02d}".format(d) for d in range(0, 16, 4)]
data = {}
for day in days:
    conv_nc = Dataset(path + 'conv_' + day + '.nc', 'r')
    fluxes_nc = Dataset(path + 'fluxes_' + day + '.nc', 'r')
    if day == days[0]:
        # read the convective condensed water content
        data[qlconv_key] = conv_nc.variables[qlconv_key][:]*1000.
        # read the vertical coordinate for theta and rho levels
        z_rho = fluxes_nc.variables['rholev_zsea_rho'][:]*1.
        z = conv_nc.variables['thlev_zsea_theta'][:]*1.
        # interpolate the rho data onto the theta levels
        data[rho_key] = interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][:], fill_value = 'extrapolate', axis = 1)(z)
        times = conv_nc.variables[u'min10'][:]*1.
    else:
        data[qlconv_key] = np.concatenate((data[qlconv_key], conv_nc.variables[qlconv_key][:]*1000.0), axis = 0)
        interpolated_rho = interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][:], fill_value = 'extrapolate', axis = 1)(z)
        data[rho_key]    = np.concatenate((data[rho_key], interpolated_rho), axis = 0)
        times = np.concatenate((times, conv_nc.variables[u'min10'][:]*1.), axis = 0)
    conv_nc.close()

# convert to liquid water path
data['lwp'] = integrate.trapz(y = data[qlconv_key]*data[rho_key][1,:,:,:], x = z, axis = 1)
times += 240.

# plot
for it in range(48,96):#range(len(times)):
    print it
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
    my_lwp = ax.contourf(data['lwp'][it,:,:], cmap = 'Greys', levels = np.linspace(0.001, 5000.0, 11), extend = 'max')
    plt.colorbar(my_lwp, ax = ax, label = u'LWP (g kg$^{-1}$)', orientation = 'horizontal')
    plt.title('Time = ' + str(int(times[it])))
    plt.show()

"""
    plt.savefig('../fake_conv_lwp/T' + "{0:04}".format(int(times[it])) + '.png', dpi = 150)
    plt.close('all')
"""

