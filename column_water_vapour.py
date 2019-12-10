import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import integrate, interpolate

hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]
base_path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'
mv_key = u'STASH_m01s00i391'
rho_key = u'STASH_m01s00i389'
for hour in hours:
    mr_nc = Dataset(base_path + 'mr_' + hour + '.nc', 'r')
    fluxes_nc = Dataset(base_path + 'fluxes_' + hour + '.nc', 'r')
    if hour == '00':
        z_theta = mr_nc.variables['thlev_zsea_theta'][:]
        z_rho = fluxes_nc.variables['rholev_zsea_rho'][:]
        CIWV = integrate.trapz(y = mr_nc.variables[mv_key][:]*interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][:], fill_value = 'extrapolate', axis = 1)(z_theta), x = z_theta, axis = 1)
        times = mr_nc.variables['min10_0'][:]
    else:
        CIWV = np.concatenate((CIWV, integrate.trapz(y = mr_nc.variables[mv_key][:]*interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][:], fill_value = 'extrapolate', axis = 1)(z_theta), x = z_theta, axis = 1)), axis = 0)
        times = np.concatenate((times, mr_nc.variables['min10_0'][:]),axis =0)
    mr_nc.close()
    fluxes_nc.close()

CIWV_mean = np.mean(CIWV, axis = (1, 2))
fig = plt.figure(tight_layout = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(times/60., CIWV_mean, 'k', lw = 2)
ax.set_ylabel('Column Integrated Water Vapour (mm)')
ax.set_xlabel('Time (hrs)')
ax.set_xticks(xrange(0, 25, 3))
ax.set_xlim([0, 24])
ax.set_title('Horizontally Averaged (Whole Domain) \nColumn Integrated Water Vapour Time Series')
plt.savefig('../CIWV.png')
plt.show()


