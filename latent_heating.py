import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate,integrate
from netCDF4 import Dataset
from SkewT_archer import Lv

Lv = 2.501e6
#hours = ["{0.02d}".format(h) for h in xrange(0, 24, 3)]
hours = ['{0:02d}'.format(h) for h in xrange(0, 24, 3)]
lhf_key = u'STASH_m01s03i222' # need for the surface latent heat flux
mcl_key = u'STASH_m01s00i392' # need for the total cloud liquid water
mr_key = u'STASH_m01s00i394' # need for the total suspended rain liquid

zth_key = 'thlev_zsea_theta' # theta levels
zrh_key = 'rholev_zsea_rho' # rho levels. ugh.

# Define a coordinate system
x = np.arange(0., 116000., 100.)
y = np.arange(0., 31900., 100.)
X, Y = np.meshgrid(x, y)

base_path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'

for hour in hours:
    print 'latent_heating.py -> starting hour: ' + hour
    # Open the netCDF
    mr_nc = Dataset(base_path + '/mr_' + hour + '.nc', 'r')
    fluxes_nc = Dataset(base_path + '/fluxes_' + hour + '.nc', 'r')
    
    # Read the data we care about
    mcl_data = mr_nc.variables[mcl_key][:]*1.
    mr_data  = mr_nc.variables[mr_key][:]*1.
    lhf_data = fluxes_nc.variables[lhf_key][:]*1.
    rho_data = fluxes_nc.variables[rho_key][:]*1.
    
    
    if hour == '00':
        z_theta = mr_nc.variables[zth_key][:]*1.
        z_rho   = fluxes_nc.variables[zrh_key][:]*1.
        
        q = np.array([interpolate.interp1d(x = z_rho, y = rho_data[it,:,:,:], axis = 0, fill_value = 'extrapolate')(z_theta) for it in rho_data.shape[0]])*(mcl_data + mr_data)
        Ldq_ts = Lv*np.nanmean(np.array([integrate.trapz(y = q[it,:,:,:], x = z_theta, axis = 0) for it in xrange(q.shape[0])]), axis = (1, 2))
        LHF_ts = Lv*np.nanmean(lhf_data[:,0,:,:], axis = (1, 2))
        times = mr_nc.variables['min10_0'][:]
    else:
        q = np.array([interpolate.interp1d(x = z_rho, y = rho_data[it,:,:,:], axis = 0, fill_value = 'extrapolate')(z_theta) for it in rho_data.shape[0]])*(mcl_data + mr_data)
        Ldq_ts = np.concatenate((Ldq_ts, Lv*np.nanmean(np.array([integrate.trapz(y = q[it,:,:,:], x = z_theta, axis = 0) for it in xrange(q.shape[0])]), axis = (1, 2))), axis = 0)
        LHF_ts = np.concatenate((LHF_ts, Lv*np.nanmean(lhf_data[:,0,:,:], axis = (1, 2))), axis = 0)
        times = np.concatenate((times, mr_nc.variables['min10_0'][:]), axis = 0)
    mr_nc.close()
    fluxes_nc.close()

times0 = 0.5*(times[2:] + times[1:-1])/60.
Ldq_ts0 = (Ldq_ts[2:] - Ldq_ts[1:-1])/600.
LHF_ts0 = 0.5*(LHF_ts[1:] + LHF_ts[:-1])

fig = plt.figure(tight_layout = True)
ax = fig.add_subplot(1, 1, 1)
ax.plot(times0, Ldq_ts0, 'g', lw = 2)
ax.plot(times0, LHF_ts0, 'b', lw = 2)
ax.set_xlabel('Time (hrs)')
ax.set_xlim([0, 24])
ax.set_xticks(xrange(0, 24, 3))
ax.set_ylabel('Latent Heating')

plt.show()




