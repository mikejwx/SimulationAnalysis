import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

hours = ["{0:02d}".format(hour) for hour in xrange(0, 24, 3)]
p_key = u'STASH_m01s00i408'
hovmoller = np.zeros((145, 1160))
for hour in hours:
    p_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control/fluxes_' + hour + '.nc', 'r')
    z = p_nc.variables['thlev_zsea_theta'][:]*1.
    my_z = 350.
    iz = np.where(np.abs(z-my_z) == np.min(np.abs(z-my_z)))[0][0]
    p_data = p_nc.variables[p_key][:,iz,:,:]*1.
    for i in xrange(p_data.shape[0]):
        hovmoller[i+18*hours.index(hour),:] = p_data[i,160,:] - np.mean(p_data[i,:,:])
    #fig = plt.figure()
    #ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
    #pp = ax.contourf(np.mean(np.transpose(np.transpose(p_data) - np.mean(p_data, axis = (1, 2))), axis = 0), cmap = 'bwr', levels = [i for i in np.arange(-10, 10.01, 0.25) if i != 0], extend = 'both')
    #fig.colorbar(pp, ax = ax)
    #plt.show()

plt.contourf(np.arange(0., 116000., 100.)/1000., np.arange(0, 1440.1, 10.)/60., hovmoller, cmap = 'bwr', levels = [i for i in np.arange(-10, 10.01, 0.25) if i != 0], extend = 'both')
plt.show()


