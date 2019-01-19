import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

hours = ["{0:02d}".format(h) for h in xrange(0, 24, 3)]

rain_key = u'STASH_m01s00i394'
rho_key = u'STASH_m01s00i389'
for hour in hours:
    mr_nc = Dataset('../mr_'+ hour + '.nc', 'r')
    fluxes_nc = Dataset('../fluxes_' + hour + '.nc', 'r')
    z = mr_nc.variables['thlev_zsea_theta'][:]*1.
    # to convert surface rain water mixing ratio to accumulated rain depth:
    # (rainwater mixing ratio)*(dry air density)*(depth of grid box)*(600 second time step)*(1000 millimeters in a meter)/(1000 kg of water per unit volume)
    # Then sum that quantity over the simulation duration
    # N.B. z[0] 0m
    if hour == '00':
        rain_accum = np.sum(z[1]*fluxes_nc.variables[rho_key][1:,0,:,:]*mr_nc.variables[rain_key][1:,0,:,:]*600., axis = 0)
    else:
        rain_accum += np.sum(z[1]*fluxes_nc.variables[rho_key][:,0,:,:]*mr_nc.variables[rain_key][:,0,:,:]*600., axis = 0)
    mr_nc.close()
    fluxes_nc.close()

plt.contourf(rain_accum)
plt.colorbar()
plt.show()


