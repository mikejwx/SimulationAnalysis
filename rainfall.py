import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

hours = ["{0:0d}".format(h) for h in xrange(0, 24, 3)]

rain_key = u'STASH_m01s00i394'
for hour in hours:
    mr_nc = Dataset('mr_'+ hour + '.nc', 'r')
    if hour == '00':
        rain_accum = np.sum(mr_nc.variables[rain_key][1:,0,:,:], axis = 0)
    else:
        rain_accum += np.sum(mr_nc.variables[rain_key][:,0,:,:], axis = 0)

plt.contourf(rain_accum)
plt.show()
