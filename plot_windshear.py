import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import u_key, v_key, zi_new_key

### Wind Speed experiments ###
paths = {'U10' : '/nerc/n02/n02/xb899100/CloudTrail/Control/',
         'U05' : '/nerc/n02/n02/xb899100/CloudTrail/U05/',
         'U06' : '/nerc/n02/n02/xb899100/CloudTrail/U06/',
         'U07' : '/nerc/n02/n02/xb899100/CloudTrail/U07/',
         'U08' : '/nerc/n02/n02/xb899100/CloudTrail/U08/',
         'U09' : '/nerc/n02/n02/xb899100/CloudTrail/U09/'}

my_keys = ['U05', 'U06', 'U07', 'U08', 'U09', 'U10']
winds = {}
# read the initial wind profiles for each simulation
for key in my_keys:
    winds[key] = []
    wind_nc = Dataset(paths[key] + 'wind_00.nc', 'r')
    zi_nc = Dataset(paths[key] + 'zi_00.nc', 'r')
    u_data = wind_nc.variables[u_key][0,:,:,:].mean(axis = (1, 2))
    v_data = wind_nc.variables[v_key][0,:,:,:].mean(axis = (1, 2))
    zi_data = zi_nc.variables[zi_new_key][0,:,:].mean()
    z = wind_nc.variables['thlev_zsea_theta'][:]*1.
    
    iz = np.where(np.abs(z - zi_data) == np.min(np.abs(z - zi_data)))[0][0]
    winds[key].append(-float(key[1:]) - u_data[:iz].mean())
    winds[key].append(0 - v_data[:iz].mean())
    
    wind_nc.close()
    zi_nc.close()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
ax.set_xlim([-2, 1])
ax.set_ylim([-1, 2])
for key in my_keys:
    # plot an arrow for the shear vector
    ax.arrow(0, 0, winds[key][0], winds[key][1], color = 'k', head_width = 0.05)
    ax.text(winds[key][0], winds[key][1], key, va = 'top', ha = 'right')

ax.plot([-2, 2],[0, 0], color = 'grey', ls = ':')
ax.plot([0, 0], [-2, 2], color = 'grey', ls = ':')
ax.set_xlabel('u (m s$^{-1}$)')
ax.set_ylabel('v (m s$^{-1}$)')
ax.set_title('Shear Vectors')
plt.savefig('../shear_vectors.png', dpi = 150, bbox_inches = 'tight')
plt.show()

