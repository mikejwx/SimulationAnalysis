"""
Attempting to compute momentum flux profiles for the cloud layer in the moisture
experiments.
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import *
from scipy import interpolate

# read in the wind data
paths = {'BLm25' : '/nerc/n02/n02/xb899100/CloudTrail/RH_BLm25/',
         'FAm25' : '/nerc/n02/n02/xb899100/CloudTrail/RH_FAm25/',
         'Control_short' : '/nerc/n02/n02/xb899100/CloudTrail/Control_short/'}

my_data = {}
time_keys = {}
for path_key in paths.keys():
    my_data[path_key] = {}
    with Dataset(paths[path_key] + 'u_00.nc', 'r') as u_nc:
        my_data[path_key][u_key] = u_nc.variables[u_key][:]*1.
        my_data[path_key]['z_rho'] = u_nc.variables['rholev_zsea_rho'][:]*1.
        time_keys[path_key] = [tkey for tkey in u_nc.variables.keys() if 'min' in tkey][0]
        my_data[path_key][time_keys[path_key]] = u_nc.variables[time_keys[path_key]][:]*1.
    
    with Dataset(paths[path_key] + 'bouy_00.nc', 'r') as w_nc:
        my_data[path_key][w_key] = w_nc.variables[w_key][:]*1.
        my_data[path_key]['z']   = w_nc.variables['thlev_zsea_theta'][:]*1.

for path_key in paths.keys():
    t_idx = [it for it in range(my_data[path_key][time_keys[path_key]].size) if my_data[path_key][time_keys[path_key]][it] <= 120.0] 
    my_data[path_key]['wp'] = np.array([np.transpose(np.transpose(my_data[path_key][w_key][it,:,:,:]) - my_data[path_key][w_key].mean(axis = (0, 2, 3))) for it in t_idx])
    my_data[path_key]['up_rho'] = np.array([np.transpose(np.transpose(my_data[path_key][u_key][it,:,:,:]) - my_data[path_key][u_key].mean(axis = (0, 2, 3))) for it in t_idx])
    my_data[path_key]['up'] = np.array([interpolate.interp1d(x = my_data[path_key]['z_rho'], y = my_data[path_key]['up_rho'][it,:,:,:], axis = 0, fill_value = 'extrapolate')(my_data[path_key]['z']) for it in t_idx])
    my_data[path_key]['wpup'] = (my_data[path_key]['up']*my_data[path_key]['wp']).mean(axis = (0, 2, 3))

my_color = {'BLm25' : 'r', 'FAm25' : 'b', 'Control_short' : 'k'}

fig = plt.figure()
axa = fig.add_subplot(1, 1, 1)

[axa.plot(my_data[path_key]['wpup'], my_data[path_key]['z']/1000., my_color[path_key], label = path_key) for path_key in paths.keys()]

axa.set_xlabel('$(\overline{w^{\prime} u^{\prime}})$ (m$^{2}$ s$^{-2}$)')
axa.set_ylabel('Height (km)')
axa.set_ylim([0, 1])
axa.legend(loc = 0, frameon = 0)
plt.savefig('../momentum_fluxes_by_experiment.png', dpi = 150, bbox_inches = 'tight')
plt.show()



