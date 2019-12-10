"""
Attempting to compute momentum flux profiles for the cloud layer in the moisture
experiments.
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import *

# read in the wind data
paths = {'BLm25' : '/nerc/n02/n02/xb899100/CloudTrail/RH_BLm25/',
         'FAm25' : '/nerc/n02/n02/xb899100/CloudTrail/RH_FAm25/',
         'Control_short' : '/nerc/n02/n02/xb899100/CloudTrail/Control_short/'}

my_data = {}
for path_key in paths.keys():
    my_data[path_key] = {}
    with Dataset(paths[path_key] + 'wind_04.nc', 'r') as wind_nc:
        my_data[path_key][u_key] = wind_nc.variables[u_key][:]*1.
        my_data[path_key][w_key] = wind_nc.variables[v_key][:]*1.
        my_data[path_key]['z']   = wind_nc.variables['thlev_zsea_theta'][:]*1.

for path_key in paths.keys():
    my_data[path_key]['wp'] = np.array([np.transpose(np.transpose(my_data[path_key][w_key][it,:,:,:]) - my_data[path_key][w_key].mean(axis = (0, 2, 3))) for it in range(my_data[path_key][w_key].shape[0])])
    my_data[path_key]['up'] = np.array([np.transpose(np.transpose(my_data[path_key][u_key][it,:,:,:]) - my_data[path_key][u_key].mean(axis = (0, 2, 3))) for it in range(my_data[path_key][u_key].shape[0])])
    my_data[path_key]['wpup'] = (my_data[path_key]['up']*my_data[path_key]['wp']).mean(axis = (0, 2, 3))

my_color = {'BLm25' : 'r', 'FAm25' : 'b', 'Control_short' : 'k'}

fig = plt.figure()
axa = fig.add_subplot(1, 1, 1)

[axa.plot(my_data[path_key]['wpup'], my_data[path_key]['z']/1000., my_color[path_key], label = path_key) for path_key in paths.keys()]

axa.set_xlabel('$(\overline{u^{\prime} w^{\prime}})$ (m$^{2}$ s$^{-2}$)')
axa.set_ylabel('Height (km)')
axa.set_ylim([0, 1])
axa.legend(loc = 0, frameon = 0)
plt.show()



