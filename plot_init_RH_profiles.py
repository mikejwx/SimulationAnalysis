import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import q_key, temp_key, pthe_key
from SkewT_archer import getQ

# create the variable storage infrastructure, using dictionaries
my_RH_data = {}

paths = ['/nerc/n02/n02/xb899100/CloudTrail/RH_BLm25/', '/nerc/n02/n02/xb899100/CloudTrail/RH_FAm25/', '/nerc/n02/n02/xb899100/CloudTrail/Control_short/']
for path in paths:
    path_key = path.split('/')[-2]
    my_RH_data[path_key] = {}
    with Dataset(path + 'bouy_00.nc', 'r') as bouy_nc:
        my_RH_data[path_key][q_key]    = bouy_nc.variables[q_key][1,:-1,:,:].mean(axis = (1,2))
        my_RH_data[path_key][temp_key] = bouy_nc.variables[temp_key][1,:-1,:,:].mean(axis = (1,2))
        z = bouy_nc.variables['thlev_zsea_theta'][:-1]*1.
    
    with Dataset(path + 'fluxes_00.nc', 'r') as fluxes_nc:
        my_RH_data[path_key][pthe_key] = fluxes_nc.variables[pthe_key][1,:-1,:,:].mean(axis = (1,2))
    
    my_RH_data[path_key]['q_sat'] = getQ(my_RH_data[path_key][temp_key], [100.], my_RH_data[path_key][pthe_key], t_units = 'K', p_units = 'Pa')
    my_RH_data[path_key]['RH']    = 100.*my_RH_data[path_key][q_key]/my_RH_data[path_key]['q_sat']
    

my_colors = {'RH_BLm25'      : 'red',
             'RH_FAm25'      : 'blue',
             'Control_short' : 'k'}
my_lw     = {'RH_BLm25'      : '2',
             'RH_FAm25'      : '2',
             'Control_short' : '1'}
my_labels = {'RH_BLm25'      : 'BLm25',
             'RH_FAm25'      : 'FAm25',
             'Control_short' : 'REF'}

fig = plt.figure()
axa = fig.add_subplot(1, 2, 1)
axa.set_ylabel('Height (km)')
axa.set_ylim([0, 4])
axa.set_xlabel('RH (%)')

axb = fig.add_subplot(1, 2, 2)
axb.set_ylim([0, 4])
axb.set_yticklabels([''])
axb.set_xlabel(u'q$_{v}$ (g kg$^{-1}$)')

for key in my_RH_data.keys():
    axa.plot(my_RH_data[key]['RH'], z/1000., color = my_colors[key], lw = my_lw[key], label = my_labels[key], ls = ['--' if 'm25' in key else '-'][0])
    axb.plot(my_RH_data[key][q_key]*1000., z/1000., color = my_colors[key], lw = my_lw[key], ls = ['--' if 'm25' in key else '-'][0])

axa.legend(loc = 0, frameon = 0)
axa.text(7.5, 3.7, 'a)')
axb.text(0.075*20., 3.7, 'b)', bbox={'fc':'w', 'ec':'None'})
plt.savefig('../RH_exp_IC.png', dpi = 150, bbox_inches = 'tight')
plt.show()

# plot the final conditions?
my_RH_data = {}

for path in paths:
    path_key = path.split('/')[-2]
    my_RH_data[path_key] = {}
    with Dataset(path + 'bouy_00.nc', 'r') as bouy_nc:
        time_key = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
        times = bouy_nc.variables[time_key][:]*1.
        it = np.where(np.abs(times - 120.) == np.min(np.abs(times - 120)))[0][0]
        my_RH_data[path_key][q_key]    = bouy_nc.variables[q_key][it,:-1,:,:].mean(axis = (1,2))
        my_RH_data[path_key][temp_key] = bouy_nc.variables[temp_key][it,:-1,:,:].mean(axis = (1,2))
        z = bouy_nc.variables['thlev_zsea_theta'][:-1]*1.
    
    with Dataset(path + 'fluxes_00.nc', 'r') as fluxes_nc:
        my_RH_data[path_key][pthe_key] = fluxes_nc.variables[pthe_key][it,:-1,:,:].mean(axis = (1,2))
    
    my_RH_data[path_key]['q_sat'] = getQ(my_RH_data[path_key][temp_key], [100.], my_RH_data[path_key][pthe_key], t_units = 'K', p_units = 'Pa')
    my_RH_data[path_key]['RH']    = 100.*my_RH_data[path_key][q_key]/my_RH_data[path_key]['q_sat']
    

fig = plt.figure()
axa = fig.add_subplot(1, 2, 1)
axa.set_ylabel('Height (km)')
axa.set_ylim([0, 4])
axa.set_xlabel('RH (%)')

axb = fig.add_subplot(1, 2, 2)
axb.set_ylim([0, 4])
axb.set_yticklabels([''])
axb.set_xlabel(u'q$_{v}$ (g kg$^{-1}$)')

for key in my_RH_data.keys():
    axa.plot(my_RH_data[key]['RH'], z/1000., color = my_colors[key], lw = my_lw[key], label = my_labels[key], ls = ['--' if 'm25' in key else '-'][0])
    axb.plot(my_RH_data[key][q_key]*1000., z/1000., color = my_colors[key], lw = my_lw[key], ls = ['--' if 'm25' in key else '-'][0])

axa.legend(loc = 0, frameon = 0)
axa.text(7.5, 3.7, 'a)')
axb.text(0.075*20., 3.7, 'b)', bbox={'fc':'w', 'ec':'None'})
plt.savefig('../RH_exp_FC.png', dpi = 150, bbox_inches = 'tight')
plt.show()

