import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import mcl_key, mr_key, w_key, theta_key, q_key, temp_key, pthe_key
from SkewT_archer import getQ
from scipy import ndimage
from datetime import datetime as dt

"""
Updraughts are where there is liquid water present and vertical velocty is positive
Cores are where there is liquid ater, positive vertical velocity and positive buoyancy
"""

### Read the data ###
path_start = '/nerc/n02/n02/xb899100/CloudTrail/Control_'
path_end   = '_HRIC_INV/'
experiments = ['0200m', '0400m', '0800m', '1600m']
exp_labels = {'0200m':'DX0200', '0400m':'DX0400', '0800m':'DX0800','1600m':'DX1600'}
hours = ["{0:02d}".format(hour) for hour in range(0, 24, 3)]
exp_cols = {'0200m':'brown', '0400m':'red', '0800m' : 'orange', '1600m':'gold'}
my_data = {}
for exp in experiments:
    print '[' + dt.now().strftime('%H:%M:%S') + '] Starting experiment: ' + exp
    my_data[exp] = {}
    dx = float(exp[:-1])/1000.
    print '[' + dt.now().strftime('%H:%M:%S') + '] Reading the cloud liquid water content data'
    for hour in hours:
        print 'hour = ' + hour
        with Dataset(path_start + exp + path_end + 'mr_' + hour + '.nc', 'r') as mr_nc:
            if hour == hours[0]:
                my_data[exp][mcl_key] = mr_nc.variables[mcl_key][1:,:,:,:]*1.
                my_data[exp][mr_key]  = mr_nc.variables[mr_key][1:,:,:,:]*1.
                z = mr_nc.variables['thlev_zsea_theta'][:]*1.
                mcl_time_key = [tkey for tkey in mr_nc.variables.keys() if 'min' in tkey][0]
                my_data[exp]['mcl_times'] = mr_nc.variables[mcl_time_key][1:]*1.
            else:
                my_data[exp][mcl_key] = np.concatenate((my_data[exp][mcl_key], mr_nc.variables[mcl_key][:]*1.), axis = 0)
                my_data[exp][mr_key]  = np.concatenate((my_data[exp][mr_key], mr_nc.variables[mr_key][:]*1.), axis = 0)
                my_data[exp]['mcl_times'] = np.concatenate((my_data[exp]['mcl_times'], mr_nc[mcl_time_key][:]*1.), axis = 0)
    
        with Dataset(path_start + exp + path_end + 'bouy_' + hour + '.nc', 'r') as bouy_nc:
            if hour == hours[0]:
                my_data[exp][theta_key] = bouy_nc.variables[theta_key][1:,:,:,:]*1.
                my_data[exp][w_key]     = bouy_nc.variables[w_key][1:,:,:,:]*1.
                my_data[exp][q_key]     = bouy_nc.variables[q_key][1:,:,:,:]*1.
                my_data[exp][temp_key]  = bouy_nc.variables[temp_key][1:,:,:,:]*1.
            else:
                my_data[exp][theta_key] = np.concatenate((my_data[exp][theta_key], bouy_nc.variables[theta_key][:]*1.), axis = 0)
                my_data[exp][w_key]     = np.concatenate((my_data[exp][w_key], bouy_nc.variables[w_key][:]*1.), axis = 0)
                my_data[exp][q_key]     = np.concatenate((my_data[exp][q_key], bouy_nc.variables[q_key][:]*1.), axis = 0)
                my_data[exp][temp_key]  = np.concatenate((my_data[exp][temp_key], bouy_nc.variables[temp_key][:]*1.), axis = 0)
        
        with Dataset(path_start + exp + path_end + 'fluxes_' + hour + '.nc', 'r') as fluxes_nc:
            if hour == hours[0]:
                my_data[exp][pthe_key]  = fluxes_nc.variables[pthe_key][1:,:,:,:]*1.
            else:
                my_data[exp][pthe_key]  = np.concatenate((my_data[exp][pthe_key], fluxes_nc.variables[pthe_key][:]*1.), axis = 0)
    
    # define theta_v
    my_data[exp]['theta_v'] = my_data[exp][theta_key]*(1. + 0.608*my_data[exp][q_key])
    my_data[exp]['theta_v_prime'] = np.array([np.transpose(np.transpose(my_data[exp]['theta_v'][it,:,:,:]) - my_data[exp]['theta_v'][it,:,:,:].mean(axis = (1, 2))) for it in range(my_data[exp]['mcl_times'].size)])
    my_data[exp]['q_total'] = my_data[exp][q_key] + my_data[exp][mcl_key] + my_data[exp][mr_key]
    my_data[exp]['q_total_prime'] = np.array([np.transpose(np.transpose(my_data[exp]['q_total'][it,:,:,:]) - my_data[exp]['q_total'][it,:,:,:].mean(axis = (1, 2))) for it in range(my_data[exp]['mcl_times'].size)])
    # define RH
    my_data[exp]['RH'] = 100.*my_data[exp][q_key]/getQ(my_data[exp][temp_key][:]*1., [100.], my_data[exp][pthe_key][:]*1., t_units = 'K', p_units = 'Pa')
    # define the updraught mask
    my_data[exp]['up_mask'] = np.where((my_data[exp][mcl_key][:] > 0)*(my_data[exp][w_key] > 0), 1.0, np.nan)
    # define the core mask
    my_data[exp]['core_mask'] = np.where((my_data[exp][mcl_key][:] > 0)*(my_data[exp][w_key] > 0)*(my_data[exp]['theta_v_prime'] > 0), 1.0, np.nan)

# calculate the profiles first
for exp in experiments:
    t_idx = [it for it in range(my_data[exp]['mcl_times'].size) if (360. <= my_data[exp]['mcl_times'][it])*(my_data[exp]['mcl_times'][it] <= 1080.)]
    my_data[exp]['theta_v_prime_u'] = np.nanmean((my_data[exp]['theta_v_prime']*my_data[exp]['up_mask'])[t_idx,:,:,:], axis = (0, 2, 3))
    my_data[exp]['theta_v_prime_c'] = np.nanmean((my_data[exp]['theta_v_prime']*my_data[exp]['core_mask'])[t_idx,:,:,:], axis = (0, 2, 3))
    my_data[exp]['q_total_prime_u'] = 1000.*np.nanmean((my_data[exp]['q_total_prime']*my_data[exp]['up_mask'])[t_idx,:-1,:,:], axis = (0, 2, 3))
    my_data[exp]['q_total_prime_c'] = 1000.*np.nanmean((my_data[exp]['q_total_prime']*my_data[exp]['core_mask'])[t_idx,:-1,:,:], axis = (0, 2, 3))

### Make the plot ###
fig = plt.figure(figsize = (10, 6), tight_layout = True)
for exp in experiments:
    # Theta_v anomaly
    axa = fig.add_subplot(1, 2, 1)
    axa.plot(my_data[exp]['theta_v_prime_u'], z/1000., color = exp_cols[exp], lw = 2, ls = '--')
    axa.plot(my_data[exp]['theta_v_prime_c'], z/1000., color = exp_cols[exp], lw = 2, label = exp_labels[exp])
    axa.plot([0,0], [0, 3.5], color = 'grey', ls = ':')
    axa.set_ylim([0, 3.5])
    axa.set_xlim([-1,1])
    axa.set_ylabel('Height (km)')
    axa.text(2*0.1-1, 3.5*0.9, 'a)')
    axa.annotate('updraught', (-0.5, 2.5), (-0.5, 1.8), va = 'center', ha = 'center', arrowprops = {'color' : 'k', 'width':1})
    axa.annotate('core', (0.5, 0.95), (0.5, 0.5), va = 'center', ha = 'center', arrowprops = {'color' : 'k', 'width':1})
    # q_total anomaly
    axb = fig.add_subplot(1, 2, 2)
    axb.plot(my_data[exp]['q_total_prime_u'], z[:-1]/1000., color = exp_cols[exp], lw = 2, ls = '--')
    axb.plot(my_data[exp]['q_total_prime_c'], z[:-1]/1000., color = exp_cols[exp], lw = 2)
    axb.set_ylim ([0, 3.5])
    axb.set_xlim([0, 15])
    axb.set_yticklabels([''])
    axb.text(15*0.1, 3.5*0.9, 'b)')
    axb.annotate('updraught', (8.9, 3.0), (4.5, 2.7), va = 'center', ha = 'center', arrowprops = {'color' : 'k', 'width':1})
    axb.annotate('core', (5.3, 1.25), (8, 1), va = 'center', ha = 'center', arrowprops = {'color' : 'k', 'width':1})

axa.legend(loc = 0, frameon = 0, fontsize = 12)
axa.set_xlabel('$<\\theta_{v}^{u,c}> - \\theta_{v}^{e}$ (K)')
axb.set_xlabel('$<q_{T}^{u,c}> - q_{T}^{e}$ (g kg$^{-1}$)')
plt.savefig('../Ch6_Figure13.png', dpi = 250, bbox_inches = 'tight')
plt.show()

#plt.plot(my_data['0200m']['RH'][t_idx,:,:,:].mean(axis = (0, 2, 3)), z)
#plt.plot(my_data['1600m']['RH'][t_idx,:,:,:].mean(axis = (0, 2, 3)), z)
#plt.xlim([0, 100])
#plt.show()



