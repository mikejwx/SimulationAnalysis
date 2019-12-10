import numpy as np
import matplotlib.pyplot as plt
from STASH_keys import theta_key, q_key, u_key, v_key
from netCDF4 import Dataset
from scipy import ndimage

# Define path to data
path = '/nerc/n02/n02/xb899100/CloudTrail/U05/'
hours = ["{:02d}".format(h) for h in range(0, 16, 4)]
for hour in hours:
    with Dataset(path + 'bouy_' + hour + '.nc', 'r') as bouy_nc:
        if hour == hours[0]:
            theta_data = bouy_nc.variables[theta_key][1:,0,:,:]*1.
            q_data = bouy_nc.variables[q_key][1:,0,:,:]*1.
            time_key = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
            times = bouy_nc.variables[time_key][1:]*1.
        else:
            theta_data = np.concatenate((theta_data, bouy_nc.variables[theta_key][:,0,:,:]*1.), axis = 0)
            q_data = np.concatenate((q_data, bouy_nc.variables[q_key][:,0,:,:]*1.), axis = 0)
            times = np.concatenate((times, bouy_nc.variables[time_key][:]*1.), axis = 0)
    
    with Dataset(path + 'wind_' + hour + '.nc', 'r') as wind_nc:
        if hour == hours[0]:
            wind_data = np.sqrt(wind_nc.variables[u_key][1:,0,:,:]**2 + wind_nc.variables[v_key][1:,0,:,:]**2)
        else:
            wind_data = np.concatenate((wind_data, np.sqrt(wind_nc.variables[u_key][:,0,:,:]**2 + wind_nc.variables[v_key][:,0,:,:]**2)), axis = 0)

thetav = theta_data*(1. + 0.608*q_data)
thetav_anom = np.array([thetav[idx,:,:] - np.nanmean(thetav[idx,:,:]) for idx in range(thetav.shape[0])])

g = 9.81
c = np.sqrt(g*np.abs(thetav_anom)*500./thetav)

# For each time, identify all of the cold pools, and find the maximum height in each cloud
cp_mask = np.where(thetav_anom <= 0, 1.0, 0.0)
mean_c = []
n_cp = []
for idx_t in [48]:#range(len(times)):
    # Use scipy.ndimage to count and label the cold pools
    coldpools, n_cps = ndimage.label(cp_mask[idx_t,:,:])
    print 'Time = ' + str(int(times[idx_t])) + ', # cp = ' + str(n_cps)
    cp_c = np.array([np.nanmean(np.where(coldpools == coldpool, c[idx_t,:,:], np.nan)) for coldpool in range(n_cps)])
    mean_c.append(cp_c)
    n_cp.append(n_cps)


