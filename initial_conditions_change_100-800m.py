import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import theta_key, q_key, u_key, v_key
from analysis_tools import fromComponents, bilinear_interpolation, get_cs_coords, lcl as get_lcl, downwind_rectangle
from scipy import integrate, interpolate
from datetime import datetime as dt
from STASH_keys import rho_key, u_key, v_key, w_key, pthe_key, q_key, theta_key
from SkewT_archer import PTtoTemp

"""
Answering the question:
Can we get away with using the 100 m initial conditions and fluxes at 800 m?
"""
hours = ["{0:02d}".format(hour) for hour in range(0, 24, 3)]
print 'The 100 m simulation'
path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'
sim100 = {}

for hour in hours:
    # Open the nc files
    bouy_nc = Dataset(path + 'bouy_' + hour + '.nc', 'r')
    u_nc = Dataset(path + 'u_' + hour + '.nc', 'r')
    v_nc = Dataset(path + 'v_' + hour + '.nc', 'r')
    
    if hour == hours[0]:
        # define the coordinate system
        z_theta = bouy_nc.variables['thlev_zsea_theta'][:]*1.
        z_rho   = u_nc.variables['rholev_zsea_rho'][:]*1.
        X, Y = np.meshgrid(np.arange(bouy_nc.variables[theta_key].shape[3])*100.0, np.arange(bouy_nc.variables[theta_key].shape[2])*100.0)
        time_key = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
        sim100['times'] = bouy_nc.variables[time_key]
        
        # define the CT region, and mask over it
        R_i = 1000.0*(50.0/np.pi)**0.5 # island radius, metres
        # centrepoint coordinates of the island
        x_c = 100000.0 + R_i
        y_c = 4*R_i
        iz = np.where(np.abs(z_rho - 650.0) == np.min(np.abs(z_rho - 650.0)))[0][0]
        speed, wind_dir = fromComponents(u = np.nanmean(u_nc.variables[u_key][0,:iz,:,:]), v = np.nanmean(v_nc.variables[v_key][0,:iz,:,:]))
        mask, y_prime, x_prime = downwind_rectangle(wind_dir = wind_dir, x_c = x_c, y_c = y_c, X = X, Y = Y, island_radius = R_i, dist_0 = 0, dist_1 = 100000.0, half_width = 7500.0)
        
        # mask over it and store to the dictionary
        sim100['theta'] = np.nanmean(mask*bouy_nc.variables[theta_key][1:,:,:,:], axis = (2, 3))
        sim100['q']     = np.nanmean(mask*bouy_nc.variables[q_key][1:,:,:,:], axis = (2, 3))
        sim100['u']     = np.nanmean(mask*u_nc.variables[u_key][1:,:,:,:], axis = (2, 3))
        sim100['v']     = np.nanmean(mask*v_nc.variables[v_key][1:,:,:,:], axis = (2, 3))
    else:
        sim100['theta'] = np.concatenate((sim100['theta'], np.nanmean(mask*bouy_nc.variables[theta_key][:], axis = (2, 3))), axis = 0)
        sim100['q']     = np.concatenate((sim100['q'], np.nanmean(mask*bouy_nc.variables[q_key][:], axis = (2, 3))), axis = 0)
        sim100['u']     = np.concatenate((sim100['u'], np.nanmean(mask*u_nc.variables[u_key][:], axis = (2, 3))), axis = 0)
        sim100['v']     = np.concatenate((sim100['v'], np.nanmean(mask*v_nc.variables[v_key][:], axis = (2, 3))), axis = 0)
        sim100['times'] = np.concatenate((sim100['times'], bouy_nc.variables[time_key][:]), axis = 0)
    bouy_nc.close()
    u_nc.close()
    v_nc.close()

### Plot the change due to the cloud trail in 100m simulation
print 'plotting'
fig = plt.figure(figsize = (18, 6))
axa = fig.add_subplot(1, 3, 1)
axb = fig.add_subplot(1, 3, 2)
axb.set_yticklabels([''])
axc = fig.add_subplot(1, 3, 3)
axc.set_yticklabels([''])

# figure out which indexes correspond to a relevant cloud trail period
dt = 100000.0/speed # time taken for the first island heating to reach the end of the rectangle, seconds
time0 = 360.0 # time island is first heated, minutes
time1 = time0 + dt/60.0 # start looking from this time onward.
# island is no longer heated at 18:00
time2 = time1 + 3*60.0
idx = [i for i in range(len(sim100['times'])) if (time1 <= sim100['times'][i]) and (sim100['times'][i] <= time2)]
axa.plot(sim100['theta'][0,:], z_theta, 'k--')
axa.plot(np.nanmean(sim100['theta'][idx,:], axis = 0), z_theta, 'k', lw = 2)
axa.set_xlabel(u'$\\theta$ (K)')

axb.plot(sim100['q'][0,:]*1000., z_theta, 'k--', label = 'initial')
axb.plot(np.nanmean(sim100['q'][idx,:], axis = 0)*1000., z_theta, 'k', lw = 2, label = 'CT region/times')
axb.set_xlabel(u'q$_{v}$ (g kg$^{-1}$)')
axb.legend(loc = 0)
axb.set_title('100 m in CT region')

axc.plot(sim100['u'][0,:], z_rho, 'b--')
axc.plot(np.nanmean(sim100['u'][idx,:], axis = 0), z_rho, 'b', lw = 2)
axc.plot(sim100['v'][0,:], z_rho, 'r--')
axc.plot(np.nanmean(sim100['v'][idx,:], axis = 0), z_rho, 'r', lw = 2)

plt.savefig('../initial_conditions_change_100m_CTbox.png', dpi = 150, bbox_inches = 'tight')
plt.show()

path = '/work/n02/n02/xb899100/cylc-run/u-bn261/share/data/history/'
sim800 = {}

for hour in hours:
    # Open the nc files
    bouy_nc = Dataset(path + 'bouy_' + hour + '.nc', 'r')
    u_nc = Dataset(path + 'u_' + hour + '.nc', 'r')
    v_nc = Dataset(path + 'v_' + hour + '.nc', 'r')
    
    if hour == hours[0]:
        # define the coordinate system
        z_theta = bouy_nc.variables['thlev_zsea_theta'][:]*1.
        z_rho   = u_nc.variables['rholev_zsea_rho'][:]*1.
        X, Y = np.meshgrid(np.arange(bouy_nc.variables[theta_key].shape[3])*800.0, np.arange(bouy_nc.variables[theta_key].shape[2])*800.0)
        time_key = [tkey for tkey in bouy_nc.variables.keys() if 'min' in tkey][0]
        sim800['times'] = bouy_nc.variables[time_key]
        # define the CT region, and mask over it
        R_i = 1000.0*(50.0/np.pi)**0.5 # island radius, metres
        # centrepoint coordinates of the island
        x_c = 100000.0 + 2*R_i
        y_c = 4*R_i
        iz = np.where(np.abs(z_rho - 650.0) == np.min(np.abs(z_rho - 650.0)))[0][0]
        speed, wind_dir = fromComponents(u = np.nanmean(u_nc.variables[u_key][0,:iz,:,:]), v = np.nanmean(v_nc.variables[v_key][0,:iz,:,:]))
        mask, y_prime, x_prime = downwind_rectangle(wind_dir = wind_dir, x_c = x_c, y_c = y_c, X = X, Y = Y, island_radius = R_i, dist_0 = 0, dist_1 = 100000.0, half_width = 7500.0)
        
        # mask over it and store to the dictionary
        sim800['theta'] = np.nanmean(mask*bouy_nc.variables[theta_key][1:,:,:,:], axis = (2, 3))
        sim800['q']     = np.nanmean(mask*bouy_nc.variables[q_key][1:,:,:,:], axis = (2, 3))
        sim800['u']     = np.nanmean(mask*u_nc.variables[u_key][1:,:,:,:], axis = (2, 3))
        sim800['v']     = np.nanmean(mask*v_nc.variables[v_key][1:,:,:,:], axis = (2, 3))
    else:
        sim800['theta'] = np.concatenate((sim800['theta'], np.nanmean(mask*bouy_nc.variables[theta_key][:], axis = (2, 3))), axis = 0)
        sim800['q']     = np.concatenate((sim800['q'], np.nanmean(mask*bouy_nc.variables[q_key][:], axis = (2, 3))), axis = 0)
        sim800['u']     = np.concatenate((sim800['u'], np.nanmean(mask*u_nc.variables[u_key][:], axis = (2, 3))), axis = 0)
        sim800['v']     = np.concatenate((sim800['v'], np.nanmean(mask*v_nc.variables[v_key][:], axis = (2, 3))), axis = 0)
        sim800['times'] = np.concatenate((sim800['times'], bouy_nc.variables[time_key][:]), axis = 0)
    bouy_nc.close()
    u_nc.close()
    v_nc.close()

### Plot the change due to the cloud trail in 800m simulation
print 'plotting'
fig = plt.figure(figsize = (18, 6))
axa = fig.add_subplot(1, 3, 1)
axb = fig.add_subplot(1, 3, 2)
axb.set_yticklabels([''])
axc = fig.add_subplot(1, 3, 3)
axc.set_yticklabels([''])
dt = 100000.0/speed # time taken for the first island heating to reach the end of the rectangle, seconds
time0 = 360.0 # time island is first heated, minutes
time1 = time0 + dt/60.0 # start looking from this time onward.
# island is no longer heated at 18:00
time2 = 18*60.0
idx = [i for i in range(len(sim800['times'])) if (time1 <= sim800['times'][i]) and (sim800['times'][i] <= time2)]
#idx = np.where(np.abs(sim800['times'] - np.nanmean(sim800['times'][idx])) == np.min(np.abs(sim800['times'] - np.nanmean(sim800['times'][idx]))))[0][0]
#idx = [idx-1, idx, idx+1]
axa.plot(sim800['theta'][0,:], z_theta, 'k--')
axa.plot(np.nanmean(sim800['theta'][idx,:], axis = 0), z_theta, 'k', lw = 2)
axa.set_xlabel(u'$\\theta$ (K)')

axb.plot(sim800['q'][0,:]*1000., z_theta, 'k--', label = 'initial')
axb.plot(np.nanmean(sim800['q'][idx,:], axis = 0)*1000., z_theta, 'k', lw = 2, label = 'CT region/times')
axb.set_xlabel(u'q$_{v}$ (g kg$^{-1}$)')
axb.legend(loc = 0)
axb.set_title('800 m')

axc.plot(sim800['u'][0,:], z_rho, 'b--')
axc.plot(np.nanmean(sim800['u'][idx,:], axis = 0), z_rho, 'b', lw = 2)
axc.plot(sim800['v'][0,:], z_rho, 'r--')
axc.plot(np.nanmean(sim800['v'][idx,:], axis = 0), z_rho, 'r', lw = 2)

plt.savefig('../initial_conditions_chang e_800m.png', dpi = 150, bbox_inches = 'tight')
plt.show()



