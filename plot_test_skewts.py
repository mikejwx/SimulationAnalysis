import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import theta_key, q_key, prho_key, pthe_key, lwp_key, u_key, v_key, w_key
from SkewT_archer import *
from scipy import interpolate

hours = ["{0:02d}".format(h) for h in range(0,13,3)]
path = '/work/n02/n02/xb899100/cylc-run/u-bn540/share/data/history/'
data = {}

for hour in hours:
    bouy_nc = Dataset(path + 'bouy_' + hour + '.nc', 'r')
    mr_nc = Dataset(path + 'mr_' + hour + '.nc', 'r')
    fluxes_nc = Dataset(path + 'fluxes_' + hour + '.nc', 'r')
    u_nc = Dataset(path + 'u_' + hour + '.nc', 'r')
    v_nc = Dataset(path + 'v_' + hour + '.nc', 'r')
    if hour == hours[0]:
        data[theta_key] = bouy_nc.variables[theta_key][1:,:,:,:]*1.
        data[q_key] = bouy_nc.variables[q_key][1:,:,:,:]*1.
        data[prho_key] = fluxes_nc.variables[prho_key][1:,:,:,:]*1.
        data[u_key] = u_nc.variables[u_key][1:,:,:,:]*1.
        data[v_key] = v_nc.variables[v_key][1:,:,:,:]*1.
        data[w_key] = bouy_nc.variables[w_key][1:,:,:,:]*1.
        data['z'] = bouy_nc.variables['thlev_zsea_theta'][:]*1.
        z_rho = fluxes_nc.variables['rholev_zsea_rho'][:]*1.
    else:
        data[theta_key] = np.concatenate((data[theta_key], bouy_nc.variables[theta_key][:]*1.), axis = 0)
        data[q_key] = np.concatenate((data[q_key], bouy_nc.variables[q_key][:]*1.), axis = 0)
        data[prho_key] = np.concatenate((data[prho_key], fluxes_nc.variables[prho_key][:]*1.), axis = 0)
        data[u_key] = np.concatenate((data[u_key], u_nc.variables[u_key][:]*1.), axis = 0)
        data[v_key] = np.concatenate((data[v_key], v_nc.variables[v_key][:]*1.), axis = 0)
        data[w_key] = np.concatenate((data[w_key], bouy_nc.variables[w_key][:]*1.), axis = 0)
    
    bouy_nc.close()
    mr_nc.close()
    fluxes_nc.close()
    u_nc.close()
    v_nc.close()

data[pthe_key] = np.array([interpolate.interp1d(z_rho, data[prho_key][idx,:,:,:], axis = 0, fill_value = 'extrapolate')(data['z']) for idx in range(data[prho_key].shape[0])])

data['temp'] = PTtoTemp(data[theta_key], data[pthe_key], t_units = 'K', p_units = 'Pa')
data['tdew'] = getDew(data[q_key], data[pthe_key], q_units = 'kg/kg', p_units ='Pa')

"""
# make a skew T at peak island heating
import os
os.system('mkdir ./scratch')
times = np.arange(10, data[w_key].shape[0]*10+1, 10)*1.
for time in times:
    idx_t = np.where(np.abs(times - time) == np.min(np.abs(times - time)))[0][0]
    plotSkewT(np.nanmean(data['temp'][idx_t,:,:,:], axis = (1,2))-273.15, np.nanmean(data['tdew'][idx_t,:,:,:], axis = (1,2))-273.15, np.nanmean(data[pthe_key][idx_t,:,:,:], axis = (1,2))/100., u = np.array([-999]), v = np.array([-999]), CAPE = True, my_title = 'T+' + str(int(time)) + 'mins')
    plt.savefig('./scratch/skewt_' + "{0:04d}".format(int(time))+'.png', dpi = 150)
    plt.close('all')

os.system('convert -delay 30 -loop 0 ./scratch/*.png ./animation.gif')
os.system('rm -rf ./scratch')
"""
data['RH'] = 100.*data[q_key]/getQ(data['temp'], [100.], data[pthe_key], t_units = 'K', p_units = 'Pa')

plt.plot(np.nanmean(data['RH'][0,:,:,:], axis = (1, 2)), np.nanmean(data[pthe_key][0,:,:,:], axis = (1, 2)))
plt.plot(np.nanmean(data['RH'][-1,:,:,:], axis = (1, 2)), np.nanmean(data[pthe_key][-1,:,:,:], axis = (1, 2)))
plt.ylim([100000., 10000.])
plt.show()

# look at the zonal winds at the start and at hour 14 (last output before crash, roughly 45 minutes before crash)
times = np.arange(10, data[w_key].shape[0]*10+1, 10)*1.
time = 130.
it = np.where(np.abs(times - time) == np.min(np.abs(times - time)))[0][0]

plt.plot(np.nanpercentile(data[u_key][0,:,:,:], 90, axis = (1, 2)), z_rho, 'k--')
plt.fill_betweenx(z_rho, np.nanpercentile(data[u_key][0,:,:,:], 25, axis = (1, 2)), np.nanpercentile(data[u_key][0,:,:,:], 75, axis = (1, 2)), facecolor = 'k', alpha = 0.5, edgecolor = 'None')
plt.plot(np.nanpercentile(data[u_key][0,:,:,:], 50, axis = (1, 2)), z_rho, 'k', lw = 2)
plt.plot(np.nanpercentile(data[u_key][0,:,:,:], 10, axis = (1, 2)), z_rho, 'k--')
plt.show()

plt.plot(np.nanpercentile(data[u_key][it,:,:,:], 90, axis = (1, 2)), z_rho, 'k--')
plt.fill_betweenx(z_rho, np.nanpercentile(data[u_key][it,:,:,:], 25, axis = (1, 2)), np.nanpercentile(data[u_key][it,:,:,:], 75, axis = (1, 2)), facecolor = 'k', alpha = 0.5, edgecolor = 'None')
plt.plot(np.nanpercentile(data[u_key][it,:,:,:], 50, axis = (1, 2)), z_rho, 'k', lw = 2)
plt.plot(np.nanpercentile(data[u_key][it,:,:,:], 10, axis = (1, 2)), z_rho, 'k--')
plt.show()

# look at the meridional winds at the start and at hour 14 (last output before crash, roughly 45 minutes before crash)
plt.plot(np.nanpercentile(data[v_key][0,:,:,:], 90, axis = (1, 2)), z_rho, 'k--')
plt.fill_betweenx(z_rho, np.nanpercentile(data[v_key][0,:,:,:], 25, axis = (1, 2)), np.nanpercentile(data[v_key][0,:,:,:], 75, axis = (1, 2)), facecolor = 'k', alpha = 0.5, edgecolor = 'None')
plt.plot(np.nanpercentile(data[v_key][0,:,:,:], 50, axis = (1, 2)), z_rho, 'k', lw = 2)
plt.plot(np.nanpercentile(data[v_key][0,:,:,:], 10, axis = (1, 2)), z_rho, 'k--')
plt.show()

plt.plot(np.nanpercentile(data[v_key][it,:,:,:], 90, axis = (1, 2)), z_rho, 'k--')
plt.fill_betweenx(z_rho, np.nanpercentile(data[v_key][it,:,:,:], 25, axis = (1, 2)), np.nanpercentile(data[v_key][it,:,:,:], 75, axis = (1, 2)), facecolor = 'k', alpha = 0.5, edgecolor = 'None')
plt.plot(np.nanpercentile(data[v_key][it,:,:,:], 50, axis = (1, 2)), z_rho, 'k', lw = 2)
plt.plot(np.nanpercentile(data[v_key][it,:,:,:], 10, axis = (1, 2)), z_rho, 'k--')
plt.show()

# look at the vertical winds at the start and at hour 14 (last output before crash, roughly 45 minutes before crash)
plt.plot(np.nanpercentile(data[w_key][0,:,:,:], 99, axis = (1, 2)), data['z'], 'k--')
plt.fill_betweenx(data['z'], np.nanpercentile(data[w_key][0,:,:,:], 25, axis = (1, 2)), np.nanpercentile(data[w_key][0,:,:,:], 75, axis = (1, 2)), facecolor = 'k', alpha = 0.5, edgecolor = 'None')
plt.plot(np.nanpercentile(data[w_key][0,:,:,:], 50, axis = (1, 2)), data['z'], 'k', lw = 2)
plt.plot(np.nanpercentile(data[w_key][0,:,:,:], 1, axis = (1, 2)), data['z'], 'k--')
plt.show()

plt.plot(np.nanmax(data[w_key][it,:,:,:], axis = (1, 2)), data['z'], 'k--')
plt.fill_betweenx(data['z'], np.nanpercentile(data[w_key][it,:,:,:], 25, axis = (1, 2)), np.nanpercentile(data[w_key][it,:,:,:], 75, axis = (1, 2)), facecolor = 'k', alpha = 0.5, edgecolor = 'None')
plt.plot(np.nanpercentile(data[w_key][it,:,:,:], 50, axis = (1, 2)), data['z'], 'k', lw = 2)
plt.plot(np.nanmin(data[w_key][it,:,:,:], axis = (1, 2)), data['z'], 'k--')
plt.show()

iz = np.where(np.abs(z_rho - 2700.0) == np.min(np.abs(z_rho - 2700.0)))[0][0]
plt.plot(np.arange(10, data[w_key].shape[0]*10+1, 10)/60., np.nanmax(data[w_key][:,iz,:,:], axis = (1, 2)))
plt.show()

