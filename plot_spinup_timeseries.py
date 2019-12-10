import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import zi_new_key, rho_key, u_key, v_key, shf_key, lhf_key
from SkewT_archer import Lv
days = ["{0:02d}".format(day) for day in range(1, 11)]

#path = '/work/n02/n02/xb899100/cylc-run/u-bn460/share/data/history/'
path = '/nerc/n02/n02/xb899100/CloudTrail/Spinup_Control/'
print path
spinup_data = {}
for day in days:
    print day
    # Read the data
    zi_nc     = Dataset(path + 'zi_' + day + '.nc', 'r') # we want the boundary layer height from here
    u_nc      = Dataset(path + 'u_' + day + '.nc', 'r') # we want the u and v wind components from here
    v_nc      = Dataset(path + 'v_' + day + '.nc', 'r')
    fluxes_nc = Dataset(path + 'fluxes_' + day + '.nc', 'r') # we want rho, and the surface fluxes from here
    
    if day == days[0]:
        print 'first day'
        print 'zi'
        spinup_data['zi'] = np.nanmean(zi_nc.variables[zi_new_key][:], axis = (1, 2))
        print 'shf'
        spinup_data['shf'] = np.nanmean(fluxes_nc.variables[shf_key][1:,0,:,:], axis = (1, 2))
        print 'lhf'
        spinup_data['lhf'] = np.nanmean(fluxes_nc.variables[lhf_key][1:,0,:,:], axis = (1, 2))*Lv
        # Get the height coordinates
        z_rho = fluxes_nc.variables['thlev_zsea_theta'][:]*1.
        
        # get the times
        time_key = 'time'
        spinup_data['time'] = zi_nc.variables[time_key][:]*1.
        # get a mass-weighted mixed layer mean wind
        u_temp = []
        v_temp = []
        for it in range(u_nc.variables[u_key].shape[0]-1):
            iz = np.where(np.abs(z_rho - spinup_data['zi'][it]) == np.min(np.abs(z_rho - spinup_data['zi'][it])))[0][0] + 1
            u_temp.append(np.nanmean([fluxes_nc.variables[rho_key][it+1,idx_z,:,:]*u_nc.variables[u_key][it+1,idx_z,:,:]/np.nanmean(fluxes_nc.variables[rho_key][it+1,:iz,:,:], axis = 0) for idx_z in range(iz)]))
            v_temp.append(np.nanmean([fluxes_nc.variables[rho_key][it+1,idx_z,:,:]*v_nc.variables[v_key][it+1,idx_z,:,:]/np.nanmean(fluxes_nc.variables[rho_key][it+1,:iz,:,:], axis = 0) for idx_z in range(iz)]))
        print 'winds'
        spinup_data['u'] = np.array(u_temp)
        spinup_data['v'] = np.array(v_temp)
    else:
        spinup_data['zi']  = np.concatenate((spinup_data['zi'], np.nanmean(zi_nc.variables[zi_new_key][:], axis = (1, 2))), axis = 0)
        spinup_data['shf'] = np.concatenate((spinup_data['shf'], np.nanmean(fluxes_nc.variables[shf_key][:,0,:,:], axis = (1, 2))), axis = 0)
        spinup_data['lhf'] = np.concatenate((spinup_data['lhf'], np.nanmean(fluxes_nc.variables[lhf_key][:,0,:,:], axis = (1, 2))*Lv), axis = 0)
        
        # get a mass-weighted mixed layer mean wind
        u_temp = []
        v_temp = []
        for it in range(u_nc.variables[u_key].shape[0]):
            iz = np.where(np.abs(z_rho - spinup_data['zi'][it]) == np.min(np.abs(z_rho - spinup_data['zi'][it])))[0][0] + 1
            u_temp.append(np.nanmean([fluxes_nc.variables[rho_key][it,idx_z,:,:]*u_nc.variables[u_key][it,idx_z,:,:]/np.nanmean(fluxes_nc.variables[rho_key][it,:iz,:,:], axis = 0) for idx_z in range(iz)]))
            v_temp.append(np.nanmean([fluxes_nc.variables[rho_key][it,idx_z,:,:]*v_nc.variables[v_key][it,idx_z,:,:]/np.nanmean(fluxes_nc.variables[rho_key][it,:iz,:,:], axis = 0) for idx_z in range(iz)]))
            
        spinup_data['u'] = np.concatenate((spinup_data['u'], np.array(u_temp)), axis = 0)
        spinup_data['v'] = np.concatenate((spinup_data['v'], np.array(v_temp)), axis = 0)
        spinup_data['time'] = np.concatenate((spinup_data['time'], zi_nc.variables['time'][:]), axis = 0)
    zi_nc.close()
    u_nc.close()
    v_nc.close()
    fluxes_nc.close()

# Plot a figure
fig = plt.figure(figsize = (12, 12/3.), tight_layout = True)
axa = fig.add_subplot(1, 3, 1)
axb = fig.add_subplot(1, 3, 2)
axc = fig.add_subplot(1, 3, 3)

axa.plot(spinup_data['time']/1440.0, spinup_data['u'], 'k', lw = 2, label = u'u-wind')
axa.plot(spinup_data['time']/1440.0, spinup_data['v'], 'k', lw = 0.5, label = u'v-wind')
axa.plot([6.0, 6.0], [-10.0, 0.0], color = 'grey', ls = ':')
axa.set_ylabel('Wind (m s$^{-1}$)')
axa.set_xlabel('Time (days)')
axa.set_ylim([-10.0, 0.0])
axa.legend(loc = 0, frameon = False)
axa.text(0.5, -1.25, 'a)')

axb.plot(spinup_data['time']/1440.0, spinup_data['zi'], 'k', lw = 2)
axb.plot([6.0, 6.0], [300, 900], color = 'grey', ls = ':')
axb.set_ylabel(u'z$_{h}$ (m)')
axb.set_xlabel(u'Time (days)')
axb.set_ylim([300, 900])
axb.text(0.5, 825, 'b)')

axc.plot(spinup_data['time']/1440.0, spinup_data['shf'], 'r', lw = 2, label = 'SHF')
axc.plot([0], [0], 'b', label = 'LHF') # just a point so that the item can be added to the legend
axc.set_ylabel(u'Sensible Heat Flux (W m$^{-2}$)', color = 'red')
axc.set_xlabel('Time (days)')
axc.text(0.5, 17.5, 'c)')
axc.set_yticklabels(range(0, 21, 5), color = 'red')

axd = axc.twinx()
axd.plot(spinup_data['time']/1440.0, spinup_data['lhf'], 'b', lw = 2)
axd.set_ylabel(u'Latent Heat Flux (W m$^{-2}$)', color = 'b')
axd.set_yticklabels(range(150, 230, 10), color = 'b')
axc.plot([6.0, 6.0], [0, 20], color = 'grey', ls = ':')
axc.legend(loc = 1, frameon = False)
axc.set_ylim([0, 20])

plt.savefig('../CT_spinup.png', dpi = 150, bbox_inches = 'tight')
plt.show()



