import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from netCDF4 import Dataset
from scipy import interpolate, integrate
import os
from analysis_tools import get_TKE, get_CC
from datetime import datetime as dt

"""
Code to read in the u, v, w, and mcl data and plot vertical profiles of resolved
scale turbulent kinetic energy and cloud cover. These plots can be made at
several times during the duration of the simulation, or for mean over given
time intervals.
"""

# Read in the data
hours = ["{0:02d}".format(hour) for hour in xrange(0, 24, 3)]
nhours = len(hours)

# initialise a tke time series variable
tke_ts = np.zeros(144)
mcl = np.zeros((144, 127))
mcf = np.zeros((144, 127))

# Define some keys
u_key   = 'STASH_m01s00i002'
v_key   = 'STASH_m01s00i003'
w_key   = 'STASH_m01s00i150'
mcl_key = 'STASH_m01s00i392'
mcf_key = 'STASH_m01s00i393'
rho_key = 'STASH_m01s00i389'

print '[' + dt.now().strftime('%H:%M:%S') + '] Starting...'
for hour in hours:
    print '[' + dt.now().strftime('%H:%M:%S') + '] Opening the netCDFs for hour ' + hour
    wind_nc = Dataset('wind_'+hour+'.nc', 'r') # contains u, v, and w
    mr_nc = Dataset('mr_'+hour+'.nc', 'r')
    rho_nc = Dataset('fluxes_'+hour+'.nc', 'r')

    if day == '01':
        print '[' + dt.now().strftime('%H:%M:%S') + '] Grabbing the height coordinates'
        z_rho = fluxes_nc.variables['rholev_zsea_rho'][:]
        z_theta = wind_nc.variables['thlev_zsea_theta'][:]
        iz_3000 = np.where(abs(z_theta - 3000) == np.min(abs(z_theta - 3000)))[0][0]+1


    # u and v need to be regridded onto theta levels
    print '[' + dt.now().strftime('%H:%M:%S') + '] Regridding rho'
    rho_regrid = regrid(wind_nc, rho_nc, rho_key)
    
    # Closing the netCDFs
    rho_nc.close()
    
    TKE = get_TKE(wind_nc.variables[u_key][:], wind_nc.variables[v_key][:], wind_nc.variables[w_key][:], rho_regrid)
    for it in xrange(TKE.shape[0]):
        tke_ts[hours.index(hour)*18 + it] = integrate.trapz(y = TKE[it,:iz_3000], x = z_theta[:iz_3000])
    
    tke_profile = np.mean(get_TKE(wind_nc.variables[u_key][:], wind_nc.variables[v_key][:], wind_nc.variables[w_key][:]), axis = 0)    
    cc_profile = get_CC(mr_nc.variables[mcl_key][:]+mr_nc.variables[mcf_key][:])
    
    wind_nc.close()

    mcl[hours.index(hour)*18:(hours.index(hour)+1)*18,:] = np.max(mr_nc.variables[mcl_key][:], axis = (2,3))
    mcf[hours.index(hour)*18:(hours.index(hour)+1)*18,:] = np.max(mr_nc.variables[mcf_key][:], axis = (2,3))
    
    mr_nc.close()
    
    print '[' + dt.now().strftime('%H:%M:%S') + '] Making the TKE-CloudCover profiles plot'
    
    fig = plt.figure()
    ax0 = fig.add_subplot(1,2,1)
    ax0.plot(tke_profile, z_theta, 'k', lw = 2, label = 'TKE Profile')
    ax0.set_ylabel('Height (m)')
    ax0.set_xlabel('TKE ($m^2 s^{-2}$)')
    ax0.set_ylim([0, 15000])
    ax0.set_xlim([0, 0.5])
    
    ax1 = fig.add_subplot(1,2,2)
    ax1.plot(cc_profile, z_theta, 'k--', lw = 2, label = 'Cloud Cover Profile')
    ax1.set_xlabel('Cloud Cover')
    ax1.set_xlim([0, 0.1])
    ax1.set_ylim([0, 15000])
    
    fig.suptitle('Profiles averages over day ' + day, fontsize = 14)
    lgd = fig.legend(loc = 'upper center', bbox_to_anchor = (-0.25, -0.1))
    plt.savefig('tke_cc_profiles_day_' + day + '.png', dpi = 100, bbox_extra_artists=(lgd,), bbox_inches = 'tight')
    plt.close('all')


print '[' + dt.now().strftime('%H:%M:%S') + '] Making the vertically integrated TKE timeseries plot'
times = np.arange(1., 1440., 10.)
plt.plot(times/60., tke_ts, 'k', lw = 2)
plt.xlabel('Time (hours)')
plt.ylabel('Vertically Integrated TKE (kg s$^{-2}$)')
plt.title('Integrated below ' + str(z_theta[iz_3000-1]) + ' m')
plt.savefig('../VerticallyIntegratedTKE.png', dpi = 100)
plt.close('all')

print '[' + dt.now().strftime('%H:%M:%S') + '] Making the cloud liquid and cloud ice hovmoller plot'
fig = plt.figure(figsize = (15, 8))
ax = fig.add_subplot(1,1,1)

CL = ax.contourf(times/60., z_theta, np.transpose(mcl), levels = np.linspace(1e-08, 5e-03, 21), cmap = 'Greys')
ax.set_yscale('log', basey = 10, nonposy = 'clip', subsy = [2, 3, 4, 5, 6, 7, 8, 9])
ax.set_ylim([100, 15000])
ax.set_xlim([0, 10])
plt.colorbar(CL, label = 'Cloud Liquid')

CI = ax.contourf(times/60., z_theta, np.transpose(mcf), levels = np.linspace(1e-08, 5e-05, 21), cmap = 'Blues')
ax.set_yscale('log', basey = 10, nonposy = 'clip', subsy = [2, 3, 4, 5, 6, 7, 8, 9])
ax.set_ylim([100, 15000])
ax.set_xlim([0, 10])
plt.colorbar(CI, label = 'Cloud Ice')

ax.set_ylabel('Height (m)')
ax.set_yticks([100, 500, 1000, 2500, 5000, 10000, 15000, 40000])
ax.set_yticklabels([100, 500, 1000, 2500, 5000, 10000, 15000, 40000])
ax.set_xlabel('Time (hours)')
plt.savefig('cloud_timeseries.png', dpi = 100)
plt.show()



