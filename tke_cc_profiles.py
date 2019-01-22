import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate, integrate
import os
from analysis_tools import regrid, get_TKE, get_CC
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

for hour in hours:
    print 'Open the netCDF for hour ' + hour
    u_nc = Dataset('u_'+hour+'.nc', 'r')
    v_nc = Dataset('v_'+hour+'.nc', 'r')
    w_nc = Dataset('bouy_'+hour+'.nc', 'r')
    mr_nc = Dataset('mr_'+hour+'.nc', 'r')
    rho_nc = Dataset('fluxes_'+hour+'.nc', 'r')

    if day == '01':
        print 'Grabbing the height coordinates'
        z_rho = u_nc.variables['rholev_zsea_rho'][:]
        z_theta = w_nc.variables['thlev_zsea_theta'][:]
        iz_3000 = np.where(abs(z_theta - 3000) == np.min(abs(z_theta - 3000)))[0][0]+1


    # u and v need to be regridded onto theta levels
    print 'Regridding the winds...'
    u_regrid = regrid(w_nc, u_nc, u_key)
    v_regrid = regrid(w_nc, v_nc, v_key)
    rho_regrid = regrid(w_nc, rho_nc, rho_key)
    
    # Closing the netCDFs
    u_nc.close()
    v_nc.close()
    rho_nc.close()
    
    TKE = get_TKE(u_regrid, v_regrid, w_nc.variables[w_key][:], rho_regrid)
    for it in xrange(TKE.shape[0]):
        tke_ts[hours.index(hour)*18 + it] = integrate.trapz(y = TKE[it,:iz_3000], x = z_theta[:iz_3000])
    
    tke_profile = np.mean(get_TKE(u_regrid, v_regrid, w_nc.variables[w_key][:]), axis = 0)    
    cc_profile = get_CC(mr_nc.variables[mcl_key][:]+mr_nc.variables[mcf_key][:])
    
    w_nc.close()

    mcl[hours.index(hour)*18:(hours.index(hour)+1)*18,:] = np.max(mr_nc.variables[mcl_key][:], axis = (2,3))
    mcf[hours.index(hour)*18:(hours.index(hour)+1)*18,:] = np.max(mr_nc.variables[mcf_key][:], axis = (2,3))
    
    mr_nc.close()
    
    print 'Making the plot...'
    fig = plt.figure()
    plt.subplot(121)
    plt.plot(tke_profile, z_theta, 'k', lw = 2, label = 'TKE Profile')
    plt.ylabel('Height (m)')
    plt.xlabel('TKE ($m^2 s^{-2}$)')
    plt.ylim([0, 15000])
    plt.xlim([0, 0.5])
    
    plt.subplot(122)
    plt.plot(cc_profile, z_theta, 'k--', lw = 2, label = 'Cloud Cover Profile')
    plt.xlabel('Cloud Cover')
    plt.xlim([0, 0.1])
    plt.ylim([0, 15000])
    plt.suptitle('Profiles averages over day ' + day, fontsize = 14)
    lgd = plt.legend(loc = 'upper center', bbox_to_anchor = (-0.25, -0.1))
    plt.savefig('tke_cc_profiles_day_' + day + '.png', dpi = 150, bbox_extra_artists=(lgd,), bbox_inches = 'tight')
    plt.close('all')
    #plt.show()

times = np.arange(1., 1440., 10.)
plt.plot(times/60., tke_ts, 'k', lw = 2)
plt.xlabel('Time (hours)')
plt.ylabel('Vertically Integrated TKE (kg s$^{-2}$)')
plt.title('Integrated below ' + str(z_theta[iz_3000-1]) + ' m')
plt.savefig('../VerticallyIntegratedTKE.png', dpi = 150)
plt.show()

fig = plt.figure()
fig.set_size_inches(15, 8)
ax = plt.subplot()
CL = ax.contourf(times/60., z_theta, np.transpose(mcl), levels = np.linspace(1e-08, 5e-03, 21), cmap = 'Greys')
ax.set_yscale('log', basey = 10, nonposy = 'clip', subsy = [2, 3, 4, 5, 6, 7, 8, 9])
ax.set_ylim([100, 15000])
ax.set_xlim([0, 10])
plt.colorbar(mappable = CL, label = 'Cloud Liquid')
CI = ax.contourf(times/60., z_theta, np.transpose(mcf), levels = np.linspace(1e-08, 5e-05, 21), cmap = 'Blues')
ax.set_yscale('log', basey = 10, nonposy = 'clip', subsy = [2, 3, 4, 5, 6, 7, 8, 9])
ax.set_ylim([100, 15000])
ax.set_xlim([0, 10])
plt.colorbar(CI, label = 'Cloud Ice')
ax.set_ylabel('Height (m)')
ax.set_yticks([100, 500, 1000, 2500, 5000, 10000, 15000, 40000])
ax.set_yticklabels([100, 500, 1000, 2500, 5000, 10000, 15000, 40000])
ax.set_xlabel('Time (hours)')
plt.savefig('cloud_timeseries.png', dpi = 150)
plt.show()



