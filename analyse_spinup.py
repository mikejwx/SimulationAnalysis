import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from analysis_tools import RDP, lcl, find_h, get_CC, get_CTZ
from scipy import integrate, interpolate

"""
Create plots of the time series of choice variables, and plots of the hodographs
for other choice variables.

This is to highlight when the spin-up simulation nears a quasi-equilibrium state.
"""

days = ["{0:02d}".format(day) for day in xrange(1, 11)]
#paths = {'control' : '/nerc/n02/n02/xb899100/CloudTrail/Control_Spinup/'}
paths = {'U05'     : '/nerc/n02/n02/xb899100/CloudTrail/U05_Spinup/'}

# Define keys for variables in our netCDF
u_key     = 'STASH_m01s00i002' # u_DD
v_key     = 'STASH_m01s00i003' # v_DD
theta_key = 'STASH_m01s00i004' # bouy_DD
mv_key    = 'STASH_m01s00i391' # fluxes_DD
lwp_key   = 'STASH_m01s30i405' # lwp_DD
mcl_key   = 'STASH_m01s00i392' # mr_DD
shf_key   = 'STASH_m01s03i217' # fluxes_DD
lhf_key   = 'STASH_m01s03i234' # fluxes_DD
p_key     = 'STASH_m01s00i407' # fluxes_DD
rho_key   = 'STASH_m01s00i389' # fluxes_DD

# Define some constants
cpd = 1005.
Rd = 287.05
p0 = 1e5

for key in paths.keys():
    for day in days:
        # Read in the netCDF files
        u_nc      = Dataset(paths[key]+'u_'+day+'.nc', 'r')
        v_nc      = Dataset(paths[key]+'v_'+day+'.nc', 'r')
        bouy_nc   = Dataset(paths[key]+'bouy_'+day+'.nc', 'r')
        lwp_nc    = Dataset(paths[key]+'lwp_'+day+'.nc', 'r')
        mr_nc     = Dataset(paths[key]+'mr_'+day+'.nc', 'r')
        fluxes_nc = Dataset(paths[key]+'fluxes_'+day+'.nc', 'r')
        
        # Time Series
        if day == '01':
            # Choose a height to do the wind time series and hodograph
            z_rho = u_nc.variables['rholev_zsea_rho'][:]*1.
            z_theta = mr_nc.variables['thlev_zsea_theta'][:]*1.
            z_theta0 = bouy_nc.variables['thlev_zsea_theta'][:]*1.
            iz_rho = np.where(np.abs(z_rho - 350.) == np.min(np.abs(z_rho - 350.)))[0][0]
            # compute the horizontal means for the time series
            print day + 'Calculating horizontal mean u'
            u_ts = np.nanmean(u_nc.variables[u_key][:,iz_rho,:,:], axis = (1, 2))
            print day + 'Calculating u hovmoller'
            u_hov = np.nanmean(u_nc.variables[u_key][:], axis = (2, 3))
            print day + 'Calculating horizontal mean v'
            v_ts = np.nanmean(v_nc.variables[v_key][:,iz_rho,:,:], axis = (1, 2))
            print day + 'Calculating v hovmoller'
            v_hov = np.nanmean(v_nc.variables[v_key][:], axis = (2, 3))
            print day + 'Calculating horizontal mean liquid water path'
            lwp_ts = np.nanmean(lwp_nc.variables[lwp_key][:], axis = (1, 2))
            print day + 'Calculating horizontal mean cloud cover'
            cc_ts = np.nanmean(np.where(lwp_nc.variables[lwp_key][:] > 1e-16, 1., 0.), axis = (1, 2))
            print day + 'Calculating mcl hovmoller'
            mcl_hov = np.nanmean(mr_nc.variables[mcl_key][:], axis = (2, 3))
            print day + 'Calculating horizontal mean cloud top height'
            ctz_ts = np.nanmean(get_CTZ(mr_nc.variables[mcl_key][:], z_theta), axis = (1, 2))
            print day + 'Calculating horizontal mean sensible heat flux'
            shf_ts = np.nanmean(fluxes_nc.variables[shf_key][:], axis = (1, 2))
            print day + 'Calculating horizontal mean latent heat flux'
            lhf_ts = np.nanmean(fluxes_nc.variables[lhf_key][:], axis = (1, 2))
            print day + 'Calculating temperature'
            temperature = np.array([interpolate.interp1d(x = z_theta0, y = bouy_nc.variables[theta_key][it,:,:,:], axis = 0, fill_value = 'extrapolate')(z_rho)/(p0/fluxes_nc.variables[p_key][it,:,:,:])**(Rd/cpd) for it in xrange(len(bouy_nc.variables['min10'][:]))])
            print day + 'Calculating q'
            q = np.array([interpolate.interp1d(x = z_theta, y = mr_nc.variables[mv_key][it,:,:,:]/(1 + mr_nc.variables[mv_key][it,:,:,:]), fill_value = 'extrapolate', axis = 0)(z_rho) for it in xrange(len(mr_nc.variables['min10'][:]))])
            print day + 'Calculating lifting condensation level'
            lcl_ts = np.nanmean(np.array([lcl(temperature[it,0,:,:], q[it,0,:,:], fluxes_nc.variables[p_key][it,:,:,:], z_rho) for it in xrange(q.shape[0])]), axis = (1, 2))
            print day + 'Calculating virtual potential temperature'
            theta_v = np.array([interpolate.interp1d(x = z_theta0, y = bouy_nc.variables[theta_key][it,:,:,:], axis = 0, fill_value = 'extrapolate')(z_rho) for it in xrange(len(bouy_nc.variables['min10'][:]))])*(1. + 0.608*q)
            print day + 'Calculating potential temperature hovmoller'
            theta_hov = np.nanmean(bouy_nc.variables[theta_key][:], axis = (2, 3))
            print day + 'Calculating boundary layer height'
            zi_ts = np.nanmean(np.array([find_h(theta_v[it,:,:,:], u_nc.variables[u_key][it,:,:,:], v_nc.variables[v_key][it,:,:,:], z_rho) for it in xrange(theta_v.shape[0])]), axis = (1, 2))
            print day + 'Calculating column integrated water vapour'
            ciwv_ts = np.nanmean(np.array([integrate.trapz(y = interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][it,:,:,:], axis = 0, fill_value = 'extrapolate')(z_theta)*mr_nc.variables[mv_key][it,:,:,:], x = z_theta, axis = 0) for it in xrange(q.shape[0])]), axis = (1, 2))
            print day + 'Calculating mv hovmoller'
            mv_hov = np.nanmean(mr_nc.variables[mv_key][:], axis = (2, 3))
            times = bouy_nc.variables['min10'][:]
        else:
            print day + 'Calculating horizontal mean u'
            u_ts = np.concatenate((u_ts, np.nanmean(u_nc.variables[u_key][:,iz_rho,:,:], axis = (1, 2))), axis = 0)
            print day + 'Calculating horizontal mean v'
            v_ts = np.concatenate((v_ts, np.nanmean(v_nc.variables[v_key][:,iz_rho,:,:], axis = (1, 2))), axis = 0)
            print day + 'Calculating horizontal mean liquid water path'
            lwp_ts = np.concatenate((lwp_ts, np.nanmean(lwp_nc.variables[lwp_key][:], axis = (1, 2))), axis = 0)
            print day + 'Calculating horizontal mean cloud cover'
            cc_ts = np.concatenate((cc_ts, np.nanmean(np.where(lwp_nc.variables[lwp_key][:] > 1e-16, 1., 0.), axis = (1, 2))), axis = 0)
            print day + 'Calculating horizontal mean cloud top height'
            ctz_ts = np.concatenate((ctz_ts, np.nanmean(get_CTZ(mr_nc.variables[mcl_key][:], z_theta), axis = (1, 2))), axis = 0)
            print day + 'Calculating horizontal mean sensible heat flux'
            shf_ts = np.concatenate((shf_ts, np.nanmean(fluxes_nc.variables[shf_key][:], axis = (1, 2))), axis = 0)
            print day + 'Calculating horizontal mean latent heat flux'
            lhf_ts = np.concatenate((lhf_ts, np.nanmean(fluxes_nc.variables[lhf_key][:], axis = (1, 2))), axis = 0)
            print day + 'Calculating temperature'
            temperature = np.array([interpolate.interp1d(x = z_theta0, y = bouy_nc.variables[theta_key][it,:,:,:], axis = 0, fill_value = 'extrapolate')(z_rho)/(p0/fluxes_nc.variables[p_key][it,:,:,:])**(Rd/cpd) for it in xrange(len(bouy_nc.variables['min10'][:]))])
            print day + 'Calculating q'
            q = np.array([interpolate.interp1d(x = z_theta, y = mr_nc.variables[mv_key][it,:,:,:]/(1 + mr_nc.variables[mv_key][it,:,:,:]), fill_value = 'extrapolate', axis = 0)(z_rho) for it in xrange(len(mr_nc.variables['min10'][:]))])
            print day + 'Calculating lifting condensation level'
            lcl_ts = np.concatenate((lcl_ts, np.nanmean(np.array([lcl(temperature[it,0,:,:], q[it,0,:,:], fluxes_nc.variables[p_key][it,:,:,:], z_rho) for it in xrange(q.shape[0])]), axis = (1, 2))), axis = 0)
            print day + 'Calculating virtual potential temperature'
            theta_v = np.array([interpolate.interp1d(x = z_theta0, y = bouy_nc.variables[theta_key][it,:,:,:], axis = 0, fill_value = 'extrapolate')(z_rho) for it in xrange(len(bouy_nc.variables['min10'][:]))])*(1. + 0.608*q)
            print day + 'Calculating boundary layer height'
            zi_ts = np.concatenate((zi_ts, np.nanmean(np.array([find_h(theta_v[it,:,:,:], u_nc.variables[u_key][it,:,:,:], v_nc.variables[v_key][it,:,:,:], z_rho) for it in xrange(theta_v.shape[0])]), axis = (1, 2))), axis = 0)
            print day + 'Calculating column integrated water vapour'
            ciwv_ts = np.concatenate((ciwv_ts, np.nanmean(np.array([integrate.trapz(y = interpolate.interp1d(x = z_rho, y = fluxes_nc.variables[rho_key][it,:,:,:], axis = 0, fill_value = 'extrapolate')(z_theta)*mr_nc.variables[mv_key][it,:,:,:], x = z_theta, axis = 0) for it in xrange(q.shape[0])]), axis = (1, 2))), axis = 0)
            print day + 'Calculating u hovmoller'
            u_hov = np.concatenate((u_hov, np.nanmean(u_nc.variables[u_key][:], axis = (2, 3))), axis = 0)
            print day + 'Calculating v hovmoller'
            v_hov = np.concatenate((v_hov, np.nanmean(v_nc.variables[v_key][:], axis = (2, 3))), axis = 0)
            print day + 'Calculating mcl hovmoller'
            mcl_hov = np.concatenate((mcl_hov, np.nanmean(mr_nc.variables[mcl_key][:], axis = (2, 3))), axis = 0)
            print day + 'Calculating mv hovmoller'
            mv_hov = np.concatenate((mv_hov, np.nanmean(mr_nc.variables[mv_key][:], axis = (2, 3))), axis = 0)
            print day + 'Calculating potential temperature hovmoller'
            theta_hov = np.concatenate((theta_hov, np.nanmean(bouy_nc.variables[theta_key][:], axis = (2, 3))), axis = 0)
            times = np.concatenate((times, bouy_nc.variables['min10'][:]), axis = 0)
        # Hovmoller t-z

    # Plot the time series
    # [[u, v, hodograph], [lwp, cc, ctz]]
    fig = plt.figure(tight_layout = True, figsize=(13,12))
    ax0 = fig.add_subplot(2, 3, 1)
    ax0.plot(times/(24.*60.), u_ts, 'b', lw = 2)
    ax0.set_xlabel('Time (days)')
    ax0.set_xlim([0, 10])
    ax0.set_ylabel('u (m/s)')
    ax0.set_title('Zonal Winds Time Series')
    
    ax1 = fig.add_subplot(2, 3, 2)
    ax1.plot(times/(24.*60.), v_ts, 'g', lw = 2)
    ax1.set_xlabel('Time (days)')
    ax1.set_xlim([0, 10])
    ax1.set_ylabel('v (m/s)')
    ax1.set_title('Meridional Winds Time Series')
    
    ax2 = fig.add_subplot(2, 3, 3)
    ax2.plot(u_ts, v_ts, 'k', lw = 2)
    ax2.set_xlabel('u (m/s)')
    ax2.set_ylabel('v (m/s)')
    ax2.set_title('Hodograph')
    
    ax3 = fig.add_subplot(2, 3, 4)
    ax3.plot(np.arange(1, 14400.1, 1.)/(24.*60.), lwp_ts*1000., 'k', lw = 2)
    ax3.set_xlabel('Time (days)')
    ax3.set_xlim([0, 10])
    ax3.set_ylabel('lwp (g/m$^{-2}$)')
    ax3.set_title('Liquid Water Path Time Series')
    
    ax4 = fig.add_subplot(2, 3, 5)
    ax4.plot(np.arange(1, 14400.1, 1.)/(24.*60.), cc_ts, 'k', lw = 2)
    ax4.set_xlabel('Time (days)')
    ax4.set_xlim([0, 10])
    ax4.set_ylabel('cloud cover)')
    ax4.set_title('Cloud Cover Time Series')
    
    ax5 = fig.add_subplot(2, 3, 6)
    ax5.plot(times/(24.*60.), ctz_ts, 'k', lw = 2)
    ax5.set_xlabel('Time (days)')
    ax5.set_xlim([0, 10])
    ax5.set_ylabel('cloud top height (m)')
    ax5.set_title('Cloud Top Height Time Series')
    plt.savefig('../'+key+'_spin-up-ic_winds_clouds.png', dpi = 150.)
    plt.show()
    
    # shf and lhf
    fig = plt.figure(tight_layout = True, figsize=(6,12))
    ax0 = fig.add_subplot(2, 1, 1)
    ax0.plot(np.arange(1, 240.1, 1.)/(24.), shf_ts, 'r', lw = 2)
    ax0.set_ylabel('Heat Flux (W/m$^{2}$)')
    ax0.set_title('Sensible Heat Flux Time Series')
    
    ax1 = fig.add_subplot(2, 1, 2, sharex = ax0)
    ax1.plot(np.arange(1, 240.1, 1.)/(24.), lhf_ts, 'b', lw = 2)
    ax1.set_ylabel('Heat Flux (W/m$^{2}$)')
    ax1.set_xlabel('Time (days)')
    ax1.set_xlim([0, 10])
    ax1.set_title('Latent Heat Flux Time Series')
    plt.savefig('../'+key+'_spin-up-ic_surface_fluxes.png', dpi = 150.)
    plt.show()
    
    # zi and lcl
    fig = plt.figure(tight_layout = True, figsize=(9,9))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(times/(24.*60.), lcl_ts, 'k--', lw = 2)
    ax.plot(times/(24.*60.), zi_ts, 'k', lw = 2)
    ax.set_xlabel('Time (days)')
    ax.set_xlim([0, 10])
    ax.set_ylabel('Height (m)')
    ax.set_title('LCL and z$_i$ Time Series')
    plt.savefig('../'+key+'_spin-up-ic_zi_lcl.png', dpi = 150.)
    plt.show()
    
    # Column integrated water vapour
    fig = plt.figure(tight_layout = True, figsize=(9,9))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(times/(24.*60.), ciwv_ts, 'k', lw = 2)
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Column Integrated Water Vapour (mm)')
    ax.set_xlim([0, 10])
    ax.set_title('Column Integrated Water Vapour Time Series')
    plt.savefig('../'+key+'_spin-up-ic_PWAT.png', dpi = 150.)
    plt.show()
    
    # Plot the hovmollers
    fig = plt.figure(tight_layout = True, figsize=(18,12))
    ax = fig.add_subplot(1, 1, 1)
    U_plt = ax.contourf(times/(24.*60.), z_rho, np.transpose(u_hov), levels = [level for level in np.arange(-10, 10.1, 1) if level != 0], cmap = 'bwr', extend = 'both')
    plt.colorbar(U_plt, ax = ax, ticks = np.arange(-10, 10.1, 2.), label = 'u (m/s)')
    ax.semilogy()
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Height (m)')
    ax.set_xlim([0, 10])
    ax.set_title('Zonal Wind Hovmoller')
    ax.set_yticklabels([1, 10, 100, 1000, 10000])
    plt.savefig('../'+key+'_spin-up-ic_u_hov.png', dpi = 150.)
    plt.show()
    
    fig = plt.figure(tight_layout = True, figsize=(18,12))
    ax = fig.add_subplot(1, 1, 1)
    V_plt = ax.contourf(times/(24.*60.), z_rho, np.transpose(v_hov), levels = [level for level in np.arange(-10, 10.1, 1) if level != 0], cmap = 'bwr', extend = 'both')
    plt.colorbar(V_plt, ax = ax, ticks = np.arange(-10, 10.1, 2.), label = 'v (m/s)')
    ax.semilogy()
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Height (m)')
    ax.set_xlim([0, 10])
    ax.set_title('Meridional Wind Hovmoller')
    ax.set_yticklabels([1, 10, 100, 1000, 10000])
    plt.savefig('../'+key+'_spin-up-ic_v_hov.png', dpi = 150.)
    plt.show()
    
    fig = plt.figure(tight_layout = True, figsize=(18,12))
    ax = fig.add_subplot(1, 1, 1)
    mv_plt = ax.contourf(times/(24.*60.), z_theta[1:], np.transpose(mv_hov[:,1:])*1000., levels = np.arange(0., 25.1, 2.5), cmap = 'BrBG', extend = 'max')
    plt.colorbar(mv_plt, ax = ax, ticks = np.arange(0., 25.1, 5.), label = '$m_{v}$ (g/kg)')
    ax.semilogy()
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Height (m)')
    ax.set_xlim([0, 10])
    ax.set_title('Water Vapour mixing ratio Hovmoller')
    ax.set_yticks([1, 10, 100, 1000, 10000])
    ax.set_yticklabels([1, 10, 100, 1000, 10000])
    plt.savefig('../'+key+'_spin-up-ic_mv_hov.png', dpi = 150.)
    plt.show()
    
    fig = plt.figure(tight_layout = True, figsize=(18,12))
    ax = fig.add_subplot(1, 1, 1)
    theta_plt = ax.contourf(times/(24.*60.), z_theta[1:], np.transpose(theta_hov[:,1:] - theta_hov[0,1:]), levels = [level for level in np.arange(-2., 2.01, 0.1) if level != 0], cmap = 'bwr', extend = 'both')
    plt.colorbar(theta_plt, ax = ax, ticks = np.arange(-2., 2.01, 0.25), label = '$\\theta^{\prime}$ (K)')
    ax.semilogy()
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Height (m)')
    ax.set_xlim([0, 10])
    ax.set_title('Potential Temperature Anomaly from IC Hovmoller')
    ax.set_yticks([1, 10, 100, 1000, 10000])
    ax.set_yticklabels([1, 10, 100, 1000, 10000])
    plt.savefig('../'+key+'_spin-up-ic_theta_hov.png', dpi = 150.)
    plt.show()

