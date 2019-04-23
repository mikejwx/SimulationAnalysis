import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, interpolate
from netCDF4 import Dataset
from SkewT_archer import cpd, Rd, Lv

days = ["{0:02d}".format(day) for day in xrange(1, 11)]
theta_key = 'STASH_m01s00i004'
q_key = 'STASH_m01s00i010'
lhf_key = u'STASH_m01s03i234'
shf_key = u'STASH_m01s03i217'

paths = {'control' : '/nerc/n02/n02/xb899100/CloudTrail/Control_Spinup/',
         'U05'     : '/nerc/n02/n02/xb899100/CloudTrail/U05_Spinup/',
         'U05v2'   : '/work/n02/n02/xb899100/cylc-run/u-bg952/share/data/history/'}
for path_key in paths.keys():
    path = paths[path_key]
    for day in days:
        bouy_nc = Dataset(path + 'bouy_' + day + '.nc', 'r')
        mr_nc = Dataset(path + 'mr_' + day + '.nc', 'r')
        fluxes_nc = Dataset(path + 'fluxes_' + day + '.nc', 'r')
        
        if day == days[0]:
            z = bouy_nc.variables['thlev_zsea_theta'][:]
            time_key = [key for key in bouy_nc.variables.keys() if 'min' in key][0]
            times = bouy_nc.variables[time_key][:]
            theta_mean = np.nanmean(bouy_nc.variables[theta_key][:], axis = (2, 3))
            q_mean = np.nanmean(mr_nc.variables[q_key][:], axis = (2, 3))
            # Deal with the fluxes differently
            time_key_f = [key for key in fluxes_nc.variables.keys() if 'hr' in key][0]
            times_f = fluxes_nc.variables[time_key_f][:]
            shf_mean = np.nanmean(fluxes_nc.variables[shf_key][:], axis = (1, 2))
            lhf_mean = np.nanmean(fluxes_nc.variables[lhf_key][:], axis = (1, 2))
        else:
            times = np.concatenate((times, bouy_nc.variables[time_key][:]), axis = 0)        
            theta_mean = np.concatenate((theta_mean, np.nanmean(bouy_nc.variables[theta_key][:], axis = (2, 3))), axis = 0)
            q_mean = np.concatenate((q_mean, np.nanmean(mr_nc.variables[q_key][:], axis = (2, 3))), axis = 0)
            
            times_f = np.concatenate((times_f, fluxes_nc.variables[time_key_f][:]), axis = 0)
            shf_mean = np.concatenate((shf_mean, np.nanmean(fluxes_nc.variables[shf_key][:], axis = (1, 2))), axis = 0)
            lhf_mean = np.concatenate((lhf_mean, np.nanmean(fluxes_nc.variables[lhf_key][:], axis = (1, 2))), axis = 0)
        bouy_nc.close()
        mr_nc.close()
        fluxes_nc.close()
            
    ### Define our subsidence profile
    w_subs = np.array([0.,-0.00156,-0.00379,-0.00588,-0.00603,-0.00498,-0.00384,-0.00337,-0.00315,-0.00311,-0.00291,-0.00300,-0.00418,-0.00208,0.,0.])
    z_subs = np.array([0.0,148.,392.,628.,857.,1119.,1381.,1644.,1906.,2168.,2430.,2693.,3217.,3867.,4326.0,40000.])
    Q_rad = np.array([-2.3148e-05, -2.3148e-05, 0.0, 0.0])
    z_rad = np.array([0., 3217., 4326., 40000.])

    # Assume surface density
    rho_star = 102000./(Rd*302.3)

    # Term 1: surface fluxes
    # Energy
    shf_ts = shf_mean
    # Moisture
    lhf_ts = lhf_mean

    # Term 2: Subsidence
    # Energy INTEGRAL(- w * dtheta/dz)Dz
    z_half = 0.5*(z[1:] + z[:-1])
    dtheta_dz = (theta_mean[:,1:] - theta_mean[:,:-1])/(z[1:] - z[:-1])
    w = interpolate.interp1d(x = z_subs, y = w_subs, fill_value = 'extrapolate')(z_half)
    theta_subs = cpd*rho_star*integrate.trapz(y = - w*dtheta_dz, x = z_half, axis = 1)
    # Moisture INTEGRAL(- w* dq/dz)Dz
    dq_dz = (q_mean[:, 1:] - q_mean[:,:-1])/(z[1:] - z[:-1])
    q_subs = Lv*rho_star*integrate.trapz(y = - w*dq_dz, x = z_half, axis = 1)

    # Term 3: Radiation
    dtheta_dt_rad = cpd*rho_star*integrate.trapz(y = Q_rad, x = z_rad)

    # Convert times from minutes to days
    times = times/1440.

    # Convert times_f from hours to days
    times_f = times_f/24.

    fig = plt.figure(tight_layout = True, figsize = (15,7))
    ax = fig.add_subplot(1, 2, 1)
    ax.plot(times, theta_subs, 'b', label = '$\\rho c_{pd} \int - w \\frac{\partial \\theta}{\partial z} dz$')
    ax.plot(times_f, shf_ts, 'k', lw = 2, label = '$\\rho c_{pd} (\overline{w^{\prime} \\theta^{\prime}})_{sfc}$')
    ax.plot([times[0], times[-1]], [dtheta_dt_rad, dtheta_dt_rad], 'r', label = '\\rho c_{pd} $\int \\frac{\partial \\theta}{\partial t}_{rad} dz$')
    ax.set_xlim([0, 10])
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Flux (W/m$^{2}$)')
    ax.legend(loc = 0)
    ax.set_title('Budget for ' + path_key)
    
    ax1 = fig.add_subplot(1, 2, 2)
    ax1.plot(times, q_subs, 'brown', label = '$\\rho L_{v} \int - w \\frac{\partial q}{\partial z} dz$')
    ax1.plot(times_f, lhf_ts, 'k', lw = 2, label = '$\\rho L_{v} (\overline{w^{\prime} q^{\prime}})_{sfc}$')
    ax1.set_xlim([0, 10])
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Flux (W/m$^{2}$)')
    ax1.legend(loc = 0)
    
    plt.savefig('../idealised_budget_for_' + path_key + '.png', dpi = 150)
    #plt.show()


