import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, interpolate
from netCDF4 import Dataset
from SkewT_archer import cpd, Rd, Lv, p0

days = ["{0:02d}".format(day) for day in xrange(1, 11)]
from STASH_keys import rho_key, lhf_key, shf_key, tinc_ideal_key, qinc_ideal_key, prho_key, tinc_rain_key, qinc_rain_key, tinc_bl_cld_key, qinc_bl_cld_key, tinc_total_key, qinc_total_key

paths = {'control' : '/nerc/n02/n02/xb899100/CloudTrail/Control_Spinup/'}#,
#         'U05'     : '/nerc/n02/n02/xb899100/CloudTrail/U05_Spinup/',
#         'U05v2'   : '/nerc/n02/n02/xb899100/CloudTrail/U05v2_Spinup/',
#         'U05v3'   : '/nerc/n02/n02/xb899100/CloudTrail/U05v3_Spinup/'
#         '1km'     : '/nerc/n02/n02/xb899100/CloudTrail/Spinup_1km/'}

for path_key in paths.keys():
    path = paths[path_key]
    for day in days:
        print day
        fluxes_nc = Dataset(path + 'fluxes_' + day + '.nc', 'r')
        qinc_nc   = Dataset(path + 'qinc_' + day + '.nc', 'r')
        tinc_nc   = Dataset(path + 'tempinc_' + day + '.nc', 'r')
        
        if day == days[0]:
            z = fluxes_nc.variables['thlev_zsea_theta'][:]
            z_rho = fluxes_nc.variables['rholev_zsea_rho'][:]
            # Deal with the fluxes differently
            #time_key = [key for key in fluxes_nc.variables.keys() if 'hr' in key][0]
            #times = fluxes_nc.variables[time_key_f][:]
            
            # Get air density and the surface fluxes
            rho_mean = np.nanmean(fluxes_nc.variables[rho_key][:], axis = (2, 3))
            shf_mean = np.nanmean(fluxes_nc.variables[shf_key][:,0,:,:], axis = (1, 2))
            lhf_mean = np.nanmean(fluxes_nc.variables[lhf_key][:,0,:,:], axis = (1, 2))
            
            # Get the idealised increments to temp and q
            tideal_mean = np.nanmean(tinc_nc.variables[tinc_ideal_key][:], axis = (2, 3))
            train_mean  = np.nanmean(tinc_nc.variables[tinc_rain_key][:], axis = (2, 3))
            tblcld_mean = np.nanmean(tinc_nc.variables[tinc_bl_cld_key][:], axis = (2, 3))
            ttotal_mean = np.nanmean(tinc_nc.variables[tinc_total_key][:], axis = (2, 3))
            qideal_mean = np.nanmean(qinc_nc.variables[qinc_ideal_key][:], axis = (2, 3))
            qrain_mean  = np.nanmean(qinc_nc.variables[qinc_rain_key][:], axis = (2, 3))
            qblcld_mean = np.nanmean(qinc_nc.variables[qinc_bl_cld_key][:], axis = (2, 3))
            qtotal_mean = np.nanmean(qinc_nc.variables[qinc_total_key][:], axis = (2, 3))
            
            # Get the pressure
            p_mean = np.nanmean(fluxes_nc.variables[prho_key][:], axis = (2, 3))
        else:
            #times = np.concatenate((times, fluxes_nc.variables[time_key][:]), axis = 0)
            
            # Get air density and surface fluxes
            rho_mean = np.concatenate((rho_mean, np.nanmean(fluxes_nc.variables[rho_key][:], axis = (2, 3))), axis = 0)
            shf_mean = np.concatenate((shf_mean, np.nanmean(fluxes_nc.variables[shf_key][:,0,:,:], axis = (1, 2))), axis = 0)
            lhf_mean = np.concatenate((lhf_mean, np.nanmean(fluxes_nc.variables[lhf_key][:,0,:,:], axis = (1, 2))), axis = 0)
            
            # Get the idealised increments to temp and q
            tideal_mean = np.concatenate((tideal_mean, np.nanmean(tinc_nc.variables[tinc_ideal_key][:], axis = (2, 3))), axis = 0)
            train_mean  = np.concatenate((train_mean, np.nanmean(tinc_nc.variables[tinc_rain_key][:], axis = (2, 3))), axis = 0)
            tblcld_mean = np.concatenate((tblcld_mean, np.nanmean(tinc_nc.variables[tinc_bl_cld_key][:], axis = (2, 3))), axis = 0)
            ttotal_mean = np.concatenate((ttotal_mean, np.nanmean(tinc_nc.variables[tinc_total_key][:], axis = (2, 3))), axis = 0)
            qideal_mean = np.concatenate((qideal_mean, np.nanmean(qinc_nc.variables[qinc_ideal_key][:], axis = (2, 3))), axis = 0)
            qrain_mean  = np.concatenate((qrain_mean, np.nanmean(qinc_nc.variables[qinc_rain_key][:], axis = (2, 3))), axis = 0)
            qblcld_mean = np.concatenate((qblcld_mean, np.nanmean(qinc_nc.variables[qinc_bl_cld_key][:], axis = (2, 3))), axis = 0)
            qtotal_mean = np.concatenate((qtotal_mean, np.nanmean(qinc_nc.variables[qinc_total_key][:], axis = (2, 3))), axis = 0)
            
            # Get the pressure
            p_mean = np.concatenate((p_mean, np.nanmean(fluxes_nc.variables[prho_key][:], axis = (2, 3))), axis = 0)
        fluxes_nc.close()
        tinc_nc.close()
        qinc_nc.close()

# Define exner pressure
pi_mean = interpolate.interp1d(x = z_rho, y = (100.*p0/p_mean)**(Rd/cpd), fill_value = 'extrapolate', axis = 1)(z)
### Define our radiative cooling profile ###
Q_rad = np.array([-2.0, -2.0, 0.0, 0.0])/86400. # K/sec
z_rad = np.array([0., 3217., 4326., 40000.]) # m
Q_rad = interpolate.interp1d(x = z_rad, y = Q_rad, fill_value = 'extrapolate')(z)/pi_mean

# Assume surface density
rho_star = 101700./(Rd*302.3)

# Assume timestep
tstep = 3.
# Term 1: Increments due to surface fluxes
# Energy = shf_mean
# Moisture = lhf_mean

# Term 2: Idealised increments (i.e. due to subsidence)
# Energy
theta_subs = cpd*integrate.trapz(y = rho_mean*(tideal_mean/tstep - Q_rad), x = z, axis = 1)
# Moisture
q_subs = Lv*integrate.trapz(y = qideal_mean/tstep, x = z, axis = 1)

# Term 3: Idealised increment due to 'radiation'
dtheta_dt_rad = cpd*integrate.trapz(y = rho_mean*Q_rad, x = z, axis = 1)

# Term 4: Increment due to large scale rain scheme
train_ts = cpd*integrate.trapz(y = rho_mean*train_mean/tstep, x = z, axis = 1)
qrain_ts = Lv*integrate.trapz(y = rho_mean*qrain_mean/tstep, x = z, axis = 1)

# Term 5: BL and ls_cld schemes
tblcld_ts = cpd*integrate.trapz(y = rho_mean*tblcld_mean/tstep, x = z, axis = 1)
qblcld_ts = Lv*integrate.trapz(y = rho_mean*qblcld_mean/tstep, x = z, axis = 1)

# Total
ttotal_ts = cpd*integrate.trapz(y = rho_mean*ttotal_mean/tstep, x = z, axis = 1)
qtotal_ts = Lv*integrate.trapz(y = rho_mean*qtotal_mean/tstep, x = z, axis = 1)

n_ts = 4*144.
my_box = np.sin(np.arange(n_ts)*np.pi/n_ts)
ttotal_smooth = np.convolve(ttotal_ts, my_box, mode = 'same')/np.sum(my_box)
qtotal_smooth = np.convolve(qtotal_ts, my_box, mode = 'same')/np.sum(my_box)

fig = plt.figure(tight_layout = True, figsize = (15,7))
ax = fig.add_subplot(1, 2, 1)
ax.plot(np.linspace(0., 14400., len(ttotal_ts))/(60.*24.), ttotal_ts, 'bo', label = 'raw points')
ax.plot(np.linspace(0., 14400., len(ttotal_ts))/(60.*24.), ttotal_smooth, 'r', label = '4-day smooth')
ax.legend(loc = 0)
ax.grid()

ax1 = fig.add_subplot(1, 2, 2)
ax1.plot(np.linspace(0., 14400., len(qtotal_ts))/(60.*24.), qtotal_ts, 'bo', label = 'raw points')
ax1.plot(np.linspace(0., 14400., len(qtotal_ts))/(60.*24.), qtotal_smooth, 'r', label = '4-day smooth')
ax1.grid()
plt.show()

"""
fig = plt.figure(tight_layout = True, figsize = (15,7))
ax = fig.add_subplot(1, 2, 1)
ax.plot(np.linspace(0., 14400., len(theta_subs))/(24.*60.), theta_subs, 'r--', label = '$ - c_{pd} \int \\rho w \\frac{\partial \\theta}{\partial z} dz$')
ax.plot(np.linspace(0., 14400., len(shf_mean))/(24.*60.), shf_mean, 'r', lw = 2, label = '$\\rho c_{pd} (\overline{w^{\prime} \\theta^{\prime}})_{sfc}$')
ax.plot(np.linspace(0., 14400., len(dtheta_dt_rad))/(24.*60.), dtheta_dt_rad, 'b', label = '$c_{pd} \int \\rho \\frac{\partial \\theta}{\partial t}_{rad} dz$')
ax.plot(np.linspace(0., 14400., len(train_ts))/(24.*60.), train_ts, 'purple', label = 'ls_rain')
#ax.plot(np.linspace(0., 14400., len(tblcld_ts))/(24.*60.), tblcld_ts, 'grey', label = 'bl + ls_cld')
#ax.plot(np.linspace(0., 14400., len(shf_mean))/(24.*60.), theta_subs + shf_mean + dtheta_dt_rad + train_ts + tblcld_ts, 'k--', lw = 2, label = 'residual')
ax.plot([0., 10.], [0., 0.], color = 'grey', ls = ':')
ax.set_xlim([0, 10])
ax.set_xlabel('Time (days)')
ax.set_ylabel('Flux (W/m$^{2}$)')
ax.legend(loc = 0)
ax.set_title('Budget for ' + path_key)

ax1 = fig.add_subplot(1, 2, 2)
ax1.plot(np.linspace(0., 14400., len(q_subs))/(24.*60.), q_subs, 'brown', label = '$- L_{v} \int w \\frac{\partial q}{\partial z} dz$')
ax1.plot(np.linspace(0., 14400., len(lhf_mean))/(24.*60.), Lv*lhf_mean, 'darkgreen', label = '$\\rho L_{v} (\overline{w^{\prime} q^{\prime}})_{sfc}$')
ax1.plot(np.linspace(0., 14400., len(qrain_ts))/(24.*60.), qrain_ts, 'purple', label = 'ls_rain')
#ax1.plot(np.linspace(0., 14400., len(qblcld_ts))/(24.*60.), qblcld_ts, 'grey', label = 'bl + ls_cld')
#ax1.plot(np.linspace(0., 14400., len(lhf_mean))/(24.*60.), q_subs + Lv*lhf_mean + qrain_ts + qblcld_ts, 'k--', lw = 2, label = 'residual')
ax1.plot([0., 10.], [0., 0.], color = 'grey', ls = ':')
ax1.set_xlim([0, 10])
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Flux (W/m$^{2}$)')
ax1.legend(loc = 0)

#plt.savefig('../idealised_budget_for_' + path_key + '.png', dpi = 150)
plt.show()
"""

