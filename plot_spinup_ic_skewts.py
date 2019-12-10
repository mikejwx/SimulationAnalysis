import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import *
from scipy import interpolate, integrate
<<<<<<< HEAD
from STASH_keys import prho_key, lhf_key, shf_key, theta_key, mv_key, mcl_key, mr_key, ls_rain_amt_key
from netCDF4 import Dataset

# Read the Spinup Data
for day in ["{0:02d}".format(d) for d in range(1, 11)]:
    bouy_nc = Dataset('/work/n02/n02/xb899100/cylc-run/u-bc341/share/data/history/bouy_' + day + '.nc', 'r')
    mr_nc = Dataset('/work/n02/n02/xb899100/cylc-run/u-bc341/share/data/history/mr_' + day + '.nc', 'r')
    if day == '01':
        theta = bouy_nc.variables[theta_key][1:,:,:,:].mean(axis = (2, 3))
        thetal = (bouy_nc.variables[theta_key][1:,:,:,:] - Lv*mr_nc.variables[mcl_key][1:,:,:,:]/cpd).mean(axis = (2, 3))
        z_theta = bouy_nc.variables['thlev_zsea_theta'][:]*1.
    else:
        theta = np.concatenate((theta, bouy_nc.variables[theta_key][:].mean(axis = (2, 3))), axis = 0)
        thetal = np.concatenate((thetal, (bouy_nc.variables[theta_key][:] - Lv*mr_nc.variables[mcl_key][:]/cpd).mean(axis = (2, 3))), axis = 0)
    
    bouy_nc.close()
    mr_nc.close()
    
    with Dataset('/work/n02/n02/xb899100/cylc-run/u-bc341/share/data/history/mr_' + day + '.nc', 'r') as mr_nc:
        if day == '01':
            mt = (mr_nc.variables[mv_key][1:,:,:,:] + mr_nc.variables[mcl_key][1:,:,:,:] + mr_nc.variables[mr_key][1:,:,:,:]).mean(axis = (2, 3))
            mv = mr_nc.variables[mv_key][1:,:,:,:].mean(axis = (2, 3))
        else:
            mt = np.concatenate((mt, (mr_nc.variables[mv_key][:] + mr_nc.variables[mcl_key][:] + mr_nc.variables[mr_key][:]).mean(axis = (2, 3))), axis = 0)
            mv = np.concatenate((mv, mr_nc.variables[mv_key][:].mean(axis = (2, 3))), axis = 0)
    with Dataset('/work/n02/n02/xb899100/cylc-run/u-bc341/share/data/history/fluxes_' + day + '.nc', 'r') as fluxes_nc:
        if day == '01':
            pressure = fluxes_nc.variables[prho_key][1:,:,:,:].mean(axis = (2, 3))
            shf = fluxes_nc.variables[shf_key][1:,0,:,:].mean(axis = (1, 2))
            lhf = fluxes_nc.variables[lhf_key][1:,0,:,:].mean(axis = (1, 2))
        else:
            pressure = np.concatenate((pressure, fluxes_nc.variables[prho_key][:].mean(axis = (2, 3))), axis = 0)
            shf = np.concatenate((shf, fluxes_nc.variables[shf_key][:,0,:,:].mean(axis = (1, 2))), axis = 0)
            lhf = np.concatenate((lhf, fluxes_nc.variables[lhf_key][:,0,:,:].mean(axis = (1, 2))), axis = 0)

qT = mt/(1 + mt)
qv = mv/(1 + mv)
#thetal = theta - (Lv/cpd)*mcl

### Compute the vertically integrated budgets ###
=======

# Read the initiatial conditions
exp_names = ['Control']

for exp_name in exp_names:
    spinup_data = {}
    with open('../InitialFields_Temperature/InitialFields_Temperature_' + exp_name + '.txt', 'r') as init_temp:
        for line in init_temp:
            if 'num_theta_init_heights' in line:
                spinup_data['num_theta_init_heights'] = int(line.strip('num_theta_init_heights: '))
            elif 'theta_init_height' in line:
                spinup_data['theta_init_height'] = np.array([float(hgt) for hgt in line.strip('theta_init_height: ').split(',')])
            elif 'theta_init_data' in line:
                spinup_data['theta_init_data'] = np.array([float(data) for data in line.strip('theta_init_data: ').split(',')])
    
    with open('../InitialFields_Moisture/InitialFields_Moisture_' + exp_name + '.txt', 'r') as init_moist:
        for line in init_moist:
            if 'num_mv_init_heights' in line:
                spinup_data['num_mv_init_heights'] = int(line.strip('num_mv_init_heights: '))
            elif 'mv_init_height' in line:
                spinup_data['mv_init_height'] = np.array([float(hgt) for hgt in line.strip('mv_init_height: ').split(',')])
            elif 'mv_init_data' in line:
                spinup_data['mv_init_data'] = np.array([float(data) for data in line.strip('mv_init_data: ').split(',')])
    
    with open('../InitialFields_Wind/InitialFields_Wind_' + exp_name + '.txt', 'r') as init_wind:
        for line in init_wind:
            if 'num_uv_init_heights' in line:
                spinup_data['num_uv_init_heights'] = int(line.strip('num_uv_init_heights: '))
            elif 'uv_init_height' in line:
                spinup_data['uv_init_height'] = np.array([float(hgt) for hgt in line.strip('uv_init_height: ').split(',')])
            elif 'u_init_data' in line:
                spinup_data['u_init_data'] = np.array([float(data) for data in line.strip('u_init_data: ').split(',')])
            elif 'v_init_data' in line:
                spinup_data['v_init_data'] = np.array([float(data) for data in line.strip('v_init_data: ').split(',')])
    
    p_sfc = 101700. # the surface pressure used to initialise the simulations
    spinup_data['temperature'] = spinup_data['theta_init_data'] - g*spinup_data['theta_init_height']/cpd
    spinup_data['pressure']    = np.zeros_like(spinup_data['temperature'])
    spinup_data['pressure'][0] = p_sfc*1. # set the surface pressure
    
    dz = 1.
    z = 0
    temperature_interp = interpolate.interp1d(x = spinup_data['theta_init_height'], y = spinup_data['temperature'])
    for k in range(1, spinup_data['num_theta_init_heights']):
        rho_0 = spinup_data['pressure'][k-1]/(Rd*temperature_interp(z)) # get the air density at the level below
        p_temp = spinup_data['pressure'][k-1] - g*rho_0*dz
        z += dz
        
        # use that air density to compute the pressure slightly above and iterate until just below the next level
        while (z + dz) < spinup_data['theta_init_height'][k]:
            rho_0 = p_temp/(Rd*temperature_interp(z))
            p_temp = p_temp - g*rho_0*dz
            z += dz
        
        # do the remaining distance to get to the next level
        rho_0 = p_temp/(Rd*temperature_interp(z))
        spinup_data['pressure'][k] = p_temp - g*rho_0*(spinup_data['theta_init_height'][k] - z)
        z += (spinup_data['theta_init_height'][k] - z)
    
    # iterate to get a better temperature estimate
    spinup_data['temperature'] = PTtoTemp(spinup_data['theta_init_data'], spinup_data['pressure'], t_units = 'K', p_units = 'Pa')
    
    # iterate to get a better pressure estimate
    spinup_data['pressure'][0] = p_sfc*1. # set the surface pressure
    
    dz = 1.
    z = 0
    temperature_interp = interpolate.interp1d(x = spinup_data['theta_init_height'], y = spinup_data['temperature'])
    for k in range(1, spinup_data['num_theta_init_heights']):
        rho_0 = spinup_data['pressure'][k-1]/(Rd*temperature_interp(z)) # get the air density at the level below
        p_temp = spinup_data['pressure'][k-1] - g*rho_0*dz
        z += dz
        
        # use that air density to compute the pressure slightly above and iterate until just below the next level
        while (z + dz) < spinup_data['theta_init_height'][k]:
            rho_0 = p_temp/(Rd*temperature_interp(z))
            p_temp = p_temp - g*rho_0*dz
            z += dz
        
        # do the remaining distance to get to the next level
        rho_0 = p_temp/(Rd*temperature_interp(z))
        spinup_data['pressure'][k] = p_temp - g*rho_0*(spinup_data['theta_init_height'][k] - z)
        z += (spinup_data['theta_init_height'][k] - z)
    
    # get the dewpoint
    mv_interp = interpolate.interp1d(x = spinup_data['mv_init_height'], y = spinup_data['mv_init_data'], fill_value = 'extrapolate')
    spinup_data['q']        = getQ(spinup_data['temperature'], mv_interp(spinup_data['theta_init_height'])*100., spinup_data['pressure'], t_units = 'K', p_units = 'Pa')
    spinup_data['dewpoint'] = getDew(spinup_data['q'], spinup_data['pressure'], q_units = 'kg/kg', p_units = 'Pa')


### Compute the vertically integrated budgets

>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
# Define the idealised forcing profiles
Q_rad   = np.array([-2.0, -2.0, 0.0, 0.0])/86400.0 # K/day -> K/s
Q_rad_z = np.array([0.0, 3217.0, 4326.0, 40000.0])

w_subs   = np.array([0.0, -0.156, -0.379, -0.588, -0.603, -0.498, -0.384, -0.337, -0.315, -0.311, -0.291, -0.300, -0.418, -0.208, 0.0, 0.0])/100.0 # cm/s -> m/s
w_subs_z = np.array([0.0, 148.0, 392.0, 628.0, 857.0, 1119.0, 1381.0, 1644.0, 1906.0, 2168.0, 2430.0, 2693.0, 3217.0, 3867.0, 4326.0, 40000.0])

<<<<<<< HEAD
dz = 1.
# Interpolate everything onto the same grid
z = np.arange(0, 5000.1, dz)
z_half = np.arange(dz/2., 5000.1, dz)
=======
# Interpolate everything onto the same grid
z = np.arange(0, 40000.1, dz)
z_half = np.arange(dz/2., 40000.1, dz)

# initial variables
theta_fun = interpolate.interp1d(spinup_data['theta_init_height'], spinup_data['theta_init_data'], fill_value = 'extrapolate')
q_fun = interpolate.interp1d(spinup_data['mv_init_height'], spinup_data['q'], fill_value = 'extrapolate')
p_fun = interpolate.interp1d(spinup_data['theta_init_height'], spinup_data['pressure'], fill_value = 'extrapolate')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

# forcings
Q_rad_fun = interpolate.interp1d(Q_rad_z, Q_rad, fill_value = 'extrapolate')
w_subs_fun = interpolate.interp1d(w_subs_z, w_subs, fill_value = 'extrapolate')

<<<<<<< HEAD
T_sfc = 302.3
p_sfc = 101700.
p0 = 1e05
rho_star = p_sfc/(Rd*T_sfc)
# initialise surface fluxes
plt.wpthp = np.zeros((theta.shape[0]))
wpqtp = np.zeros_like(wpthp)
thetav = theta*(1+0.608*qv)
for it in range(theta.shape[0]):
    # initial variables
    theta_fun = interpolate.interp1d(z_theta, theta[it,:], fill_value = 'extrapolate')
    q_fun = interpolate.interp1d(z_theta, qv[it,:], fill_value = 'extrapolate')
    p_fun = interpolate.interp1d(0.5*(z_theta[1:] + z_theta[:-1]), pressure[it,:], fill_value = 'extrapolate')
    
    rho = p_fun(z_half)/(Rd*theta_fun(z_half)/((p0/p_fun(z_half))**(Rd/cpd)))
    exner = (p_fun(z_half)/p0)**(Rd/cpd)
    
    # Compute the gradients
    dt_dz = (theta_fun(z[1:]) - theta_fun(z[:-1]))/(z[1:] - z[:-1])
    dq_dz = (q_fun(z[1:]) - q_fun(z[:-1]))/(z[1:] - z[:-1])
    
    # subsidence warming
    Q_subs = - rho*exner*w_subs_fun(z_half)*dt_dz
    
    # subsidence drying
    drying = - rho*w_subs_fun(z_half)*dq_dz
    
    # Integrate to compute the surface fluxes which balance those tendencies
    wpthp[it] = cpd*(-integrate.trapz(x = z_half, y = rho*exner*Q_rad_fun(z_half) + Q_subs))
    wpqtp[it] = - Lv*integrate.trapz(x = z_half, y = drying)

fig = plt.figure()
axa = fig.add_subplot(1, 2, 1)
axa.plot(np.arange(20, 14400.1, 10.), wpthp, lw = 2, color = 'r')
axa.plot(np.arange(20., 14400.1, 10.), shf, 'r--')

axb = fig.add_subplot(1, 2, 2)
axb.plot(np.arange(20, 14400.1, 10.), wpqtp, lw = 2, color = 'b')
axb.plot(np.arange(20, 14400.1, 10.), Lv*lhf, 'b--')
print wpthp[-1]
print wpqtp[-1]
plt.show()

# check the precip over the last four days
for day in ['07', '08', '09', '10']:
    with Dataset('/work/n02/n02/xb899100/cylc-run/u-bc341/share/data/history/lwp_' + day + '.nc', 'r') as lwp_nc:
        if day == '07':
            rain = lwp_nc.variables[ls_rain_amt_key][:]*1.
        else:
            rain = np.concatenate((rain, lwp_nc.variables[ls_rain_amt_key][:]*1.), axis = 0)

print Lv*np.nanmean(rain - rain[0,:,:])/(4*86400.)
=======
# Compute the gradients
dt_dz = (theta_fun(z[1:]) - theta_fun(z[:-1]))/(z[1:] - z[:-1])
dq_dz = (q_fun(z[1:]) - q_fun(z[:-1]))/(z[1:] - z[:-1])

# subsidence warming
Q_subs = - w_subs_fun(z_half)*dt_dz

# subsidence drying
drying = - w_subs_fun(z_half)*dq_dz

# Integrate to compute the surface fluxes which balance those tendencies
T_sfc = 302.3
rho_star = p_sfc/(Rd*T_sfc)
exner = (p_fun(0)/(p0*100.0))**(Rd/cpd)
wpthp = exner*rho_star*cpd*(-integrate.cumtrapz(x = z_half, y = (Q_rad_fun(z_half) + Q_subs)))
wpqp  = -rho_star*Lv*integrate.cumtrapz(x = z_half, y = drying)

fig = plt.figure()
axa = fig.add_subplot(1, 2, 1)
axa.plot(-wpqp + wpqp[-1], z[1:-1])
axa.set_ylim([0, 4000])

axb = fig.add_subplot(1, 2, 2)
axb.plot(-wpthp + wpthp[-1], z[1:-1])
axb.set_ylim([0, 4000])
plt.show()

print wpthp[-1]
print wpqp[-1]
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

