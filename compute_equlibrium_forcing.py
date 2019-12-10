import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
from SkewT_archer import getQ, g, cpd, Rd, Lv, p0

# Define the simplified initial conditions
theta_init = np.array([299.0, 299.0, 299.01, 302.47, 311.38, 318.05, 323.86, 328.95, 335.98, 339.10, 342.42, 354.57, 402.01, 1200.0])
theta_init_z = np.array([0.0, 148.0, 509.0, 857.0, 3217.0, 3840.0, 5931.0, 7635.0, 9714.0, 10963.0, 12408.0, 14189.0, 16649.0, 40000.0])

RH_init = np.array([0.775, 0.814, 0.890, 0.554, 0.386, 0.200, 0.200])
RH_init_z = np.array([0.0, 148.0, 509.0, 857.0, 3217.0, 3509.0, 40000.0])

p_sfc = 101700.0 # Pa
pressure_init = [p_sfc]
# Use the ideal gas law for dry air and the hydrostatic equation to compute p
temperature_init = theta_init - (g/cpd)*theta_init_z
temperature_fun = interpolate.interp1d(x = theta_init_z, y = temperature_init, fill_value = 'extrapolate')
dz = 1
for z in range(1, 40001, dz):
    rho = pressure_init[-1]/(Rd*temperature_fun(z-0.5*dz)) # ideal gas law
    pressure_init.append(pressure_init[-1] - rho*g*dz) # hydrostatic balance

pressure_init_theta = np.array([pressure_init[i] for i in range(len(pressure_init)) if i in theta_init_z])
pressure_init_RH    = np.array([pressure_init[i] for i in range(len(pressure_init)) if i in RH_init_z])

# Use getQ to convert from RH to q
q_init = getQ(temperature_fun(RH_init_z), RH_init*100.0, pressure_init_RH, t_units = 'K', p_units = 'Pa')

# Define the idealised forcing profiles
Q_rad   = np.array([-2.0, -2.0, 0.0, 0.0])/86400.0 # K/day -> K/s
Q_rad_z = np.array([0.0, 3217.0, 4326.0, 40000.0])

w_subs   = np.array([0.0, -0.156, -0.379, -0.588, -0.603, -0.498, -0.384, -0.337, -0.315, -0.311, -0.291, -0.300, -0.418, -0.208, 0.0, 0.0])/100.0 # cm/s -> m/s
w_subs_z = np.array([0.0, 148.0, 392.0, 628.0, 857.0, 1119.0, 1381.0, 1644.0, 1906.0, 2168.0, 2430.0, 2693.0, 3217.0, 3867.0, 4326.0, 40000.0])

# Interpolate everything onto the same grid
z = np.arange(0, 40000.1, dz)
z_half = np.arange(dz/2., 40000.1, dz)

# initial variables
theta_fun = interpolate.interp1d(theta_init_z, theta_init, fill_value = 'extrapolate')
q_fun = interpolate.interp1d(RH_init_z, q_init, fill_value = 'extrapolate')
p_fun = interpolate.interp1d(np.arange(0.0, 40000.1, 1.0), pressure_init, fill_value = 'extrapolate')

# forcings
Q_rad_fun = interpolate.interp1d(Q_rad_z, Q_rad, fill_value = 'extrapolate')
w_subs_fun = interpolate.interp1d(w_subs_z, w_subs, fill_value = 'extrapolate')

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
wpthp = exner*rho_star*cpd*(-integrate.trapz(x = z_half, y = (Q_rad_fun(z_half) + Q_subs)))
wpqp  = rho_star*Lv*integrate.trapz(z_half, drying)

print wpthp
print wpqp

fig = plt.figure(figsize = (12, 6), tight_layout = True)
axa = fig.add_subplot(1, 2, 1)
axa.plot(86400*Q_rad_fun(z_half), z_half/1000.0, 'k', lw = 0.5, label = u'$(\\frac{\partial \\theta}{\partial t})_{rad}$')
axa.plot(86400*Q_subs, z_half/1000.0, 'k', lw = 2, label = u'$ - w_{subs} \\frac{\partial \\theta}{\partial z}$')
axa.set_title(str(round(wpthp, 4)) + u' W m$^{-2}$')
axa.set_ylim([0, 5])
axa.set_xlabel(u'$\\frac{\partial \\theta}{\partial t}$ (K day$^{-1}$)')
axa.set_ylabel('Height (km)')
axa.legend(loc = 0, frameon = 0)
axa.set_xlim([-2.5, 6])

axb = fig.add_subplot(1, 2, 2)
axb.plot(86400*drying, z_half/1000.0, 'k', lw = 2, label = u'$ - w_{subs} \\frac{\partial q_{v}}{\partial z}$')
axb.set_title(str(round(wpqp, 4)) + u' W m$^{-2}$')
axb.set_ylim([0, 5])
axb.set_yticklabels([''])
axb.set_xlabel(u'$\\frac{\partial q}{\partial t}$ (kg kg$^{-1}$ day$^{-1}$)')
axb.legend(loc = 0, frameon = 0)
axb.set_xlim([-0.008, 0.001])

plt.savefig('../equilibrium_forcing.png', dpi = 150, bbox_inches = 'tight')
plt.show()

