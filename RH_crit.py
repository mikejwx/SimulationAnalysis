"""
Code to estimate the RH_crit profile that would be diagnosed in our simulations.
"""
import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import getQ as q_sat
from SkewT_archer import Lv, cpd, Rv, Rd, temp2theta, EPS, g
from netCDF4 import Dataset
from STASH_keys import theta_key, pthe_key, temp_key, u_key, v_key, mv_key, mcl_key, w_key

# Read in the initial data from the Control simulation
bouy_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control/bouy_00.nc', 'r')
fluxes_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control/fluxes_00.nc', 'r')
wind_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control/wind_00.nc', 'r')
mr_nc = Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control/mr_00.nc', 'r')

# Read ni the variables needed
theta = bouy_nc.variables[theta_key][-1,:-1,:,:]
T     = bouy_nc.variables[temp_key][-1,:-1,:,:]
p     = fluxes_nc.variables[pthe_key][-1,:-1,:,:]
mv    = mr_nc.variables[mv_key][-1,:-1,:,:]
q     = mv/(1. + mv)
mcl   = mr_nc.variables[mcl_key][-1,:-1,:,:]
qcl   = mcl /(1. + mcl)
u     = wind_nc.variables[u_key][-1,:-1,:,:]
v     = wind_nc.variables[v_key][-1,:-1,:,:]
up    = np.transpose(np.transpose(wind_nc.variables[u_key][-1,:-1,:,:]) - np.nanmean(wind_nc.variables[u_key][-1,:-1,:,:], axis = (1, 2)))
vp    = np.transpose(np.transpose(wind_nc.variables[v_key][-1,:-1,:,:]) - np.nanmean(wind_nc.variables[v_key][-1,:-1,:,:], axis = (1, 2)))
wp    = np.transpose(np.transpose(wind_nc.variables[w_key][-1,:-1,:,:]) - np.nanmean(wind_nc.variables[w_key][-1,:-1,:,:], axis = (1, 2)))
E     = 0.5*np.sqrt(up**2 + vp**2 + wp**2)

z_the = bouy_nc.variables['thlev_zsea_theta'][:-1]

# Close the netCDF
bouy_nc.close()
fluxes_nc.close()
wind_nc.close()
mr_nc.close()

# Define an approximate clausius clapeyron equation
def CC(temperature, pressure):
    """
    Clausius-Clapeyron...
    """
    qs = q_sat(temperature, [100.], pressure)
    CC = EPS*Lv*qs/(Rd*temperature*temperature)
    
    return CC

# Define an approximate gradient function
def gradient(var, z):
    return np.transpose(np.transpose(var[1:,:,:] - var[:-1,:,:])/(z[1:] - z[:-1]))

# Define an approximate interpolation
def interp(var):
    return 0.5*(var[1:] + var[:-1])

# Define the stability function
def Ri(temperature, theta_l, qT):
    beta_T = interp(1./temperature)
    cv = (1./EPS) - 1.
    V_fac = interp(1. + cv*q - qcl)
    beta_q = cv/V_fac
    deltaB = g*(beta_T*np.transpose(np.transpose(gradient(theta_l, z_the))*(z_the[1:] - z_the[:-1])) + beta_q*np.transpose(np.transpose(gradient(qT, z_the))*(z_the[1:] - z_the[:-1])))
    return (np.transpose(np.transpose(deltaB)/(z_the[1:] - z_the[:-1])))/(Shear**2)

def fx(Ri, z):
    Pr_N = 0.7
    g0 = 10.
    Dm = g0/4.
    Dh = g0/25.
    k = 0.4
    z0m = 0.0002
    fx = np.where(Ri < 0, (1./Pr_N)*(1. - g0*Ri/(1. + Dh*(np.abs(Ri)**(0.5)))), 1./(Pr_N*(1. + g0*Ri)))
    return fx

def stability(Ri, z):
    Kh = np.transpose((l**2)*np.transpose(Shear*fx(Ri, z_the)))
    np.where(Kh < W1D*75., W1D*75., Kh)
    return Kh/np.transpose(l*np.transpose(interp(np.sqrt(2.*E))))

# Compute any additionals
T_L     = T - Lv*qcl/cpd
theta_l = np.transpose(np.transpose(T_L) + g*z_the/cpd)
qT      = q + qcl

a = (1. + (Lv/cpd)*CC(T_L, p))**(-1.)
b = a*(T/theta)*CC(T_L, p)

zh = 700.
W1D = 0.29
z0m = 0.0002
cs = 0.2
dx = 100.
l_bl = ((1./(0.4*(z_the + z0m))) + (1./np.max([40., 0.15*zh])))**(-1)
l_smag = np.sqrt((1./((cs*dx)**2) + 1./((0.4*(z_the + z0m))**2))**(-1))
l = interp(W1D*l_bl + (1. - W1D)*l_smag)

# Compute the vertical gradients
Shear     = np.sqrt(gradient(u, z_the)**2 + gradient(v, z_the)**2)

Sh        = stability(Ri(T, theta_l, qT), z_the)
dthetaldz = gradient(theta_l, z_the)
dqTdz     = gradient(qT, z_the)

thetalp2   = np.transpose(np.transpose(15.*Sh*(dthetaldz**2))*(l**2))
qTp2       = np.transpose(np.transpose(15.*Sh*(dqTdz**2))*(l**2))
qTpthetalp = np.transpose(np.transpose(15.*Sh*(dthetaldz*dqTdz))*(l**2))

# Interpolate a and b
a = interp(a)
b = interp(b)

# Calculate the half-width
sigma_s2 = (a**2.)*qTp2 - 2*a*b*qTpthetalp + (b**2.)*thetalp2
sigma_s = np.sqrt(sigma_s2)

plt.plot(np.nanmean(thetalp2, axis = (1, 2)), interp(z_the), label = '$\overline{\\theta_{l}^{\prime 2}}$')
plt.plot(np.nanmean(qTp2, axis = (1, 2)), interp(z_the), label = '$\overline{q_{T}^{\prime 2}}$')
plt.plot(np.nanmean(qTpthetalp, axis = (1, 2)), interp(z_the), label = '$\overline{q_{T}^{\prime} \\theta_{l}^{\prime}}$')
plt.legend(loc = 0)
plt.show()


plt.plot(np.nanmean((a**2.)*qTp2, axis = (1, 2)), interp(z_the), label = 'Term 1')
plt.plot(np.nanmean(- 2*a*b*qTpthetalp, axis = (1, 2)), interp(z_the),label = 'Term 2')
plt.plot(np.nanmean((b**2.)*thetalp2, axis = (1, 2)), interp(z_the), label = 'Term 3')
plt.plot(np.nanmean(sigma_s2, axis = (1, 2)), interp(z_the), 'k--', lw = 2, label = '$\sigma ^{2}$')
plt.legend(loc = 0)
plt.show()

plt.plot(np.nanmean(a, axis = (1, 2)), interp(z_the), label = 'a')
plt.plot(np.nanmean(b, axis = (1, 2)), interp(z_the), label = 'b')
plt.legend(loc = 0)
plt.show()

RH_crit = 1. - np.sqrt(6.)*sigma_s/(a*interp(q_sat(T_L, [100.], p, t_units = 'K', p_units = 'Pa')))

RH_C = np.nanmin(RH_crit*100., axis = (1, 2))
plt.plot(RH_C, interp(z_the))
plt.xlim([-5., 105.])
plt.ylim([0, 15000])
plt.show()

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, adjustable = 'box', aspect = 1)
my_plt = ax.contourf(RH_crit[iz,:,:]*100., levels = np.arange(0., 100.1, 5.), extend = 'max', cmap = 'Reds')
fig.colorbar(my_plt, ax = ax)
for iz in xrange(1, len(interp(z_the))):
    ax.cla()
    ax.contourf(RH_crit[iz,:,:]*100., levels = np.arange(0., 100.1, 5.), extend = 'max', cmap = 'Reds')
    ax.set_title(str(int(interp(z_the)[iz])))
    plt.pause(0.1)


