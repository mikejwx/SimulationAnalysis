"""
For use in extending the BOMEX initial conditions up to 40 km. The original
BOMEX initial conditions end at 3 km. We use the original conditions below 3 km
and then assume a pseudoadiabatic ascent above that to 12 km, followed by an 
isothermal layer to 15 km and increasing theta to 1400 K at 40 km.
"""

import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import getQ, getGM, Rd, cpd, p0, g, EPS
from scipy import interpolate, integrate

# BOMEX ends at z = 3000 m with theta = 311.85 K. With a moist hydrostatic
# balance the pressure at this height is about 705 hPa.
z = np.array([0., 520., 1480., 2000., 3000.]) # heights
theta = np.array([298.7, 298.7, 302.4, 308.2, 311.85]) # Potential Temperature profile
mv = np.array([17.3, 16.6, 10.8, 4.2, 3.0])*1e-03 # Mixing ratio profile
p_sfc = 101500. # Pa

p = np.zeros_like(z)
p[0] = p_sfc*1.
for k in xrange(1, len(z)):
    dz = z[k] - z[k-1]
    Tv = (theta[k-1]/((p0*100./p[k-1])**(Rd/cpd)))*(1. + 0.608*mv[k-1])
    rho = p[k-1]/(Rd*Tv)
    p[k] = p[k-1] - g*rho*dz

### Extended heights ###
z_ext = np.array([4000., 5500., 7000., 9500., 12000., 15000., 18000., 21000., 27000., 33000., 40000.])

# Assume a constant relative humidity above 3 km
# Get the first point
# Calculate the distance between levels
dz = 1.
# Calculate the temperature at the lower level
T = (theta[-1]/((p0*100./p[-1])**(Rd/cpd)))
# Calculate the virtual temperature at the lower level
Tv = T*(1. + 0.608*mv[-1])
# Calculate the air density at the lower level
rho = p[-1]/(Rd*Tv)
# Use hydrostatic balance and air density at the lower level to calculate the pressure at the upper level
p_ext = [p[-1] - g*rho*dz]

# Calculate the temperature at the upper level
T_new = T - getGM(T, 0.5*(p[-1]+p_ext[0]), t_units = 'K', p_units = 'Pa')*dz
# Calculate the specific humidity at the upper level
RH_const = 100.*mv[-1]/getQ(T, 100., p[-1], t_units = 'K', p_units = 'Pa')
q = getQ(T_new, RH_const, p_ext[0], t_units = 'K', p_units = 'Pa')[0]
# Convert specific humidity to mixing ratio
mv_ext = [q/(1. - q)]
# Calculate the new potential temperature at the upper level
theta_ext = [T_new*(p0*100./p_ext[0])**(Rd/cpd)]
z_ext_integration = [z.max() + dz]

### repeat for the remaining levels
for k in xrange(1, int((z_ext.max()-(z.max()+dz))/dz)):
    T = (theta_ext[k-1]/((p0*100./p_ext[k-1])**(Rd/cpd)))
    Tv = T*(1. + 0.608*mv_ext[k-1])
    rho = p_ext[k-1]/(Rd*Tv)
    p_ext.append(p_ext[k-1] - g*rho*dz)
    
    T_new = T - getGM(T, 0.5*(p_ext[k-1]+p_ext[k]), t_units = 'K', p_units = 'Pa')*dz
    q = getQ(T_new, RH_const, 0.5*(p_ext[k-1]+p_ext[k]), t_units = 'K', p_units = 'Pa')[0]
    mv_ext.append(q/(1. - q))
    theta_ext.append(T_new*(p0*100./p_ext[k])**(Rd/cpd))
    z_ext_integration.append(z_ext_integration[-1] + dz)

theta_ext = np.array(theta_ext)
mv_ext = np.array(mv_ext)
p_ext = np.array(p_ext)
z_ext_integration = np.array(z_ext_integration)

print 'p = '
print interpolate.interp1d(x = z_ext_integration, y = p_ext, fill_value = 'extrapolate')(z_ext)

print 'Theta = '
extended_theta = interpolate.interp1d(x = z_ext_integration, y = theta_ext, fill_value = 'extrapolate')(z_ext)

# fill in isothermal from 12000 to 15000 m then increasing gradient by 10% to model top
for k in range(5, len(z_ext)):
    extended_theta[k] = extended_theta[k-1] + (1. + 0.1*(k-5))*g*(z_ext[k] - z_ext[k-1])/cpd
print extended_theta

print 'm_v = '
print interpolate.interp1d(x = z_ext_integration, y = mv_ext, fill_value = 'extrapolate')(z_ext)





