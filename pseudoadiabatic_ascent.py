"""
For use in extending the BOMEX initial conditions up to 40 km. The original
BOMEX initial conditions end at 3 km. We use the original conditions below 3 km
and then assume a pseudoadiabatic ascent above that to 12 km, followed by an 
isothermal layer to 15 km and increasing theta to 1400 K at 40 km.
"""

import numpy as np
import matplotlib.pyplot as plt
from SkewT_archer import getGM, Rd, cpd, p0, g, EPS
from scipy import interpolate, integrate

# BOMEX ends at z = 3000 m with theta = 311.85 K. With a moist hydrostatic
# balance the pressure at this height is about 707 hPa.
z = np.array([0., 520., 1480., 2000., 3000.]) # heights
theta = np.array([298.7, 298.7, 302.4, 308.2, 311.85]) # Potential Temperature profile
mv = np.array([17.3, 16.6, 10.8, 4.2, 3.0])*1e-03 # Mixing ratio profile
p_sfc = 101500. # Pa

p = np.zeros_like(z)
p[0] = p_sfc*1.
for k in xrange(1, len(z)):
    dz = z[k] - z[k-1]
    T = (theta[k-1]/((p0*100./p[k-1])**(Rd/cpd)))*(1. + 0.608*mv[k-1])
    rho = p[k-1]/(Rd*T)
    p[k] = p[k-1] - g*rho*dz

print p
### Extended heights ###
z_ext = np.array([4000., 5500., 7000., 9500., 12000., 15000., 18000., 21000., 27000., 33000., 40000.])
theta_ext = np.zeros_like(z_ext)
p_ext     = np.zeros_like(z_ext)
mv_ext    = np.zeros_like(z_ext)

# Assume a constant relative humidity above 3 km
# Get the first point
(theta[-1]/((p0*100./p[-1])**(Rd/cpd)))*(1. + 0.608*mv[-1])
