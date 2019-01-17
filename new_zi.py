"""
Attempting to find the boundary layer height using method from 
Vogelezang, D. H. P. and Holtslag A. A. M. 1996, equation 2

Ri_g = (g/theta_vs)(theta_vh - theta_vs)(h-z_s)/((u_h - u_s)^2 + (v_h - v_s)^2)
Looking for where Ri_g = 0.3
h = boundary layer top, and variable with subscript h denote variables observed at height h
z_s = depth of surface layer (?), variables with subscript s are at observed at this height z_s
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

def find_h(theta_v, u, v, z):
    """
    Iteratively solves equation 2 from Vogelezang and Holtslag (1996) to
    estimate the boundary layer depth.
    
    This method uses a modified Richardson number approach and vertical profiles
    of virtual potential temperature (theta_v), u- and v- wind components to
    estimate the boundary layer height.
    """
    # Requires numpy arrays
    import numpy as np
    # Requires interpolation routines from scipy
    from scipy import interpolate
    
    # Constants
    g = 9.81
    
    # Need a first guess for the height of the boundary layer.
    h   = np.arange(1.0, 6000.1, 1.0) # m
    z_s = 0.1*h
    
    theta_v = interpolate.interp1d(x = z, y = theta_v, fill_value = 'extrapolate')
    u       = interpolate.interp1d(x = z, y = u, fill_value = 'extrapolate')
    v       = interpolate.interp1d(x = z, y = v, fill_value = 'extrapolate')
    
    Ri_g = (g/theta_v(z_s))*(theta_v(h) - theta_v(z_s))*(h - z_s)/((u(h) - u(z_s))**2. + (v(h) - v(z_s))**2.)
    
    z_i = h[np.where(np.min(np.abs(Ri_g - 0.30)) == np.abs(Ri_g - 0.30))[0][0]]
    
    return z_i

# test it out

u     = np.array([5.91, 5.62, 5.61, 6.75, 8.28, 8, 7.86, 6.94, 7.09, 7.37, 7.27, 7.27, 7.88, 8.32, 9.16, 9.63, 10.08, 10.15, 9.17, 7.91, 7.81, 7.49])
v     = np.array([4.96, 6.70, 8.01, 7.76, 6.94, 8, 8.14, 8.28, 8.15, 7.90, 7.27, 7.27, 6.61, 6.05, 5.72, 4.90, 3.87, 3.69, 4.67, 5.74, 5.88, 6.28])
z     = np.array([37, 184, 344, 515, 791, 866, 904, 1190, 1307, 1396, 1575, 1595, 1778, 1964, 2110, 2313, 2542, 2597, 2821, 3085, 3120, 3225])

theta_v = np.array([302.7, 302.1, 302.1, 302.1, 303.1, 303.4, 303.5, 304.3, 304.6, 305.1, 305.8, 305.9, 306.7, 307.6, 307.9, 308.6, 310.0, 310.0, 310.0, 311.1, 311.5, 311.8])

zi = find_h(theta_v, u, v, z)

print zi
fig = plt.figure()
ax = fig.add_subplot(1, 2, 1)
ax.plot(theta_v, z, 'r', lw = 2)
ax.plot([np.min(theta_v)-1., np.max(theta_v)+1.], [zi, zi], 'k--')
ax.set_xlim([np.min(theta_v)-1., np.max(theta_v)+1.])
ax = fig.add_subplot(1, 2, 2)
ax.plot(u, z, 'k', lw = 2)
ax.plot(v, z, 'k:', lw = 2)
ax.plot([np.min([u, v])-1., np.max([u, v])+1.], [zi, zi], 'k--')
ax.set_xlim([np.min([u, v])-1., np.max([u, v])+1.])
plt.show()

