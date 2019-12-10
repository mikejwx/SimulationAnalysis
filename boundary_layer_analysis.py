import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import *
from scipy import interpolate
from scipy.spatial.distance import pdist, squareform
from datetime import datetime as dt
from random import sample
"""
Tackling some interesting points about how the representation of boundary layer 
eddies changes with increasing grid spacing.

i)   the separation lengthscale in boundary layer eddies.
ii)  the mean flux profiles (w'qT') and (w'thetal')
iii) the number of updraughts as a function of height...
iv)  turbulence spectra? 
"""

################################################################################
### Read data ###
with Dataset('/nerc/n02/n02/xb899100/CloudTrail/Control2/bouy_09.nc', 'r') as wind_nc:
    w_zi = wind_nc.variables[theta_key][-1,:,:,:]*1.
    z = wind_nc.variables['thlev_zsea_theta'][:]*1.

iz = np.where(np.abs(z - 350) == np.min(np.abs(z - 350)))[0][0]
w_zi = w_zi[0,:,:]

################################################################################
### Spatial autocorrelation lengthscale between boundary layer eddies ###

# cartesian coordinate system
x, y = np.meshgrid(np.arange(w_zi.shape[1])*100., np.arange(w_zi.shape[0])*100.)

def spatial_corr_lengthscale(x, y, in_field, max_dist, my_size = 2000):
    """
    Compute the spatial autocorrelation lengthscale as the smallest spatial 
    distance at which the spatial autocorrelation drops below e^-1, averaged
    for all rows and columns.
    Requires the x and y coordinates (for cartesian data) and the 2D field of 
    the data itself.
    """
    # assume isotropic horizontal grid, i.e. dx = dy
    dx = (x[0,1] - x[0,0])
    # range of distances starting at one grid length away
    distances = np.arange(0, max_dist + dx/2., dx)
    # will take literally 12 hours to do the whole scene, so let's break it down
    # initialize an array to contain the data
    my_data = np.zeros((distances.size, my_size))
    f_flat = in_field.flatten()
    # take our random samples from f_flat
    random_samples = sample(list(enumerate(f_flat)),my_size)
    # the above should return a list of tuples [(idx, sample), (idx, sample), ...]
    my_idx = [my_sample[0] for my_sample in random_samples]
    my_data[0,:] = np.array([my_sample[1] for my_sample in random_samples])
    x_flat = x.flatten()
    y_flat = y.flatten()
    for i in my_idx:
        R = np.sqrt((x - x_flat[i])**2 + (y - y_flat[i])**2).flatten()
        my_data[:,my_idx.index(i)] = np.array([f_flat[((dist-dx/2) < R)*(R < (dist+dx/2))].mean() for dist in distances])
    
    rho_c = np.array([np.corrcoef(my_data[0,:], my_data[j,:])[0,1] for j in range(distances.size)])
    
    return distances[np.min(np.where(rho_c < np.exp(-1))[0])]

### Flux profiles ###

### Updraught numbers ###

### Turbulence Spectrum ###

# rename numpy's 'fast fourier transform shift' function to 'shift'
shift = np.fft.fftshift
# get domain dimentions
nN, nE = w_zi.shape
dN, dE = 100., 100.
# get the spectrum by performing a 2D fft and then the fftshift on that
spec = shift(np.fft.fft2(w_zi))

# get the dimensions in inverse distance space
kE = np.fft.fftfreq(nE, dE)
kN = np.fft.fftfreq(nN, dN)
k = kN if kN.size < kE.size else kE
k = k[k>0]
# compute the distance from the corners of the domain in inverse distance space
k_rad = np.sqrt(kN[:,np.newaxis]**2 + kE[np.newaxis,:]**2)
# get the power spectrum as the square of the spectrum
pspec = np.abs(spec)**2
# set to zero in the corners of the domain because it causes badness(?)
#pspec[k_rad == 0.] = 0.

# shift the dimensions in inverse distance space to get back to the original coordinates
kE = shift(kE)
kN = shift(kN)
power_interp = interpolate.RectBivariateSpline(kN, kE, pspec)

def power2DMean(k, N=1024):
    """
    Computes the 2D mean power spectrum
    Takes the wavenumbers for which we want to compute the power, k, and
    considers all the points at that radius to compute the amplitude.
    
    N is the azimuthal resolution, 360 is ~two steps every degree
    """
    theta = np.linspace(-np.pi, np.pi, N, False)
    power = np.empty_like(k)
    for i in xrange(k.size):
        kE = np.sin(theta) * k[i]
        kN = np.cos(theta) * k[i]
        power[i] = np.mean(power_interp.ev(kN, kE))
    
    return 4*np.pi*power/pspec.size

plt.loglog(k, power2DMean(k))
plt.show()

