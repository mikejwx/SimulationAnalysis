import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from STASH_keys import lwp_key
from scipy import ndimage
"""
Assess organisation in our simulations using an organisation index:
SCAI (Tobin et al., 2012)

This index depends on the number of convective objects and the degree of 
clumping into clusters of convective objects.
"""
path = '/nerc/n02/n02/xb899100/CloudTrail/Control2/'
with Dataset(path + 'lwp_00.nc', 'r') as lwp_nc:
    lwp_data = lwp_nc.variables[lwp_key][:]*1.

# manufacture coordinate system
x, y = np.meshgrid(np.arange(lwp_data.shape[2])*0.8, np.arange(lwp_data.shape[1])*0.8)
L = np.sqrt(x.max()**2 + y.max()**2)
a = 0.8
Nmax = (L/a)**2

cloud_mask = np.where(lwp_data > 0, 1.0, 0.0)
SCAI = []
time = []
ns = []
d0s = []
for it in range(0, cloud_mask.shape[0], 15):
    # Number of convective clusters
    conv_obj, N = ndimage.label(cloud_mask[it,:,:])
    # N depends on domain size, the spatial resolution of the data, and the threshold used
    # N can be normalised by a potential maximum of number 'Nmax' proportional to the ratio (L/a)^2
    # where L is the length of the domain and a is the grid length
    
    # clumping of the clusters
    ds = []
    for obj_noi in range(1, N+1, 30):
        # compute the centre of mass of the object based on lwp
        # for object i
        x_ci = np.nansum(np.where(conv_obj == obj_noi, x*lwp_data[it,:,:], 0.0))/np.nansum(np.where(conv_obj == obj_noi, lwp_data[it,:,:], 0.0))
        y_ci = np.nansum(np.where(conv_obj == obj_noi, y*lwp_data[it,:,:], 0.0))/np.nansum(np.where(conv_obj == obj_noi, lwp_data[it,:,:], 0.0))
        for obj_noj in range(1, N+1, 30):
            if obj_noj != obj_noi:
                # i.e. we don't want to include the distance between the object in itself as this is by definition = 0
                # compute the centre of mass of the object based on lwp
                # for object j
                x_cj = np.nansum(np.where(conv_obj == obj_noj, x*lwp_data[it,:,:], 0.0))/np.nansum(np.where(conv_obj == obj_noj, lwp_data[it,:,:], 0.0))
                y_cj = np.nansum(np.where(conv_obj == obj_noj, y*lwp_data[it,:,:], 0.0))/np.nansum(np.where(conv_obj == obj_noj, lwp_data[it,:,:], 0.0))
                # distance between those two centres of mass
                d = np.sqrt((x_ci - x_cj)**2 + (y_ci - y_cj)**2)
                ds.append(d)
    
    ds = np.array(ds)
    N -= 1
    n = N*(N-1.)/2.
    D0 = np.exp(np.log(ds).sum()/len(ds))
    D1 = np.mean(ds)
    SCAI.append((N/Nmax)*(D0/L)*1000)
    ns.append(N)
    d0s.append(D0)
    time.append(it)

# D0 is the geometric mean of the distance between convective objects
# D1 is the arithmetic mean of the distance between convective objects
# The centres of mass of convective objects are used to compute their distances of separation

plt.scatter(ns, d0s)
plt.show()

plt.plot(time, SCAI)
plt.show()

