import numpy as np
from multiprocessing import Pool

def bilinear_interpolation_new(x_in, y_in, z_in, x_out, y_out, kind = 0, d = 2000.0, p = 0.5, operation = np.nanmean):
    """
    Performs bilinear interpolation on a given 2D array, z_in.
    ----------------------------------------------------------------------------
    INPUT:
    x_in = 2D array of x-coordinates, input
    y_in = 2D array of y-coordinates, input
    z_in = 3D array of data points, input z_in[nz, ny, nx]

    x_out = one x-coordinate to be interpolated onto
    y_out = one y-coordinate to be interpolated onto
    
    kind = option for how to interpolate
    option 0: nearest neighbor
    option 1: inverse distance weighted within radius d to power p
    option 2: bilinear interpolation
    option 3: mean within radius d
    
    OUTPUT:
    z_out = 2D array, output z_out[nz, nr] where nr is the length of x_out
    e.g. nz could be the number of height levels, and nr could be the number of 
    downwind points
    
    interpolation scheme finds the nearest points and uses scipy's bisplrep to 
    do the interpolation.
    """
    import numpy as np
    # make use of the bi-periodic boundary conditions to not truncate at boundaries
    x_in_p = np.concatenate((x_in-np.max(x_in), x_in, x_in+np.max(x_in)), axis = 1)
    x_in_p = np.concatenate((x_in_p, x_in_p, x_in_p), axis = 0)
    y_in_p = np.concatenate((y_in-np.max(y_in), y_in, y_in+np.max(y_in)), axis = 0)
    y_in_p = np.concatenate((y_in_p, y_in_p, y_in_p), axis = 1)
    # Initialise our output array
    z_out = np.zeros((z_in.shape[0], len(x_out)))
    
    # For each point to be interpolated onto
    # Find the distance to the input data coordinates
    r = np.sqrt((x_in_p - x_out)**2 + (y_in_p - y_out)**2)
    
    if kind == 0:
        # Nearest neighbor approach
        iy, ix = np.where(r == np.min(r))
        z_out = z_in[:,iy[0]%z_in.shape[1],ix[0]%z_in.shape[2]]
    elif kind == 1:
        # idw approach
        # Use all points within 'd' m of x_out, y_out for weighted mean
        iy, ix = np.where(r < d)
        w = 0
        
        for j in xrange(len(iy)):
            w += (1./r[iy[j], ix[j]]**p)
            z_out += (1./r[iy[j], ix[j]]**p)*z_in[:,iy[j]%z_in.shape[1], ix[j]%z_in.shape[2]]
        z_out /= w
    elif kind == 2:
        # Find the nearest point
        iy, ix = np.where(r == np.min(r))
        iy = iy[0]%z_in.shape[1]
        ix = ix[0]%z_in.shape[2]
        # Determine which quadrant (i.e. up and left, up and right, down and
        # right, or down and left) this point is with respect to the output 
        # point
        dx = x_in[iy, ix] - x_out# if +ve, input point is to the right
        dy = y_in[iy, ix] - y_out # if +ve, input point is up
        if dx > 0:
            # nearest input point is to the right
            if dy > 0:
                # nearest input point is up
                z_out_up   = (z_in[:,iy,ix] - z_in[:,iy,ix-1])*(x_out - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + z_in[:,iy,ix-1]
                z_out_down = (z_in[:,iy-1,ix] - z_in[:,iy-1,ix-1])*(x_out - x_in[iy-1,ix-1])/(x_in[iy-1,ix] - x_in[iy-1,ix-1]) + z_in[:,iy-1,ix-1]
                y_in_up    = (y_in[iy,ix] - y_in[iy,ix-1])*(x_out - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + y_in[iy,ix-1]
                y_in_down  = (y_in[iy-1,ix] - y_in[iy-1,ix-1])*(x_out - x_in[iy-1,ix-1])/(x_in[iy-1,ix] - x_in[iy-1,ix-1]) + y_in[iy-1,ix-1]
                z_out      = (z_out_up - z_out_down)*(y_out - y_in_down)/(y_in_up - y_in_down) + z_out_down
            elif dy <= 0:
                # nearest input point is down
                z_out_up   = (z_in[:,iy+1,ix] - z_in[:,iy+1,ix-1])*(x_out - x_in[iy+1,ix-1])/(x_in[iy+1,ix] - x_in[iy+1,ix-1]) + z_in[:,iy+1,ix-1]
                z_out_down = (z_in[:,iy,ix] - z_in[:,iy,ix-1])*(x_out - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + z_in[:,iy,ix-1]
                y_in_up    = (y_in[iy+1,ix] - y_in[iy+1,ix-1])*(x_out - x_in[iy+1,ix-1])/(x_in[iy+1,ix] - x_in[iy+1,ix-1]) + y_in[iy+1,ix-1]
                y_in_down  = (y_in[iy,ix] - y_in[iy,ix-1])*(x_out - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + y_in[iy,ix-1]
                z_out      = (z_out_up - z_out_down)*(y_out - y_in_down)/(y_in_up - y_in_down) + z_out_down
        elif dx <= 0:
            # nearest input point is to the left
            if dy > 0:
                # nearest input point is up
                z_out_up   = (z_in[:,iy,ix+1] - z_in[:,iy,ix])*(x_out - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + z_in[:,iy,ix]
                z_out_down = (z_in[:,iy-1,ix+1] - z_in[:,iy-1,ix])*(x_out - x_in[iy-1,ix])/(x_in[iy-1,ix+1] - x_in[iy-1,ix]) + z_in[:,iy-1,ix]
                y_in_up    = (y_in[iy,ix+1] - y_in[iy,ix])*(x_out - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + y_in[iy,ix]
                y_in_down  = (y_in[iy-1,ix+1] - y_in[iy-1,ix])*(x_out - x_in[iy-1,ix])/(x_in[iy-1,ix+1] - x_in[iy-1,ix]) + y_in[iy-1,ix]
                z_out      = (z_out_up - z_out_down)*(y_out - y_in_down)/(y_in_up - y_in_down) + z_out_down
            elif dy <= 0:
                # nearest input point is down
                z_out_up   = (z_in[:,iy+1,ix+1] - z_in[:,iy+1,ix])*(x_out - x_in[iy+1,ix])/(x_in[iy+1,ix+1] - x_in[iy+1,ix]) + z_in[:,iy+1,ix]
                z_out_down = (z_in[:,iy,ix+1] - z_in[:,iy,ix])*(x_out - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + z_in[:,iy,ix]
                y_in_up    = (y_in[iy+1,ix+1] - y_in[iy+1,ix])*(x_out - x_in[iy+1,ix])/(x_in[iy+1,ix+1] - x_in[iy+1,ix]) + y_in[iy+1,ix]
                y_in_down  = (y_in[iy,ix+1] - y_in[iy,ix])*(x_out - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + y_in[iy,ix]
                z_out      = (z_out_up - z_out_down)*(y_out - y_in_down)/(y_in_up - y_in_down) + z_out_down
    elif kind == 3:
        iy, ix = np.where(r < d)
        z_out = operation(z_in[:,iy%z_in.shape[1],ix%z_in.shape[2]], axis = 1)
    
    return z_out

from datetime import datetime as dt
from netCDF4 import Dataset
from analysis_tools import bilinear_interpolation, get_cs_coords

theta_key = u'STASH_m01s00i004'
bouy_nc = Dataset('../bouy_09.nc', 'r')
z = bouy_nc.variables['thlev_zsea_theta'][:]
X = np.arange(0., 116000., 100.)
Y = np.arange(0., 31900., 100.)
X, Y = np.meshgrid(X, Y)

x_cs, y_cs = get_cs_coords(108000., 15950., 80., X, Y, max_r = 20000.)

output = np.zeros((bouy_nc.variables[theta_key].shape[1], len(x_cs)))
start = dt.now()
output = bilinear_interpolation(x_in = X, y_in = Y, z_in = bouy_nc.variables[theta_key][-1,:,:,:], x_out = x_cs, y_out = y_cs, kind = 2)

end = dt.now()
dt_old = end - start
print 'Bilinear Interpolation Original took = ' + str(dt_old)

def do_work(i):
    output = bilinear_interpolation_new(x_in = X, y_in = Y, z_in = bouy_nc.variables[theta_key][-1,:,:,:], x_out = x_cs[i], y_out = y_cs[i], kind = 2)

start = dt.now()
p = Pool()
output2 = np.array(p.map(do_work, xrange(len(x_cs))))
p.close()
end = dt.now()
dt_new = end - start
print 'Bilinear Interpolation New took = ' + str(dt_new)

R = - np.sign(x_cs)*np.sqrt((x_cs - 108000.)**2. + (y_cs - 15950.)**2.)

import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(1, 2, 1)
o1 = ax.contourf(z, R, output, cmap = 'viridis')
fig.colorbar(o1, ax = ax)
ax.set_title('Old = ' + str(dt_old))

ax = fig.add_subplot(1,2,2)
o2 = ax.contourf(z, R, output2, cmap = 'viridis')
fig.colorbar(o2, ax = ax)
ax.set_title('New = ' + str(dt_new))
plt.show()

