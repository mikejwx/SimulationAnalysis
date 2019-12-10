"""
Functions to coarse grain simulations in whole-number multiples of the input
grid-spacing. Only horizontal coarse-graining is supported.
"""
import numpy as np
def coarse_grain(input_data, reduce_y, reduce_x, operation = np.nanmean, l_periodic = True):
    """
    input_data = input data array, must have the following dimensions
                 f[(t), (z), y, x]
    reduce_x   = number of grid spaces to coarse grain to
                 i.e. if 2, then steps of 2 x-coordinates are averaged together
    reduce_y   = number of grid spaces to coarse grain to
    l_periodic = logical flag to make use of periodic boundary conditions, in 
                 case there is a remainder of point in a given direction
    ----------------------------------------------------------------------------
    output_data = coarse grained array to be output
    """
    # Determine the shape of the input data and behave differently for different
    # data shapes
    in_shape = input_data.shape
    
    if len(in_shape) == 2:
        # This is a 2D array, assume f[y,x]
        output_data = np.zeros_like(input_data[::reduce_y, ::reduce_x])
        for j in xrange(output_data.shape[0]):
            for i in xrange(output_data.shape[1]):
                j_in = [(j*reduce_y+J)%input_data.shape[0] for J in xrange(reduce_y)]
                i_in = [(i*reduce_x+I)%input_data.shape[1] for I in xrange(reduce_x)]
                j_in, i_in = np.meshgrid(j_in, i_in)
                output_data[j,i] = operation(input_data[j_in,i_in])
    elif len(in_shape) == 3:
        # This is a 3D array, assume f[(t/z), y, x]
        output_data = np.zeros_like(input_data[:,::reduce_y, ::reduce_x])
        for j in xrange(output_data.shape[1]):
            for i in xrange(output_data.shape[2]):
                j_in = [(j*reduce_y+J)%input_data.shape[1] for J in xrange(reduce_y)]
                i_in = [(i*reduce_x+I)%input_data.shape[2] for I in xrange(reduce_x)]
                j_in, i_in = np.meshgrid(j_in, i_in)
                output_data[:,j,i] = operation(input_data[:,j_in,i_in], axis = (1, 2))
    elif len(in_shape) == 4:
        # This is a 4d array, assume f[t, z, y, x]
        output_data = np.zeros_like(input_data[:,:,::reduce_y, ::reduce_x])
        for j in xrange(output_data.shape[2]):
            for i in xrange(output_data.shape[3]):
                j_in = [(j*reduce_y+J)%input_data.shape[2] for J in xrange(reduce_y)]
                i_in = [(i*reduce_x+I)%input_data.shape[3] for I in xrange(reduce_x)]
                j_in, i_in = np.meshgrid(j_in, i_in)
                output_data[:,:,j,i] = operation(input_data[:,:,j_in,i_in], axis = (2, 3))
    return output_data

def fix_x(array, axis = 1):
    """
    Expects a 2D array
    expects you to want to even the x-axis, i.e. axis = 3
    """
    if array.shape[axis]%2 == 1:
        # if there's an odd number of x points, increase that dimension by 1
        if axis == 0:
            fixed_array = np.concatenate((array, np.zeros((1,array.shape[1]))), axis = 0)
            fixed_array[-1,:] = 0.5*(fixed_array[0,:] + fixed_array[-2,:])
        elif axis == 1:
            fixed_array = np.concatenate((array, np.zeros((array.shape[0], 1))), axis = 1)
            fixed_array[:,-1] = 0.5*(fixed_array[:,0] + fixed_array[:,-2])
    else:
        fixed_array = array*1.
    return fixed_array

"""
# Do a test
lwp_path   = '/nerc/n02/n02/xb899100/BOMEX/DX050_128/water.nc'
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
lwp_key = u'STASH_m01s30i405'
lwp_nc = Dataset(lwp_path, 'r')

lwp_data = lwp_nc.variables[lwp_key][-1,:,:]*1000.
X, Y = np.meshgrid(np.arange(0, 6400., 50.), np.arange(0., 6400., 50.))
lwp_nc.close()

fig = plt.figure(tight_layout = True)
ax0 = fig.add_subplot(2, 3, 1)
ax0.contourf(lwp_data, cmap = 'Greys_r', levels = np.arange(0, 250.1, 25.))
ax0.set_title('50m')

ax1 = fig.add_subplot(2, 3, 2)
ax1.contourf(coarse_grain(lwp_data, reduce_y = 2, reduce_x = 2), cmap = 'Greys_r', levels = np.arange(0, 250.1, 25.))
ax1.set_title('100m')

ax2 = fig.add_subplot(2, 3, 3)
ax2.contourf(coarse_grain(lwp_data, reduce_y = 4, reduce_x = 4), cmap = 'Greys_r', levels = np.arange(0, 250.1, 25.))
ax2.set_title('200m')

ax3 = fig.add_subplot(2, 3, 4)
ax3.contourf(coarse_grain(lwp_data, reduce_y = 8, reduce_x = 8), cmap = 'Greys_r', levels = np.arange(0, 250.1, 25.))
ax3.set_title('400m')

ax4 = fig.add_subplot(2, 3, 5)
ax4.contourf(coarse_grain(lwp_data, reduce_y = 16, reduce_x = 16), cmap = 'Greys_r', levels = np.arange(0, 250.1, 25.))
ax4.set_title('800m')

plt.show()
"""
