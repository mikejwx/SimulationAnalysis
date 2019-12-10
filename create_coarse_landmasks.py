"""
Code to create new land masks for the coarser resolution experiments.
"""
import numpy as np
import matplotlib.pyplot as plt
from coarse_graining import coarse_grain as cg
from netCDF4 import Dataset
from math import ceil

def fix_x(array, axis = 3):
    """
    Expects a 4D array
    expects you to want to even the x-axis, i.e. axis = 3
    """
    if array.shape[axis]%2 == 1:
        # if there's an odd number of x points, increase that dimension by 1
        fixed_array = np.zeros((array.shape[axis-3], array.shape[axis-2], array.shape[axis-1], array.shape[axis]+1))
        fixed_array[:,:,:,:array.shape[axis]] += array
    else:
        fixed_array = array*1.
    return fixed_array

def make_veg(lsm, my_layer = 6):
    """
    Function to make a vegetation array
    """
    veg_out = np.zeros_like(lsm)
    for layer in range(9):
        if layer != 0:
            if layer == (my_layer - 1):
                veg_out = np.concatenate((veg_out, np.zeros_like(lsm) + 1), axis = 1)
            else:
                veg_out = np.concatenate((veg_out, np.zeros_like(lsm)), axis = 1)
    
    return veg_out

def create_nc(filename, new_array, d, in_path, original_fn, var):
    # Changing the longitude and latitude dimensions, also changing the lsm variable
    dim_exclude = ['longitude', 'latitude']
    var_exclude = ['longitude', 'latitude', var]
    
    with Dataset(in_path + original_fn, 'r') as old_nc, Dataset(in_path + filename, 'w') as new_nc:
        # Create some dimensions for the new netCDF
        D = old_nc.variables['longitude'][1] - old_nc.variables['longitude'][0]
        X = np.arange(new_array.shape[3])*d
        Y = np.arange(new_array.shape[2])*d
        
        # Copy global attributes all at once via dictionary
        for global_attr in old_nc.__dict__:
            print global_attr
         
        new_nc.setncatts(old_nc.__dict__)
        
        # Copy dimensions
        for name, dimension in old_nc.dimensions.items():
            # Copy all of the dimensions except latitude and longitude
            if name not in dim_exclude:
                new_nc.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))
            else:
                # If latitude and longitude, copy the name but put the new dimension sizes in from the new array
                new_nc.createDimension(
                    name, (new_array.shape[2] if name == 'latitude' else new_array.shape[3]))
        
        # Copy all file data except for the excluded
        for name, variable in old_nc.variables.items():
            x = new_nc.createVariable(name, variable.datatype, variable.dimensions)
            # Copy variable attributes all at once via dictionary
            new_nc[name].setncatts(old_nc[name].__dict__)
            if name not in var_exclude:
                new_nc[name][:] = old_nc[name][:]
            else:
                if name == var:
                    new_nc[name][:] = new_array[:]
                elif name == 'longitude':
                    new_nc[name][:] = X
                elif name == 'latitude':
                    new_nc[name][:] = Y
        new_nc.set_fill_off()

### User inputs ###
l_testing = True
new_resolution = 1600.0
original_resolution = 100.0
###################

# Read the existing land mask
lsm_path = '/work/n02/n02/xb899100/island_masks/'
original_lsm = Dataset(lsm_path + 'lsm50.nc', 'r')
original_veg = Dataset(lsm_path + 'island_050_vegfrac.nc', 'r')
lsm_key = 'lsm'
veg_key = 'field1391'

# Coarse grain to the other resolutions
lsm_0100 = original_lsm.variables[lsm_key][:]
new_step = int(new_resolution/original_resolution)
lsm_new0 = cg(lsm_0100, new_step, new_step, operation = np.nanmax)

# There is a quirk in the UM such that the x-coordinate cannot have an odd length.
# make it even here:
lsm_new = fix_x(lsm_new0)
veg_new = make_veg(lsm_new)

original_lsm.close()
original_veg.close()

# Create the new lsm netCDF
create_nc(filename = 'lsm50_' + "{0:04d}".format(int(new_resolution)) + '.nc', new_array = lsm_new, d = new_resolution, in_path = lsm_path, original_fn = 'lsm50.nc', var = lsm_key)
create_nc(filename = 'island_050_vegfrac_' + "{0:04d}".format(int(new_resolution)) + '.nc', new_array = veg_new, d = new_resolution, in_path = lsm_path, original_fn = 'island_050_vegfrac.nc', var = veg_key)

# If testing
if l_testing:
    fig = plt.figure()
    ax0 = fig.add_subplot(2, 1, 1, adjustable = 'box', aspect = 1)
    ax0.contourf(lsm_0100[0,0,:,:])
    
    ax1 = fig.add_subplot(2, 1, 2, adjustable = 'box', aspect = 1)
    ax1.contourf(lsm_new[0,0,:,:])
    
    plt.show()


