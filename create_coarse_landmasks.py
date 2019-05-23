"""
Code to create new land masks for the coarser resolution experiments.
"""
import numpy as np
import matplotlib.pyplot as plt
from coarse_graining import coarse_grain as cg
from netCDF4 import Dataset
from math import ceil

l_testing = True

# Read the existing land mask
lsm_path = '/work/n02/n02/xb899100/island_masks/'
original_lsm = Dataset(lsm_path + 'lsm50.nc', 'r')
lsm_key = 'lsm'

# Coarse grain to the other resolutions
lsm_0100 = original_lsm.variables[lsm_key][:]
lsm_0800 = cg(lsm_0100, 8, 8, operation = np.nanmax)
lsm_1600 = cg(lsm_0100, 16, 16, operation = np.nanmax)

original_lsm.close()

def createlsm_nc(filename, new_lsm_array, d):
    # Changing the longitude and latitude dimensions, also changing the lsm variable
    dim_exclude = ['longitude', 'latitude']
    var_exclude = ['longitude', 'latitude', 'lsm']
    
    with Dataset(lsm_path + 'lsm50.nc', 'r') as lsm_nc, Dataset(lsm_path + filename, 'w') as new_lsm_nc:
        # Create some dimensions for the new netCDF
        D = lsm_nc.variables['longitude'][1] - lsm_nc.variables['longitude'][0]
        X = np.arange(0., lsm_nc.variables['longitude'][:].max() + D, d)
        Y = np.arange(0., lsm_nc.variables['latitude'][:].max() + D, d)
        
        # Copy global attributes all at once via dictionary
        new_lsm_nc.setncatts(lsm_nc.__dict__)
        
        # Copy dimensions
        for name, dimension in lsm_nc.dimensions.items():
            # Copy all of the dimensions except latitude and longitude
            if name not in dim_exclude:
                new_lsm_nc.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))
            else:
                # If latitude and longitude, copy the name but put the new dimension sizes in from the new array
                new_lsm_nc.createDimension(
                    name, (new_lsm_array.shape[2] if name == 'latitude' else new_lsm_array.shape[3]))
        
        # Copy all file data except for the excluded
        for name, variable in lsm_nc.variables.items():
            x = new_lsm_nc.createVariable(name, variable.datatype, variable.dimensions)
            # Copy variable attributes all at once via dictionary
            new_lsm_nc[name].setncatts(lsm_nc[name].__dict__)
            if name not in var_exclude:
                new_lsm_nc[name][:] = lsm_nc[name][:]
            else:
                if name == 'lsm':
                    new_lsm_nc[name][:] = new_lsm_array[:]
                elif name == 'longitude':
                    new_lsm_nc[name][:] = X
                elif name == 'latitude':
                    new_lsm_nc[name][:] = Y

# Create the 800 m lsm netCDF
createlsm_nc('lsm50_0800.nc', lsm_0800, 800.)

# Create the 1600 m lsm netCDF
createlsm_nc('lsm50_1600.nc', lsm_1600, 1600.)

# If testing
if l_testing:
    fig = plt.figure()
    ax0 = fig.add_subplot(3, 1, 1, adjustable = 'box', aspect = 1)
    ax0.contourf(lsm_0100[0,0,:,:])
    
    ax1 = fig.add_subplot(3, 1, 2, adjustable = 'box', aspect = 1)
    ax1.contourf(lsm_0800[0,0,:,:])
    
    ax2 = fig.add_subplot(3, 1, 3, adjustable = 'box', aspect = 1)
    ax2.contourf(lsm_1600[0,0,:,:])
    
    plt.show()


