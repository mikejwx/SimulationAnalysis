"""
Code to regrid data from the original cartesian coordinates to coordinates
which follow the horizontal mean wind direction.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
<<<<<<< HEAD
from analysis_tools import fromComponents, get_cs_coords, bilinear_interpolation, in_plane_winds, getML_mean, send_email
=======
from analysis_tools import fromComponents, get_cs_coords, bilinear_interpolation, in_plane_winds, getML_mean
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
from netCDF4 import Dataset
from datetime import datetime as dt
from multiprocessing import Pool
import os

def get_swath_coord(iC):
    swath_x, swath_y = get_cs_coords(x_cs[iC], y_cs[iC], orientation_in + 90., X, Y, h = res, isPeriodic = True, max_r = swath_width/2.)
    return swath_x, swath_y

################################################################################
#                                                                              #
# Section One: Define the direction to cut a swath and compute swath coords    #
#                                                                              #
################################################################################
"""
Define the swath to be along the mean mixed layer wind direction, and we want 
our swath to go across the island
"""

# Define the keys for the wind components
print '[' + dt.now().strftime("%H:%M:%S") + '] Reading keys'
from STASH_keys import u_key, v_key, zi_new_key, mcl_key

# Get a list of all of the wind nc files
print '[' + dt.now().strftime("%H:%M:%S") + '] Determining the swath orientation'
<<<<<<< HEAD
my_path = '/nerc/n02/n02/xb899100/CloudTrail/Control_1600m/'
=======
my_path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
my_files = os.listdir(my_path)
wind_files = [my_file for my_file in my_files if 'wind' in my_file]
zi_files = [my_file for my_file in my_files if 'zi' in my_file]
wind_files.sort()
zi_files.sort()

if len(wind_files) == 8:
    print '[' + dt.now().strftime("%H:%M:%S") + '] This is a long experiment'
    l_short = False
else:
    print '[' + dt.now().strftime("%H:%M:%S") + '] This is a short experiment'
    l_short = True

# Get the horizontally averaged winds
for ncfile in wind_files:
    print '[' + dt.now().strftime("%H:%M:%S") + '] ' + ncfile
    wind_nc = Dataset(my_path + ncfile, 'r')
    if ncfile == wind_files[0]:
        z = wind_nc.variables['thlev_zsea_theta'][:]
        u = np.nanmean(wind_nc.variables[u_key][:], axis = (2, 3))
        v = np.nanmean(wind_nc.variables[v_key][:], axis = (2, 3))
    else:
        u = np.concatenate((u, np.nanmean(wind_nc.variables[u_key][:], axis = (2, 3))), axis = 0)
        v = np.concatenate((v, np.nanmean(wind_nc.variables[v_key][:], axis = (2, 3))), axis = 0)
    wind_nc.close()

print '[' + dt.now().strftime("%H:%M:%S") + '] Get the time averaged winds'
u = np.nanmean(u, axis = 0)
v = np.nanmean(v, axis = 0)

print '[' + dt.now().strftime("%H:%M:%S") + '] Get the horizontally averaged mixed layer depth'
for ncfile in zi_files:
    print '[' + dt.now().strftime("%H:%M:%S") + '] ' + ncfile
    zi_nc = Dataset(my_path + ncfile, 'r')
    if ncfile == zi_files[0]:
        zi = np.nanmean(zi_nc.variables[zi_new_key][:], axis = (1, 2))
    else:
        zi = np.concatenate((zi, np.nanmean(zi_nc.variables[zi_new_key][:], axis = (1, 2))), axis = 0)
    zi_nc.close()

# Get the time averaged mixed layer depth
zi = np.nanmean(zi, axis = 0)

# Get the mixed layer mean wind direction
u_ML = getML_mean(u, z, zi)
v_ML = getML_mean(v, z, zi)
orientation_in = fromComponents(u_ML, v_ML)[1]

print '[' + dt.now().strftime("%H:%M:%S") + '] ' + str(orientation_in)

print '[' + dt.now().strftime("%H:%M:%S") + '] Defining the swath coordinates to interpolate onto'
# We know from our domain definition where the centre of the island should be
R_i = 1000.0*(50.0/np.pi)**0.5 # Island radius m
x_c = 100000.0 + R_i           # x-coordinate for the island centre
y_c = 4*R_i                    # y-coordinate for the island centre

# Create the coordinate system for this data
<<<<<<< HEAD
dx = 1600.0
X, Y = np.meshgrid(np.arange(74)*dx, np.arange(20)*dx)

# Get coordinates of the cross section along the swath
h = dx*1. # Resolution in the along-swath direction
=======
X, Y = np.meshgrid(np.arange(0., 116000., 100.), np.arange(0., 31900., 100.))

# Get coordinates of the cross section along the swath
h = 100. # Resolution in the along-swath direction
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
x_cs, y_cs = get_cs_coords(x_c, y_c, orientation_in, X, Y, h = h)

# Truncate to one domain - i.e. remove points that exit a boundary
removal = []
for i in xrange(len(x_cs)):
    if (np.min(X) > x_cs[i]) or (np.max(X) < x_cs[i]) or (np.min(Y) > y_cs[i]) or (np.max(Y) < y_cs[i]):
        removal.append(i)

x_cs = np.delete(x_cs, removal)
y_cs = np.delete(y_cs, removal)

R_i /= 1000. # Convert the island radius to km

# Calculate the distance from the island centre in the along swath direction
R = -np.array([np.round(x, 0) for x in np.sign(x_cs - x_c)*np.sqrt((x_cs - x_c)**2 + (y_cs - y_c)**2)]) - R_i*1000.

# Define how wide and long we want the swath to be
swath_width  = 10000.
<<<<<<< HEAD
res = dx*1. # Resolution in the across-swath direction
=======
res = 100. # Resolution in the across-swath direction
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

# Use multiprocessing to quickly compute all the coordinates within the swath
p = Pool()
tempSwath = p.map(get_swath_coord, xrange(len(R)))
p.close()
p.join()

# Decompose the returned tempSwath into x and y 2D coordinate arrays
# Initialise the coordinate arrays
swath_x = np.zeros((int(swath_width/res + 1), len(R)))
swath_y = np.zeros((int(swath_width/res + 1), len(R)))
for iC in xrange(len(R)):
    swath_x[:,iC], swath_y[:,iC] = tempSwath[iC]

print '[' + dt.now().strftime("%H:%M:%S") + '] Creating the swath'
"""
Compute the swath along the cloud trail region. This is along the 
orientation_in direction.

----------------------------------------------------------------------------
INPUT:
"""
<<<<<<< HEAD
path = my_path
=======
path = '/nerc/n02/n02/xb899100/CloudTrail/Control/'
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
nc_in = 'mr'
var_in = mcl_key
nc_out = 'mcl_swath'

"""
OUTPUT:
Creates a netCDF in path that contains nc_in regridded along the swath, includes
a land mask layer in the new coordinate system.
"""

# Use the naming convention for the input files
if l_short:
    hours = ["{0:02d}".format(h) for h in xrange(0, 16, 4)]
else:
    hours = ["{0:02}".format(h) for h in xrange(0, 24, 3)]

def createSwathNC(hour):
    ############################################################################
    #                                                                          #
    # Section Two: Read in the netCDF, do any processing and concatenate       #
    #                                                                          #
    ############################################################################
    
    print '[' + dt.now().strftime("%H:%M:%S") + '] Reading input variable for ' + hour
    input_nc = Dataset(path + nc_in + '_' + hour + '.nc', 'r')
    
    time_key = [i for i in input_nc.variables.keys() if 'min' in i][0]
    # Initialise the any timeseries arrays
    input_times = input_nc.variables[time_key][:]
    output_var = np.zeros((input_nc.variables[var_in].shape[0], input_nc.variables[var_in].shape[1], swath_x.shape[0], swath_x.shape[1]))
    for it in xrange(input_nc.variables[var_in][:].shape[0]):
        print '[' + dt.now().strftime("%H:%M:%S") + '] Interpolating time index ' + "{0:02d}".format(it)
        output_var[it,:,:,:] = bilinear_interpolation(X, Y, input_nc.variables[var_in][it,:,:,:], swath_x.flatten(), swath_y.flatten(), kind = 2).reshape((len(z), int(swath_width/res + 1), len(R)))
    
    ############################################################################
    #                                                                          #
    # Section Three: Save new netCDF                                           #
    #                                                                          #
    ############################################################################
    
    print '[' + dt.now().strftime("%H:%M:%S") + '] Creating new netCDF'
    # Create a new netcdf for the regridded wind components
    output_nc = Dataset(path + nc_out + '_' + hour + '.nc', 'w')
    
    print '[' + dt.now().strftime("%H:%M:%S") + '] Creating dimensions for the netCDF'
    # Create dimensions for that netcdf
    time_dim = output_nc.createDimension(time_key, input_nc.variables[var_in][:].shape[0])
    z_dim    = output_nc.createDimension('thlev_zsea_theta', input_nc.variables[var_in][:].shape[1])
    y_dim    = output_nc.createDimension('y_prime', swath_y.shape[0])
    x_dim    = output_nc.createDimension('x_prime', swath_x.shape[1])
    
    print '[' + dt.now().strftime("%H:%M:%S") + '] Creating variables for the netCDF'
    # Create variables for each dimension
    times_var = output_nc.createVariable(time_key, np.float32, (time_key,))
    z_var     = output_nc.createVariable('thlev_zsea_theta', np.float32, ('thlev_zsea_theta',))
    ys_var    = output_nc.createVariable('y_prime', np.float32, ('y_prime','x_prime'))
    xs_var    = output_nc.createVariable('x_prime', np.float32, ('y_prime','x_prime'))
    
    # Create the variables to store the wind components
    out_var   = output_nc.createVariable(var_in, np.float64, (time_key, 'thlev_zsea_theta', 'y_prime', 'x_prime'))
    lsm_var   = output_nc.createVariable('lsm', np.float32, ('y_prime', 'x_prime'))
    
    print '[' + dt.now().strftime("%H:%M:%S") + '] Populating the dimension variables'
    # populate the dimension variables
    times_var[:] = input_times
    z_var[:]     = input_nc.variables['thlev_zsea_theta'][:]
    ys_var[:]    = swath_y[:]
    xs_var[:]    = swath_x[:]
    
    print '[' + dt.now().strftime("%H:%M:%S") + '] Populating the variable'
    out_var[:] = output_var[:]*1.
    lsm_var[:] = lsm_var_interp[:]*1.
    
    print '[' + dt.now().strftime("%H:%M:%S") + '] Closing the input ncs'
    input_nc.close()
    
    print '[' + dt.now().strftime("%H:%M:%S") + '] Closing the output_nc'
    output_nc.close()
    
    print '[' + dt.now().strftime("%H:%M:%S") + '] Complete.'

# Create the land-sea mask interpolated field
<<<<<<< HEAD
lsm_nc   = Dataset('/work/n02/n02/xb899100/island_masks/lsm50_1600m.nc', 'r')
=======
lsm_nc   = Dataset('/work/n02/n02/xb899100/island_masks/lsm50.nc', 'r')
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
lsm_var_interp = bilinear_interpolation(X, Y, lsm_nc.variables['lsm'][0,:,:,:], swath_x.flatten(), swath_y.flatten(), kind = 2).reshape((int(swath_width/res + 1), len(R)))
lsm_nc.close()

p = Pool(processes = len(hours))
p.map(createSwathNC, hours)
p.close()
p.join()

<<<<<<< HEAD
send_email(message = 'Completed mcl swath for ' + my_path.split('/')[-2] + ' experiment.', subject = 'regrid_mcl_along_flow.py', attachments = [''], isAttach = False)
=======
#in_plane_winds(u, v, orientation = 90.)



>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696

