### Find the cloud weighted wind direction ###
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate, integrate

def toComponents(speed, direction):
    "takes the wind speed and wind direction and breaks the wind speed into"
    "u- and v- components to be returned."
    "speed = wind speed"
    "direction = wind direction"
    u = - speed*np.sin(direction*np.pi/180)
    v = - speed*np.cos(direction*np.pi/180)
    
    return u, v

def fromComponents(u, v, isList = False):
    "Takes the two horizontal wind components and combines them into a wind"
    "speed and direction to be returned."
    "u = the zonal wind component"
    "v = the meridional wind component"
    speed = np.sqrt(np.array(u)**2 + np.array(v)**2)
    if isList:
        direction = np.arcsin(v/speed)*180/np.pi
        for iV in range(len(u)):
            if speed[iV] > 0:
                if (u[iV] >= 0) and (v[iV] > 0):
                    #if u > 0 and v > 0 SW quadrant
                    direction[iV] = 270.0 - direction[iV]
                elif (u[iV] > 0) and (v[iV] <= 0):
                    #if u > 0 and v < 0 NW quadrant
                    direction[iV] = 270.0 - direction[iV]
                elif (u[iV] <= 0) and (v[iV] > 0):
                    #if u < 0 and v > 0 SE quadrant
                    direction[iV] = 90.0 + direction[iV]
                elif (u[iV] < 0) and (v[iV] <= 0):
                    #if u < 0 and v < 0 NE quadrant
                    direction[iV] = 90.0 + direction[iV]
            else:
                direction[iV] = 0.
    else:
        direction = np.arcsin(v/speed)*180/np.pi
        if speed > 0:
            if (u >= 0) and (v > 0):
                #if u > 0 and v > 0 SW quadrant
                direction = 270.0 - direction
            elif (u > 0) and (v <= 0):
                #if u > 0 and v < 0 NW quadrant
                direction = 270.0 - direction
            elif (u <= 0) and (v > 0):
                #if u < 0 and v > 0 SE quadrant
                direction = 90.0 + direction
            elif (u < 0) and (v <= 0):
                #if u < 0 and v < 0 NE quadrant
                direction = 90.0 + direction
        else:
            direction = 0.
    return speed, direction

def regrid(target_grid, current_grid, current_data_key):
    """
    Function to attempt to regrid an already read-in netCDF B (current_grid) to 
    the same grid as already read-in netCDF A (target_grid).
    Should only work if there is some overlap between the range of the functions
    otherwise the interpolation won't make sense.
    
    Regrids in all four dimensions, although going from e.g. 1 minute data to
    10 or 15 minute data shouldn't really require regridding. But going from
    10 minute to 1 minute data definitely requires regridding and the linear 
    interpolation approach used here is probably not appropriate.
    """
    # Read the data that contains the target grid
    target_x_key     = [key for key in target_grid.variables.keys() if ('longitude' in key) and ('_l' not in key)][0]
    target_x = target_grid.variables[target_x_key][:]
    target_y_key     = [key for key in target_grid.variables.keys() if ('latitude' in key) and ('_l' not in key)][0]
    target_y = target_grid.variables[target_y_key][:]
    target_z_key     = [key for key in target_grid.variables.keys() if ('zsea' in key) and ('bounds' not in key)][0]
    target_z = target_grid.variables[target_z_key][:]
    target_t_key     = [key for key in target_grid.variables.keys() if 'min' in key][0]
    target_t = target_grid.variables[target_t_key][:]
    
    # Read the netCDF that needs to be transformed
    current_x_key     = [key for key in current_grid.variables.keys() if ('longitude' in key) and ('_l' not in key)][0]
    current_x = current_grid.variables[current_x_key][:]
    current_y_key     = [key for key in current_grid.variables.keys() if ('latitude' in key) and ('_l' not in key)][0]
    current_y = current_grid.variables[current_y_key][:]
    current_z_key     = [key for key in current_grid.variables.keys() if ('zsea' in key) and ('bounds' not in key)][0]
    current_z = current_grid.variables[current_z_key][:]
    current_t_key     = [key for key in current_grid.variables.keys() if 'min' in key][0]
    current_t = current_grid.variables[current_t_key][:]
    current_data = current_grid.variables[current_data_key][:]
     
    ### Procedure ###
    # Will usually only need to interpolate in z, sometimes in x or y, probably
    # never required to interpolate in more than two dimensions
    # 1. interpolate in time
    if target_t_key != current_t_key:
        print 'Regridding in time...'
        # Interpolate in time
        current_data = interpolate.interp1d(current_t, current_data, axis = 0, fill_value = 'extrapolate')(target_t)
        print 'Complete.'
    # 2. interpolate in z
    if target_z_key != current_z_key:
        print 'Regridding in height...'
        # Interpolate in z
        current_data = interpolate.interp1d(current_z, current_data, axis = 1, fill_value = 'extrapolate')(target_z)
        print 'Complete.'
    # 3. interpolate in y
    if target_y_key != current_y_key:
        print 'Regridding in latitude...'
        # Interpolate in y
        for i in xrange(current_data.shape[3]):
            current_data[:,:,:,i] = interpolate.interp1d(current_y[:,i], current_data[:,:,:,i], axis = 2, fill_value = 'extrapolate')(target_y[:,i])
        print 'Complete.'
    # 4. interpolate in x
    if target_x_key != current_x_key:
        print 'Regridding in longitude...'
        # Interpolate in x
        for j in xrange(current_data.shape[2]):
            current_data[:,:,j,:] = interpolate.interp1d(current_x[j,:], current_data[:,:,j,:], axis = 2, fill_value = 'extrapolate')(target_x[j,:])
        print 'Complete.'
    
    return current_data

days = ["{0:02d}".format(x) for x in xrange(1, 11)]
ndays = len(days)
u_key = u'STASH_m01s00i002'
v_key = u'STASH_m01s00i003'
q_key = u'STASH_m01s00i392'
for day in days:
    print 'Starting day ' + day
    u = Dataset('u_'+day+'.nc', 'r')
    v = Dataset('v_'+day+'.nc', 'r')
    mr = Dataset('mr_'+day+'.nc', 'r')
    
    if day == '01':
        dt_i = u.variables[u_key].shape[0]
        u_cld = np.zeros(dt_i*ndays)
        v_cld = np.zeros_like(u_cld)
        z = mr.variables['thlev_zsea_theta'][:]*1.
    # get the cloud liquid water mixing ratio, u- and v- wind components
    cld = mr.variables[q_key][:]*1.
    u_rg = regrid(mr, u, u_key)
    v_rg = regrid(mr, v, v_key)
    ## weight wind components by the cloud liquid water mixing ratio
    # 1. multiply by the cld liquid water mixing ratio
    u_weight = u_rg*cld
    v_weight = v_rg*cld
    # 2. vertically integrate the weighted profiles
    u_weight = integrate.trapz(u_weight, z, axis = 1)
    v_weight = integrate.trapz(v_weight, z, axis = 1)
    # 3. divide by the vertically integrated cloud liquid water mixing ratio
    total_cld = integrate.trapz(cld, z, axis = 1)
    u_weight /= total_cld
    v_weight /= total_cld
    # 4. take the horizontal mean
    u_weight = np.nanmean(np.where((total_cld > 0), u_weight, np.nan), axis = (1, 2))
    v_weight = np.nanmean(np.where((total_cld > 0), v_weight, np.nan), axis = (1, 2))
    # 5. store to u_cld and v_cld
    u_cld[dt_i*days.index(day):dt_i*(days.index(day)+1)] = u_weight
    v_cld[dt_i*days.index(day):dt_i*(days.index(day)+1)] = v_weight
    u.close()
    v.close()
    mr.close()
# convert from components to wind direction in degrees
wind_direction = fromComponents(u_cld, v_cld, 1)[1]

plt.plot(wind_direction)
plt.show()
print 'The mean wind direction over the last 4 days is: ' + str(round(np.mean(wind_direction[-4*dt_i:]), 1))

