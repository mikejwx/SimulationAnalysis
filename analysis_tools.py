import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
from netCDF4 import Dataset

################################################################################
#                                                                              #
# Boundary Layer diagnostics                                                   #
#                                                                              #
################################################################################

def zi(theta_v, z):
    """
    Function to calculate the mixed layer depth as a function of the virtual
    potential temeprature.
    
    Input:
    theta_v = virtual potential temperature (K), theta_v = theta_v(z, y, x)
    z       = the heights (m) at which theta_v are computed
    
    Output:
    zi      = the mixed layer depth
    
    Assumes that the mixed layer depth is defined as the level to which a
    surface based parcel of air (i.e. theta_v_parcel = theta_v(z = 0)) must be 
    lifted in order for the parcel to become neutrally buoyant (i.e. 
    theta_v_parcel = theta_v_env(z = zi).
    """
    
    # Use a surface based parcel
    theta_v_parcel = theta_v[0,:,:]
    # Find the height of the level of neutral buoyancy
    z_i = np.zeros_like(theta_v_parcel)
    for j in xrange(theta_v.shape[1]):
        for i in xrange(theta_v.shape[2]):
            ik = np.max(np.where(np.abs(theta_v_parcel[j,i] - theta_v[2:,j,i]) == np.min(np.abs(theta_v_parcel[j,i] - theta_v[2:,j,i])))[0])
            if theta_v[ik,j,i] > theta_v_parcel[j,i]:
                # Means the nearest point is higher height than the boundary layer depth
                # Linearly interpolate to get z_i
                m = 1
                while theta_v[ik-m,j,i] == theta_v[ik,j,i]:
                    m += 1
                n = m - 1
                z_i[j,i] = (z[ik-n] - z[ik-m])*(theta_v_parcel[j,i] - theta_v[ik-m,j,i])/(theta_v[ik-n,j,i] - theta_v[ik-m,j,i]) + z[ik-m]
            elif theta_v[ik,j,i] < theta_v_parcel[j,i]:
                # Means the nearest point is at a lower height than the boundary layer depth
                m = 1
                while theta_v[ik+m,j,i] == theta_v[ik,j,i]:
                    m += 1
                n = m - 1
                z_i[j,i] = (z[ik+m] - z[ik+n])*(theta_v_parcel[j,i] - theta_v[ik+n,j,i])/(theta_v[ik+m,j,i] - theta_v[ik+n,j,i]) + z[ik+n]
            else:
                z_i[j,i] = z[ik]
            #z_i[j, i] = interpolate.interp1d(theta_v[t2:, j, i], z[2:], fill_value = "extrapolate")(theta_v_parcel[j,i])
            #if z_i[j, i] < 0.:
            #    z_i[j, i] = 0.
            
            #k = 2
            #while theta_v_parcel[j,i] > theta_v[k,j,i]:
            #    k += 1
            #z_i[j,i] = z[k]
    return z_i

def lcl(temp, q, pres, z):
    """
    lcl calculates the surface-based lifting condensation level as a function of
    the temperature (temp), specific humidity (q), and pressure (pres).
    These input are used to find the temperature and dew point which are used in
    Bolton (1980)'s equation for the LCL.
    
    Input:
    temp = air temperature (K), temp = temp[z=0, y, x]
    q    = specific humidity (kg/kg), q = q[z=0, y, x]
    pres = pressure (Pa), pres = pres[z, y, x]
    z    = heights (m)
    
    Output:
    lcl  = the lifting condensation level (m)
    """
    # Define some constants
    T0 = 273.15 # freezing point of water
    Lv = 2.501e6 # Latent heat of vapourisation
    Rv = 461. # Gas constant for water vapour
    e0 = 611.2 # vapour pressure of water vapour at T0
    Rd = 287.05
    E = Rd/Rv
    cp = 1005.
    
    # calculate the dew point temperature
    e = q*pres[0,:,:]/(E - q*E + q)
    dewpoint = (1./T0 - (Rv/Lv)*np.log(e/e0))**(-1.)
    
    term_A = 1./(dewpoint - 56)
    term_B = np.log(temp/dewpoint)/800.
    T_LCL = 1/(term_A + term_B) + 56
    
    p_LCL = pres[0,:,:]*(T_LCL/temp)**(cp/Rd)
    z_LCL = np.zeros_like(p_LCL)
    for j in xrange(pres.shape[1]):
        for i in xrange(pres.shape[2]):
            z_LCL[j,i] = interpolate.interp1d(pres[:,j,i], z)(p_LCL[j,i])
    
    return z_LCL

def find_h(theta_v, u, v, z):
    """
    Iteratively solves equation 2 from Vogelezang and Holtslag (1996) to
    estimate the boundary layer depth.
    
    This method uses a modified Richardson number approach and vertical profiles
    of virtual potential temperature (theta_v), u- and v- wind components to
    estimate the boundary layer height.
    
    theta_v = f(z)
    u       = f(z)
    v       = f(z)
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

def get_TKE(U, V, W, rho = 1., start = 0, end = None):
    """
    Function to compute the vertical profile of turbulent kinetic energy from 
    given 4D U, V, and W fields.
    
    First, calculates the perturbations from horizontal mean profiles of U, V, 
    and W. Then, calculates the turbulent kinetic energy, and then the 
    horizontal mean so that there is a time series of profiles.
    
    U, V, and W must all be on the same grid.
    """
    U_mean = horizontalMean(U)
    V_mean = horizontalMean(V)
    W_mean = horizontalMean(W)
    U_p = np.zeros_like(U)
    V_p = np.zeros_like(V)
    W_p = np.zeros_like(W)
    for t in xrange(U.shape[0]):
        for k in xrange(U.shape[1]):
            U_p[t,k,:,:] = U[t,k,:,:] - U_mean[t,k]
            V_p[t,k,:,:] = V[t,k,:,:] - V_mean[t,k]
            W_p[t,k,:,:] = W[t,k,:,:] - W_mean[t,k]
    if type(rho) != float:
        TKE_profile = 0.5*horizontalMean(rho)*(horizontalMean(U_p**2)+ horizontalMean(V_p**2) + horizontalMean(W_p**2))
    else:
        TKE_profile = 0.5*rho*(horizontalMean(U_p**2)+ horizontalMean(V_p**2) + horizontalMean(W_p**2))
    
    return TKE_profile

################################################################################
#                                                                              #
# Cloud diagnostics                                                            #
#                                                                              #
################################################################################

def get_CC(m_cl, start = 0, end = None, threshold = 0.):
    """
    Function to compute the vertical profile of cloud cover from given 4D m_cl 
    field.
    
    First, finds the number of cloudy pixels, i.e. those that contain non-zero 
    m_cl. This is a binary, 0 = no cloud, 1 = cloud. Then calculates the cloud 
    cover on each vertical level as the mean of the binary field at each level. 
    Finally, does a time mean over a given period. Returns a vertical profile of
    this cloud cover estimate.
    """
    cloud_mask = m_cl[:]
    cloud_mask[m_cl > threshold] = 1.
    cloud_mask[m_cl <= threshold] = 0.
    cc_profile = np.mean(horizontalMean(cloud_mask)[start:end,:], axis = 0)
    
    return cc_profile

################################################################################
#                                                                              #
# Regridding, interpolation and spatial smoothing routines                     #
#                                                                              #
################################################################################

def cartesian2polar(x, y, x_o = 0., y_o = 0.):
    """
    Convert from cartesian to polar coordinates
    x, y = the cartesian coordinates
    x_o, y_o = the coordinates of the origin
    """
    dx = x - x_o
    dy = y - y_o
    r = np.sqrt(dx**2 + dy**2)
    
    """
    Four cases:
    1. if dx and dy are positive
    2. if dx is positive and dy is negative
    3. if dx and dy are negative
    4. if dx is negative and dy is positive
    """
    theta1 = np.arcsin(dx/r)
    theta2 = np.arcsin(-dy/r) + np.pi/2.
    theta3 = np.arcsin(-dx/r) + np.pi
    theta4 = np.arcsin(dy/r)  + 3.*np.pi/2.
    
    theta1[dx < 0.] = 0.
    theta1[dy < 0.] = 0.
    
    theta2[dx < 0.] = 0.
    theta2[dy > 0.] = 0.

    theta3[dx > 0.] = 0.
    theta3[dy > 0.] = 0.
    
    theta4[dx > 0.] = 0.
    theta4[dy < 0.] = 0.
    
    # Add together and convert from radians to degrees
    theta = (theta1 + theta2 + theta3 + theta4)*180./np.pi
    
    return r, theta

def polar2cartesian(r, theta, x_o, y_o):
    """
    Converts from polar coordinates to cartesian coordinates
    r = distance from the origin
    theta = angle from x = 0 in degrees
    x_o, y_o = the origin in cartesian space
    """
    # Convert from degrees to radians
    theta = theta*np.pi/180.
    x = r*np.sin(theta)
    y = r*np.cos(theta)
    
    # return the origin to 0,0
    x += x_o
    y += y_o
    
    return x, y

def bilinear_interpolation(x_in, y_in, z_in, x_out, y_out, kind = 0, d = 2000.0, p = 0.5, fv = np.nan):
    """
    Performs bilinear interpolation on a given 2D array, z_in.
    x_in = 2D array of x-coordinates, input
    y_in = 2D array of y-coordinates, input
    z_in = 3D array of data points, input z_in[nz, ny, nx]

    x_out = 1D array of x-coordinates to be interpolated onto
    y_out = 1D array of y-coordinates to be interpolated onto

    interpolation scheme finds the nearest points and uses scipy's bisplrep to 
    do the interpolation.
    """
    # Initialise our output array
    z_out = np.zeros((z_in.shape[0], len(x_out)))
    
    # For each point to be interpolated onto
    for i in xrange(len(x_out)):
        # Find the distance to the input data coordinates
        r = np.sqrt((x_in - x_out[i])**2 + (y_in - y_out[i])**2)
        
        if kind == 0:
            # Nearest neighbor approach
            iy, ix = np.where(r == np.min(r))
            z_out[:,i] = z_in[:,iy[0],ix[0]]
        elif kind == 1:
            # idw approach
            # Use all points within 'd' m of x_out, y_out for weighted mean
            iy, ix = np.where(r < d)
            w = 0
            
            for j in xrange(len(iy)):
                w += (1./r[iy[j], ix[j]]**p)
                z_out[:,i] += (1./r[iy[j], ix[j]]**p)*z_in[:,iy[j], ix[j]]
            
            z_out[:,i] /= w
        elif kind == 2:
            # Find the nearest point
            iy, ix = np.where(r == np.min(r))
            iy = iy[0]
            ix = ix[0]
            # Determine which quadrant (i.e. up and left, up and right, down and
            # right, or down and left) this point is with respect to the output 
            # point
            dx = x_in[iy, ix] - x_out[i] # if +ve, input point is to the right
            dy = y_in[iy, ix] - y_out[i] # if +ve, input point is up
            if np.min([j in xrange(z_in.shape[1]) for j in [iy-1, iy, iy+1]]) and np.min([j in xrange(z_in.shape[2]) for j in [ix-1, ix, ix+1]]):
                if dx > 0:
                    # nearest input point is to the right
                    if dy > 0:
                        # nearest input point is up
                        z_out_up = (z_in[:,iy,ix] - z_in[:,iy,ix-1])*(x_out[i] - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + z_in[:,iy,ix-1]
                        z_out_down = (z_in[:,iy-1,ix] - z_in[:,iy-1,ix-1])*(x_out[i] - x_in[iy-1,ix-1])/(x_in[iy-1,ix] - x_in[iy-1,ix-1]) + z_in[:,iy-1,ix-1]
                        y_in_up = (y_in[iy,ix] - y_in[iy,ix-1])*(x_out[i] - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + y_in[iy,ix-1]
                        y_in_down = (y_in[iy-1,ix] - y_in[iy-1,ix-1])*(x_out[i] - x_in[iy-1,ix-1])/(x_in[iy-1,ix] - x_in[iy-1,ix-1]) + y_in[iy-1,ix-1]
                        z_out[:,i] = (z_out_up - z_out_down)*(y_out[i] - y_in_down)/(y_in_up - y_in_down) + z_out_down
                    elif dy <= 0:
                        # nearest input point is down
                        z_out_up = (z_in[:,iy+1,ix] - z_in[:,iy+1,ix-1])*(x_out[i] - x_in[iy+1,ix-1])/(x_in[iy+1,ix] - x_in[iy+1,ix-1]) + z_in[:,iy+1,ix-1]
                        z_out_down = (z_in[:,iy,ix] - z_in[:,iy,ix-1])*(x_out[i] - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + z_in[:,iy,ix-1]
                        y_in_up = (y_in[iy+1,ix] - y_in[iy+1,ix-1])*(x_out[i] - x_in[iy+1,ix-1])/(x_in[iy+1,ix] - x_in[iy+1,ix-1]) + y_in[iy+1,ix-1]
                        y_in_down = (y_in[iy,ix] - y_in[iy,ix-1])*(x_out[i] - x_in[iy,ix-1])/(x_in[iy,ix] - x_in[iy,ix-1]) + y_in[iy,ix-1]
                        z_out[:,i] = (z_out_up - z_out_down)*(y_out[i] - y_in_down)/(y_in_up - y_in_down) + z_out_down
                elif dx <= 0:
                    # nearest input point is to the left
                    if dy > 0:
                        # nearest input point is up
                        z_out_up = (z_in[:,iy,ix+1] - z_in[:,iy,ix])*(x_out[i] - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + z_in[:,iy,ix]
                        z_out_down = (z_in[:,iy-1,ix+1] - z_in[:,iy-1,ix])*(x_out[i] - x_in[iy-1,ix])/(x_in[iy-1,ix+1] - x_in[iy-1,ix]) + z_in[:,iy-1,ix]
                        y_in_up = (y_in[iy,ix+1] - y_in[iy,ix])*(x_out[i] - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + y_in[iy,ix]
                        y_in_down = (y_in[iy-1,ix+1] - y_in[iy-1,ix])*(x_out[i] - x_in[iy-1,ix])/(x_in[iy-1,ix+1] - x_in[iy-1,ix]) + y_in[iy-1,ix]
                        z_out[:,i] = (z_out_up - z_out_down)*(y_out[i] - y_in_down)/(y_in_up - y_in_down) + z_out_down
                    elif dy <= 0:
                        # nearest input point is down
                        z_out_up = (z_in[:,iy+1,ix+1] - z_in[:,iy+1,ix])*(x_out[i] - x_in[iy+1,ix])/(x_in[iy+1,ix+1] - x_in[iy+1,ix]) + z_in[:,iy+1,ix]
                        z_out_down = (z_in[:,iy,ix+1] - z_in[:,iy,ix])*(x_out[i] - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + z_in[:,iy,ix]
                        y_in_up = (y_in[iy+1,ix+1] - y_in[iy+1,ix])*(x_out[i] - x_in[iy+1,ix])/(x_in[iy+1,ix+1] - x_in[iy+1,ix]) + y_in[iy+1,ix]
                        y_in_down = (y_in[iy,ix+1] - y_in[iy,ix])*(x_out[i] - x_in[iy,ix])/(x_in[iy,ix+1] - x_in[iy,ix]) + y_in[iy,ix]
                        z_out[:,i] = (z_out_up - z_out_down)*(y_out[i] - y_in_down)/(y_in_up - y_in_down) + z_out_down
            else:
                z_out[:,i] = fv
        elif kind == 3:
            iy, ix = np.where(r < d)
            z_out[:,i] = np.nanmean(z_in[:,iy,ix], axis = 1)
        
    return z_out

def get_cs_coords(x_c, y_c, direction, h = 100.):
    """
    Uses an input coordinate (x_c, y_c) and a wind direction (direction) to 
    calculate the coordinates of a line that goes through the input coordinate
    and would be parallel to the wind direction. These coordinates are given
    at a horizontal resolution (h) in metres.
    ---------------------------------------------------------------------------
    x_c: the x-coordinate of the input point (e.g. island centre)
    y_c: the y-coordinate of the input point
    direction: the direction of the cross section
    h: the resolution at which to take points along the cross section
    """
    
    dx = (h*np.sin(np.pi*direction/180.0))
    dy = (h*np.cos(np.pi*direction/180.0))

    x_cs = [x_c]
    y_cs = [y_c]

    # Populate my list of x and y coordinates along the cross section
    while (x_cs[-1] < np.max(x))*(y_cs[-1] < np.max(y)):
        x_cs.append(x_cs[-1] + dx)
        y_cs.append(y_cs[-1] + dy)

    x_cs = x_cs[::-1]
    y_cs = y_cs[::-1]

    while (x_cs[-1] > np.min(x))*(y_cs[-1] > np.min(y)):
        x_cs.append(x_cs[-1] - dx)
        y_cs.append(y_cs[-1] - dy)

    # check that all the coordinates along the cross section are within the domain
    i = 0
    while i < len(x_cs):
        if (x_cs[i] > np.max(x)) or (x_cs[i] < np.min(x)) or (y_cs[i] > np.max(y)) or (y_cs[i] < np.min(y)):
            del x_cs[i]
            del y_cs[i]
        else:
            i += 1

    # store to arrays
    x_cs = np.array(x_cs)
    y_cs = np.array(y_cs)
    
    return x_cs, y_cs

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

def cross_section_x(data, x_rot, x_pos, t_pos = -1):
    """
    data = data[time, z, y, x]
    x_rot = rotated x-coordinate
    x_pos = x_coordinate to do cross section at
    """
    data_xs = np.zeros_like(data[t_pos,:,:,x_pos/100])
    x_target = x_rot - x_pos
    
    my_ixs = []
    for iy in range(w.shape[2]):
        ix = np.where(abs(x_target[iy,:]) == np.min(abs(x_target[iy,:])))[0][0]
        my_ixs.append(ix)
        if x_target[iy, ix] > 0.:
            # interpolate ix and ix-1
            data_xs[:,iy] = (data[t_pos,:,iy,ix] - data[t_pos,:,iy,ix-1])*(x_pos - x_rot[iy,ix-1])/(x_rot[iy,ix] - x_rot[iy,ix-1]) + data[t_pos,:,iy,ix-1]
        elif x_target[iy, ix] > 0:
            # interpolate ix+1 and ix
            data_xs[:,iy] = (data[t_pos,:,iy,ix+1] - data[t_pos,:,iy,ix])*(x_pos - x_rot[iy,ix])/(x_rot[iy,ix+1] - x_rot[iy,ix]) + data[t_pos,:,iy,ix]
        else:
            # no need to interpolate as x_target == 0
            data_xs[:,iy] = data[t_pos,:,iy,ix]
    return data_xs, my_ixs

def ddx(data, X, periodic = True):
    """
    Calculate the horizontal gradient in the x-direction.

    data     = 2, 3, or 4 dimensional array
    data     = f([t],[z],y,x)
    X        = 2D array of x-coordinates
    X        = f(y,x)
    periodic = Boolean, 
               .True. means periodic boundary conditions
               .False. means zero gradient boundary conditions.
    """
    print 'WARNING: gradients are returned at staggered coordinates'
    d_dx = np.zeros_like(data)
    if len(data.shape) == 2:
        d_dx[:,:-1] = (data[:,1:] - data[:,:-1])/(X[:,1:] - X[:,:-1])
        d_dx[:,-1] = (data[:,0] - data[:,-1])/(X[:,0] - X[:,-1])
    elif len(data.shape) == 3:
        d_dx[:,:,:-1] = (data[:,:,1:] - data[:,:,:-1])/(X[:,1:] - X[:,:-1])
        d_dx[:,:,-1] = (data[:,:,0] - data[:,:,-1])/(X[:,0] - X[:,-1])
    elif len(data.shape) == 4:
        d_dx[:,:,:,:-1] = (data[:,:,:,1:] - data[:,:,:,:-1])/(X[:,1:] - X[:,:-1])
        d_dx[:,:,:,-1] = (data[:,:,:,0] - data[:,:,:,-1])/(X[:,0] - X[:,-1])
    
    return d_dx

def ddy(data, Y, periodic = True):
    """
    Calculate the horizontal gradient in the x-direction.

    data     = 2, 3, or 4 dimensional array
    data     = f([t],[z],y,x)
    Y        = 2D array of y-coordinates
    Y        = f(y,x)
    periodic = Boolean, 
               .True. means periodic boundary conditions
               .False. means zero gradient boundary conditions.
    """
    print 'WARNING: gradients are returned at staggered coordinates'
    d_dy = np.zeros_like(data)
    if len(data.shape) == 2:
        d_dy[:-1,:] = (data[1:,:] - data[:-1,:])/(Y[1:,:] - Y[:-1,:])
        d_dy[-1,:] = (data[0,:] - data[-1,:])/(Y[0,:] - Y[-1,:])
    elif len(data.shape) == 3:
        d_dy[:,:-1,:] = (data[:,1:,:] - data[:,:-1,:])/(Y[1:,:] - Y[:-1,:])
        d_dy[:,-1,:] = (data[:,0,:] - data[:,-1,:])/(Y[0,:] - Y[-1,:])
    elif len(data.shape) == 4:
        d_dy[:,:,:-1,:] = (data[:,:,1:,:] - data[:,:,:-1,:])/(Y[1:,:] - Y[:-1,:])
        d_dy[:,:,-1,:] = (data[:,:,0,:] - data[:,:,-1,:])/(Y[0,:] - Y[-1,:])
    
    return d_dy

def circular_smoothing(X, Y, data, r):
    """
    Smooth the data by taking the mean of all points within radius r.
    Requires the following input:
    X, a 2D array f(y,x) of x-coordinates
    Y, a 2D array f(y,x) of y-coordinates
    data a 2, 3, or 4D array f([t],[z],y,x) of data to be smoothed
    r, a single float in the same units as the X and Y coordinates to define a 
    radius with which to smooth the data.
    
    N.B. a future version of this function may attempt to do some form of
    weighting to make the points closest to the actual point weigh more than
    the points further away.
    """
    x_periodic = np.concatenate((X-np.max(X), X, X+np.max(X)), axis = 1)
    x_periodic = np.concatenate((x_periodic, x_periodic, x_periodic), axis = 0)
    y_periodic = np.concatenate((Y-np.max(Y), Y, Y+np.max(Y)), axis = 0)
    y_periodic = np.concatenate((y_periodic, y_periodic, y_periodic), axis = 1)
    smoothed_data = np.zeros_like(data)
    if len(data.shape) == 2:
        for i in xrange(data.shape[0]):
            for j in xrange(data.shape[1]):
                iy, ix = np.where(np.sqrt((X[i,j] - x_periodic)**2 + (Y[i,j] - y_periodic)**2) <= r)
                smoothed_data[i,j] = np.nanmean(data[iy%X.shape[0],ix%X.shape[1]])
            plt.plot(X[i,:], data[i,:])
            plt.plot(X[i,:], smoothed_data[i,:])
            plt.show()
    elif len(data.shape) == 3:
        for tz0 in xrange(data.shape[0]):
            for i in xrange(data.shape[1]):
                for j in xrange(data.shape[2]):
                    iy, ix = np.where(np.sqrt((X[i,j] - x_periodic)**2 + (Y[i,j] - y_periodic)**2) <= r)
                    smoothed_data[tz0,i,j] = np.nanmean(data[tz0,iy%X.shape[0],ix%X.shape[1]])
    elif len(data.shape) == 4:
        for tz0 in xrange(data.shape[0]):
            for tz1 in xrange(data.shape[1]):
                for i in xrange(data.shape[2]):
                    for j in xrange(data.shape[3]):
                        iy, ix = np.where(np.sqrt((X[i,j] - x_periodic)**2 + (Y[i,j] - y_periodic)**2) <= r)
                        smoothed_data[tz0,tz1,i,j] = np.nanmean(data[tz0,tz1,iy%X.shape[0],ix%X.shape[1]])
    
    return smoothed_data

################################################################################
#                                                                              #
# Minimise number of points in a profile                                       #
#                                                                              #
################################################################################

def RDP(x, y, e):
    """
    Algorithm to choose the minimum number of points that retains the shape of 
    y with maximum error, "error". Returns the values of y, and the coordinates 
    x for those values.
    """
    # find the end of the data
    end = len(x)
    index0 = 0
    index1 = end -1
    # The first point of our data is the same as the first point in our simplification
    x_simplified = [x[index0]]
    y_simplified = [y[index0]]
    #plt.plot(x, y, 'k-', marker = 'o', lw = 2)
    #plt.show()
    while (index0 < (end - 1)):
        #first estimate is a straight line between the first and last point of the data
        y_est = (y[index1] - y[index0])*(x - x[index0])/(x[index1] - x[index0]) + y[index0]
        d = np.abs(y - y_est)[index0:index1]
        dmax = np.max(d)
        while (dmax > e)*(index0 < index1):
            #plt.plot(x, y, 'k-', marker = 'o', lw = 2)
            #plt.plot(x[index0:(index1+1)], y_est[index0:(index1+1)], 'b--')
            index1 = index0 + np.where(d == dmax)[0][0] # find the index of distance_max > error
            #plt.plot(x[index1], y[index1], 'ro')
            #plt.show()
            y_est = (y[index1] - y[index0])*(x - x[index0])/(x[index1] - x[index0]) + y[index0] # draw straight line to that point
            d = np.abs(y[index0:index1] - y_est[index0:index1]) # recalculate the differences
            dmax = np.max(d) # find the maximum differences
        x_simplified.append(x[index1])
        y_simplified.append(y[index1])
        index0 = index1
        index1 = end - 1
    return np.array(x_simplified), np.array(y_simplified)

################################################################################
#                                                                              #
# Convert winds from u, v, to speed, direction (and reverse), and along flow   #
#                                                                              #
################################################################################

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

def transform_winds(u, v):
    """
    Computes the flow-relative wind anomalies.
    At each vertical level, computes the along and across flow anomalies (s & n)
    
    These anomalies are  with respect to the domain total horizontal-mean winds
    at that level.
    
    N.B.
        u and v must be regridded to the same horizontal grid points BEFORE they
        are passed to this function!
    
    IN:
        u = u[nt, nz, ny, nx] -> wind component in x-direction
        v = v[nt, nz, ny, nx] -> wind component in y-direction
    OUT:
        s = s[nt, nz, ny, nx] -> wind component in direction of mean flow
        n = n[nt, nz, ny, nx] -> wind component across direction of mean flow
    """
    # Init arrays for s and n
    s = np.zeros_like(u) # along flow
    n = np.zeros_like(u) # cross flow
    
    for it in xrange(u.shape[0]):
        
        for k in xrange(u.shape[1]):
            u_slice = u[it,k,:,:]
            v_slice = v[it,k,:,:]
            # Find the mean wind direction
            wind_speed = np.sqrt(u_slice**2 + v_slice**2)
            wind_direction = np.where((v_slice >= 0), (np.arctan(u_slice/v_slice)*180./np.pi-180.)%360., (np.arctan(u_slice/v_slice)*180./np.pi)%360.)
            mean_dir = np.where((v_slice.mean() >= 0), (np.arctan(u_slice.mean()/v_slice.mean())*180./np.pi-180.)%360., (np.arctan(u_slice.mean()/v_slice.mean())*180./np.pi)%360.)
            # Rotate the coordinates by subtracting mean_dir
            rotated_dir = wind_direction - mean_dir
            n[it,k,:,:] = wind_speed*np.sin(rotated_dir*np.pi/180.) # cross-wind
            s[it,k,:,:] = wind_speed*np.cos(rotated_dir*np.pi/180.) # along-wind
            
    return s, n
################################################################################
#                                                                              #
# Other functions                                                              #
#                                                                              #
################################################################################
def send_email(message, subject, attachments, isAttach = True):
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.image import MIMEImage
    from email.mime.multipart import MIMEMultipart
    
    # Create the container (outer) email message.
    msg = MIMEMultipart()
    msg['Subject'] = subject
    
    me = 'm.c.johnston@pgr.reading.ac.uk'
    
    msg['From'] = me
    msg['To'] = me
    msg.preamble = '\n'
    
    if isAttach:
        # Assume we know that the image files are all in PNG format
        if type(attachments) != list:
            attachments = [attachments]
        for my_file in attachments:
            # Open the files in binary mode.  Let the MIMEImage class automatically
            # guess the specific image type.
            fp = open(my_file, 'rb')
            img = MIMEImage(fp.read())
            fp.close()
            msg.attach(img)
    
    # Create a text/plain message
    body = MIMEText(message) # convert the body to a MIME compatible string
    msg.attach(body) # attach it to your main message
    
    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('localhost')
    s.sendmail(me, me, msg.as_string())
    s.quit()
    
    return None

def summary(aIN):
    """
    Function to print some basic statistics for an input array.
    Requires numpy, prints max, 95th, 75th, 50th, mean 25th, 5th percentiles, and min.
    """
    
    print "\n===================="
    print "\n  MAX: "+str(np.nanmax(aIN))
    print "\n  95th: " + str(np.percentile(aIN, 95))
    print "\n  75th: " + str(np.percentile(aIN, 75))
    print "\n  50th: " + str(np.percentile(aIN, 50))
    print "\n  MEAN: " + str(np.nanmean(aIN))
    print "\n  25th: " + str(np.percentile(aIN, 25))
    print "\n  5th: " + str(np.percentile(aIN, 5))
    print "\n  MIN: " + str(np.nanmin(aIN))
    print "\n===================="

